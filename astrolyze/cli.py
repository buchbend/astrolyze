"""The astrolyze command-line interface (typer + rich) — ADR-0011/0012.

astrolyze is AI-independent: a human, or a coding agent *identically*, drives the same
``load -> moment0 -> convert -> plot`` spine from the shell, producing ordinary, reviewable
commands rather than reaching for a second AI-only interface (ADR-0012, no MCP).

Two commands cover the tracer-bullet acceptance path (PRD 0001):

- ``astrolyze info PATH`` — parse and show the metadata schema of a FITS file and whether it
  is *complete* (read-only; an incomplete header still opens — ADR-0006 ii — and the gaps are
  named, never guessed).
- ``astrolyze moment0 PATH`` — velocity-integrate a cube, convert the result to a target unit,
  and save a house-style figure (beam on WCS axes, cividis).

The CLI layer **owns output**: it uses ``rich`` for tables/messages and ``typer.Exit`` for
exit codes. The no-``print``/no-``SystemExit`` rule applies to the *library*, not here — the
CLI is exactly where user-facing output belongs. Heavy imports (matplotlib via viz,
spectral-cube via core) are deferred into the command bodies so ``astrolyze --help`` stays
fast and ``import astrolyze`` remains side-effect free.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import astropy.units as u
import typer
from rich.console import Console
from rich.table import Table
from rich.tree import Tree

import astrolyze

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="astrolyze — radio/sub-mm PPV cubes and spectra, the house way.",
)
# `astrolyze manifest …` — read-only inspection of an experiment's dataset registry. The
# registry is *generated* by `ingest`; this group only queries it (no hand-editing path).
manifest_app = typer.Typer(
    no_args_is_help=True, help="Inspect an experiment's dataset manifest."
)
app.add_typer(manifest_app, name="manifest")
console = Console()
err_console = Console(stderr=True)

# Fields shown by ``info``, in schema-declaration order, each with a display formatter.
_SCHEMA_ROWS = (
    ("object", str),
    ("telescope", str),
    ("species", str),
    ("rest_frequency", lambda q: f"{q.to(u.GHz):.6g}"),
    ("velocity_convention", lambda c: c.value),
    ("beam", lambda b: _format_beam(b)),
    ("bunit", str),
    ("distance", str),
    ("calibration_error", str),
    ("name_tag", str),
)


def _format_beam(beam) -> str:
    """``1.20" x 1.00" @ 30.0 deg`` from a radio_beam.Beam."""
    maj = beam.major.to(u.arcsec).value
    min_ = beam.minor.to(u.arcsec).value
    pa = beam.pa.to(u.deg).value
    return f'{maj:.2f}" x {min_:.2f}" @ {pa:.1f} deg'


def _version_callback(value: bool) -> None:
    if value:
        # highlight=False: keep the version a plain, parseable string (rich would otherwise
        # syntax-colour the numbers).
        console.print(f"astrolyze {astrolyze.__version__}", highlight=False)
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the astrolyze version and exit.",
    ),
) -> None:
    """astrolyze command-line interface."""


@app.command()
def init(
    directory: Path = typer.Argument(
        ...,
        file_okay=False,
        help="Experiment directory to scaffold (created if it does not exist).",
    ),
) -> None:
    """Scaffold the fixed experiment skeleton at DIRECTORY (ADR-0009).

    Creates ``data/{raw,interim,processed}``, ``outputs/{figures,tables}``, ``logs/`` and a
    default ``config.toml``. Idempotent: re-running on an existing experiment adds nothing and
    never overwrites a ``config.toml`` you have edited. ``raw/`` is sacred — astrolyze never
    writes to or renames anything inside it.
    """
    from astrolyze.experiment import Experiment  # deferred (pulls dynaconf)

    existed = directory.exists()
    experiment = Experiment.init(directory)

    tree = Tree(f"[bold]{experiment.root}[/bold]")
    data = tree.add("data/")
    data.add("raw/        [dim]immutable inputs — sacred[/dim]")
    data.add("interim/    [dim]derived intermediates[/dim]")
    data.add("processed/  [dim]analysis-ready products[/dim]")
    outputs = tree.add("outputs/")
    outputs.add("figures/    [dim]house-style plots[/dim]")
    outputs.add("tables/     [dim]ECSV/CSV results[/dim]")
    tree.add("logs/         [dim]experiment run log(s)[/dim]")
    tree.add(f"{experiment.config.name}   [dim]dynaconf settings[/dim]")
    console.print(tree)

    verb = "ready" if existed else "created"
    console.print(f"[green]{verb}[/green] experiment at {experiment.root}")


@app.command()
def ingest(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="Experiment directory to ingest (its data/raw/ is scanned).",
    ),
) -> None:
    """Validate every dataset in DIRECTORY's ``data/raw/`` and register what is complete.

    The merciless gate (ADR-0009): each raw file's header is checked against the metadata
    schema. Files carrying the mandatory physical context (rest frequency, velocity
    convention) are *accepted* and recorded in the manifest; files missing it are *rejected*,
    with the exact missing fields named so you know what to fix. Rejected files are never
    registered, and ``raw/`` is never modified. Fix a header and re-run — re-ingest updates the
    same row, it does not duplicate.
    """
    # deferred (pulls dynaconf / SQLAlchemy)
    from astrolyze.experiment import Experiment
    from astrolyze.experiment import ingest as run_ingest

    experiment = Experiment(directory)
    if not experiment.raw.is_dir():
        err_console.print(
            f"[red]error:[/red] {directory} is not an astrolyze experiment "
            "(no data/raw/). Run `astrolyze init` first."
        )
        raise typer.Exit(code=1)

    report = run_ingest(experiment)

    if report.accepted:
        accepted = Table(title="accepted — registered in the manifest")
        accepted.add_column("dataset", style="bold")
        accepted.add_column("object")
        accepted.add_column("species")
        for item in report.accepted:
            m = item.record.metadata
            accepted.add_row(item.source_path, m.object or "—", m.species or "—")
        console.print(accepted)

    if report.rejected:
        rejected = Table(title="rejected — not registered")
        rejected.add_column("dataset", style="bold")
        rejected.add_column("reason")
        for item in report.rejected:
            reason = (
                f"unreadable: {item.error}"
                if item.error is not None
                else "missing mandatory context: " + ", ".join(item.missing)
            )
            rejected.add_row(item.source_path, reason)
        console.print(rejected)

    colour = "green" if report.n_rejected == 0 else "yellow"
    console.print(
        f"[{colour}]ingested[/{colour}] {experiment.root}: "
        f"{report.n_accepted} accepted, {report.n_rejected} rejected"
    )


@app.command()
def info(
    path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        readable=True,
        help="FITS cube or map to inspect.",
    ),
) -> None:
    """Show the astrolyze metadata schema parsed from PATH and whether it is complete.

    Read-only and never raises on a partial header: a real archival file missing its rest
    frequency or velocity convention still opens, and the missing mandatory context is named
    so you know what a later unit/velocity operation will require (ADR-0003/0006).
    """
    from astrolyze import io  # deferred: keeps --help fast

    metadata = io.load(path).metadata

    table = Table(title=f"astrolyze metadata — {path.name}")
    table.add_column("field", style="bold")
    table.add_column("value")
    for name, formatter in _SCHEMA_ROWS:
        value = getattr(metadata, name)
        table.add_row(name, "—" if value is None else formatter(value))
    console.print(table)

    if metadata.is_complete:
        console.print("[green]complete[/green] — has the mandatory physical context.")
    else:
        console.print(
            "[yellow]incomplete[/yellow] — missing mandatory context: "
            f"{', '.join(metadata.missing)} (the file still opens; an operation needing it "
            "will raise rather than guess)."
        )


@app.command()
def moment0(
    path: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        readable=True,
        help="FITS cube to velocity-integrate.",
    ),
    unit: str = typer.Option(
        "K km/s", "--unit", "-u", help="Target unit for the integrated map."
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="PNG to write (default: <input>_moment0.png next to the input).",
    ),
    cmap: str = typer.Option("cividis", "--cmap", help="Matplotlib colormap."),
    temperature_scale: Optional[str] = typer.Option(
        None,
        "--temperature-scale",
        help="rayleigh_jeans | planck — required only for a brightness-temperature "
        "conversion; never defaulted (ADR-0003).",
    ),
    beam: bool = typer.Option(
        True, "--beam/--no-beam", help="Draw the beam ellipse on the figure."
    ),
) -> None:
    """Velocity-integrate PATH (moment 0), convert to UNIT, and save a house-style figure.

    This is the tracer-bullet spine end to end: ``load -> Cube -> moment0() -> .to(unit) ->
    plot()``. The object carries its own beam / rest frequency / velocity convention, so the
    conversion needs no extra context from you — except the genuinely-ambiguous
    Rayleigh-Jeans-vs-Planck scale, which is passed through only if you give it and is never
    assumed.
    """
    from astrolyze import io  # deferred (pulls spectral-cube / matplotlib)
    from astrolyze.core import Cube
    from astrolyze.units import MissingContextError

    out_path = output or path.with_name(f"{path.stem}_moment0.png")
    scale_kwargs = {"temperature_scale": temperature_scale} if temperature_scale else {}
    try:
        cube = Cube.from_loaded(io.load(path))
        converted = cube.moment0().to(unit, **scale_kwargs)
        fig, _ = converted.plot(cmap=cmap, add_beam=beam)
    except (MissingContextError, u.UnitConversionError) as exc:
        err_console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1)

    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    console.print(f"[green]wrote[/green] {out_path}  ({converted.unit})")


@manifest_app.command("list")
def manifest_list(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="Experiment directory whose manifest to list.",
    ),
    obj: Optional[str] = typer.Option(
        None, "--object", "-O", help="Filter by object (e.g. NGC0628)."
    ),
    species: Optional[str] = typer.Option(
        None, "--species", "-s", help="Filter by species (e.g. CO21)."
    ),
) -> None:
    """List the datasets registered in DIRECTORY's manifest (optionally filtered).

    Read-only: the manifest is generated and kept in sync by ``astrolyze ingest`` — this just
    shows what is registered, so you can find datasets by object or species instead of grepping
    filenames.
    """
    # deferred (dynaconf / SQLAlchemy)
    from astrolyze.experiment import Experiment, Manifest

    manifest = Manifest.for_experiment(Experiment(directory))
    filters = {}
    if obj is not None:
        filters["object"] = obj
    if species is not None:
        filters["species"] = species
    records = manifest.query(**filters) if filters else manifest.all()

    if not records:
        scope = " matching the filter" if filters else ""
        console.print(f"[yellow]no datasets registered{scope}[/yellow] in {directory}")
        return

    table = Table(title=f"dataset manifest — {directory}")
    table.add_column("id", justify="right")
    table.add_column("dataset", style="bold")
    table.add_column("object")
    table.add_column("telescope")
    table.add_column("species")
    table.add_column("rest freq")
    table.add_column("doi")
    for record in sorted(records, key=lambda r: r.id):
        m = record.metadata
        # `is not None`, not truthiness: bool() of an astropy Quantity is ambiguous and raises.
        rest = (
            f"{m.rest_frequency.to(u.GHz):.6g}" if m.rest_frequency is not None else "—"
        )
        table.add_row(
            str(record.id),
            record.source_path,
            m.object or "—",
            m.telescope or "—",
            m.species or "—",
            rest,
            record.doi or "—",
        )
    console.print(table)


@app.command()
def narrate(
    directory: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=False,
        help="Experiment directory to write a narrative note for.",
    ),
    run: Optional[str] = typer.Option(
        None,
        "--run",
        help="Run id to document (default: the most recent run).",
    ),
) -> None:
    """Offer a human-readable narrative note for a run (ADR-0010).

    Scaffolds (or re-opens) a markdown note in DIRECTORY's ``logs/``, pre-filled with a
    reference to the run's machine record and the operations it performed, plus empty prose
    sections for the scientific *why / what I found / what it means*. This is an **offer**, never
    a requirement: re-running it never overwrites a note you have started editing, and astrolyze
    never blocks on it. Keep any account you write checkable against the referenced run log.
    """
    from astrolyze.experiment import Experiment  # deferred (pulls dynaconf)
    from astrolyze.experiment import narrate as run_narrate

    experiment = Experiment(directory)
    if not experiment.logs.is_dir():
        err_console.print(
            f"[red]error:[/red] {directory} is not an astrolyze experiment "
            "(no logs/). Run `astrolyze init` first."
        )
        raise typer.Exit(code=1)

    before = set(experiment.logs.glob("*.md"))
    note = run_narrate(experiment, run_id=run)
    verb = "opened existing" if note in before else "scaffolded"

    console.print(f"[green]{verb}[/green] narrative note: {note}")
    console.print(
        "[dim]optional, never required — write the why/what-you-found in your editor; "
        "keep it checkable against the run log it references.[/dim]"
    )


if __name__ == "__main__":  # pragma: no cover
    app()
