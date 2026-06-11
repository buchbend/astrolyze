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
# `astrolyze collection …` — read-only inspection of a published corpus (issue #57). Mirrors the
# `manifest` group's shape but over a *collection* (a catalog.parquet over Zarr stores) rather
# than an experiment's registry; both are query-only views, never an editing path.
collection_app = typer.Typer(
    no_args_is_help=True, help="Browse a published corpus (collection)."
)
app.add_typer(collection_app, name="collection")
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


@collection_app.command("list")
def collection_list(
    path: str = typer.Argument(
        ...,
        help="Corpus root to browse — a local directory today, an s3:// URL later "
        "(same call shape; fsspec resolves both).",
    ),
) -> None:
    """List a published corpus object-first: one row per source, its lines/surveys/beam range.

    Opens the corpus at PATH through the read-only :class:`~astrolyze.collection.Collection`
    (fsspec, so a local directory and an object-store URL share this command) and prints the same
    object-first overview the Python ``Collection.list()`` returns — surveys, species, store
    count, and beam range per source. A catalog-less directory is transparently scanned (#61), so
    a bare directory of Zarr stores lists too; an unknown ``catalog_schema_version`` exits non-zero
    with a clear message, the CLI surfacing the library's refusal rather than papering over it
    (ADR-0003).
    """
    # deferred (pulls pyarrow / fsspec) so `astrolyze --help` stays fast.
    from astrolyze.collection import Collection
    from astrolyze.collection.catalog import CatalogSchemaError

    try:
        collection = Collection.open(path)
        summaries = collection.list()
    except FileNotFoundError:
        err_console.print(
            f"[red]error:[/red] no catalog at {path} "
            "(a published corpus carries a catalog.parquet at its root)."
        )
        raise typer.Exit(code=1)
    except CatalogSchemaError as exc:
        err_console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1)

    if not summaries:
        console.print(f"[yellow]no datasets[/yellow] in {path}")
        return

    table = Table(title=f"collection — {path}  (catalog {collection.catalog_version})")
    table.add_column("object", style="bold")
    table.add_column("surveys")
    table.add_column("species")
    table.add_column("stores", justify="right")
    table.add_column("beam range")
    for summary in summaries:
        table.add_row(
            summary.object or "—",
            ", ".join(summary.surveys) or "—",
            ", ".join(summary.species) or "—",
            str(summary.n_stores),
            _format_beam_range(summary.beam_range_arcsec),
        )
    console.print(table)


def _format_beam_range(beam_range) -> str:
    """``1.50"–13.00"`` from a (min, max) major-axis arcsec pair; a single value collapses."""
    lo, hi = beam_range
    if lo is None or hi is None:
        return "—"
    if lo == hi:
        return f'{lo:.2f}"'
    return f'{lo:.2f}"–{hi:.2f}"'


@collection_app.command("describe")
def collection_describe(
    object: str = typer.Argument(
        ..., help="Source name to expand (e.g. NGC3521), as it appears in the catalog."
    ),
    path: str = typer.Argument(
        ...,
        help="Corpus root to browse — a local directory today, an s3:// URL later "
        "(same call shape; fsspec resolves both).",
    ),
    deep: bool = typer.Option(
        False,
        "--deep",
        help="Open each store to add native channel width + velocity coverage "
        "(the catalog does not carry these). Off by default — a shallow describe reads "
        "the catalog only and issues no store reads, staying cheap on a remote corpus.",
    ),
) -> None:
    """Expand OBJECT into per-store physical detail: beam, bunit, rest frequency, provenance.

    The drill-down behind ``collection list`` (PRD #56 user story 3): one row per store of the
    source, with its beam / bunit / rest frequency / transition / survey-telescope-checksum. The
    answer comes from the catalog **row first** — a plain ``describe`` opens nothing — so it stays
    cheap on a remote corpus; ``--deep`` opens each store's spectral axis to add native channel
    width and velocity coverage (the fields the catalog does not carry). An unknown OBJECT exits
    non-zero with the known objects named, the CLI surfacing the library's refusal (ADR-0003).
    """
    # deferred (pulls pyarrow / fsspec) so `astrolyze --help` stays fast.
    from astrolyze.collection import Collection
    from astrolyze.collection.catalog import CatalogSchemaError

    try:
        collection = Collection.open(path)
        details = collection.describe(object, deep=deep)
    except FileNotFoundError:
        err_console.print(
            f"[red]error:[/red] no catalog at {path} "
            "(a published corpus carries a catalog.parquet at its root)."
        )
        raise typer.Exit(code=1)
    except CatalogSchemaError as exc:
        err_console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1)
    except KeyError as exc:
        # describe() raises KeyError(message); str() of it carries the quotes, .args[0] is clean.
        err_console.print(f"[red]error:[/red] {exc.args[0]}")
        raise typer.Exit(code=1)

    table = Table(title=f"{object} — {path}  (catalog {collection.catalog_version})")
    table.add_column("survey", style="bold")
    table.add_column("telescope")
    table.add_column("species")
    table.add_column("transition")
    table.add_column("beam")
    table.add_column("bunit")
    if deep:
        table.add_column("chan width")
        table.add_column("velocity range")
    for detail in details:
        row = [
            detail.survey or "—",
            detail.telescope or "—",
            detail.species or "—",
            detail.transition or "—",
            _format_detail_beam(detail),
            detail.bunit or "—",
        ]
        if deep:
            row.append(_format_kms(detail.channel_width_kms))
            row.append(_format_velocity_range(detail))
        table.add_row(*row)
    console.print(table)


def _format_detail_beam(detail) -> str:
    """``13.00" x 11.00" @ 30.0 deg`` from a StoreDetail's beam columns; ``—`` if unstated."""
    maj, min_, pa = (
        detail.beam_major_arcsec,
        detail.beam_minor_arcsec,
        detail.beam_pa_deg,
    )
    if maj is None or min_ is None or pa is None:
        return "—"
    return f'{maj:.2f}" x {min_:.2f}" @ {pa:.1f} deg'


def _format_kms(value) -> str:
    """``2.00 km/s`` from a velocity Quantity value in km/s; ``—`` when not read (shallow)."""
    return "—" if value is None else f"{value:.2f} km/s"


def _format_velocity_range(detail) -> str:
    """``0.00–4.00 km/s`` from a StoreDetail's (min, max) velocity; ``—`` when not read."""
    lo, hi = detail.velocity_min_kms, detail.velocity_max_kms
    if lo is None or hi is None:
        return "—"
    return f"{lo:.2f}–{hi:.2f} km/s"


@collection_app.command("query")
def collection_query(
    path: str = typer.Argument(
        ...,
        help="Corpus root to browse — a local directory today, an s3:// URL later "
        "(same call shape; fsspec resolves both).",
    ),
    obj: Optional[str] = typer.Option(
        None, "--object", "-O", help="Filter by object (e.g. NGC3521)."
    ),
    survey: Optional[str] = typer.Option(
        None, "--survey", help="Filter by survey (e.g. HERACLES)."
    ),
    telescope: Optional[str] = typer.Option(
        None, "--telescope", help="Filter by telescope (e.g. IRAM30M)."
    ),
    species: Optional[str] = typer.Option(
        None, "--species", "-s", help="Filter by species (e.g. CO)."
    ),
    transition: Optional[str] = typer.Option(
        None, "--transition", help="Filter by transition (e.g. 2-1)."
    ),
) -> None:
    """List the corpus stores matching the given metadata filters (one row per store).

    Slices the corpus along the catalog's own axes (object / survey / telescope / species /
    transition, PRD #56 user story 4) without parsing filenames; multiple filters are conjunctive
    (a store must match all). Mirrors the Python ``Collection.query(**filters)``, which returns a
    composable sub-collection — here the flat per-store result is rendered as a rich table. An
    unknown ``catalog_schema_version`` exits non-zero (ADR-0003).
    """
    # deferred (pulls pyarrow / fsspec) so `astrolyze --help` stays fast.
    from astrolyze.collection import Collection
    from astrolyze.collection.catalog import CatalogSchemaError

    filters = {}
    if obj is not None:
        filters["object"] = obj
    if survey is not None:
        filters["survey"] = survey
    if telescope is not None:
        filters["telescope"] = telescope
    if species is not None:
        filters["species"] = species
    if transition is not None:
        filters["transition"] = transition

    try:
        collection = Collection.open(path)
        records = collection.query(**filters).records
    except FileNotFoundError:
        err_console.print(
            f"[red]error:[/red] no catalog at {path} "
            "(a published corpus carries a catalog.parquet at its root)."
        )
        raise typer.Exit(code=1)
    except CatalogSchemaError as exc:
        err_console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1)

    if not records:
        scope = " matching the filter" if filters else ""
        console.print(f"[yellow]no records{scope}[/yellow] in {path}")
        return

    table = Table(
        title=f"collection query — {path}  (catalog {collection.catalog_version})"
    )
    table.add_column("object", style="bold")
    table.add_column("survey")
    table.add_column("telescope")
    table.add_column("species")
    table.add_column("transition")
    table.add_column("store")
    for record in records:
        r = record.row
        table.add_row(
            r.object or "—",
            r.survey or "—",
            r.telescope or "—",
            r.species or "—",
            r.transition or "—",
            record.store_path,
        )
    console.print(table)


@collection_app.command("covering")
def collection_covering(
    path: str = typer.Argument(
        ...,
        help="Corpus root to browse — a local directory today, an s3:// URL later "
        "(same call shape; fsspec resolves both).",
    ),
    ra: Optional[float] = typer.Option(
        None, "--ra", help="Right ascension in degrees (ICRS). Use with --dec."
    ),
    dec: Optional[float] = typer.Option(
        None, "--dec", help="Declination in degrees (ICRS). Use with --ra."
    ),
    coord: Optional[str] = typer.Option(
        None,
        "--coord",
        "-c",
        help="The position as an ICRS coordinate string (e.g. '10h00m -0d02m' or "
        "'170.0 -0.04'), instead of --ra/--dec.",
    ),
) -> None:
    """List the corpus cubes whose footprint COVERS a sky position (one row per covering store).

    The shell mirror of ``Collection.covering(SkyCoord)`` (PRD #56 user stories 5/6/7): give a
    position by ``--ra``/``--dec`` (degrees, ICRS) or a single ``--coord`` string, and it lists
    every cube that covers it, across surveys. The semantics are the library's: the catalog
    footprint (center+radius, plus an exact MOC when mocpy is installed) only *prefilters*
    candidates; the final containment is decided by each candidate store's **own** WCS at open time
    (the catalog is an index, never the geometry authority). A position covered by nothing reports
    "no records" and exits cleanly. An unknown ``catalog_schema_version`` exits non-zero (ADR-0003).
    """
    # deferred (pulls pyarrow / fsspec / astropy.coordinates) so `astrolyze --help` stays fast.
    from astrolyze.collection import Collection
    from astrolyze.collection.catalog import CatalogSchemaError

    position = _parse_position(ra, dec, coord)

    try:
        collection = Collection.open(path)
        records = collection.covering(position).records
    except FileNotFoundError:
        err_console.print(
            f"[red]error:[/red] no catalog at {path} "
            "(a published corpus carries a catalog.parquet at its root)."
        )
        raise typer.Exit(code=1)
    except CatalogSchemaError as exc:
        err_console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=1)

    pos_str = position.to_string("hmsdms")
    if not records:
        console.print(f"[yellow]no records covering[/yellow] {pos_str} in {path}")
        return

    table = Table(
        title=f"collection covering {pos_str} — {path}  "
        f"(catalog {collection.catalog_version})"
    )
    table.add_column("object", style="bold")
    table.add_column("survey")
    table.add_column("telescope")
    table.add_column("species")
    table.add_column("transition")
    table.add_column("store")
    for record in records:
        r = record.row
        table.add_row(
            r.object or "—",
            r.survey or "—",
            r.telescope or "—",
            r.species or "—",
            r.transition or "—",
            record.store_path,
        )
    console.print(table)


def _parse_position(ra, dec, coord):
    """Build the ICRS :class:`~astropy.coordinates.SkyCoord` covering target from the CLI options.

    Accepts either ``--ra``/``--dec`` (degrees) or a single ``--coord`` string, but **not** both
    and **not** neither — a position must be stated explicitly (astrolyze never guesses one,
    ADR-0003). A non-parseable coordinate string exits non-zero with a clear message rather than
    a traceback."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    has_radec = ra is not None or dec is not None
    if has_radec and coord is not None:
        err_console.print(
            "[red]error:[/red] give a position with EITHER --ra/--dec OR --coord, not both."
        )
        raise typer.Exit(code=2)
    if coord is not None:
        try:
            return SkyCoord(coord, unit=(u.deg, u.deg), frame="icrs")
        except Exception as exc:  # noqa: BLE001 — surface a bad coord as a clean CLI error
            err_console.print(
                f"[red]error:[/red] could not parse --coord {coord!r} as a sky position ({exc})."
            )
            raise typer.Exit(code=2)
    if ra is None or dec is None:
        err_console.print(
            "[red]error:[/red] a position is required — give --ra and --dec (degrees) "
            "or --coord (a coordinate string). astrolyze never guesses a position (ADR-0003)."
        )
        raise typer.Exit(code=2)
    return SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")


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


@app.command()
def explore(
    collection: str = typer.Argument(
        ...,
        help="Corpus root to browse in the web explorer — a local directory or an s3:// URL "
        "(same call shape; fsspec resolves both).",
    ),
    host: str = typer.Option(
        "127.0.0.1",
        "--host",
        help="Interface to bind. Defaults to loopback (a local analysis tool); pass 0.0.0.0 "
        "to expose it on the network.",
    ),
    port: int = typer.Option(8000, "--port", help="Port to serve on."),
) -> None:
    """Serve the published corpus at COLLECTION as a browsable web explorer.

    The GUI sibling of the ``astrolyze collection …`` commands (issue #66): starts a local FastAPI
    server whose endpoints are thin wrappers over the same read-only
    :class:`~astrolyze.collection.Collection` API, and serves a single-page app with the
    object-first **list** view and the per-store **detail** view. Open the printed URL in a
    browser. Works on a local corpus and an ``s3://`` URL alike (it rides ``Collection.open``'s
    fsspec layer).

    The web stack (fastapi + uvicorn) is the optional ``astrolyze[web]`` extra — opt-in, not a hard
    dependency. On a bare install this command exits non-zero with a single actionable line
    (``pip install 'astrolyze[web]'``) rather than an opaque traceback (the same opt-in-extra
    contract as ``astrolyze[s3]`` / ``astrolyze[coverage]``).
    """
    # Imported here, not at module top, so `astrolyze --help` and a bare install never touch the
    # web extra: astrolyze.web's own imports of fastapi/uvicorn stay lazy behind require_web_extra.
    from astrolyze.web import WebExtraNotInstalled, serve

    try:
        console.print(
            f"[green]serving[/green] {collection} at "
            f"[bold]http://{host}:{port}[/bold]  (Ctrl-C to stop)"
        )
        serve(collection, host=host, port=port)
    except WebExtraNotInstalled as exc:
        # The hint contains "astrolyze[web]". rich would (a) parse "[web]" as a markup tag and
        # drop it, and (b) syntax-highlight the brackets, splicing ANSI codes mid-string. Disable
        # BOTH markup and highlight so the install command survives verbatim and copy-pastes clean.
        err_console.print("[red]error:[/red]", end=" ")
        err_console.print(str(exc), markup=False, highlight=False)
        raise typer.Exit(code=1)


if __name__ == "__main__":  # pragma: no cover
    app()
