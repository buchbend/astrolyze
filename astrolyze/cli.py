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

import astrolyze

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="astrolyze — radio/sub-mm PPV cubes and spectra, the house way.",
)
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


if __name__ == "__main__":  # pragma: no cover
    app()
