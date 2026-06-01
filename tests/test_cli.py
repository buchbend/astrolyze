"""Tests for astrolyze.cli — the typer/rich CLI (issue #7, ADR-0011/0012).

The CLI is a *first-class deliverable*: a human (or an agent, identically) drives the same
load -> moment0 -> convert -> plot spine from the shell. These are behavioral/contract tests
over the public command surface, run against synthetic in-memory cubes (no large real files):

- ``info`` reads the metadata schema and reports completeness; it never raises on an
  incomplete header (that is the io contract — ADR-0006 ii), it shows the gaps;
- ``moment0`` runs the tracer (load -> moment0 -> ``.to(unit)`` -> plot) and writes a figure;
- a conversion needing context the header does not carry makes the command exit non-zero with
  a clear message — the CLI surfaces the library's refusal, it does not paper over it
  (no silent physics, ADR-0003);
- the CLI layer owns its output (rich + ``typer.Exit``); the *library* still contains no
  ``print`` / ``SystemExit`` (asserted in the per-module suites).
"""

import warnings

import matplotlib

matplotlib.use("Agg")  # headless: the CLI saves figures, never opens a window

import numpy as np
import pytest
from astropy.io import fits
from typer.testing import CliRunner

from astrolyze.cli import app

warnings.filterwarnings("ignore", module="spectral_cube")

runner = CliRunner()

REST_HZ = 230.538e9  # CO(2-1)


def _all_output(result) -> str:
    """stdout plus stderr (click 8.2+ captures them separately) — robust for assertions."""
    out = result.stdout or ""
    try:
        out += result.stderr or ""
    except ValueError:  # stderr was not captured separately
        pass
    return out


def _base_header(bunit):
    """A 3D cube header on the radio velocity convention with a real beam."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = bunit
    h["BMAJ"], h["BMIN"], h["BPA"] = 3.33e-4, 2.78e-4, 30.0  # ~1.2"x1.0"
    return h


def _write(path, header, shape=(5, 6, 6)):
    """A small synthetic cube with a gradient so the moment-0 image has real range."""
    data = np.arange(np.prod(shape), dtype="float32").reshape(shape)
    fits.writeto(path, data, header)
    return path


@pytest.fixture
def complete_cube(tmp_path):
    """A schema-complete brightness-temperature cube (K, rest freq + radio convention)."""
    h = _base_header("K")
    h["RESTFRQ"] = (REST_HZ, "Hz")
    return _write(tmp_path / "ngc0628_co21.fits", h)


@pytest.fixture
def incomplete_cube(tmp_path):
    """A Jy/beam cube missing the rest frequency: holding it is fine (lazy), but a
    Jy/beam -> K conversion needs the rest frequency it does not carry."""
    h = _base_header("Jy/beam")  # no RESTFRQ written
    return _write(tmp_path / "archival.fits", h)


# --------------------------------------------------------------------------------------
# info: read the schema, report completeness (read-only, never raises)
# --------------------------------------------------------------------------------------
def test_info_reports_object_and_completeness(complete_cube):
    result = runner.invoke(app, ["info", str(complete_cube)])
    assert result.exit_code == 0
    out = _all_output(result)
    assert "NGC0628" in out
    assert "ALMA" in out
    # the completeness verdict is surfaced (the header carries the mandatory context).
    assert "complete" in out.lower()


def test_info_flags_incomplete_header_without_raising(incomplete_cube):
    result = runner.invoke(app, ["info", str(incomplete_cube)])
    assert result.exit_code == 0  # incomplete loads (ADR-0006 ii); info never raises
    out = _all_output(result)
    assert "incomplete" in out.lower()  # the file is flagged, not silently accepted
    assert "rest_frequency" in out  # and the missing mandatory field is named


# --------------------------------------------------------------------------------------
# moment0: the tracer — load -> moment0 -> .to(unit) -> plot -> save figure
# --------------------------------------------------------------------------------------
def test_moment0_writes_figure(complete_cube, tmp_path):
    out_png = tmp_path / "mom0.png"
    result = runner.invoke(
        app, ["moment0", str(complete_cube), "-u", "K km/s", "-o", str(out_png)]
    )
    assert result.exit_code == 0, _all_output(result)
    assert out_png.exists() and out_png.stat().st_size > 0


def test_moment0_default_output_path_next_to_input(complete_cube):
    result = runner.invoke(app, ["moment0", str(complete_cube), "-u", "K km/s"])
    assert result.exit_code == 0, _all_output(result)
    default = complete_cube.with_name(complete_cube.stem + "_moment0.png")
    assert default.exists()


def test_moment0_missing_context_exits_nonzero_with_clear_message(incomplete_cube):
    # Jy/beam -> K km/s needs the rest frequency the header lacks: the library refuses and
    # the CLI surfaces that as a non-zero exit, not a traceback or a silent guess.
    result = runner.invoke(app, ["moment0", str(incomplete_cube), "-u", "K km/s"])
    assert result.exit_code != 0
    out = _all_output(result).lower()
    assert "rest_frequency" in out or "context" in out


# --------------------------------------------------------------------------------------
# version
# --------------------------------------------------------------------------------------
def test_version_flag():
    import astrolyze

    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert astrolyze.__version__ in _all_output(result)
