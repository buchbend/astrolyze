"""The tracer-bullet path — PRD 0001 acceptance (issue #7).

One thin vertical slice exercised end-to-end through the public seams:

    Cube.from_loaded(io.load(path)).moment0().to("K km/s").plot()

Three flavours, in increasing realism:

- ``test_tracer_spine_synthetic`` — the whole spine on a tiny in-memory cube (fast, CI).
- ``test_tracer_real_cutout`` — the same spine on a committed 128x128x50 cutout of the real
  PHANGS-ALMA NGC 628 CO(2-1) cube (the always-on *real-data* smoke; see
  ``tests/data/PROVENANCE.md``).
- ``test_tracer_full_cube`` — the spine on the full ~958 MB cube, run only when
  ``$ASTROLYZE_TRACER_CUBE`` points at it (the literal PRD acceptance, off by default).

For a brightness-temperature (K) cube, ``moment0`` yields K·(m/s) and ``.to("K km/s")`` is a
pure metric rescale — no rest frequency or RJ/Planck scale is needed, so the spine reads
cleanly. The figure carries the house display: a beam ellipse on WCS axes.
"""

import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless

import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from matplotlib.patches import Ellipse

from astrolyze.core import Cube
from astrolyze.io import load

warnings.filterwarnings("ignore", module="spectral_cube")

REST_HZ = 230.538e9  # CO(2-1)
CUTOUT = Path(__file__).parent / "data" / "ngc0628_co21_cutout.fits.gz"


def _synthetic_k_cube(path):
    """A small schema-complete brightness-temperature cube (K, radio convention)."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = "NGC0628", "ALMA", "K"
    h["RESTFRQ"] = (REST_HZ, "Hz")
    h["BMAJ"], h["BMIN"], h["BPA"] = 3.33e-4, 2.78e-4, 30.0
    data = np.arange(5 * 6 * 6, dtype="float32").reshape((5, 6, 6))
    fits.writeto(path, data, h)
    return path


def _run_spine(path):
    """load -> Cube -> moment0 -> .to('K km/s') -> plot. Returns (map, fig, ax)."""
    cube = Cube.from_loaded(load(path))
    mom0 = cube.moment0().to("K km/s")
    fig, ax = mom0.plot()
    return mom0, fig, ax


def _assert_house_figure(mom0, ax):
    assert mom0.unit.is_equivalent(u.K * u.km / u.s)
    assert any(isinstance(p, Ellipse) for p in ax.patches)  # the beam is drawn


def test_tracer_spine_synthetic(tmp_path):
    mom0, fig, ax = _run_spine(_synthetic_k_cube(tmp_path / "ngc0628_co21.fits"))
    _assert_house_figure(mom0, ax)


@pytest.mark.skipif(not CUTOUT.exists(), reason="committed cutout fixture missing")
def test_tracer_real_cutout():
    mom0, fig, ax = _run_spine(CUTOUT)
    _assert_house_figure(mom0, ax)
    # real CO emission: the moment-0 map has genuine, finite signal.
    values = u.Quantity(mom0.data).value
    assert np.isfinite(values).all()
    assert np.nanmax(values) > 1.0  # K km/s


@pytest.mark.skipif(
    not os.environ.get("ASTROLYZE_TRACER_CUBE"),
    reason="set ASTROLYZE_TRACER_CUBE to the full NGC 628 cube to run the acceptance",
)
def test_tracer_full_cube():
    mom0, fig, ax = _run_spine(os.environ["ASTROLYZE_TRACER_CUBE"])
    _assert_house_figure(mom0, ax)
    assert np.nanmax(u.Quantity(mom0.data).value) > 1.0
