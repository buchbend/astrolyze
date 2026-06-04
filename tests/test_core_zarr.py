"""Tests for ``Cube.to_zarr`` / ``Cube.from_zarr`` — the Zarr path on the object model
(issue #23, ADR-0004/0006).

Written first (red/green TDD). The obligation: the Zarr backend (issue #23) reaches the object
model through the *same* io seam as FITS, so a ``Cube`` saved to a Zarr store and reloaded
carries the **same physical context** (beam + rest frequency + velocity convention) as the FITS
path — astrolyze adds the context carry and the dispatch, not a storage layer.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.io import load
from astrolyze.units import VelocityConvention

# spectral-cube emits cosmetic warnings (no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


def _cube_header() -> fits.Header:
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


@pytest.fixture
def fits_cube(tmp_path):
    """A complete :class:`Cube` built through the FITS io seam (the eager reference path)."""
    path = tmp_path / "ngc0628_co21.fits"
    data = np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5)
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


# --------------------------------------------------------------------------------------
# AC: Cube.to_zarr / from_zarr produce a Cube carrying the same context as the FITS path
# --------------------------------------------------------------------------------------
def test_to_zarr_then_from_zarr_carries_the_same_context(fits_cube, tmp_path):
    store = fits_cube.to_zarr(tmp_path / "z")
    assert store.exists()

    roundtrip = Cube.from_zarr(store)
    assert isinstance(roundtrip, Cube)
    # The same physical context the FITS path carried (ADR-0004) ...
    assert u.isclose(roundtrip.rest_frequency, fits_cube.rest_frequency, rtol=1e-12)
    assert roundtrip.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(roundtrip.beam.major, fits_cube.beam.major, rtol=1e-12)
    assert u.isclose(roundtrip.beam.minor, fits_cube.beam.minor, rtol=1e-12)
    assert u.isclose(roundtrip.beam.pa, fits_cube.beam.pa, rtol=1e-12)
    assert roundtrip.metadata.object == fits_cube.metadata.object
    # ... and the same spectral axis (the WCS travelled verbatim through the store).
    assert u.allclose(
        roundtrip.spectral_axis.to(u.m / u.s),
        fits_cube.spectral_axis.to(u.m / u.s),
        rtol=1e-9,
    )


def test_from_zarr_cube_supports_context_carrying_transitions(fits_cube, tmp_path):
    store = fits_cube.to_zarr(tmp_path / "z")
    roundtrip = Cube.from_zarr(store)
    # A type transition still carries context off the Zarr-sourced cube (ADR-0004).
    mom0 = roundtrip.moment0()
    assert u.isclose(mom0.beam.major, BEAM.major, rtol=1e-12)
    assert mom0.metadata.velocity_convention is VelocityConvention.RADIO


def test_from_zarr_subcube_slice_returns_a_cube(fits_cube, tmp_path):
    store = fits_cube.to_zarr(tmp_path / "z")
    roundtrip = Cube.from_zarr(store)
    sub = roundtrip[:2]
    assert isinstance(sub, Cube)
    assert sub.shape[0] == 2
    # Context is carried through the slice.
    assert sub.metadata.velocity_convention is VelocityConvention.RADIO
