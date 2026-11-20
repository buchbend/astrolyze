"""Tests for ``Cube.resample_to_beam_sampling`` — a fixed pixels-per-beam sky grid (issue #94).

Written first (red/green TDD). Contract: the new celestial pixel is the beam MINOR-axis FWHM
divided by ``pixels_per_beam``; the spectral axis is unchanged; a finer-than-native pixel raises
:class:`LossyDirectionError`; integrated flux density is conserved; a passed ``NoiseModel`` comes
back on the new grid with σ per pixel unchanged (point-resampling a beam-correlated field).
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_scales

from astrolyze.core import Cube, LossyDirectionError, NoiseQuality
from astrolyze.core import noise as noise_mod
from astrolyze.io import load

warnings.filterwarnings("ignore", module="spectral_cube")
warnings.filterwarnings("ignore", module="astropy")

BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
PIXEL = (
    0.72 * u.arcsec
)  # native pixel: ~14 px per minor-axis FWHM (strongly oversampled)
NZ = 6


def _write_cube(path, data, bunit="MJy/sr"):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -PIXEL.to_value(u.deg), 20.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = PIXEL.to_value(u.deg), 20.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = "NGC0628", "ALMA", bunit
    h["RESTFRQ"] = (230.538e9, "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    fits.writeto(path, np.asarray(data, dtype="float32"), h)
    return Cube.from_loaded(load(path))


def _pixel_scales(cube) -> u.Quantity:
    return proj_plane_pixel_scales(cube._sc.wcs.celestial) * u.deg


def _flux(cube) -> np.ndarray:
    """Per-channel integrated flux density (Jy): Σ I_ν × pixel solid angle."""
    sx, sy = _pixel_scales(cube)
    omega = (sx * sy).to(u.sr)
    data = cube._data_quantity.to(u.MJy / u.sr)
    return (np.nansum(data, axis=(1, 2)) * omega).to_value(u.Jy)


@pytest.fixture
def gaussian_cube(tmp_path):
    """A beam-sized Gaussian source of known flux, centred on a 96x96 native grid."""
    ny = nx = 96
    y, x = np.mgrid[:ny, :nx]
    s_major = (BEAM.major / PIXEL).decompose().value / 2.3548
    s_minor = (BEAM.minor / PIXEL).decompose().value / 2.3548
    blob = np.exp(-0.5 * (((x - 47.3) / s_minor) ** 2 + ((y - 48.6) / s_major) ** 2))
    data = np.broadcast_to(blob, (NZ, ny, nx)) * np.arange(1, NZ + 1)[:, None, None]
    return _write_cube(tmp_path / "gauss.fits", data)


def test_pixel_size_is_minor_fwhm_over_n_and_spectral_axis_unchanged(gaussian_cube):
    out = gaussian_cube.resample_to_beam_sampling(3.0)
    expected = (BEAM.minor / 3.0).to(u.deg)
    np.testing.assert_allclose(
        _pixel_scales(out).to_value(u.deg), expected.value, rtol=1e-6
    )
    np.testing.assert_array_equal(
        out._sc.spectral_axis.value, gaussian_cube._sc.spectral_axis.value
    )
    assert out.shape[0] == gaussian_cube.shape[0]
    assert out.shape[1] < gaussian_cube.shape[1]  # coarser grid, same footprint
    assert out.metadata.beam == gaussian_cube.metadata.beam


def test_finer_than_native_raises(gaussian_cube):
    # 10"/20 = 0.5" < 0.72" native pixel: would invent sampling (ADR-0003).
    with pytest.raises(LossyDirectionError):
        gaussian_cube.resample_to_beam_sampling(20.0)


def test_gaussian_flux_conserved_within_1_percent(gaussian_cube):
    out = gaussian_cube.resample_to_beam_sampling(3.0)
    ratio = _flux(out) / _flux(gaussian_cube)
    print("flux ratio per channel:", ratio)
    np.testing.assert_allclose(ratio, 1.0, rtol=0.01)


def test_noise_model_regridded_sigma_unchanged(tmp_path):
    # White noise convolved with the beam, normalised to σ = 1 (beam-correlated field).
    field = noise_mod.synthesize_correlated_noise(
        (64, 128, 128), beam=BEAM, pixel_scale=PIXEL, rng=np.random.default_rng(94)
    )
    cube = _write_cube(tmp_path / "noise.fits", field)
    model = cube.estimate_noise()

    out, out_model = cube.resample_to_beam_sampling(3.0, noise=model)

    assert out_model.quality is NoiseQuality.PROPAGATED
    assert out_model.sigma_map.shape == out.shape[1:]
    sigma_in = np.nanmedian(model.sigma_map.data.value)
    sigma_prop = np.nanmedian(out_model.sigma_map.data.value)
    # The physics oracle: σ re-measured on the resampled noise itself.
    sigma_meas = np.nanmedian(out.estimate_noise().sigma_map.data.value)
    print(
        f"σ in {sigma_in:.4f}, propagated {sigma_prop:.4f}, re-measured {sigma_meas:.4f}"
    )
    assert sigma_prop == pytest.approx(sigma_in, rel=0.10)
    assert sigma_meas == pytest.approx(sigma_in, rel=0.10)
