"""Tests for the pluggable noise-estimator suite behind ``Cube.estimate_noise(method=…)`` (#28).

Written first (red/green TDD). Issue #27 shipped the :class:`~astrolyze.core.NoiseModel` type,
the ``NoiseQuality`` flag, the separability test, and a single default estimator (``mad_std``).
This file is the correctness obligation for issue #28: the *pluggable estimator suite* that sits
behind the same entry point. The contracts pinned here:

1. **Signal-free-channel RMS (``method="rms"``).** RMS taken over the channels that carry no line
   signal produces a :class:`~astrolyze.core.NoiseModel` flagged
   :data:`~astrolyze.core.NoiseQuality.MEASURED`, and its scalar σ recovers the cube's known σ.

2. **Robust MAD / σ-clip (``method="mad"``).** A σ-clipping robust estimator built on
   :mod:`astropy.stats` likewise produces a ``MEASURED`` model recovering the known σ.

3. **Survey weight / RMS-map ingest.** :meth:`NoiseModel.from_weight_map` builds a model from a
   published weight map (σ = 1/√w), and :meth:`NoiseModel.from_rms_map` from a published RMS map
   (σ = the map), instead of *estimating* σ from the data. Both carry the cube's context.

4. **Pluggable via a registry / protocol.** A caller can register a brand-new estimator with
   :func:`~astrolyze.core.register_estimator` and select it through ``estimate_noise(method=…)``
   *without editing core* — the routing is open, not a closed ``if/elif`` ladder.

5. **No silent physics (ADR-0003).** An estimator handed a cube with no usable signal-free data
   yields :data:`~astrolyze.core.NoiseQuality.UNRELIABLE` and a ``NaN`` σ — never a fabricated
   number. An unknown method is refused, not silently substituted.

6. **Separability (reused from #27).** Each estimator's model selects the separable vs full-3D
   representation via the same separability test, and both reconstruct a full ``.sigma_cube``.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import (
    Cube,
    Map,
    NoiseModel,
    NoiseQuality,
    Spectrum,
    available_estimators,
    register_estimator,
)
from astrolyze.io import load

# spectral-cube emits cosmetic warnings (no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
SIGMA = 0.2  # the true per-voxel noise level (K) of the synthetic cube


def _cube_header(bunit="K"):
    """A 3D cube header on the radio velocity convention with the full astrolyze schema."""
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
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


def _write_cube(path, data, bunit="K"):
    fits.writeto(path, np.asarray(data, dtype="float32"), _cube_header(bunit))
    return Cube.from_loaded(load(path))


@pytest.fixture
def noise_cube(tmp_path):
    """A signal-free Gaussian-noise cube (K): every estimator should recover ~SIGMA."""
    rng = np.random.default_rng(20280101)
    data = rng.normal(0.0, SIGMA, size=(120, 16, 16))
    return _write_cube(tmp_path / "noise_cube.fits", data)


@pytest.fixture
def line_cube(tmp_path):
    """A noise cube with a bright emission line over the central channels.

    The line lives in a contiguous block of channels (a Gaussian profile in v, a blob in the
    plane); the rest are signal-free. A signal-free-channel RMS estimator must measure σ from
    the line-free channels only, so it recovers ~SIGMA despite the bright contaminating line.
    """
    rng = np.random.default_rng(20280102)
    nz, ny, nx = 120, 16, 16
    data = rng.normal(0.0, SIGMA, size=(nz, ny, nx))
    channels = np.arange(nz)
    profile = 5.0 * np.exp(
        -0.5 * ((channels - 60) / 6.0) ** 2
    )  # bright line ~ channel 60
    yy, xx = np.mgrid[0:ny, 0:nx]
    blob = np.exp(-0.5 * (((xx - 8) / 3.0) ** 2 + ((yy - 8) / 3.0) ** 2))
    data += profile[:, None, None] * blob[None, :, :]
    return _write_cube(tmp_path / "line_cube.fits", data)


# --------------------------------------------------------------------------------------
# AC: method="rms" (signal-free channels) -> a MEASURED NoiseModel recovering the known σ
# --------------------------------------------------------------------------------------
def test_rms_estimator_is_measured_and_recovers_sigma(noise_cube):
    model = noise_cube.estimate_noise(method="rms")
    assert isinstance(model, NoiseModel)
    assert model.quality is NoiseQuality.MEASURED
    assert model.method == "rms"
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=0.15)


def test_rms_estimator_ignores_line_channels(line_cube):
    # The RMS estimator must measure σ over the *signal-free* channels only, so the bright line
    # does not inflate the recovered σ above the true noise level.
    model = line_cube.estimate_noise(method="rms")
    assert model.quality is NoiseQuality.MEASURED
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=0.2)


# --------------------------------------------------------------------------------------
# AC: a robust mad/σ-clip method -> a MEASURED NoiseModel recovering the known σ
# --------------------------------------------------------------------------------------
def test_mad_sigma_clip_estimator_is_measured_and_recovers_sigma(noise_cube):
    model = noise_cube.estimate_noise(method="mad")
    assert isinstance(model, NoiseModel)
    assert model.quality is NoiseQuality.MEASURED
    assert model.method == "mad"
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=0.15)


def test_mad_sigma_clip_is_robust_to_a_line(line_cube):
    # σ-clipping rejects the bright line voxels, so the robust σ stays near the true noise level
    # even though the line is present in the data the estimator sees.
    model = line_cube.estimate_noise(method="mad")
    assert model.quality is NoiseQuality.MEASURED
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=0.25)


# --------------------------------------------------------------------------------------
# AC: NoiseModel.from_weight_map / from_rms_map build a model from a survey map
# --------------------------------------------------------------------------------------
def test_from_rms_map_builds_a_model_from_a_survey_rms_map(noise_cube):
    # A published RMS map: σ is the map itself. Build a model on the cube's context from it.
    sigma_map = noise_cube[0]._with_data(np.full((16, 16), SIGMA) * u.K)
    model = NoiseModel.from_rms_map(noise_cube, sigma_map)
    assert isinstance(model, NoiseModel)
    assert model.quality is NoiseQuality.MEASURED
    # σ comes straight from the map, not estimated from the data.
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=1e-6)
    # The spatial σ product is the supplied map (in the cube's unit), context carried.
    np.testing.assert_allclose(
        np.asarray(model.sigma_map.data.to_value(u.K)), SIGMA, rtol=1e-6
    )
    assert u.isclose(model.sigma_map.beam.major, BEAM.major, rtol=1e-12)


def test_from_weight_map_builds_a_model_from_a_survey_weight_map(noise_cube):
    # A published weight map: σ = 1/√w. A uniform weight of 1/SIGMA**2 gives σ = SIGMA.
    weight = noise_cube[0]._with_data(np.full((16, 16), 1.0 / SIGMA**2) / u.K**2)
    model = NoiseModel.from_weight_map(noise_cube, weight)
    assert isinstance(model, NoiseModel)
    assert model.quality is NoiseQuality.MEASURED
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=1e-6)


def test_from_weight_map_zero_weight_is_unreliable_not_fabricated(noise_cube):
    # A zero-weight (un-observed) map carries no σ: refuse to invent one (ADR-0003).
    weight = noise_cube[0]._with_data(np.zeros((16, 16)) / u.K**2)
    model = NoiseModel.from_weight_map(noise_cube, weight)
    assert model.quality is NoiseQuality.UNRELIABLE
    assert not np.isfinite(model.scalar.to_value(u.K))


def test_map_ingest_products_are_first_class(noise_cube):
    sigma_map = noise_cube[0]._with_data(np.full((16, 16), SIGMA) * u.K)
    model = NoiseModel.from_rms_map(noise_cube, sigma_map)
    assert isinstance(model.sigma_cube, Cube)
    assert isinstance(model.sigma_map, Map)
    assert isinstance(model.sigma_spectrum, Spectrum)
    assert model.sigma_cube.shape == noise_cube.shape


# --------------------------------------------------------------------------------------
# AC: estimators are pluggable via a registry/protocol (no core edit needed)
# --------------------------------------------------------------------------------------
def test_a_caller_supplied_estimator_is_usable_without_editing_core(noise_cube):
    # A caller registers its own estimator and selects it by name — the routing is open.
    sentinel = 0.123

    def constant_sigma(data):
        """A toy estimator: report a known constant σ everywhere (MEASURED)."""
        from astrolyze.core.noise import SeparableNoise

        nz = data.shape[0]
        ny, nx = data.shape[1:]
        rep = SeparableNoise(
            sigma_xy=np.full((ny, nx), sentinel),
            sigma_v=np.full(nz, sentinel),
            scalar=sentinel,
            acf=np.array([1.0]),
        )
        return rep, NoiseQuality.MEASURED

    register_estimator("constant-test", constant_sigma)
    assert "constant-test" in available_estimators()

    model = noise_cube.estimate_noise(method="constant-test")
    assert model.method == "constant-test"
    assert model.quality is NoiseQuality.MEASURED
    assert model.scalar.to_value(u.K) == pytest.approx(sentinel, rel=1e-6)


def test_builtin_methods_are_registered(noise_cube):
    names = available_estimators()
    for method in ("mad_std", "rms", "mad"):
        assert method in names


# --------------------------------------------------------------------------------------
# AC: separability test selects separable vs full-3D for the new estimators too
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("method", ["rms", "mad", "mad_std"])
def test_each_estimator_uses_separable_form_on_stationary_noise(noise_cube, method):
    model = noise_cube.estimate_noise(method=method)
    assert model.is_separable is True
    assert model.sigma_cube.shape == noise_cube.shape


# --------------------------------------------------------------------------------------
# AC: no usable signal-free data -> UNRELIABLE, never fabricated (every estimator)
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("method", ["rms", "mad", "mad_std"])
def test_no_usable_signal_free_data_is_unreliable(tmp_path, method):
    blank = np.full((10, 4, 4), np.nan)
    cube = _write_cube(tmp_path / "blank.fits", blank)
    model = cube.estimate_noise(method=method)
    assert model.quality is NoiseQuality.UNRELIABLE
    assert not np.isfinite(model.scalar.to_value(u.K))


# --------------------------------------------------------------------------------------
# House rule: unknown method refused, not silently substituted (ADR-0003)
# --------------------------------------------------------------------------------------
def test_unknown_method_is_refused(noise_cube):
    with pytest.raises(ValueError, match="method"):
        noise_cube.estimate_noise(method="does-not-exist")
