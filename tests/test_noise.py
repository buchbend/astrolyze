"""Tests for the context-carrying ``NoiseModel`` + ``Cube.estimate_noise()`` (issue #27).

Written first (red/green TDD). This file *is* the correctness obligation for adding a
first-class noise representation to astrolyze (ADR-0003/0004/0005). The contracts pinned here:

1. **Products are first-class astrolyze objects (ADR-0004).** ``estimate_noise()`` returns a
   :class:`~astrolyze.core.NoiseModel` whose products are real wrappers: ``.sigma_cube`` is a
   ``Cube``, ``.sigma_map`` a ``Map``, ``.sigma_spectrum`` a ``Spectrum`` — each in the cube's
   intensity unit — and ``.scalar`` a ``Quantity``. The type transitions **reuse** the existing
   machinery (the model holds a parent ``Cube`` and slices/moments through it); astrolyze does
   not re-derive a cube-building path here.

2. **Carries context.** The model carries the cube's beam + rest frequency + velocity
   convention (so the products it builds are themselves context-carrying), and the per-cube
   spectral noise autocorrelation ``.spectral_acf``.

3. **No silent physics (ADR-0003).** A normal estimate (real signal-free data present) is
   flagged ``MEASURED``; a cube with *no usable signal-free data* is flagged ``UNRELIABLE`` —
   the model never fabricates a σ to fill the gap.

4. **Separable vs full-3D representation.** When a separability test passes the model stores
   the compact ``σ_xy · σ_v`` form (+ ACF + scalar + beam); otherwise it falls back to a full
   3D σ-cube. **Both** reconstruct ``.sigma_cube`` to the same per-voxel σ.

5. **Unit hub + viz seam.** ``.to(unit)`` follows the data unit (it converts the noise like any
   intensity, supplying the cube's context); ``.plot()`` routes through :mod:`astrolyze.viz`.

6. **Zarr companion group.** The model round-trips as a *companion group* alongside the #23
   Zarr cube store, preserving the estimator method + a schema version as provenance.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention — here a noise-only cube (pure
Gaussian scatter, no source) so a robust estimator has clean signal-free data to measure.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, Map, NoiseModel, NoiseQuality, Spectrum
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
    """A signal-free Gaussian-noise cube (K): a robust estimator should recover ~SIGMA.

    Stationary scatter (one σ everywhere) is the textbook *separable* case: σ(x,y,v) is a
    constant, which factors trivially into σ_xy · σ_v. Sized large enough that the separability
    test sees the real (flat) spatial shape, not estimation scatter."""
    rng = np.random.default_rng(20270127)
    data = rng.normal(0.0, SIGMA, size=(120, 16, 16))
    return _write_cube(tmp_path / "noise_cube.fits", data)


@pytest.fixture
def structured_noise_cube(tmp_path):
    """A *non-separable* noise cube: per-voxel σ whose spatial shape changes with frequency.

    The lower half of the band has an x-gradient noise pattern, the upper half a y-gradient: no
    single spatial shape times a single spectral shape (a rank-1 σ_xy·σ_v) can reproduce it, so
    this must trip the separability test into the full-3D σ-cube fallback."""
    rng = np.random.default_rng(20270128)
    nz, ny, nx = 120, 16, 16
    yy, xx = np.mgrid[0:ny, 0:nx]
    pattern_a = 0.1 + 0.9 * (xx / nx)  # gradient in x (lower half of the band)
    pattern_b = 0.1 + 0.9 * (yy / ny)  # gradient in y (upper half of the band)
    sigma_field = np.empty((nz, ny, nx))
    sigma_field[: nz // 2] = pattern_a
    sigma_field[nz // 2 :] = pattern_b
    data = rng.normal(0.0, 1.0, size=(nz, ny, nx)) * sigma_field
    return _write_cube(tmp_path / "structured_cube.fits", data)


# --------------------------------------------------------------------------------------
# AC: estimate_noise() returns a NoiseModel carrying context (beam + rest freq + convention)
# --------------------------------------------------------------------------------------
def test_estimate_noise_returns_a_noisemodel_carrying_context(noise_cube):
    from astrolyze.units import VelocityConvention

    model = noise_cube.estimate_noise()
    assert isinstance(model, NoiseModel)
    # The physical context travels onto the companion (ADR-0004): beam + rest freq + convention.
    assert u.isclose(model.beam.major, BEAM.major, rtol=1e-12)
    assert u.isclose(model.rest_frequency, REST, rtol=1e-12)
    assert model.velocity_convention is VelocityConvention.RADIO
    # It records *how* it was produced (the default estimator) for provenance.
    assert isinstance(model.method, str) and model.method


# --------------------------------------------------------------------------------------
# AC: products are first-class objects in the cube's unit; transitions reuse the machinery
# --------------------------------------------------------------------------------------
def test_products_are_first_class_objects_in_the_cubes_unit(noise_cube):
    model = noise_cube.estimate_noise()

    assert isinstance(model.sigma_cube, Cube)
    assert isinstance(model.sigma_map, Map)
    assert isinstance(model.sigma_spectrum, Spectrum)
    assert isinstance(model.scalar, u.Quantity)

    # each product carries the cube's intensity unit (K here) ...
    assert model.sigma_cube.unit.is_equivalent(noise_cube.unit)
    assert model.sigma_map.unit.is_equivalent(noise_cube.unit)
    assert model.sigma_spectrum.unit.is_equivalent(noise_cube.unit)
    assert model.scalar.unit.is_equivalent(noise_cube.unit)

    # ... and the σ-cube matches the parent cube's shape (a full PPV σ field).
    assert model.sigma_cube.shape == noise_cube.shape


def test_products_carry_context_through_the_reused_transitions(noise_cube):
    # The σ map / spectrum are produced by reusing the Cube->Map/Spectrum machinery, so they
    # carry the same context (ADR-0004) rather than being re-derived from scratch.
    model = noise_cube.estimate_noise()
    assert u.isclose(model.sigma_map.beam.major, BEAM.major, rtol=1e-12)
    assert model.sigma_cube.metadata.velocity_convention is (
        noise_cube.metadata.velocity_convention
    )


def test_scalar_recovers_the_true_noise_level(noise_cube):
    # The single robust σ is close to the cube's true scatter (a sanity check on the estimator).
    model = noise_cube.estimate_noise()
    assert model.scalar.to_value(u.K) == pytest.approx(SIGMA, rel=0.15)


# --------------------------------------------------------------------------------------
# AC: .spectral_acf returns the per-cube spectral autocorrelation
# --------------------------------------------------------------------------------------
def test_spectral_acf_is_returned(noise_cube):
    model = noise_cube.estimate_noise()
    acf = model.spectral_acf
    # The ACF is exposed as a Spectrum/array; pull its values either way.
    values = np.asarray(getattr(acf, "flux", acf))
    values = np.asarray(getattr(values, "value", values), dtype=float)
    assert values.ndim == 1
    # White noise: the autocorrelation peaks at zero lag and is normalised there.
    assert np.nanargmax(values) == 0
    assert values[0] == pytest.approx(1.0, rel=1e-6)


# --------------------------------------------------------------------------------------
# AC: the NoiseQuality flag — MEASURED normally, UNRELIABLE with no usable signal-free data
# --------------------------------------------------------------------------------------
def test_normal_estimate_is_measured(noise_cube):
    model = noise_cube.estimate_noise()
    assert model.quality is NoiseQuality.MEASURED


def test_no_usable_signal_free_data_is_unreliable_not_fabricated(tmp_path):
    # A cube with no finite data at all: there is nothing to measure a σ from. The model must
    # flag UNRELIABLE and refuse to invent a number (ADR-0003), not silently return e.g. 0.
    blank = np.full((10, 4, 4), np.nan)
    cube = _write_cube(tmp_path / "blank.fits", blank)
    model = cube.estimate_noise()
    assert model.quality is NoiseQuality.UNRELIABLE
    # No fabricated σ: the scalar is not a finite invented value.
    assert not np.isfinite(model.scalar.to_value(u.K))


# --------------------------------------------------------------------------------------
# AC: separable vs full-3D representation; both reconstruct .sigma_cube
# --------------------------------------------------------------------------------------
def test_stationary_noise_uses_the_separable_representation(noise_cube):
    model = noise_cube.estimate_noise()
    assert model.is_separable is True
    # The compact form stores σ_xy, σ_v, the scalar, the ACF and the beam (not a full σ-cube).
    assert isinstance(model.sigma_map, Map)
    assert isinstance(model.sigma_spectrum, Spectrum)


def test_structured_noise_falls_back_to_full_sigma_cube(structured_noise_cube):
    model = structured_noise_cube.estimate_noise()
    assert model.is_separable is False
    # The fallback still reconstructs a full σ-cube of the parent shape.
    assert model.sigma_cube.shape == structured_noise_cube.shape


def test_both_representations_reconstruct_sigma_cube(noise_cube, structured_noise_cube):
    # Whichever branch is taken, .sigma_cube is a full per-voxel σ field of the parent shape.
    for cube in (noise_cube, structured_noise_cube):
        model = cube.estimate_noise()
        sc = model.sigma_cube
        assert isinstance(sc, Cube)
        assert sc.shape == cube.shape
        # σ is non-negative everywhere it is finite (a scatter, never negative).
        values = np.asarray(sc._data_quantity.value)
        finite = np.isfinite(values)
        assert np.all(values[finite] >= 0.0)


def test_separable_sigma_cube_is_the_outer_product(noise_cube):
    # For the separable case, the reconstructed σ-cube equals σ_xy(x,y)·σ_v(v)/scalar — i.e. the
    # rank-1 outer product the compact form stores (not stored densely).
    model = noise_cube.estimate_noise()
    sigma_cube = np.asarray(model.sigma_cube._data_quantity.to_value(u.K))
    sxy = np.asarray(model.sigma_map.data.to_value(u.K))
    sv = np.asarray(model.sigma_spectrum.flux.to_value(u.K))
    s0 = model.scalar.to_value(u.K)
    expected = sv[:, None, None] * sxy[None, :, :] / s0
    np.testing.assert_allclose(sigma_cube, expected, rtol=1e-6, atol=1e-8)


# --------------------------------------------------------------------------------------
# AC: .to(unit) follows the data unit
# --------------------------------------------------------------------------------------
def test_to_follows_the_data_unit(tmp_path):
    # A Jy/beam noise cube converts its σ to Jy/sr through the unit hub, supplying the cube's
    # beam context — .to() follows the data exactly as it does on Cube/Map/Spectrum (ADR-0003c).
    rng = np.random.default_rng(7)
    data = rng.normal(0.0, 0.05, size=(20, 6, 6))
    cube = _write_cube(tmp_path / "jybeam.fits", data, bunit="Jy/beam")
    model = cube.estimate_noise()

    converted = model.to(u.Jy / u.sr)
    assert isinstance(converted, NoiseModel)
    assert converted.sigma_cube.unit.is_equivalent(u.Jy / u.sr)
    assert converted.scalar.unit.is_equivalent(u.Jy / u.sr)


# --------------------------------------------------------------------------------------
# AC: .plot() routes through the viz seam (ADR-0005)
# --------------------------------------------------------------------------------------
def test_plot_routes_through_the_viz_seam(noise_cube, monkeypatch):
    import astrolyze.viz as viz

    called = {}

    def fake_plot_noise(model, **kwargs):
        called["model"] = model
        return ("fig", "ax")

    monkeypatch.setattr(viz, "plot_noise", fake_plot_noise, raising=False)
    model = noise_cube.estimate_noise()
    result = model.plot()
    assert result == ("fig", "ax")
    assert called["model"] is model


def test_real_plot_returns_fig_ax(noise_cube):
    import matplotlib

    matplotlib.use("Agg")
    model = noise_cube.estimate_noise()
    fig, ax = model.plot()
    assert fig is not None and ax is not None


# --------------------------------------------------------------------------------------
# AC: Zarr companion-group round-trip preserves method + version provenance
# --------------------------------------------------------------------------------------
def test_zarr_companion_group_roundtrip_preserves_provenance(noise_cube, tmp_path):
    store = noise_cube.to_zarr(tmp_path / "z")
    model = noise_cube.estimate_noise()

    group = model.to_zarr_companion(store)
    assert group.exists()

    back = NoiseModel.from_zarr_companion(store)
    assert isinstance(back, NoiseModel)
    # Provenance survives: the estimator method + a schema version are carried.
    assert back.method == model.method
    assert back.version == model.version
    # The reconstructed σ-cube matches the original to numerical precision.
    np.testing.assert_allclose(
        np.asarray(back.sigma_cube._data_quantity.to_value(u.K)),
        np.asarray(model.sigma_cube._data_quantity.to_value(u.K)),
        rtol=1e-6,
        atol=1e-8,
    )
    # The quality flag survives the round-trip too.
    assert back.quality is model.quality


def test_zarr_companion_roundtrip_carries_context(noise_cube, tmp_path):
    from astrolyze.units import VelocityConvention

    store = noise_cube.to_zarr(tmp_path / "z")
    model = noise_cube.estimate_noise()
    model.to_zarr_companion(store)

    back = NoiseModel.from_zarr_companion(store)
    assert u.isclose(back.beam.major, BEAM.major, rtol=1e-9)
    assert back.velocity_convention is VelocityConvention.RADIO


# --------------------------------------------------------------------------------------
# House rule: an unknown estimator method is refused, not silently substituted (ADR-0003)
# --------------------------------------------------------------------------------------
def test_unknown_method_is_refused(noise_cube):
    with pytest.raises(ValueError, match="method"):
        noise_cube.estimate_noise(method="does-not-exist")
