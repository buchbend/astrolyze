"""Tests for noise propagation through matching + correlated-noise synthesis (issue #32).

Written first (red/green TDD). These tests *are* the correctness obligation for propagating a
:class:`~astrolyze.core.NoiseModel` (#27) through beam/channel matching (#31), plus the
correlated-noise synthesis utility. The contracts pinned here:

1. **Analytic propagation is the fast path (the propagation laws are the value-add).**
   - *spatial* — convolving to a larger beam scales the per-pixel σ by
     ``√(Ω_corr,in / Ω_out)``, where the correlation area is the **beam** solid angle
     (instrument-limited), **not** the pixel area;
   - *spectral* — binning/smoothing reduces σ by ``√M_eff`` with the **effective** number of
     independent samples ``M_eff = M² / Σρ`` read off the stored autocorrelation (NOT the naive
     ``σ/√M`` that assumes white noise).

2. **No silent physics (ADR-0003).** A propagated model carries
   :data:`~astrolyze.core.NoiseQuality.PROPAGATED`; when the stationarity / single-beam
   assumption breaks (a per-channel beam, a non-separable/non-stationary source model) the flag
   tells the truth and degrades to :data:`~astrolyze.core.NoiseQuality.APPROXIMATE`.

3. **Re-estimation is the validation oracle, never the hot path.** ``mad_std`` on signal-free
   voxels (the #27 estimator) *agrees* with the analytic propagation on a synthetic stationary
   cube within tolerance — but it is **never** called inside the matching op (it would fail when
   bright extended emission fills the field of view). The matching op is analytic only.

4. **Correlated-noise synthesis.** White noise convolved by the beam (spatial) and shaped by the
   ACF via FFT (spectral) reproduces the *target* spatial + spectral correlation — the recovered
   beam / ACF come back within tolerance.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, NoiseModel, NoiseQuality
from astrolyze.core import noise as noise_mod
from astrolyze.io import load

# spectral-cube emits cosmetic warnings (no-beam, stokes, big-cube) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")
warnings.filterwarnings("ignore", module="astropy")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
LARGER_BEAM = radio_beam.Beam(major=24 * u.arcsec, minor=20 * u.arcsec, pa=30 * u.deg)
SIGMA = 0.2  # the true per-voxel noise level (K) of the synthetic cube
CHANNEL_WIDTH = 2000.0 * u.m / u.s


def _cube_header(bunit="K"):
    """A 3D cube header on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = (
        CHANNEL_WIDTH.to_value(u.m / u.s),
        1.0,
        "m/s",
    )
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
    """A signal-free, stationary white-noise cube (K): one σ everywhere, uncorrelated.

    The textbook separable case (σ(x,y,v) constant) and the clean validation target for both the
    analytic propagation and the re-estimation oracle."""
    rng = np.random.default_rng(20320127)
    data = rng.normal(0.0, SIGMA, size=(160, 32, 32))
    return _write_cube(tmp_path / "noise_cube.fits", data)


# --------------------------------------------------------------------------------------
# AC: spatial propagation — RMS scales by √(Ω_corr,in / Ω_out) using the BEAM area
# --------------------------------------------------------------------------------------
def test_spatial_propagation_scales_by_beam_solid_angle_ratio():
    # The correlation area is the BEAM solid angle (instrument-limited), not the pixel area.
    # Convolving from beam_in to a larger beam_out lowers the per-pixel σ by
    # √(Ω_in / Ω_out): the noise is averaged over a larger correlated footprint.
    factor = noise_mod._spatial_rms_factor(BEAM, LARGER_BEAM)
    expected = float(
        np.sqrt((BEAM.sr / LARGER_BEAM.sr).to_value(u.dimensionless_unscaled))
    )
    assert factor == pytest.approx(expected, rel=1e-9)
    # a larger output beam means a SMALLER σ (averaging down), so the factor is < 1.
    assert 0.0 < factor < 1.0


def test_beam_area_fallback_uses_beam_not_pixel_area(noise_cube):
    # Without a measured spatial ACF (#96) the scaling depends ONLY on the beam ratio, not on
    # the pixel grid (beam-, not pixel-, limited correlation area).
    model = NoiseModel.from_rms_map(noise_cube, np.full(noise_cube.shape[1:], SIGMA))
    out_cube = noise_cube.convolve_to_beam(LARGER_BEAM)
    propagated = noise_mod.propagate(model, beam_out=LARGER_BEAM)
    ratio = propagated.scalar.to_value(u.K) / model.scalar.to_value(u.K)
    expected = float(
        np.sqrt((BEAM.sr / LARGER_BEAM.sr).to_value(u.dimensionless_unscaled))
    )
    assert ratio == pytest.approx(expected, rel=1e-6)
    # the propagated model carries the new beam as context (ADR-0004).
    assert u.isclose(propagated.beam.major, LARGER_BEAM.major, rtol=1e-9)
    assert out_cube.shape == noise_cube.shape  # spatial smoothing keeps the grid


# --------------------------------------------------------------------------------------
# AC: spectral propagation — σ/√M_eff with M_eff = M²/Σρ from the ACF (not naive σ/√M)
# --------------------------------------------------------------------------------------
def test_m_eff_from_acf_white_noise_equals_m():
    # White noise: ρ is a spike at lag 0 (Σρ = 1 over the kept window), so M_eff == M and the
    # propagation reduces to the naive √M for the uncorrelated case.
    m = 8
    acf = np.zeros(32)
    acf[0] = 1.0
    assert noise_mod._m_eff(acf, m) == pytest.approx(float(m), rel=1e-9)


def test_m_eff_from_acf_correlated_noise_is_less_than_m():
    # Correlated noise (a broad ACF) has FEWER independent samples than channels: M_eff < M, so
    # binning averages down LESS than the naive √M would claim (M_eff = M² / Σρ).
    m = 8
    rho = np.exp(-np.arange(32) / 3.0)  # a slowly-decaying (correlated) ACF
    m_eff = noise_mod._m_eff(rho, m)
    assert m_eff < m
    # the law: M_eff = M² / Σ_{|lag|<M} ρ(lag) with the two-sided sum.
    lags = np.arange(m)
    two_sided = rho[0] + 2.0 * np.sum(rho[1:m] * (1.0 - lags[1:] / m))
    assert m_eff == pytest.approx((m * m) / (m * two_sided), rel=1e-6)


def test_spectral_propagation_uses_m_eff_not_naive_sqrt_m():
    # Build a model whose stored ACF is correlated; binning by M must reduce σ by √M_eff, which
    # is strictly larger (less reduction) than the naive √M for correlated noise.
    rng = np.random.default_rng(7)
    data = rng.normal(0.0, SIGMA, size=(64, 8, 8))
    # impose spectral correlation by a small running mean along the spectral axis.
    kernel = np.array([0.25, 0.5, 0.25])
    data = np.apply_along_axis(lambda c: np.convolve(c, kernel, mode="same"), 0, data)
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        cube = _write_cube(f"{d}/corr.fits", data)
        model = cube.estimate_noise()
        m = 4
        propagated = noise_mod.propagate(model, spectral_factor=m)
        acf = np.asarray(model._rep.acf, dtype=float)
        m_eff = noise_mod._m_eff(acf, m)
        expected = model.scalar.to_value(u.K) / np.sqrt(m_eff)
        naive = model.scalar.to_value(u.K) / np.sqrt(m)
        assert propagated.scalar.to_value(u.K) == pytest.approx(expected, rel=1e-6)
        # correlated noise: M_eff < M, so the propagated σ is LARGER than the naive estimate.
        assert propagated.scalar.to_value(u.K) > naive


# --------------------------------------------------------------------------------------
# AC: spectral binning rebins the noise spectral axis to match the data (#87)
# --------------------------------------------------------------------------------------
def test_spectral_bin_rebins_noise_axis_to_match_data(noise_cube):
    # The σ companion must be paired voxel-for-voxel with the binned data — the bug is that
    # propagate() kept N channels (variance-only rescale), so after spectral_bin(k) the data had
    # N/k channels but the noise still had N. Rebinning the noise spectral axis fixes the pairing.
    model = noise_cube.estimate_noise()
    out_cube, out_model = noise_cube.spectral_bin(2, noise=model)
    assert out_model.sigma_spectrum.flux.size == out_cube.shape[0]
    assert out_model.sigma_cube.shape == out_cube.shape


def test_spectral_bin_pairing_holds_for_non_divisor_factor(noise_cube):
    # 160 is not divisible by 3: the rebin must match spectral-cube's downsample channel count
    # EXACTLY (remainder included), not a naive N // k, or the pairing breaks for odd factors.
    model = noise_cube.estimate_noise()
    out_cube, out_model = noise_cube.spectral_bin(3, noise=model)
    assert out_model.sigma_spectrum.flux.size == out_cube.shape[0]


def test_propagate_spectral_bin_returns_coarse_sigma_v(noise_cube):
    # The direct propagate() API: the new keyword-only spectral_bin rebins the stored σ_v profile
    # onto the coarse grid (160 -> 80). Fails today with TypeError: unexpected keyword 'spectral_bin'.
    model = noise_cube.estimate_noise()
    propagated = noise_mod.propagate(model, spectral_bin=2)
    assert (
        np.asarray(propagated._rep.sigma_v).size
        == np.asarray(model._rep.sigma_v).size // 2
    )
    assert isinstance(propagated, NoiseModel)


def test_spectral_bin_white_noise_recovers_sigma_over_sqrt_k(noise_cube):
    # Physics oracle on white noise: block-averaging k channels lowers σ by √k. Every coarse
    # channel must land at SIGMA/√k.
    model = noise_cube.estimate_noise()
    k = 4
    out_cube, out_model = noise_cube.spectral_bin(k, noise=model)
    sv = out_model.sigma_spectrum.flux.to_value(u.K)
    # the σ_v profile is on the BINNED axis (one σ per coarse channel), not the original 160.
    assert sv.size == out_cube.shape[0]
    assert np.allclose(sv, SIGMA / np.sqrt(k), rtol=0.15)
    # cross-check: re-estimating σ from the binned data (the mad_std oracle) agrees.
    oracle = noise_mod.reestimate(out_cube).sigma_spectrum.flux.to_value(u.K)
    assert np.nanmedian(sv) == pytest.approx(np.nanmedian(oracle), rel=0.15)


def test_spectral_bin_white_acf_stays_white(noise_cube):
    # The coarse ACF is re-derived on the binned grid: binning white noise leaves it ~white, so
    # the ACF is a spike at lag 0 with near-zero side lobes.
    model = noise_cube.estimate_noise()
    propagated = noise_mod.propagate(model, spectral_bin=4)
    acf = np.asarray(propagated._rep.acf, dtype=float)
    assert acf[0] == pytest.approx(1.0, rel=1e-9)
    assert np.all(np.abs(acf[1:5]) < 0.2)


def test_spectral_factor_stays_variance_only_no_rebin(noise_cube):
    # BACKWARD-COMPAT GUARD (passes today, must keep passing): spectral_factor is the
    # smoothing/variance-only path — spectral_smooth_to keeps the channel grid; only the new
    # spectral_bin rebins. The two must not be conflated, so the channel count is UNCHANGED.
    model = noise_cube.estimate_noise()
    propagated = noise_mod.propagate(model, spectral_factor=4)
    assert (
        np.asarray(propagated._rep.sigma_v).size == np.asarray(model._rep.sigma_v).size
    )


# --------------------------------------------------------------------------------------
# AC: the propagated model carries PROPAGATED (→ APPROXIMATE when assumptions break)
# --------------------------------------------------------------------------------------
def test_propagated_model_is_flagged_propagated(noise_cube):
    model = noise_cube.estimate_noise()
    propagated = noise_mod.propagate(model, beam_out=LARGER_BEAM)
    assert isinstance(propagated, NoiseModel)
    assert propagated.quality is NoiseQuality.PROPAGATED


def test_propagation_degrades_to_approximate_for_per_channel_beam(noise_cube):
    # A per-channel beam breaks the single spatial-correlation-kernel assumption: the result is
    # an order-of-magnitude estimate, flagged APPROXIMATE (no silent physics, ADR-0003).
    model = noise_cube.estimate_noise()
    propagated = noise_mod.propagate(model, beam_out=LARGER_BEAM, per_channel_beam=True)
    assert propagated.quality is NoiseQuality.APPROXIMATE


def test_propagation_degrades_to_approximate_for_non_stationary(
    structured_propagated_model,
):
    # Non-separable / non-stationary correlation breaks the stationary-noise assumption: the
    # analytic propagation can still produce a number but must flag it APPROXIMATE.
    assert structured_propagated_model.quality is NoiseQuality.APPROXIMATE


@pytest.fixture
def structured_propagated_model(tmp_path):
    rng = np.random.default_rng(20320128)
    nz, ny, nx = 120, 16, 16
    yy, xx = np.mgrid[0:ny, 0:nx]
    pattern_a = 0.1 + 0.9 * (xx / nx)
    pattern_b = 0.1 + 0.9 * (yy / ny)
    sigma_field = np.empty((nz, ny, nx))
    sigma_field[: nz // 2] = pattern_a
    sigma_field[nz // 2 :] = pattern_b
    data = rng.normal(0.0, 1.0, size=(nz, ny, nx)) * sigma_field
    cube = _write_cube(tmp_path / "structured.fits", data)
    model = cube.estimate_noise()
    assert model.is_separable is False  # the non-stationary case (#27 full-3D fallback)
    return noise_mod.propagate(model, beam_out=LARGER_BEAM)


# --------------------------------------------------------------------------------------
# AC: matching ops (#31) propagate a NoiseModel analytically (opt-in, hot path is analytic)
# --------------------------------------------------------------------------------------
def test_convolve_to_beam_propagates_noise_model(noise_cube):
    model = noise_cube.estimate_noise()
    out_cube, out_model = noise_cube.convolve_to_beam(LARGER_BEAM, noise=model)
    assert isinstance(out_cube, Cube)
    assert isinstance(out_model, NoiseModel)
    assert out_model.quality is NoiseQuality.PROPAGATED
    # the op routes through propagate (the ACF law, #96: white noise here, not beam-correlated)
    expected = noise_mod.propagate(model, beam_out=LARGER_BEAM).scalar.to_value(u.K)
    assert out_model.scalar.to_value(u.K) == pytest.approx(expected, rel=1e-6)
    # the propagated model carries the new, larger beam.
    assert u.isclose(out_model.beam.major, LARGER_BEAM.major, rtol=1e-9)


def test_convolve_to_beam_without_noise_is_unchanged(noise_cube):
    # Opt-in: with no noise model the op returns a plain Cube exactly as before (#31).
    out = noise_cube.convolve_to_beam(LARGER_BEAM)
    assert isinstance(out, Cube)


def test_spectral_bin_propagates_noise_model(noise_cube):
    model = noise_cube.estimate_noise()
    out_cube, out_model = noise_cube.spectral_bin(2, noise=model)
    assert isinstance(out_cube, Cube)
    assert isinstance(out_model, NoiseModel)
    acf = np.asarray(model._rep.acf, dtype=float)
    m_eff = noise_mod._m_eff(acf, 2)
    expected = model.scalar.to_value(u.K) / np.sqrt(m_eff)
    assert out_model.scalar.to_value(u.K) == pytest.approx(expected, rel=1e-6)


# --------------------------------------------------------------------------------------
# AC (#82): velocity-grid interpolation propagates σ (variance-correct linear interp)
# --------------------------------------------------------------------------------------
SIGMA_CONST = (
    0.3  # a known, CONSTANT per-voxel σ (K) for exact interpolation arithmetic
)


def _constant_noise_member(tmp_path, name, nz=3):
    """A small cube on the velocity axis [0,2,4,...] km/s plus a constant-σ (SIGMA_CONST) companion.

    The constant-σ companion is the published-noise case the interpolation propagation is exact for,
    and makes the post-regrid σ hand-computable from the linear-interpolation weights alone."""
    cube = _write_cube(tmp_path / name, np.full((nz, 4, 4), 1.0))
    model = NoiseModel.from_rms_map(cube, np.full(cube.shape[1:], SIGMA_CONST))
    return cube, model


def test_to_velocity_grid_preserves_sigma_on_coincident_grid(tmp_path):
    # Regridding onto a grid that coincides with the native channels is an identity on σ: each
    # output channel draws its full weight from one input channel (no variance redistribution).
    cube, model = _constant_noise_member(tmp_path, "coincident.fits")
    native = cube.velocity_axis()  # [0, 2, 4] km/s
    out_cube, out_model = cube.to_velocity_grid(native, noise=model)
    assert isinstance(out_cube, Cube)
    assert isinstance(out_model, NoiseModel)
    assert out_model.quality is NoiseQuality.PROPAGATED
    sigma = out_model.sigma_cube._data_quantity.to_value(u.K)
    assert np.allclose(sigma, SIGMA_CONST)


def test_to_velocity_grid_propagates_variance_through_midpoint_interpolation(tmp_path):
    # A target channel at the midpoint of two native channels is their linear interpolation (weights
    # 0.5, 0.5); under channel-independent noise its variance is 0.5²σ² + 0.5²σ² = 0.5σ², so the
    # propagated σ is √0.5·σ — NOT the input σ (interpolating σ directly would be wrong, ADR-0003).
    cube, model = _constant_noise_member(tmp_path, "midpoint.fits")
    midpoints = np.array([1.0, 3.0]) * u.km / u.s  # between native [0,2] and [2,4]
    _, out_model = cube.to_velocity_grid(midpoints, noise=model)
    sigma = out_model.sigma_cube._data_quantity.to_value(u.K)
    assert np.allclose(sigma, np.sqrt(0.5) * SIGMA_CONST)


def test_to_velocity_grid_off_coverage_sigma_is_nan(tmp_path):
    # A target channel outside the native coverage is "not observed": its σ is NaN (zero weight in
    # the coadd), exactly as the data regrid fills it with NaN — never an extrapolated value.
    cube, model = _constant_noise_member(tmp_path, "offcov.fits")
    grid = (
        np.array([2.0, 10.0]) * u.km / u.s
    )  # native covers [0,4]; 10 km/s is beyond it
    out_cube, out_model = cube.to_velocity_grid(grid, noise=model)
    sigma = out_model.sigma_cube._data_quantity.to_value(u.K)
    assert np.all(np.isfinite(sigma[0]))  # 2 km/s is in coverage
    assert np.all(np.isnan(sigma[1]))  # 10 km/s is off coverage -> NaN
    assert np.all(np.isnan(out_cube._data_quantity.to_value(u.K)[1]))  # data agrees


def test_to_velocity_grid_without_noise_is_unchanged(tmp_path):
    # Opt-in: with no noise model the regrid returns a plain Cube exactly as before (#65).
    cube, _ = _constant_noise_member(tmp_path, "nonoise.fits")
    out = cube.to_velocity_grid(cube.velocity_axis())
    assert isinstance(out, Cube)


def test_match_to_propagates_both_noise_models(tmp_path):
    rng = np.random.default_rng(11)
    a = _write_cube(tmp_path / "a.fits", rng.normal(0.0, SIGMA, size=(60, 20, 20)))
    b = _write_cube(tmp_path / "b.fits", rng.normal(0.0, SIGMA, size=(60, 20, 20)))
    b = b.convolve_to_beam(LARGER_BEAM)
    model_a = a.estimate_noise()
    model_b = b.estimate_noise()
    (ma, mb), (na, nb) = a.match_to(b, noise=(model_a, model_b))
    assert isinstance(na, NoiseModel) and isinstance(nb, NoiseModel)
    # both models land on the common beam.
    assert u.isclose(na.beam.major, nb.beam.major, rtol=1e-6)
    assert na.quality is NoiseQuality.PROPAGATED


def test_matching_op_does_not_call_reestimation(noise_cube, monkeypatch):
    # The hot path is ANALYTIC only: the matching op must NEVER call the mad_std re-estimator
    # (it fails when bright extended emission fills the FoV). Trip a tripwire if it does.
    called = {"n": 0}
    real_estimate = noise_mod.estimate

    def tripwire(*args, **kwargs):
        called["n"] += 1
        return real_estimate(*args, **kwargs)

    model = noise_cube.estimate_noise()  # one legitimate estimate up front
    monkeypatch.setattr(noise_mod, "estimate", tripwire)
    noise_cube.convolve_to_beam(LARGER_BEAM, noise=model)
    noise_cube.spectral_bin(2, noise=model)
    assert called["n"] == 0  # no re-estimation inside the matching ops


# --------------------------------------------------------------------------------------
# AC: re-estimation is the validation oracle and AGREES with analytic propagation
# --------------------------------------------------------------------------------------
def test_reestimation_oracle_agrees_with_analytic_spectral_propagation(noise_cube):
    # On a synthetic STATIONARY white-noise cube, binning the data and re-estimating σ (the
    # mad_std oracle on the now-binned, still signal-free cube) agrees with the analytic
    # σ/√M_eff propagation. White noise -> M_eff == M, so both predict σ/√M.
    model = noise_cube.estimate_noise()
    factor = 2
    binned = noise_cube.spectral_bin(factor)
    oracle = noise_mod.reestimate(binned)
    analytic = noise_mod.propagate(model, spectral_factor=factor)
    assert oracle.scalar.to_value(u.K) == pytest.approx(
        analytic.scalar.to_value(u.K), rel=0.1
    )
    # the oracle is a real NoiseModel measured from data (MEASURED), not propagated.
    assert oracle.quality is NoiseQuality.MEASURED


def test_reestimation_oracle_agrees_with_analytic_spatial_propagation(tmp_path):
    # Spatial: the √(Ω_in/Ω_out) law is the BEAM-correlation law, so the oracle validates it on a
    # cube whose noise is already correlated on the input beam (the realistic radio case — noise
    # in a Jy/beam map is beam-correlated, not per-pixel white). We synthesize exactly that, then
    # convolve to a larger beam and re-estimate σ; the oracle agrees with the analytic prediction.
    # The oracle is measured on an INTERIOR cutout — convolution NaN-pads the field edges, so a
    # real re-estimation always avoids the edge-contaminated border (as one would on sky data).
    pixel_scale = 2e-4 * u.deg
    field = noise_mod.synthesize_correlated_noise(
        (60, 160, 160),
        beam=BEAM,
        pixel_scale=pixel_scale,
        sigma=SIGMA,
        rng=np.random.default_rng(20320201),
    )
    cube = _write_cube(tmp_path / "beam_corr.fits", field)
    model = cube.estimate_noise()
    smoothed = cube.convolve_to_beam(LARGER_BEAM)
    interior = smoothed[:, 40:120, 40:120]  # away from the convolution boundary
    oracle = noise_mod.reestimate(interior)
    analytic = noise_mod.propagate(model, beam_out=LARGER_BEAM)
    assert oracle.scalar.to_value(u.K) == pytest.approx(
        analytic.scalar.to_value(u.K), rel=0.15
    )


def test_reestimate_is_the_mad_std_estimator(noise_cube):
    # The oracle is exactly the #27 mad_std estimator surfaced under a propagation-validation
    # name; it returns a context-carrying NoiseModel.
    oracle = noise_mod.reestimate(noise_cube)
    assert isinstance(oracle, NoiseModel)
    assert oracle.method == noise_mod.DEFAULT_METHOD


# --------------------------------------------------------------------------------------
# AC: correlated-noise synthesis reproduces the target spatial + spectral correlation
# --------------------------------------------------------------------------------------
def test_synthesis_reproduces_spatial_beam_correlation(noise_cube):
    # White noise convolved by the beam should come back with a spatial autocorrelation whose
    # width matches the target beam: a wider beam => a wider spatial ACF.
    shape = (40, 64, 64)
    pixel_scale = 2e-4 * u.deg  # |CDELT| of the reference header
    realization = noise_mod.synthesize_correlated_noise(
        shape,
        beam=LARGER_BEAM,
        pixel_scale=pixel_scale,
        rng=np.random.default_rng(3),
    )
    assert realization.shape == shape
    # measure the spatial ACF FWHM of a representative channel and compare to the beam.
    recovered = noise_mod._spatial_acf_fwhm_pixels(realization)
    target_fwhm_pix = float(
        (LARGER_BEAM.major / pixel_scale).to_value(u.dimensionless_unscaled)
    )
    # recovered correlation width tracks the beam major axis (within a generous tolerance:
    # convolution of white noise gives an ACF whose width is the beam width).
    assert recovered == pytest.approx(target_fwhm_pix, rel=0.4)


def test_synthesis_reproduces_spectral_acf(noise_cube):
    # White noise shaped by a target ACF via FFT should come back with that spectral ACF.
    shape = (256, 16, 16)
    target_acf = np.exp(-(np.arange(shape[0]) ** 2) / (2.0 * 2.5**2))  # a Gaussian ACF
    realization = noise_mod.synthesize_correlated_noise(
        shape,
        beam=BEAM,
        pixel_scale=2e-4 * u.deg,
        spectral_acf=target_acf,
        rng=np.random.default_rng(5),
    )
    recovered = noise_mod._measure_spectral_acf(realization)
    # the first few lags of the recovered ACF track the target (normalised at lag 0).
    n = 6
    np.testing.assert_allclose(recovered[:n], target_acf[:n], atol=0.15)


def test_synthesis_recovers_target_rms(noise_cube):
    # A synthesis with an explicit target σ comes back at that σ (within tolerance): the utility
    # produces realistic correlated noise AT a target level, not an arbitrary one.
    shape = (60, 48, 48)
    realization = noise_mod.synthesize_correlated_noise(
        shape,
        beam=BEAM,
        pixel_scale=2e-4 * u.deg,
        sigma=SIGMA,
        rng=np.random.default_rng(9),
    )
    from astropy.stats import mad_std

    assert float(mad_std(realization)) == pytest.approx(SIGMA, rel=0.15)


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code (ADR-0004) — propagation/synthesis API.
# --------------------------------------------------------------------------------------
def test_propagation_requires_a_geometry(noise_cube):
    # propagate() with neither a beam nor a spectral factor is a programming error (nothing to
    # propagate through), refused clearly rather than silently returning the input.
    model = noise_cube.estimate_noise()
    with pytest.raises(ValueError):
        noise_mod.propagate(model)


# --------------------------------------------------------------------------------------
# AC (#96): propagate(beam_out=...) through the measured spatial noise ACF
# --------------------------------------------------------------------------------------
# The beam-area law √(Ω_in/Ω_out) assumes the noise is correlated like the beam. Receiver noise
# usually is not (OTF gridding kernel narrower than the beam; interferometer noise missing short
# spacings), so propagate uses σ_out²/σ_in² = Σ_d (k⋆k)(d)·ρ(d) with the measured spatial ACF ρ.
PIX = 2e-4 * u.deg  # the _cube_header pixel scale (0.72")
MID_BEAM = radio_beam.Beam(major=18 * u.arcsec, minor=16 * u.arcsec, pa=30 * u.deg)


def _interior_std_ratio(before, after, pad):
    """σ(after)/σ(before) on the interior (``pad`` px from each edge: convolution edges bias σ)."""
    s = (slice(None), slice(pad, -pad), slice(pad, -pad))
    return float(np.std(after[s]) / np.std(before[s]))


def test_estimate_fills_a_centred_normalised_spatial_acf(noise_cube):
    acf = noise_cube.estimate_noise()._rep.spatial_acf
    assert acf is not None
    assert acf.shape[0] % 2 == 1 and acf.shape[1] % 2 == 1
    assert acf.shape == (63, 63)  # 32x32 map -> lags ±31
    assert acf[31, 31] == pytest.approx(1.0)
    # white noise: no correlation off lag 0
    off = acf.copy()
    off[31, 31] = 0.0
    assert np.max(np.abs(off)) < 0.1


def test_propagate_narrow_correlation_matches_measured_sigma(tmp_path):
    # OTF-like: white noise smoothed by a kernel (FWHM ~3.5 px = 2.5") much narrower than the
    # declared 12"x10" beam. The beam-area law over-predicts σ; the ACF law must hit it.
    from astropy.convolution import Gaussian2DKernel, convolve_fft

    rng = np.random.default_rng(96)
    white = rng.normal(0.0, 1.0, size=(48, 384, 384))
    kernel = Gaussian2DKernel(1.5)
    data = np.stack([convolve_fft(p, kernel, boundary="wrap") for p in white])
    cube = _write_cube(tmp_path / "narrow.fits", data * SIGMA)
    model = cube.estimate_noise()
    out_cube, out_model = cube.convolve_to_beam(MID_BEAM, noise=model)
    measured = _interior_std_ratio(
        cube._data_quantity.value, out_cube._data_quantity.value, pad=40
    )
    predicted = out_model.scalar.to_value(u.K) / model.scalar.to_value(u.K)
    beam_law = noise_mod._spatial_rms_factor(BEAM, MID_BEAM)
    print(
        f"narrow: measured {measured:.4f} ACF law {predicted:.4f} beam law {beam_law:.4f}"
    )
    assert out_model.quality is NoiseQuality.PROPAGATED
    assert predicted == pytest.approx(measured, rel=0.03)
    assert beam_law > 1.5 * measured  # the old law is badly off here


def test_acf_ignores_faint_emission_in_signal_free_channels(tmp_path):
    # Real "signal-free" channels still hold faint extended emission (below the channel-level
    # threshold, so _signal_free_channels keeps them); it inflates ρ at long lags and the ACF law
    # over-predicts σ. Half the channels here carry a faint broad Gaussian that passes the test.
    from astropy.convolution import Gaussian2DKernel, convolve_fft

    rng = np.random.default_rng(962)
    nz, n = 48, 384
    white = rng.normal(0.0, 1.0, size=(nz, n, n))
    noise = np.stack(
        [convolve_fft(p, Gaussian2DKernel(1.5), boundary="wrap") for p in white]
    )
    noise /= noise.std()
    yy, xx = np.mgrid[0:n, 0:n]
    blob = np.exp(-((yy - n / 2) ** 2 + (xx - n / 2) ** 2) / (2 * 60.0**2))
    line = np.arange(nz) % 2 == 1
    amp = np.where(line, rng.uniform(0.3, 0.8, nz), 0.0)
    data = (noise + amp[:, None, None] * blob) * SIGMA
    assert noise_mod._signal_free_channels(data).sum() > nz // 2  # emission passes
    cube = _write_cube(tmp_path / "faint.fits", data)
    model = cube.estimate_noise()
    out_cube, out_model = cube.convolve_to_beam(MID_BEAM, noise=model)
    quiet = ~line  # σ measured on the emission-free channels only
    measured = _interior_std_ratio(
        cube._data_quantity.value[quiet], out_cube._data_quantity.value[quiet], pad=40
    )
    predicted = out_model.scalar.to_value(u.K) / model.scalar.to_value(u.K)
    print(f"faint emission: measured {measured:.4f} ACF law {predicted:.4f}")
    assert predicted == pytest.approx(measured, rel=0.03)


def test_acf_law_matches_beam_law_for_beam_correlated_noise(tmp_path):
    # White noise convolved with the beam: the textbook case where the beam-area law is exact.
    # A 6"x5" beam keeps the map many correlation areas wide: the per-plane mean subtraction
    # biases ρ by ~ -(correlation area / map area), which a 12" beam on 192 px would not hide.
    small, small_out = (
        radio_beam.Beam(major=6 * u.arcsec, minor=5 * u.arcsec, pa=30 * u.deg),
        radio_beam.Beam(major=9 * u.arcsec, minor=8 * u.arcsec, pa=30 * u.deg),
    )
    data = noise_mod.synthesize_correlated_noise(
        (64, 192, 192),
        beam=small,
        pixel_scale=PIX,
        sigma=SIGMA,
        rng=np.random.default_rng(961),
    )
    header = _cube_header()
    header.update(small.to_header_keywords())
    fits.writeto(tmp_path / "beamcorr.fits", data.astype("float32"), header)
    model = Cube.from_loaded(load(tmp_path / "beamcorr.fits")).estimate_noise()
    acf_law = float(
        noise_mod.propagate(model, beam_out=small_out).scalar / model.scalar
    )
    beam_law = noise_mod._spatial_rms_factor(small, small_out)
    print(f"beam-correlated: ACF law {acf_law:.4f} beam law {beam_law:.4f}")
    assert acf_law == pytest.approx(beam_law, rel=0.03)


def test_propagate_without_spatial_acf_is_beam_law_and_approximate(noise_cube):
    model = NoiseModel.from_rms_map(noise_cube, np.full(noise_cube.shape[1:], SIGMA))
    assert model._rep.spatial_acf is None
    out = noise_mod.propagate(model, beam_out=LARGER_BEAM)
    assert out.quality is NoiseQuality.APPROXIMATE
    ratio = out.scalar.to_value(u.K) / model.scalar.to_value(u.K)
    assert ratio == pytest.approx(
        noise_mod._spatial_rms_factor(BEAM, LARGER_BEAM), rel=1e-9
    )


def test_with_spatial_acf_attaches_measured_acf(noise_cube):
    model = NoiseModel.from_rms_map(noise_cube, np.full(noise_cube.shape[1:], SIGMA))
    attached = noise_mod.with_spatial_acf(model, noise_cube)
    # same noise, flat vs measured σ_xy normalisation: the ACFs agree to the estimate scatter
    np.testing.assert_allclose(
        attached._rep.spatial_acf,
        noise_cube.estimate_noise()._rep.spatial_acf,
        atol=0.03,
    )
    assert attached.quality is model.quality
    out = noise_mod.propagate(attached, beam_out=LARGER_BEAM)
    assert out.quality is NoiseQuality.PROPAGATED


def test_spatial_acf_survives_cutout_and_unit_change(noise_cube):
    model = noise_cube.estimate_noise()
    acf = model._rep.spatial_acf
    np.testing.assert_array_equal(model.cutout(0, 8, 0, 8)._rep.spatial_acf, acf)
    np.testing.assert_array_equal(model.to(u.mK)._rep.spatial_acf, acf)


def test_spatial_acf_zarr_roundtrip(noise_cube, tmp_path):
    store = noise_cube.to_zarr(tmp_path / "z")
    model = noise_cube.estimate_noise()
    model.to_zarr_companion(store)
    back = NoiseModel.from_zarr_companion(store)
    np.testing.assert_array_equal(back._rep.spatial_acf, model._rep.spatial_acf)


def test_companion_without_spatial_acf_still_loads(noise_cube, tmp_path):
    # Existing L1 companions were written before #96 and carry no spatial ACF.
    store = noise_cube.to_zarr(tmp_path / "z")
    NoiseModel.from_rms_map(
        noise_cube, np.full(noise_cube.shape[1:], SIGMA)
    ).to_zarr_companion(store)
    back = NoiseModel.from_zarr_companion(store)
    assert back._rep.spatial_acf is None


def test_rewriting_companion_without_spatial_acf_drops_stale_acf(noise_cube, tmp_path):
    # Overwriting a companion with a model lacking the ACF must not leave the old one behind.
    from dataclasses import replace

    store = noise_cube.to_zarr(tmp_path / "z")
    model = noise_cube.estimate_noise()
    model.to_zarr_companion(store)
    NoiseModel(
        model._cube,
        replace(model._rep, spatial_acf=None),
        method=model.method,
        quality=model.quality,
    ).to_zarr_companion(store)
    assert NoiseModel.from_zarr_companion(store)._rep.spatial_acf is None
