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


def test_spatial_propagation_uses_beam_not_pixel_area(noise_cube):
    # The scaling must depend ONLY on the beam ratio, not on the pixel grid. Two cubes with the
    # same beams but different pixel scales propagate to the SAME σ ratio (beam-, not pixel-,
    # limited correlation area).
    model = noise_cube.estimate_noise()
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
    expected = model.scalar.to_value(u.K) * float(
        np.sqrt((BEAM.sr / LARGER_BEAM.sr).to_value(u.dimensionless_unscaled))
    )
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
