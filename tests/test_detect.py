"""Tests for the masked, ACF+beam-corrected detection statistic with SoFiA-style reliability.

Written first (red/green TDD). These tests *are* the correctness obligation for
:mod:`astrolyze.core.detect`. The contracts pinned here:

1. **Band-width invariance (the headline property).** The masked integrated SNR of a fixed line
   is set by the *line* it detects (its ``L`` channels), NOT by the cube's channel count ``M``.
   Padding a cube with empty channels (M: 60 -> 300) leaves the statistic unchanged — whereas a
   full-band integral dilutes by ``√(M/L)``. This is *why* the module exists (issue #22).

2. **The two noise corrections are exact reductions, not fudge factors.**
   - spectral ``κ_spec = M_det / M_eff(ρ, M_det)`` (≥1; =1 white) — reuses the noise module's
     ``_m_eff`` over the *detected* channels;
   - spatial ``κ_spat = pixels-per-beam = Ω_beam / Ω_pix`` (≥1) — a sum over a beam-correlated
     footprint does NOT average down as ``√N_pix``; it is the integration-side dual of the noise
     module's ``_spatial_rms_factor``.

3. **Smooth-and-clip in PHYSICAL units.** Kernels are km/s (spectral) and beams (spatial),
   converted per-cube — so the same physical line at 2 vs 5 km/s channels detects consistently
   (a channel-count kernel would not). Single-pixel / single-channel spikes are rejected by the
   connected-area floor + the consecutive-channel rule.

4. **SoFiA reliability from the negative field.** A real line lands where the sign-flipped noise
   null cannot reach -> reliability ≈ 1; pure noise -> low reliability / no detection.

5. **No silent physics (ADR-0003).** An ``UNRELIABLE`` noise model or an all-NaN tile yields
   ``UNRELIABLE`` + NaN, never a fabricated 0; a line-free but *reliable* tile is the common,
   non-pathological ``NO_DETECTION`` case and returns finite zeros.

The reference dataset mirrors the noise tests: NGC 0628, CO(2-1) at 230.538 GHz, 12"x10" (PA 30)
beam, radio convention. ``synthesize_correlated_noise`` is the injection harness (beam-correlated,
ACF-shaped noise at a stated σ).
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, NoiseModel, NoiseQuality
from astrolyze.core import noise as noise_mod
from astrolyze.core.detect import (
    DEFAULT_LADDER,
    ComponentDetection,
    DetectionQuality,
    DetectionResult,
    ScaleProvenance,
    ScaleStep,
    detect,
    detect_components,
    detection_params,
    kappa_spatial,
    kappa_spectral,
    pixels_per_beam,
    reliability,
)
from astrolyze.io import load

warnings.filterwarnings("ignore", module="spectral_cube")
warnings.filterwarnings("ignore", module="astropy")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
SIGMA = 0.2  # the true per-voxel noise level (K)
PIXEL = 2e-4 * u.deg  # |CDELT|; 0.72 arcsec/pixel
CHANNEL_WIDTH = 2000.0 * u.m / u.s  # 2 km/s


def _cube_header(channel_width=CHANNEL_WIDTH, bunit="K"):
    """A 3D cube header on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = (
        channel_width.to_value(u.m / u.s),
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


def _write_cube(path, data, channel_width=CHANNEL_WIDTH, bunit="K"):
    fits.writeto(
        path, np.asarray(data, dtype="float32"), _cube_header(channel_width, bunit)
    )
    return Cube.from_loaded(load(path))


def _gaussian_line(shape, *, c0, fwhm_chan, y, x, fwhm_pix, peak, beam, pixel_scale):
    """A beam-shaped, Gaussian-in-velocity emission line array (no noise), peak amplitude *peak*.

    Spatial profile = a 2-D Gaussian of ``fwhm_pix`` pixels at ``(y,x)``; spectral profile = a
    1-D Gaussian of ``fwhm_chan`` channels centred on ``c0``. The outer product, scaled to *peak*.
    """
    nz, ny, nx = shape
    to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    zz = np.arange(nz)
    spec = np.exp(-0.5 * ((zz - c0) / (fwhm_chan * to_sigma)) ** 2)
    yy, xx = np.mgrid[0:ny, 0:nx]
    sky = np.exp(-0.5 * (((yy - y) ** 2 + (xx - x) ** 2) / (fwhm_pix * to_sigma) ** 2))
    return peak * spec[:, None, None] * sky[None, :, :]


def _noise(shape, *, rng, sigma=SIGMA, spectral_acf=None):
    return noise_mod.synthesize_correlated_noise(
        shape,
        beam=BEAM,
        pixel_scale=PIXEL,
        sigma=sigma,
        spectral_acf=spectral_acf,
        rng=rng,
    )


# ======================================================================================
# AC 1: the two noise corrections are exact reductions (unit laws, no cube needed)
# ======================================================================================
def test_pixels_per_beam_is_beam_over_pixel_area():
    ppb = pixels_per_beam(BEAM, PIXEL)
    omega_pix = (PIXEL.to_value(u.deg) ** 2) * (u.deg**2).to(u.sr)
    expected = float(BEAM.sr.to_value(u.sr) / omega_pix)
    assert ppb == pytest.approx(expected, rel=1e-9)
    assert ppb > 1.0  # a 12" beam is many 0.72" pixels


def test_kappa_spatial_equals_pixels_per_beam_and_is_geometry_only():
    # κ_spat is a property of the BEAM and pixel grid, not of how big the mask footprint is.
    assert kappa_spatial(BEAM, PIXEL) == pytest.approx(
        pixels_per_beam(BEAM, PIXEL), rel=1e-12
    )
    assert kappa_spatial(BEAM, PIXEL) >= 1.0


def test_kappa_spectral_white_noise_is_one():
    acf = np.zeros(64)
    acf[0] = 1.0  # white: a spike at lag 0
    assert kappa_spectral(acf, 12) == pytest.approx(1.0, rel=1e-9)


def test_kappa_spectral_correlated_is_greater_than_one():
    # A broad ACF means fewer independent channels (M_eff < M) -> the integrated noise is LARGER
    # than the naive √N -> κ_spec > 1 = n / _m_eff(acf, n).
    rho = np.exp(-np.arange(64) / 3.0)
    n = 10
    assert kappa_spectral(rho, n) == pytest.approx(
        n / noise_mod._m_eff(rho, n), rel=1e-9
    )
    assert kappa_spectral(rho, n) > 1.0


# ======================================================================================
# AC 2: band-width invariance (the headline) + the full-band contrast
# ======================================================================================
def _line_cube(path, n_total, *, c0=30, L=5, peak=8 * SIGMA, seed=1):
    """A cube with a fixed beam-shaped line in a 60-channel region, padded with noise to n_total.

    Returns ``(cube, noise_model)`` where the model is estimated from the line-FREE field: the
    band-width-invariance contract is a property of the DETECTOR given a noise model, so the test
    must not be confounded by mad_std's mild line bias (worse when the line is a larger fraction of
    a short band). The line region (first 60 channels) is byte-identical across n_total."""
    ny = nx = 48
    base_m = 60
    assert n_total >= base_m
    rng = np.random.default_rng(seed)
    field = _noise((base_m, ny, nx), rng=rng)
    if n_total > base_m:
        extra = _noise(
            (n_total - base_m, ny, nx), rng=np.random.default_rng(seed + 100)
        )
        field = np.concatenate([field, extra], axis=0)
    line = _gaussian_line(
        field.shape,
        c0=c0,
        fwhm_chan=L,
        y=ny // 2,
        x=nx // 2,
        fwhm_pix=18.0,
        peak=peak,
        beam=BEAM,
        pixel_scale=PIXEL,
    )
    cube = _write_cube(path, field + line)
    # The cubes are spectrally white (no spectral_acf in synthesis), so the exact noise model is a
    # constant σ map. Using it (vs estimate_noise) removes the σ-estimator's per-pixel scatter so
    # the band-width-invariance contract is tested on the detector alone.
    model = NoiseModel.from_rms_map(cube, np.full((ny, nx), SIGMA))
    return cube, model


def test_band_width_invariance(tmp_path):
    """THE acceptance gate: same line, M=60 vs M=300 -> the same masked integrated SNR."""
    narrow, m_narrow = _line_cube(tmp_path / "narrow.fits", 60)
    wide, m_wide = _line_cube(tmp_path / "wide.fits", 300)
    r_narrow = detect(narrow, m_narrow, rng=np.random.default_rng(0))
    r_wide = detect(wide, m_wide, rng=np.random.default_rng(0))
    assert r_narrow.quality is DetectionQuality.MEASURED
    assert r_wide.quality is DetectionQuality.MEASURED
    assert r_wide.integrated_snr == pytest.approx(r_narrow.integrated_snr, rel=0.2)


def test_full_band_integral_is_diluted_but_masked_is_not(tmp_path):
    """The contrast that justifies the module: the naive full-band SNR collapses with M, the
    masked statistic does not."""
    narrow, m_narrow = _line_cube(tmp_path / "n.fits", 60)
    wide, m_wide = _line_cube(tmp_path / "w.fits", 300)

    def full_band_snr(cube, model):
        data = np.asarray(cube.validity.data.value)
        sigma = np.asarray(model.sigma_cube.validity.data.value)
        integ = np.nansum(data, axis=0)
        var = np.nansum(sigma**2, axis=0)
        return float(np.nanmax(integ / np.sqrt(var)))

    fb_ratio = full_band_snr(wide, m_wide) / full_band_snr(narrow, m_narrow)
    masked_ratio = (
        detect(wide, m_wide, rng=np.random.default_rng(0)).integrated_snr
        / detect(narrow, m_narrow, rng=np.random.default_rng(0)).integrated_snr
    )
    assert fb_ratio < 0.7  # full band: diluted ~√(60/300) ≈ 0.45
    assert masked_ratio == pytest.approx(1.0, abs=0.2)  # masked: invariant


# ======================================================================================
# AC 3: line recovery + reliability
# ======================================================================================
def test_strong_line_is_recovered_with_high_reliability(tmp_path):
    cube, model = _line_cube(tmp_path / "line.fits", 60, peak=15 * SIGMA)
    result = detect(cube, model, rng=np.random.default_rng(0))
    assert result.quality is DetectionQuality.MEASURED
    assert result.integrated_snr > 5.0
    assert result.n_signal_voxels > 0
    assert result.area_beams >= 1.0
    # a strong, clean line sits well above the sign-flipped noise null -> high reliability.
    assert result.reliability > 0.8


def test_pure_noise_is_not_a_confident_detection(tmp_path):
    rng = np.random.default_rng(42)
    field = _noise((80, 48, 48), rng=rng)
    cube = _write_cube(tmp_path / "noise.fits", field)
    result = detect(cube, cube.estimate_noise(), rng=np.random.default_rng(0))
    # either nothing clears the mask, or what does has low reliability — never a confident line.
    if result.quality is DetectionQuality.NO_DETECTION:
        assert result.integrated_snr == 0.0
    else:
        assert not (result.reliability > 0.5)


# ======================================================================================
# AC 4: physical-unit kernels + morphological rejection
# ======================================================================================
def test_single_pixel_single_channel_spike_is_not_a_strong_detection(tmp_path):
    rng = np.random.default_rng(7)
    field = _noise((80, 48, 48), rng=rng)
    field[40, 24, 24] += 20 * SIGMA  # one voxel, no spatial/spectral neighbours
    cube = _write_cube(tmp_path / "spike.fits", field)
    result = detect(cube, cube.estimate_noise(), rng=np.random.default_rng(0))
    # A point-like spike has its flux concentrated: the κ-corrected masked integral over any
    # blob the smoothing kernels spread it into is WEAK (matched-filter mismatch). It never
    # becomes a strong, keep-worthy detection. (A sign-flip null cannot reject a positive-only
    # *artifact* — that is an upstream deglitch concern, not the reliability statistic's job.)
    if result.quality is DetectionQuality.MEASURED:
        assert result.integrated_snr < 5.0


def test_km_s_kernel_consistent_across_channel_widths(tmp_path):
    """The same PHYSICAL line written at 2 km/s and 5 km/s channels detects consistently — the
    kernels are km/s, not channels."""
    ny = nx = 48
    # 2 km/s cube: line spans ~5 channels (10 km/s); 5 km/s cube: same 10 km/s -> ~2 channels.
    rng = np.random.default_rng(3)
    f2 = _noise((60, ny, nx), rng=rng) + _gaussian_line(
        (60, ny, nx),
        c0=30,
        fwhm_chan=5,
        y=ny // 2,
        x=nx // 2,
        fwhm_pix=14.0,
        peak=8 * SIGMA,
        beam=BEAM,
        pixel_scale=PIXEL,
    )
    f5 = _noise((30, ny, nx), rng=np.random.default_rng(4)) + _gaussian_line(
        (30, ny, nx),
        c0=15,
        fwhm_chan=2,
        y=ny // 2,
        x=nx // 2,
        fwhm_pix=14.0,
        peak=8 * SIGMA,
        beam=BEAM,
        pixel_scale=PIXEL,
    )
    c2 = _write_cube(tmp_path / "w2.fits", f2, channel_width=2000.0 * u.m / u.s)
    c5 = _write_cube(tmp_path / "w5.fits", f5, channel_width=5000.0 * u.m / u.s)
    r2 = detect(c2, c2.estimate_noise(), rng=np.random.default_rng(0))
    r5 = detect(c5, c5.estimate_noise(), rng=np.random.default_rng(0))
    assert r2.quality is DetectionQuality.MEASURED
    assert r5.quality is DetectionQuality.MEASURED
    # both detect the same physical line at comparable SNR (generous tolerance: different gridding)
    assert r5.integrated_snr == pytest.approx(r2.integrated_snr, rel=0.4)


def test_continuum_m1_has_no_spectral_penalty(tmp_path):
    rng = np.random.default_rng(11)
    field = _noise((1, 48, 48), rng=rng)
    field[0] += _gaussian_line(
        (1, 48, 48),
        c0=0,
        fwhm_chan=1,
        y=24,
        x=24,
        fwhm_pix=14.0,
        peak=10 * SIGMA,
        beam=BEAM,
        pixel_scale=PIXEL,
    )[0]
    cube = _write_cube(tmp_path / "cont.fits", field)
    # continuum noise realistically comes from a published RMS map (a 2-D σ field), the intended
    # path for a single-channel map — estimate_noise needs ≥2 channels for its separability test.
    model = NoiseModel.from_rms_map(cube, np.full((48, 48), SIGMA))
    result = detect(cube, model, rng=np.random.default_rng(0))
    assert result.kappa_spec == pytest.approx(
        1.0, rel=1e-9
    )  # M_det=1 -> no spectral penalty
    assert np.isfinite(result.integrated_snr)


# ======================================================================================
# AC 5: no silent physics (ADR-0003)
# ======================================================================================
def test_unreliable_noise_model_refuses_to_fabricate(tmp_path):
    rng = np.random.default_rng(5)
    cube = _write_cube(tmp_path / "n.fits", _noise((40, 32, 32), rng=rng))
    model = cube.estimate_noise()
    # force an UNRELIABLE model (the ADR-0003 no-usable-data case)
    bad = NoiseModel(
        model._cube,
        noise_mod._unreliable_representation(cube.shape),
        method=model.method,
        quality=NoiseQuality.UNRELIABLE,
    )
    result = detect(cube, bad, rng=np.random.default_rng(0))
    assert result.quality is DetectionQuality.UNRELIABLE
    assert np.isnan(result.integrated_snr)


def test_all_nan_tile_is_unreliable_not_zero(tmp_path):
    field = np.full((40, 32, 32), np.nan)
    cube = _write_cube(tmp_path / "nan.fits", field)
    # an all-NaN cube has no usable data; estimate_noise flags UNRELIABLE, and so must detect().
    result = detect(cube, cube.estimate_noise(), rng=np.random.default_rng(0))
    assert result.quality is DetectionQuality.UNRELIABLE
    assert np.isnan(result.integrated_snr)


def test_no_detection_returns_finite_zero_not_error(tmp_path):
    rng = np.random.default_rng(99)
    cube = _write_cube(tmp_path / "empty.fits", _noise((80, 48, 48), rng=rng))
    result = detect(cube, cube.estimate_noise(), rng=np.random.default_rng(0))
    # line-free but reliable: the common SSL case. Finite zeros, NaN reliability, never raises.
    if result.quality is DetectionQuality.NO_DETECTION:
        assert result.integrated_snr == 0.0
        assert result.n_signal_voxels == 0


def test_result_carries_full_provenance(tmp_path):
    cube, model = _line_cube(tmp_path / "line.fits", 60, peak=8 * SIGMA)
    result = detect(cube, model, rng=np.random.default_rng(0))
    assert isinstance(result, DetectionResult)
    for field in (
        "integrated_snr",
        "reliability",
        "n_signal_voxels",
        "area_beams",
        "n_channels_detected",
        "kappa_spec",
        "kappa_spat",
        "pixels_per_beam",
        "bowl_fraction",
        "min_neg_sigma",
        "null_kind",
        "scales",
    ):
        assert hasattr(result, field)
    assert len(result.scales) >= 1
    assert isinstance(result.scales[0], ScaleProvenance)


def test_default_ladder_is_physical():
    # the ladder steps are physical (beams + km/s), the band-width-invariance choice.
    assert all(isinstance(s, ScaleStep) for s in DEFAULT_LADDER)
    assert any(s.spectral_fwhm_kms > 0 for s in DEFAULT_LADDER)
    assert any(s.spatial_fwhm_beams > 0 for s in DEFAULT_LADDER)


# ======================================================================================
# AC 6: the array-level surface a tiling caller (the SSL sharder) composes
# ======================================================================================
def test_detect_components_returns_footprint_bbox_around_the_line(tmp_path):
    """detect_components yields ComponentDetections whose projected bbox brackets the source — so a
    tiler can map a component to the tiles it overlaps (the cube-wide-detect-then-assign pattern)."""
    cube, model = _line_cube(tmp_path / "line.fits", 60, peak=12 * SIGMA)
    data = np.asarray(cube.validity.data.value)
    sigma = np.asarray(model.sigma_cube.validity.data.value)
    acf = np.asarray(model.spectral_acf.flux.value)
    comps, scales = detect_components(
        data,
        sigma,
        acf=acf,
        beam=BEAM,
        pixel_scale=cube.coordinates.pixel_scale,
        channel_width_kms=2.0,
    )
    assert comps and all(isinstance(c, ComponentDetection) for c in comps)
    best = max(comps, key=lambda c: c.integrated_snr)
    assert (
        best.y0 <= 24 < best.y1 and best.x0 <= 24 < best.x1
    )  # brackets the line at (24,24)
    assert best.integrated_snr > 5.0


def test_cube_wide_null_assigns_reliability_to_a_component(tmp_path):
    """The sharder pattern: build a cube-wide null once, score a positive component against it."""
    cube, model = _line_cube(tmp_path / "line.fits", 60, peak=15 * SIGMA)
    data = np.asarray(cube.validity.data.value)
    sigma = np.asarray(model.sigma_cube.validity.data.value)
    acf = np.asarray(model.spectral_acf.flux.value)
    kw = dict(
        acf=acf,
        beam=BEAM,
        pixel_scale=cube.coordinates.pixel_scale,
        channel_width_kms=2.0,
    )
    pos, _ = detect_components(data, sigma, **kw)
    neg, _ = detect_components(-data, sigma, **kw)
    pos = [c for c in pos if c.integrated_snr > 0]
    best_idx = int(np.argmax([c.integrated_snr for c in pos]))
    r = reliability(detection_params(pos), detection_params(neg))[best_idx]
    assert (
        0.0 <= r <= 1.0 and r > 0.8
    )  # a strong line is reliable against the negative null


def test_detection_params_keeps_only_positive_in_log_space():
    comps = [
        ComponentDetection(10.0, 4.0, 100, 5, 1.0, 0, 5, 0, 5),
        ComponentDetection(
            -3.0, 2.0, 20, 2, 1.0, 0, 3, 0, 3
        ),  # net-negative -> dropped
    ]
    params = detection_params(comps)
    assert params.shape == (1, 2)
    assert params[0, 0] == pytest.approx(np.log10(10.0))
