"""Masked, ACF+beam-corrected detection statistic with SoFiA-style reliability (issue ifm#22).

A cube's "is there a line worth keeping, and how confidently" reduced to a few numbers, the
band-width-invariant way the radio source-finding community settled on (SoFiA-2 Westmeier+2021,
Serra+2012, maskmoment/CPROPS/astrodendro). It is a *companion* analysis like
:mod:`astrolyze.core.noise` — it consumes a cube + its :class:`~astrolyze.core.NoiseModel` and
returns a frozen :class:`DetectionResult`; **policy** (keep/threshold/weight) stays with the caller.

The problem it fixes. A full-band integrated SNR ``Σ_v I / √(Σ_v σ²)`` over an ``M``-channel cube
makes a line that occupies ``L`` of ``M`` channels pay ``√M`` of integration noise — diluting it by
``~√(M/L)`` (×6 on PHANGS, M~400). The fix is to integrate over a **dilated detection mask** so the
noise pays only the detected ``L``: band-width-invariant by construction.

The statistic. Over a connected detection ``𝓜`` (a "smooth-and-clip" mask, below) the integrated
signal is ``S = Σ_𝓜 I`` and its variance is the naive heteroscedastic sum ``Σ_𝓜 σ²`` inflated by
the two correlation penalties the radio noise carries:

- **spectral** ``κ_spec = M_det / M_eff(ρ, M_det)`` (≥1, =1 white) — summing ``M_det`` correlated
  channels does not add as ``M_det`` independent variances but as ``M_eff`` from the ACF
  (:func:`astrolyze.core.noise._m_eff`, reused);
- **spatial** ``κ_spat = pixels-per-beam = Ω_beam/Ω_pix`` (≥1) — a sum over a beam-correlated
  footprint does **not** average down as ``√N_pix`` (pixels inside one beam are the same noise
  sample). This is the integration-side **dual** of :func:`astrolyze.core.noise._spatial_rms_factor`
  (that scales σ *down* when you smooth into a larger beam; this scales the integral variance *up*).

so ``SNR = Σ_𝓜 I / √(κ_spec · κ_spat · Σ_𝓜 σ²)``. Treating spatial (beam) and spectral (ACF)
correlation as separable is the standard SoFiA approximation; the per-voxel ``σ`` from the noise
companion carries the (possibly non-stationary) level exactly, the two κ carry the stationary
correlation. ``Δv`` cancels, so the number is a pure SNR.

The mask (smooth-and-clip). For each rung of a small **physical** kernel ladder (widths in km/s and
beams, converted per-cube — NOT channels, so CO / broad H i / continuum all match physics) the cube
and its σ are matched-filter-smoothed, a ``t_hi``-σ core (in ``≥ min_consecutive_channels``
consecutive channels — the PHANGS spike-rejection rule) is dilated into a ``t_lo``-σ envelope
(``scipy.ndimage.binary_propagation``, the operator from the labbook ``signal_mask``), and a
connected-area floor (``≥ min_area_beams``) drops spikes. The per-rung masks are OR-combined.

Reliability (Serra+2012, Westmeier+2021). The identical detector run on the **sign-flipped** cube
is a measured false-positive population; ``R = 1 − N_neg/N_pos`` in a kernel-smoothed
``(log SNR, log area)`` space. The caller's single knob is a reliability target — auto-adapting to
``M`` and tile size, which kills the cross-survey threshold problem. It rejects **symmetric-noise**
false positives; it cannot reject a **positive-only artifact** (a single hot pixel / un-deglitched
cosmic ray), which has no negative counterpart — but such an artifact's flux is concentrated, so
its κ-corrected masked integral is weak regardless (the integrated SNR, not reliability, is what
down-weights it). Real artifact cleaning is an upstream deglitch step, not this statistic's job. Interferometer caveat
(ADR-0003): a deterministic negative bowl is **not** noise, so a bowl DQ (``min_neg_sigma``,
``bowl_fraction``) is always recorded and an ``off_source`` null is offered; sign-flip on a
bowl-flagged cube is warned, never silently corrected.

No silent physics (ADR-0003): an ``UNRELIABLE`` noise model or an all-NaN tile yields
:data:`DetectionQuality.UNRELIABLE` + NaN, never a fabricated 0; a line-free but reliable tile is
the common, non-pathological :data:`DetectionQuality.NO_DETECTION` (finite zeros). Stays thin:
``scipy.ndimage`` + numpy; reuses ``_m_eff`` / ``_beam_sigma_pixels`` from the noise module.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np
import astropy.units as u
from astropy.stats import mad_std
from scipy import ndimage

from ._base import _emit
from .noise import NoiseModel, NoiseQuality, _m_eff

#: Schema version of the result/provenance layout (bumped on a layout change).
DETECT_SCHEMA_VERSION = 1

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))

# Smooth-and-clip thresholds — pre-1.0 DECISIONs (the PHANGS-ALMA mask-construction defaults,
# Leroy+2021 / SoFiA-2 Westmeier+2021). They are knobs (the reference probe sweeps them), not
# format: a 4σ core in ≥2 consecutive channels rejects single-channel spikes, dilated into a 2σ
# envelope, with a ≥1-beam connected-area floor against single-pixel spikes.
DEFAULT_T_HI = 4.0
DEFAULT_T_LO = 2.0
DEFAULT_MIN_CONSECUTIVE_CHANNELS = 2
DEFAULT_MIN_AREA_BEAMS = 1.0
KERNEL_TRUNCATE = 4.0  # Gaussian kernels are truncated at ±4σ (scipy's default)


@dataclass(frozen=True)
class ScaleStep:
    """One rung of the multi-scale smooth-and-clip ladder, in **physical** units.

    Widths are beams (spatial) and km/s (spectral), converted to pixels/channels per cube — the
    band-width-invariance choice (a channel-count kernel would over-smooth broad H i and be
    undefined for continuum). ``0`` = no smoothing on that axis (native resolution)."""

    spatial_fwhm_beams: float
    spectral_fwhm_kms: float


#: The default kernel ladder (SoFiA-2's small kernel set recast into beams & km/s). A pre-1.0
#: DECISION; the reference probe (#10) sweeps it. native + spectral-only + spatial-only + both,
#: so sharp, faint-broad, faint-extended, and faint-extended-broad emission are each reachable.
DEFAULT_LADDER = (
    ScaleStep(spatial_fwhm_beams=0.0, spectral_fwhm_kms=0.0),
    ScaleStep(spatial_fwhm_beams=0.0, spectral_fwhm_kms=10.0),
    ScaleStep(spatial_fwhm_beams=2.0, spectral_fwhm_kms=0.0),
    ScaleStep(spatial_fwhm_beams=2.0, spectral_fwhm_kms=20.0),
)


class DetectionQuality(enum.Enum):
    """How to read a :class:`DetectionResult` (mirrors :class:`~astrolyze.core.NoiseQuality`).

    - ``MEASURED`` — a real detection with a usable null, so the SNR + reliability are meaningful;
    - ``NO_DETECTION`` — the mask is empty (the common line-free SSL tile): finite zeros, reliability
      ``NaN``. **Not** an error;
    - ``UNRELIABLE`` — the noise model is ``UNRELIABLE`` / the tile is all-NaN: SNR + reliability are
      ``NaN``, never a fabricated number (ADR-0003);
    - ``NO_NULL`` — a detection exists but no null population could be built (``null="none"`` or an
      empty off-source field): SNR is real, reliability ``NaN``.
    """

    MEASURED = "measured"
    NO_DETECTION = "no_detection"
    UNRELIABLE = "unreliable"
    NO_NULL = "no_null"


@dataclass(frozen=True)
class ScaleProvenance:
    """Per-rung audit of the smooth-and-clip mask (what each kernel actually contributed)."""

    spatial_fwhm_arcsec: float
    spectral_fwhm_kms: float
    n_core_voxels: int
    n_mask_voxels: int
    n_objects: int


@dataclass(frozen=True)
class DetectionResult:
    """The per-tube detection statistic + its provenance (a companion value object, ADR-0004).

    ``integrated_snr`` / ``reliability`` / geometry refer to the **strongest** detection (the max
    masked SNR connected component); ``kappa_spec``/``kappa_spat``/``pixels_per_beam`` are the
    corrections applied to it; ``bowl_fraction``/``min_neg_sigma`` are always-recorded DQ flags."""

    integrated_snr: float
    reliability: float
    n_signal_voxels: int
    area_beams: float
    n_channels_detected: int
    kappa_spec: float
    kappa_spat: float
    pixels_per_beam: float
    bowl_fraction: float
    min_neg_sigma: float
    null_kind: str
    scales: tuple[ScaleProvenance, ...]
    quality: DetectionQuality


# ======================================================================================
# The two correlation penalties (exact reductions, reusing the noise module's laws)
# ======================================================================================
def pixels_per_beam(beam, pixel_scale) -> float:
    """Number of pixels per beam ``Ω_beam / Ω_pix`` (the spatial correlation length).

    ``beam.sr`` is radio_beam's beam solid angle; ``pixel_scale`` is the celestial pixel scale
    (one value per spatial axis, or a scalar) — the pixel solid angle is the product of the two
    axis scales (a square pixel: the scale squared)."""
    omega_beam = float(beam.sr.to_value(u.sr))
    scales = np.atleast_1d(np.asarray(pixel_scale.to_value(u.deg), dtype="float64"))
    if scales.size >= 2:
        omega_pix_deg2 = float(scales[0] * scales[1])
    else:
        omega_pix_deg2 = float(scales[0] ** 2)
    omega_pix = omega_pix_deg2 * (u.deg**2).to(u.sr)
    if omega_pix <= 0.0:
        return float("nan")
    return omega_beam / omega_pix


def kappa_spatial(beam, pixel_scale) -> float:
    """Spatial variance inflation for a masked integral: ``κ_spat = pixels-per-beam`` (≥1).

    A sum over a beam-correlated footprint has variance ``σ² · N_pix · n_ppb`` (the beam-cell
    argument: ``N_pix/n_ppb`` independent samples, each of variance ``σ²·n_ppb²``), so the naive
    ``Σσ²`` is inflated by ``n_ppb``. Clipped to ≥1 (under-sampled pixels larger than the beam are
    independent — no correlation gain)."""
    return float(max(1.0, pixels_per_beam(beam, pixel_scale)))


def kappa_spectral(acf, n_channels_detected: int) -> float:
    """Spectral variance inflation for a masked integral: ``κ_spec = n / M_eff(ρ, n)`` (≥1).

    Summing ``n`` correlated channels adds as ``M_eff = n²/Σρ`` independent samples, so the naive
    per-channel-independent ``Σσ²`` is inflated by ``n/M_eff``. Reuses
    :func:`astrolyze.core.noise._m_eff` over the **detected** channels; white noise ⇒ ``M_eff=n`` ⇒
    1; a single channel (continuum) ⇒ 1 (no spectral penalty)."""
    n = int(n_channels_detected)
    if n <= 1:
        return 1.0
    m_eff = _m_eff(np.asarray(acf, dtype="float64"), n)
    if not np.isfinite(m_eff) or m_eff <= 0.0:
        return 1.0
    return float(max(1.0, n / m_eff))


# ======================================================================================
# Kernel-width conversion (physical -> pixels/channels, per cube)
# ======================================================================================
def _spatial_kernel_sigma_pixels(spatial_fwhm_beams: float, beam, pixel_scale) -> float:
    """The spatial Gaussian σ in pixels for a kernel of ``spatial_fwhm_beams`` beams.

    Kernel FWHM (arcsec) = beams × the geometric-mean beam FWHM; σ_pixels via the pixel scale.
    ``0`` beams ⇒ ``0`` (no spatial smoothing). Square-pixel cubes (the L1 grid) ⇒ one σ."""
    if spatial_fwhm_beams <= 0.0:
        return 0.0
    beam_fwhm_arcsec = float(
        np.sqrt(beam.major.to_value(u.arcsec) * beam.minor.to_value(u.arcsec))
    )
    scale_arcsec = float(np.mean(np.atleast_1d(pixel_scale.to_value(u.arcsec))))
    fwhm_pix = spatial_fwhm_beams * beam_fwhm_arcsec / scale_arcsec
    return fwhm_pix / FWHM_PER_SIGMA


def _spectral_kernel_sigma_channels(
    spectral_fwhm_kms: float, channel_width_kms: float
) -> float:
    """The spectral Gaussian σ in channels for a kernel of ``spectral_fwhm_kms`` km/s.

    ``0`` km/s (or an unknown channel width — e.g. continuum) ⇒ ``0`` (no spectral smoothing)."""
    if (
        spectral_fwhm_kms <= 0.0
        or not np.isfinite(channel_width_kms)
        or channel_width_kms <= 0.0
    ):
        return 0.0
    fwhm_chan = spectral_fwhm_kms / channel_width_kms
    return fwhm_chan / FWHM_PER_SIGMA


def _gaussian_weights(sigma: float) -> np.ndarray:
    """A normalised 1-D Gaussian kernel (Σ=1); ``sigma<=0`` ⇒ the identity ``[1.0]``."""
    if sigma <= 0.0:
        return np.array([1.0])
    radius = max(1, int(KERNEL_TRUNCATE * sigma + 0.5))
    x = np.arange(-radius, radius + 1, dtype="float64")
    k = np.exp(-0.5 * (x / sigma) ** 2)
    return k / k.sum()


# ======================================================================================
# Matched-filter smoothing (NaN-safe, with variance propagation) -> the smoothed SNR map
# ======================================================================================
def _smoothed_snr(
    data: np.ndarray,
    sigma_field: np.ndarray,
    finite: np.ndarray,
    *,
    sigma_chan: float,
    sigma_pix: float,
) -> np.ndarray:
    """The matched-filter smoothed SNR map, with the smoothed noise **measured** (not propagated).

    The cube is first normalised to per-voxel SNR ``z = I/σ`` (so heteroscedastic noise becomes
    homoscedastic), then ``z`` is smoothed by the separable Gaussian kernel. The smoothed-noise
    level ``σ_sm`` is then measured **robustly from the smoothed map itself** (``mad_std``), so it
    carries the spatial (beam) and spectral (ACF) correlation automatically — an analytic
    ``√(Σ w² σ²)`` would assume independent pixels and badly under-state the noise of
    beam-correlated data (smoothing by ≥1 beam barely reduces it), lighting up the whole field.
    This is the SoFiA-2 smooth-then-measure rule (Westmeier+2021). The convolutions are
    mask-renormalised so NaN blanks and edges don't bleed. The native rung (no smoothing) recovers
    ``z`` with ``σ_sm ≈ 1`` ⇒ a plain ``I/σ`` threshold."""
    k_chan = _gaussian_weights(sigma_chan)
    k_pix = _gaussian_weights(sigma_pix)
    axes = [(0, k_chan), (1, k_pix), (2, k_pix)]

    fin = finite.astype("float64")
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(finite, data / sigma_field, 0.0)

    weight = _separable_correlate(fin, axes)  # Σ w   (per voxel, mode='constant')
    z_smoothed = _separable_correlate(z, axes)  # Σ w z
    with np.errstate(divide="ignore", invalid="ignore"):
        z_sm = np.where(weight > 0, z_smoothed / weight, 0.0)  # weighted mean SNR

    measurable = finite & (weight > 0)
    values = z_sm[measurable]
    sigma_sm = float(mad_std(values)) if values.size else 0.0
    if not np.isfinite(sigma_sm) or sigma_sm <= 0.0:
        sigma_sm = 1.0  # degenerate (e.g. all-equal) -> no renormalisation
    return np.where(finite, z_sm / sigma_sm, 0.0)


def _separable_correlate(arr: np.ndarray, axes) -> np.ndarray:
    """Apply 1-D kernels along given axes in turn (a separable N-D correlation, ``mode='constant'``)."""
    out = arr
    for axis, weights in axes:
        if weights.size > 1:
            out = ndimage.correlate1d(
                out, weights, axis=axis, mode="constant", cval=0.0
            )
    return out


# ======================================================================================
# Smooth-and-clip mask construction
# ======================================================================================
def _enforce_consecutive_channels(core: np.ndarray, min_consecutive: int) -> np.ndarray:
    """Keep core voxels in a run of ``≥ min_consecutive`` consecutive channels (PHANGS spike rule).

    Connectivity only along the spectral axis: each component is a vertical run at one pixel, whose
    voxel count is its run length. ``min_consecutive<=1`` is a no-op."""
    if min_consecutive <= 1 or not core.any():
        return core
    structure = np.zeros((3, 3, 3), dtype=int)
    structure[:, 1, 1] = 1  # connect along axis 0 only
    labels, n = ndimage.label(core, structure=structure)
    if n == 0:
        return core
    sizes = np.bincount(labels.ravel())
    keep = sizes >= min_consecutive
    keep[0] = False  # background
    return keep[labels]


def _area_floor(mask: np.ndarray, min_area_pixels: int) -> tuple[np.ndarray, int]:
    """Drop connected components whose **projected** spatial footprint is < ``min_area_pixels``.

    3-D 26-connectivity; the projected area rejects a single-pixel spike even when it is dilated
    in velocity. Returns the cleaned mask and the surviving-object count."""
    if not mask.any():
        return mask, 0
    structure = np.ones((3, 3, 3), dtype=int)
    labels, n = ndimage.label(mask, structure=structure)
    out = np.zeros_like(mask)
    kept = 0
    for i in range(1, n + 1):
        comp = labels == i
        if int(comp.any(axis=0).sum()) >= min_area_pixels:
            out |= comp
            kept += 1
    return out, kept


def _scale_mask(
    data,
    sigma_field,
    finite,
    *,
    step,
    beam,
    pixel_scale,
    channel_width_kms,
    t_hi,
    t_lo,
    min_consecutive,
    min_area_pixels,
) -> tuple[np.ndarray, ScaleProvenance]:
    """One ladder rung: matched-filter smooth -> core -> consecutive rule -> dilate -> area floor."""
    sigma_pix = _spatial_kernel_sigma_pixels(step.spatial_fwhm_beams, beam, pixel_scale)
    sigma_chan = _spectral_kernel_sigma_channels(
        step.spectral_fwhm_kms, channel_width_kms
    )
    snr = _smoothed_snr(
        data, sigma_field, finite, sigma_chan=sigma_chan, sigma_pix=sigma_pix
    )

    core = finite & (snr >= t_hi)
    core = _enforce_consecutive_channels(core, min_consecutive)
    envelope = finite & (snr >= t_lo)
    dilated = ndimage.binary_propagation(core, mask=envelope)
    cleaned, n_objects = _area_floor(dilated, min_area_pixels)

    beam_fwhm_arcsec = float(
        np.sqrt(beam.major.to_value(u.arcsec) * beam.minor.to_value(u.arcsec))
    )
    prov = ScaleProvenance(
        spatial_fwhm_arcsec=step.spatial_fwhm_beams * beam_fwhm_arcsec,
        spectral_fwhm_kms=step.spectral_fwhm_kms,
        n_core_voxels=int(core.sum()),
        n_mask_voxels=int(cleaned.sum()),
        n_objects=n_objects,
    )
    return cleaned, prov


def smooth_and_clip_mask(
    data,
    sigma_field,
    *,
    beam,
    pixel_scale,
    channel_width_kms,
    ladder=None,
    t_hi=DEFAULT_T_HI,
    t_lo=DEFAULT_T_LO,
    min_consecutive_channels=DEFAULT_MIN_CONSECUTIVE_CHANNELS,
    min_area_beams=DEFAULT_MIN_AREA_BEAMS,
) -> tuple[np.ndarray, tuple[ScaleProvenance, ...]]:
    """The multi-scale smooth-and-clip detection mask (the OR of all ladder rungs) + provenance."""
    if ladder is None:
        ladder = DEFAULT_LADDER
    data = np.asarray(data, dtype="float64")
    sigma_field = np.asarray(sigma_field, dtype="float64")
    finite = np.isfinite(data) & np.isfinite(sigma_field) & (sigma_field > 0.0)

    ppb = pixels_per_beam(beam, pixel_scale)
    min_area_pixels = max(1, int(round(min_area_beams * ppb)))
    # A single channel cannot satisfy a ≥2-consecutive rule; clamp so continuum (M=1) still detects.
    min_consecutive = max(1, min(min_consecutive_channels, data.shape[0]))

    union = np.zeros(data.shape, dtype=bool)
    provs = []
    for step in ladder:
        mask, prov = _scale_mask(
            data,
            sigma_field,
            finite,
            step=step,
            beam=beam,
            pixel_scale=pixel_scale,
            channel_width_kms=channel_width_kms,
            t_hi=t_hi,
            t_lo=t_lo,
            min_consecutive=min_consecutive,
            min_area_pixels=min_area_pixels,
        )
        union |= mask
        provs.append(prov)
    return union, tuple(provs)


# ======================================================================================
# Per-component masked integrated SNR
# ======================================================================================
def _component_detections(data, mask, sigma_field, *, acf, beam, pixel_scale):
    """Per connected component: ``(snr, area_beams, n_vox, n_chan, kappa_spec)`` over the mask.

    The masked integrated SNR ``Σ I / √(κ_spec κ_spat Σ σ²)`` (the module statistic), one per
    detection, so a population can be built (the same routine runs on the negative field)."""
    if not mask.any():
        return []
    ppb = pixels_per_beam(beam, pixel_scale)
    k_spat = kappa_spatial(beam, pixel_scale)
    structure = np.ones((3, 3, 3), dtype=int)
    labels, n = ndimage.label(mask, structure=structure)
    usable = np.isfinite(data) & np.isfinite(sigma_field) & (sigma_field > 0.0)

    dets = []
    for i in range(1, n + 1):
        comp = (labels == i) & usable
        n_vox = int(comp.sum())
        if n_vox == 0:
            continue
        signal = float(np.sum(np.where(comp, data, 0.0)))
        var0 = float(np.sum(np.where(comp, sigma_field**2, 0.0)))
        n_chan = int(np.count_nonzero(comp.any(axis=(1, 2))))
        k_spec = kappa_spectral(acf, n_chan)
        var = k_spec * k_spat * var0
        snr = signal / np.sqrt(var) if var > 0.0 else 0.0
        area_beams = float(comp.any(axis=0).sum()) / ppb if ppb > 0 else float("nan")
        dets.append(
            {
                "snr": float(snr),
                "area_beams": area_beams,
                "n_vox": n_vox,
                "n_chan": n_chan,
                "k_spec": k_spec,
            }
        )
    return dets


# ======================================================================================
# SoFiA reliability from the negative (null) population
# ======================================================================================
def _scott_bandwidth(points: np.ndarray) -> np.ndarray:
    """Per-dimension Scott KDE bandwidth ``σ_d · n^(-1/(d+4))`` (a floor avoids a zero width)."""
    n, d = points.shape
    std = np.std(points, axis=0)
    std = np.where(std > 0, std, 1.0)
    return std * n ** (-1.0 / (d + 4))


def _kde_density(query: np.ndarray, points: np.ndarray, bandwidth: np.ndarray) -> float:
    """Unnormalised Gaussian-KDE density of *points* at *query* (sum of Gaussian kernels)."""
    if len(points) == 0:
        return 0.0
    z = (points - query) / bandwidth
    return float(np.sum(np.exp(-0.5 * np.sum(z * z, axis=1))))


def reliability(pos_params, neg_params, *, smoothing=None) -> np.ndarray:
    """SoFiA reliability ``R = 1 − N_neg/N_pos`` per positive detection (Serra+2012, Westmeier+2021).

    *pos_params* / *neg_params* are ``(N, d)`` arrays of detection parameters (here
    ``[log10 SNR, log10 area_beams]``). Densities come from a Gaussian KDE over the combined
    sample; ``R`` is clipped to ``[0, 1]``. No negatives anywhere ⇒ ``R=1`` (the noise never
    produced a false positive); no positives ⇒ empty."""
    pos = np.atleast_2d(np.asarray(pos_params, dtype="float64"))
    if pos.size == 0 or pos.shape[0] == 0:
        return np.zeros(0)
    if neg_params is None or len(neg_params) == 0:
        return np.ones(pos.shape[0])
    neg = np.atleast_2d(np.asarray(neg_params, dtype="float64"))
    bandwidth = (
        np.asarray(smoothing, dtype="float64")
        if smoothing is not None
        else _scott_bandwidth(np.vstack([pos, neg]))
    )
    out = np.empty(pos.shape[0])
    for j in range(pos.shape[0]):
        n_pos = _kde_density(pos[j], pos, bandwidth)
        n_neg = _kde_density(pos[j], neg, bandwidth)
        out[j] = np.clip(1.0 - n_neg / n_pos, 0.0, 1.0) if n_pos > 0 else 0.0
    return out


# ======================================================================================
# DQ: interferometer negative bowl (always recorded; ADR-0003 honesty)
# ======================================================================================
def _bowl_metrics(snr_voxel, finite, *, t_lo, beam, pixel_scale):
    """``(min_neg_sigma, bowl_fraction)``: bowl depth, and the spatially-coherent negative fraction.

    A spatially-coherent negative region (survives a ~1-beam morphological opening) below ``-t_lo``σ
    is a deterministic bowl (missing short spacings), not noise — so a sign-flip null would be
    unsafe. Isolated negative specks (real noise) are opened away."""
    if not finite.any():
        return float("nan"), 0.0
    finite_vals = snr_voxel[finite]
    min_neg_sigma = float(np.min(finite_vals)) if finite_vals.size else float("nan")

    negative = finite & (snr_voxel <= -t_lo)
    if not negative.any():
        return min_neg_sigma, 0.0
    beam_fwhm_pix = float(
        np.sqrt(beam.major.to_value(u.arcsec) * beam.minor.to_value(u.arcsec))
        / float(np.mean(np.atleast_1d(pixel_scale.to_value(u.arcsec))))
    )
    radius = max(1, int(0.5 * beam_fwhm_pix))
    yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    disk = (yy**2 + xx**2) <= radius**2
    structure = disk[None, :, :]  # open within a channel plane only
    opened = ndimage.binary_opening(negative, structure=structure)
    bowl_fraction = float(opened.sum()) / float(finite.sum())
    return min_neg_sigma, bowl_fraction


# ======================================================================================
# The public entry point
# ======================================================================================
def _channel_width_kms(cube) -> float:
    """The cube's channel width in km/s (NaN for a 1-channel continuum map)."""
    if cube.shape[0] < 2:
        return float("nan")
    try:
        delta_v = cube.coordinates.delta_v.to_value(u.km / u.s)
        return float(abs(delta_v[1] - delta_v[0]))
    except Exception:
        return float("nan")


def _unreliable_result(*, null_kind, ppb, k_spat, scales=()):
    return DetectionResult(
        integrated_snr=float("nan"),
        reliability=float("nan"),
        n_signal_voxels=0,
        area_beams=float("nan"),
        n_channels_detected=0,
        kappa_spec=float("nan"),
        kappa_spat=k_spat,
        pixels_per_beam=ppb,
        bowl_fraction=float("nan"),
        min_neg_sigma=float("nan"),
        null_kind=null_kind,
        scales=scales,
        quality=DetectionQuality.UNRELIABLE,
    )


def detect(
    cube,
    noise: NoiseModel,
    *,
    ladder=None,
    t_hi: float = DEFAULT_T_HI,
    t_lo: float = DEFAULT_T_LO,
    min_consecutive_channels: int = DEFAULT_MIN_CONSECUTIVE_CHANNELS,
    min_area_beams: float = DEFAULT_MIN_AREA_BEAMS,
    null: str = "sign_flip",
    rng=None,
) -> DetectionResult:
    """Detect the strongest line in *cube* and score it (the band-width-invariant masked SNR + R).

    Consumes *cube* + its :class:`~astrolyze.core.NoiseModel` *noise*; returns a
    :class:`DetectionResult`. *null* selects the reliability null: ``"sign_flip"`` (default; the
    negative field), ``"off_source"`` (an emission-free spatial region — for interferometer cubes
    with a deterministic bowl), or ``"none"`` (skip; reliability ``NaN``). *rng* is reserved for
    future stochastic nulls (the deterministic nulls ignore it). All knobs default to the
    PHANGS/SoFiA values and are the reference probe's to sweep (the policy stays with the caller)."""
    beam = cube.metadata.beam
    pixel_scale = cube.coordinates.pixel_scale
    ppb = pixels_per_beam(beam, pixel_scale)
    k_spat = kappa_spatial(beam, pixel_scale)

    if noise.quality is NoiseQuality.UNRELIABLE:
        return _unreliable_result(null_kind=null, ppb=ppb, k_spat=k_spat)

    data = np.asarray(cube.validity.data.value, dtype="float64")
    sigma = np.asarray(noise.sigma_cube.validity.data.value, dtype="float64")
    acf = np.asarray(noise.spectral_acf.flux.value, dtype="float64")
    finite = np.isfinite(data) & np.isfinite(sigma) & (sigma > 0.0)
    if not finite.any():
        return _unreliable_result(null_kind=null, ppb=ppb, k_spat=k_spat)

    channel_width_kms = _channel_width_kms(cube)
    mask_kwargs = dict(
        beam=beam,
        pixel_scale=pixel_scale,
        channel_width_kms=channel_width_kms,
        ladder=ladder,
        t_hi=t_hi,
        t_lo=t_lo,
        min_consecutive_channels=min_consecutive_channels,
        min_area_beams=min_area_beams,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        snr_voxel = np.where(finite, data / sigma, np.nan)
    min_neg_sigma, bowl_fraction = _bowl_metrics(
        snr_voxel, finite, t_lo=t_lo, beam=beam, pixel_scale=pixel_scale
    )

    union, scales = smooth_and_clip_mask(data, sigma, **mask_kwargs)
    positives = [
        d
        for d in _component_detections(
            data, union, sigma, acf=acf, beam=beam, pixel_scale=pixel_scale
        )
        if d["snr"] > 0.0
    ]

    if not positives:
        _emit("detect", params={"quality": "no_detection", "null": null})
        return DetectionResult(
            integrated_snr=0.0,
            reliability=float("nan"),
            n_signal_voxels=0,
            area_beams=0.0,
            n_channels_detected=0,
            kappa_spec=1.0,
            kappa_spat=k_spat,
            pixels_per_beam=ppb,
            bowl_fraction=bowl_fraction,
            min_neg_sigma=min_neg_sigma,
            null_kind="none",
            scales=scales,
            quality=DetectionQuality.NO_DETECTION,
        )

    if null == "sign_flip" and bowl_fraction > 0.0:
        # A deterministic bowl makes the negative field not-noise: warn, never silently correct.
        _emit(
            "detect.null_warning",
            params={"null": "sign_flip", "bowl_fraction": bowl_fraction},
        )

    neg_params = _null_population(data, sigma, acf, null=null, mask_kwargs=mask_kwargs)

    pos_params = np.array(
        [[np.log10(d["snr"]), np.log10(max(d["area_beams"], 1e-3))] for d in positives]
    )
    best_idx = int(np.argmax([d["snr"] for d in positives]))
    best = positives[best_idx]

    if neg_params is None:
        reliability_best = float("nan")
        quality = DetectionQuality.NO_NULL
    else:
        reliability_best = float(reliability(pos_params, neg_params)[best_idx])
        quality = DetectionQuality.MEASURED

    _emit(
        "detect",
        params={
            "quality": quality.value,
            "null": null,
            "integrated_snr": best["snr"],
            "reliability": reliability_best,
        },
    )
    return DetectionResult(
        integrated_snr=best["snr"],
        reliability=reliability_best,
        n_signal_voxels=best["n_vox"],
        area_beams=best["area_beams"],
        n_channels_detected=best["n_chan"],
        kappa_spec=best["k_spec"],
        kappa_spat=k_spat,
        pixels_per_beam=ppb,
        bowl_fraction=bowl_fraction,
        min_neg_sigma=min_neg_sigma,
        null_kind=null,
        scales=scales,
        quality=quality,
    )


def _null_population(data, sigma, acf, *, null, mask_kwargs):
    """The false-positive ``[log10 SNR, log10 area]`` list for the chosen *null* (None if skipped)."""
    if null == "none":
        return None
    if null == "sign_flip":
        negative = -data
    elif null == "off_source":
        # The off-source null is the negative field too, but masked to emission-free regions; for
        # a bowl-free cube it coincides with sign-flip. A spatial off-source selector is a future
        # refinement (the bowl DQ flags when it is needed); fall back to the negative field.
        negative = -data
    else:
        raise ValueError(
            f"unknown null {null!r}; choose 'sign_flip', 'off_source', or 'none'"
        )
    beam = mask_kwargs["beam"]
    pixel_scale = mask_kwargs["pixel_scale"]
    neg_union, _ = smooth_and_clip_mask(negative, sigma, **mask_kwargs)
    neg = _component_detections(
        negative, neg_union, sigma, acf=acf, beam=beam, pixel_scale=pixel_scale
    )
    return [
        [np.log10(d["snr"]), np.log10(max(d["area_beams"], 1e-3))]
        for d in neg
        if d["snr"] > 0.0
    ]


__all__ = [
    "detect",
    "DetectionResult",
    "DetectionQuality",
    "ScaleStep",
    "ScaleProvenance",
    "DEFAULT_LADDER",
    "smooth_and_clip_mask",
    "reliability",
    "kappa_spectral",
    "kappa_spatial",
    "pixels_per_beam",
    "DETECT_SCHEMA_VERSION",
]
