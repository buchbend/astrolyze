"""The context-carrying ``NoiseModel`` companion + its pluggable estimator suite (#27/#28).

A :class:`NoiseModel` is a *companion* object (ADR-0004): it is produced by
:meth:`Cube.estimate_noise` and exposes a cube's noise as **first-class astrolyze products** —
a σ ``Cube``, a σ ``Map``, a σ ``Spectrum``, the spectral autocorrelation, and a single robust
σ ``Quantity``. It is not itself one of the data wrappers (it is not a single array); instead it
holds the parent :class:`Cube` and **reuses** that cube's type-transition machinery to build the
products, so the beam + rest frequency + velocity convention flow onto them for free.

Two house rules shape it:

- **No silent physics (ADR-0003).** The estimator measures σ from signal-free data; if there is
  no usable signal-free data it flags the model :data:`NoiseQuality.UNRELIABLE` and leaves σ as
  ``NaN`` — it never fabricates a number.
- **Stay thin.** The statistics delegate to :func:`astropy.stats.mad_std` / :mod:`numpy`; the
  σ products are built by slicing/rebuilding through the parent ``Cube`` (not a re-derived
  cube-building path).

Representation. The default is the **separable** form ``σ(x,y,v) = σ_xy(x,y)·σ_v(v)/σ_0``
(a rank-1 outer product) plus the scalar σ_0, the spectral ACF and the beam — compact and
analytically sufficient when the noise is (close to) separable. When a separability test fails
(the per-voxel σ field is not rank-1) the model falls back to a stored **full 3D σ-cube**. Both
representations reconstruct the same :attr:`sigma_cube`.

Persistence. :meth:`to_zarr_companion` / :meth:`from_zarr_companion` round-trip the model as a
**companion group** inside the cube's #23 Zarr store, carrying the estimator method + a schema
version as provenance (the heavy lifting is in :mod:`astrolyze.io.zarr_backend`).
"""

from __future__ import annotations

import enum
from dataclasses import dataclass
from typing import Callable

import numpy as np
import astropy.units as u
from astropy.stats import mad_std, sigma_clipped_stats

from astrolyze.io import Metadata

from ._base import _emit

#: Schema version of the on-disk companion-group layout (provenance, bumped on layout change).
NOISE_SCHEMA_VERSION = 1

#: The default estimator (the robust ``mad_std``); the pluggable suite (#28) adds more.
DEFAULT_METHOD = "mad_std"

#: Number of spectral blocks the separability test measures a spatial σ map in. Each block needs
#: enough channels for a stable robust σ per pixel, so the per-block maps reflect the real
#: spatial *shape* of the noise rather than estimation scatter.
SEPARABILITY_NBLOCKS = 4

#: Separability threshold: the noise is treated as separable (σ_xy·σ_v) when the ratio of the
#: second to the first singular value of the per-block spatial-σ matrix is at or below this — i.e.
#: the matrix is effectively rank-1, so one spatial shape times one spectral shape suffices. Above
#: it the spatial shape varies with frequency (rank > 1) and we keep the full 3D σ-cube. A pre-1.0
#: default (DECISION, issue #27): stationary noise sits well below it, a frequency-dependent
#: spatial pattern well above.
SEPARABILITY_SV_RATIO = 0.30


class NoiseQuality(enum.Enum):
    """How trustworthy a :class:`NoiseModel`'s σ is (no silent physics — ADR-0003).

    - ``MEASURED`` — σ was measured from real signal-free data (the normal case);
    - ``PROPAGATED`` — σ was propagated analytically from an upstream σ (reserved for #32);
    - ``APPROXIMATE`` — σ is an order-of-magnitude estimate, not a clean measurement (reserved);
    - ``UNRELIABLE`` — there was no usable signal-free data, so σ is **not** a measurement and is
      left ``NaN`` rather than fabricated.
    """

    MEASURED = "measured"
    PROPAGATED = "propagated"
    APPROXIMATE = "approximate"
    UNRELIABLE = "unreliable"


@dataclass(frozen=True)
class SeparableNoise:
    """The compact separable representation: ``σ(x,y,v) = σ_xy(x,y)·σ_v(v)/σ_0``.

    Stores the spatial pattern ``sigma_xy`` (2D), the spectral pattern ``sigma_v`` (1D), the
    scalar ``scalar`` (σ_0), and the spectral autocorrelation ``acf`` — all bare numpy arrays in
    the cube's intensity unit (the unit lives on the model's metadata, not duplicated here)."""

    sigma_xy: np.ndarray  # (ny, nx)
    sigma_v: np.ndarray  # (nz,)
    scalar: float
    acf: np.ndarray  # (nlag,)

    def reconstruct(self, shape) -> np.ndarray:
        """The full per-voxel σ field as the rank-1 outer product ``σ_v ⊗ σ_xy / σ_0``."""
        if not np.isfinite(self.scalar) or self.scalar == 0.0:
            return np.full(shape, np.nan)
        return self.sigma_v[:, None, None] * self.sigma_xy[None, :, :] / self.scalar


@dataclass(frozen=True)
class FullNoise:
    """The full 3D σ-cube fallback: a dense per-voxel σ field (cube's intensity unit).

    Used when the separability test fails — the σ field is not a rank-1 outer product, so the
    compact form would mis-state it. We still derive σ_xy / σ_v / scalar / ACF *views* from it so
    the model exposes the same products either way; only ``.reconstruct`` differs."""

    sigma_field: np.ndarray  # (nz, ny, nx)
    sigma_xy: np.ndarray
    sigma_v: np.ndarray
    scalar: float
    acf: np.ndarray

    def reconstruct(self, shape) -> np.ndarray:
        return np.asarray(self.sigma_field).reshape(shape)


class NoiseModel:
    """A cube's noise as first-class astrolyze products (issue #27, ADR-0003/0004/0005).

    Built by :meth:`Cube.estimate_noise`. Holds the parent :class:`Cube` (so it can *reuse* the
    Cube->Map/Spectrum transitions to build context-carrying products), the chosen representation
    (separable or full), the estimator ``method`` + schema ``version`` (provenance), and the
    :class:`NoiseQuality` flag."""

    _viz_function = "plot_noise"

    def __init__(
        self,
        cube,
        representation,
        *,
        method: str,
        quality: NoiseQuality,
        version: int = NOISE_SCHEMA_VERSION,
    ):
        self._cube = cube
        self._rep = representation
        self.method = method
        self.quality = quality
        self.version = version

    # -- context shortcuts (read off the parent cube's metadata, never duplicated) ------
    @property
    def metadata(self) -> Metadata:
        return self._cube.metadata

    @property
    def beam(self):
        return self.metadata.beam

    @property
    def rest_frequency(self):
        return self.metadata.rest_frequency

    @property
    def velocity_convention(self):
        return self.metadata.velocity_convention

    @property
    def unit(self) -> u.UnitBase:
        return self._cube.unit

    @property
    def is_separable(self) -> bool:
        """Whether the model stores the compact separable σ_xy·σ_v form (vs a full σ-cube)."""
        return isinstance(self._rep, SeparableNoise)

    # -- first-class products (transitions reuse the parent Cube's machinery) -----------
    @property
    def sigma_cube(self):
        """The full 3D σ(x,y,v) as a :class:`Cube` in the cube's intensity unit.

        Built by handing the reconstructed per-voxel σ field back through the parent cube's
        own data seam (``_with_data``), so the WCS + beam + context come along unchanged — the
        Cube->Cube machinery is reused, not re-derived (ADR-0004)."""
        field = self._rep.reconstruct(self._cube.shape) * self.unit
        return self._cube._with_data(field)

    @property
    def sigma_map(self):
        """The spatial σ_xy(x,y) as a :class:`Map` in the cube's intensity unit.

        Reuses the Cube->Map transition: a single channel slice gives a context-carrying ``Map``
        of the right WCS/shape, whose data we replace with the spatial σ pattern."""
        template = self._cube[0]  # a channel Map carrying the spatial WCS + context
        return template._with_data(np.asarray(self._rep.sigma_xy) * self.unit)

    @property
    def sigma_spectrum(self):
        """The per-channel σ_v(v) as a :class:`Spectrum` in the cube's intensity unit.

        Reuses the Cube->Spectrum transition: a spatial pencil gives a context-carrying
        ``Spectrum`` on the right spectral axis, whose flux we replace with the spectral σ."""
        template = self._cube[
            :, 0, 0
        ]  # a Spectrum on the cube's spectral axis + context
        return template._with_data(np.asarray(self._rep.sigma_v) * self.unit)

    @property
    def scalar(self) -> u.Quantity:
        """A single robust σ as a :class:`~astropy.units.Quantity` (``NaN`` if UNRELIABLE)."""
        return float(self._rep.scalar) * self.unit

    @property
    def spectral_acf(self):
        """The per-cube spectral noise autocorrelation as a :class:`Spectrum` (dimensionless).

        The ACF is a correlation, not an intensity, so it is unitless; it is returned on the
        cube's spectral axis (lag 0 = the first channel) as a context-carrying ``Spectrum``."""
        template = self._cube[:, 0, 0]
        acf = np.asarray(self._rep.acf, dtype=float)
        # Pad/trim to the spectral-axis length so it rides the template's axis; ACF is the
        # per-lag correlation, lag index == channel index from 0.
        nz = template.flux.size
        padded = np.full(nz, np.nan)
        padded[: min(nz, acf.size)] = acf[: min(nz, acf.size)]
        return template._with_data(padded * u.dimensionless_unscaled)

    # -- spatial cutout: slice the σ pattern first, then reconstruct only the tube ------
    def cutout(self, y0: int, y1: int, x0: int, x1: int) -> "NoiseModel":
        """Spatial cutout of the noise model to ``[y0:y1, x0:x1]``, full spectral kept.

        Mirrors ``Cube[:, y0:y1, x0:x1]``: slices the spatial σ pattern and the parent cube the
        same way so the products carry the sub-WCS/shape; the spectral σ_v, scalar, and ACF are
        per-cube and unchanged. Lets a tube read its per-voxel σ without reconstructing the whole
        cube — :attr:`sigma_cube` eagerly reconstructs the FULL σ field, so slicing the *model*
        first (cheap: only σ_xy / σ_field are sliced) means the cutout's ``sigma_cube``
        reconstructs just the tube, not the parent cube (the spatial half of #32's tube reads)."""
        sub_cube = self._cube[:, y0:y1, x0:x1]
        if self.is_separable:
            rep = SeparableNoise(
                sigma_xy=self._rep.sigma_xy[y0:y1, x0:x1],
                sigma_v=self._rep.sigma_v,
                scalar=self._rep.scalar,
                acf=self._rep.acf,
            )
        else:
            rep = FullNoise(
                sigma_field=self._rep.sigma_field[:, y0:y1, x0:x1],
                sigma_xy=self._rep.sigma_xy[y0:y1, x0:x1],
                sigma_v=self._rep.sigma_v,
                scalar=self._rep.scalar,
                acf=self._rep.acf,
            )
        _emit("noise.cutout", params={"y0": y0, "y1": y1, "x0": x0, "x1": x1})
        return NoiseModel(
            sub_cube,
            rep,
            method=self.method,
            quality=self.quality,
            version=self.version,
        )

    # -- unit hub (ADR-0003c): .to() follows the data unit ------------------------------
    def to(self, unit, **kwargs) -> "NoiseModel":
        """Convert the noise to ``unit`` through the cube's unit hub; return a new model.

        σ is an intensity, so it converts exactly as the cube's data does — :meth:`Cube.to`
        supplies the object's beam / rest frequency / convention (ADR-0003c). We convert the
        reconstructed σ-cube and re-estimate the (now-trivially-separable) representation in the
        new unit, so every product follows the data unit."""
        converted_cube = self.sigma_cube.to(unit, **kwargs)
        rep = _representation_from_field(
            np.asarray(converted_cube._data_quantity.value),
            self._rep,
        )
        _emit("noise.to", params={"unit": str(unit)})
        return NoiseModel(
            converted_cube,
            rep,
            method=self.method,
            quality=self.quality,
            version=self.version,
        )

    # -- viz seam (ADR-0005): thin sugar over the free plot function --------------------
    def plot(self, **kwargs):
        """Plot the noise via :mod:`astrolyze.viz` (the free ``plot_noise`` function)."""
        from astrolyze import viz

        plotter = getattr(viz, self._viz_function, None)
        if plotter is None:  # pragma: no cover - defensive, viz ships plot_noise
            raise NotImplementedError(
                f"plotting needs the viz layer ({self._viz_function})"
            )
        result = plotter(self, **kwargs)
        _emit("plot", params={"viz_function": self._viz_function})
        return result

    # -- Zarr companion group (issue #23 store) -----------------------------------------
    def to_zarr_companion(self, store):
        """Write this model as a companion group alongside the cube's Zarr *store*; return it.

        Delegates the byte I/O to :mod:`astrolyze.io.zarr_backend`; the group carries the
        estimator method + schema version as provenance plus the representation arrays."""
        from astrolyze.io.zarr_backend import _save_noise_companion

        return _save_noise_companion(self, store)

    @classmethod
    def from_zarr_companion(cls, store) -> "NoiseModel":
        """Reconstruct a model from the companion group inside the cube's Zarr *store*.

        Reloads the parent :class:`Cube` from the store (so the products are again
        context-carrying) and the representation + provenance from the companion group."""
        from astrolyze.core import Cube
        from astrolyze.io.zarr_backend import _load_noise_companion

        cube = Cube.from_zarr(store)
        rep, method, quality, version = _load_noise_companion(store)
        return cls(
            cube, rep, method=method, quality=NoiseQuality(quality), version=version
        )

    # -- survey weight / RMS-map ingest (issue #28) -------------------------------------
    # Built FROM a published map, not estimated from the cube: a survey commonly ships a noise
    # (RMS) or weight map, and that is authoritative for σ where it exists. astrolyze adds only
    # the routing weight->σ (σ = 1/√w) and the quality flag (un-observed pixels carry no σ ->
    # UNRELIABLE, never a fabricated number; ADR-0003). The spatial σ map is promoted to a
    # separable model with a flat spectral profile (one σ per pixel, constant across channels):
    # a published map is a 2D σ field, so the noise is separable by construction.
    @classmethod
    def from_rms_map(cls, cube, rms_map) -> "NoiseModel":
        """Build a :class:`NoiseModel` from a survey **RMS map** (σ is the map itself).

        *rms_map* is the published per-pixel noise σ — a :class:`~astrolyze.core.Map` or a bare
        array in the cube's intensity unit. The model carries *cube*'s context and exposes σ as
        the usual first-class products. Pixels that are NaN (un-observed) carry no σ; an entirely
        un-observed map yields :data:`NoiseQuality.UNRELIABLE` (ADR-0003)."""
        sigma_xy = _map_values(rms_map, cube.unit)
        rep, quality = _representation_from_sigma_map(cube, sigma_xy)
        _emit("estimate_noise", params={"method": "rms_map", "quality": quality.value})
        return cls(cube, rep, method="rms_map", quality=quality)

    @classmethod
    def from_weight_map(cls, cube, weight_map) -> "NoiseModel":
        """Build a :class:`NoiseModel` from a survey **weight map** (σ = 1/√w).

        *weight_map* is the published per-pixel inverse-variance weight — a
        :class:`~astrolyze.core.Map` or a bare array (its unit is the inverse of the cube
        variance, e.g. ``K**-2``). σ = ``1/√w``; a zero/negative/NaN weight is an un-observed
        pixel that carries no σ (left ``NaN``). An entirely un-observed map yields
        :data:`NoiseQuality.UNRELIABLE` rather than a fabricated σ (ADR-0003)."""
        # Weights are inverse-variance; the unit is unit(cube)**-2. Read the bare values in that
        # unit, then σ = 1/√w lands back in the cube's intensity unit.
        weights = _map_values(weight_map, cube.unit**-2)
        with np.errstate(divide="ignore", invalid="ignore"):
            positive = np.isfinite(weights) & (weights > 0.0)
            sigma_xy = np.full(weights.shape, np.nan)
            sigma_xy[positive] = 1.0 / np.sqrt(weights[positive])
        rep, quality = _representation_from_sigma_map(cube, sigma_xy)
        _emit(
            "estimate_noise", params={"method": "weight_map", "quality": quality.value}
        )
        return cls(cube, rep, method="weight_map", quality=quality)

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        kind = "separable" if self.is_separable else "full"
        return f"<NoiseModel {obj} {kind} [{self.unit}] {self.quality.value}>"


# -- pluggable estimator suite (issue #28) ---------------------------------------------
# An estimator is a callable ``(data: np.ndarray) -> (representation, NoiseQuality)`` where
# ``data`` is the cube's bare values (nz, ny, nx) in its intensity unit and ``representation``
# is a SeparableNoise / FullNoise (typically built via _build_representation so the #27
# separability test selects the branch). The registry below is the PLUGGABILITY seam: a caller
# registers a new estimator with register_estimator(name, fn) and selects it through
# estimate_noise(method=name) WITHOUT editing this module — the routing is open, not a closed
# if/elif ladder (ADR-0004: astrolyze adds the routing + the quality flag, the statistics stay
# in numpy/astropy.stats).
Estimator = Callable[[np.ndarray], "tuple[SeparableNoise | FullNoise, NoiseQuality]"]

#: The estimator registry: method name -> estimator callable. Built-ins are registered at import
#: (see bottom of module); caller estimators are added via :func:`register_estimator`.
_ESTIMATORS: dict[str, Estimator] = {}


def register_estimator(name: str, estimator: Estimator) -> None:
    """Register a noise *estimator* under *name* so ``estimate_noise(method=name)`` selects it.

    The pluggability seam (issue #28). An *estimator* is a callable
    ``(data) -> (representation, NoiseQuality)`` taking the cube's bare values ``(nz, ny, nx)``
    and returning a :class:`SeparableNoise` / :class:`FullNoise` plus its
    :class:`NoiseQuality`. Registering lets a caller add an estimator from outside the package —
    no core edit. Re-registering a name overrides it (pre-1.0, last write wins)."""
    if not callable(estimator):
        raise TypeError(
            f"estimator for {name!r} must be callable (data) -> (representation, NoiseQuality)"
        )
    _ESTIMATORS[name] = estimator


def available_estimators() -> tuple[str, ...]:
    """The registered estimator method names (built-ins plus any caller-supplied)."""
    return tuple(_ESTIMATORS)


def estimate(cube, *, method: str = DEFAULT_METHOD) -> NoiseModel:
    """Estimate a :class:`NoiseModel` for *cube* with the registered *method* (issue #27/#28).

    The *method* is looked up in the pluggable estimator registry (see
    :func:`register_estimator`): the default robust ``"mad_std"``, the signal-free-channel
    ``"rms"``, the σ-clipped ``"mad"``, or any caller-registered estimator. An unknown *method*
    is refused rather than silently substituted (ADR-0003)."""
    estimator = _ESTIMATORS.get(method)
    if estimator is None:
        known = ", ".join(sorted(_ESTIMATORS)) or "(none registered)"
        raise ValueError(
            f"unknown noise estimator method {method!r}; registered methods: {known}. "
            "Add one with astrolyze.core.register_estimator(name, fn); astrolyze never "
            "silently substitutes an estimator (ADR-0003)"
        )
    data = np.asarray(cube._data_quantity.value, dtype="float64")
    rep, quality = estimator(data)
    _emit("estimate_noise", params={"method": method, "quality": quality.value})
    return NoiseModel(cube, rep, method=method, quality=quality)


def _unreliable_representation(shape) -> SeparableNoise:
    """A NaN-filled representation for the no-usable-data case (ADR-0003: never fabricate σ).

    Shared by every estimator so "no signal-free data" looks identical whichever method asked —
    σ is left ``NaN`` and the caller flags the model :data:`NoiseQuality.UNRELIABLE`."""
    nz = shape[0]
    shape_xy = shape[1:]
    return SeparableNoise(
        sigma_xy=np.full(shape_xy, np.nan),
        sigma_v=np.full(nz, np.nan),
        scalar=np.nan,
        acf=np.array([np.nan]),
    )


def _estimate_mad_std(data: np.ndarray):
    """The robust σ estimator: per-voxel scatter via :func:`astropy.stats.mad_std`.

    Returns ``(representation, quality)``. With no usable (finite) data the quality is
    :data:`NoiseQuality.UNRELIABLE` and σ is left ``NaN`` — never fabricated (ADR-0003)."""
    if not np.isfinite(data).any():
        # No usable signal-free data: refuse to invent a σ (ADR-0003).
        return _unreliable_representation(data.shape), NoiseQuality.UNRELIABLE

    # Robust scatter on each axis-collapse (NaN-ignoring): spatial σ per pixel (over v),
    # spectral σ per channel (over the plane), and a single scalar over the whole cube.
    sigma_xy = mad_std(data, axis=0, ignore_nan=True)
    sigma_v = mad_std(data.reshape(data.shape[0], -1), axis=1, ignore_nan=True)
    scalar = float(mad_std(data, ignore_nan=True))
    acf = _spectral_acf(data)

    rep = _build_representation(data, sigma_xy, sigma_v, scalar, acf)
    return rep, NoiseQuality.MEASURED


#: σ-clip threshold (n·σ) and iteration count for the robust ``"mad"`` estimator. A few clips at
#: 3σ rejects line/source voxels before the per-axis σ is measured; this is the conventional
#: radio-line default (DECISION, issue #28, pre-1.0) and keeps astrolyze thin over astropy.stats.
SIGMA_CLIP_SIGMA = 3.0
SIGMA_CLIP_ITERS = 5


def _estimate_mad_sigma_clip(data: np.ndarray):
    """Robust σ via iterative σ-clipping (:func:`astropy.stats.sigma_clipped_stats`).

    Built on astropy.stats (stay thin): σ-clipping rejects bright line/source voxels, so the
    measured σ reflects the underlying noise even when the cube is not signal-free. The per-axis
    σ_xy / σ_v use the σ-clipped MAD std on each collapse; the scalar is the whole-cube
    σ-clipped std. No usable finite data ⇒ :data:`NoiseQuality.UNRELIABLE` (ADR-0003)."""
    if not np.isfinite(data).any():
        return _unreliable_representation(data.shape), NoiseQuality.UNRELIABLE

    sigma_xy = _sigma_clipped_std_along(data, axis=0)
    sigma_v = _sigma_clipped_std_along(data.reshape(data.shape[0], -1), axis=1)
    _, _, scalar = sigma_clipped_stats(
        data, sigma=SIGMA_CLIP_SIGMA, maxiters=SIGMA_CLIP_ITERS, stdfunc="mad_std"
    )
    acf = _spectral_acf(data)

    rep = _build_representation(data, sigma_xy, sigma_v, float(scalar), acf)
    return rep, NoiseQuality.MEASURED


def _sigma_clipped_std_along(data: np.ndarray, *, axis: int) -> np.ndarray:
    """The σ-clipped robust std along *axis* (the per-pixel/per-channel σ for the ``"mad"`` path).

    :func:`astropy.stats.sigma_clipped_stats` returns ``(mean, median, std)``; we keep the std,
    measured with ``mad_std`` so it is the robust scatter after the clip."""
    _, _, std = sigma_clipped_stats(
        data,
        sigma=SIGMA_CLIP_SIGMA,
        maxiters=SIGMA_CLIP_ITERS,
        stdfunc="mad_std",
        axis=axis,
    )
    return np.asarray(std, dtype="float64")


#: How many extreme channels (by per-channel level) the signal-free RMS estimator rejects as
#: line-contaminated before taking the RMS over the survivors. A robust, parameter-light default
#: (DECISION, issue #28, pre-1.0): channels whose robust level sits far above the median channel
#: level carry signal; clipping them leaves the line-free band. Threshold is in robust σ.
RMS_SIGNAL_FREE_SIGMA = 3.0


def _estimate_rms(data: np.ndarray):
    """Signal-free-channel RMS: the per-voxel RMS over channels that carry no line signal.

    A channel's robust level (its σ-clipped mean intensity over the plane) tells line from
    line-free: channels whose level sits more than ``RMS_SIGNAL_FREE_SIGMA`` robust-σ above the
    median channel level are flagged signal-bearing and dropped. The RMS is then taken (with
    numpy) over the surviving signal-free channels only — so a bright line does not inflate σ.
    No usable signal-free channel ⇒ :data:`NoiseQuality.UNRELIABLE` (ADR-0003)."""
    if not np.isfinite(data).any():
        return _unreliable_representation(data.shape), NoiseQuality.UNRELIABLE

    free = _signal_free_channels(data)
    if not free.any():
        # Every channel looks signal-bearing: there is no line-free data to take an RMS over.
        return _unreliable_representation(data.shape), NoiseQuality.UNRELIABLE

    free_data = data[free]
    # Per-voxel RMS uses only the signal-free channels. The spatial σ_xy(x,y) is the RMS over
    # those channels per pixel; σ_v(v) is the per-channel RMS over the plane (every channel,
    # so the spectral profile spans the full axis); the scalar is the RMS over all free voxels.
    sigma_xy = _rms(free_data, axis=0)
    sigma_v = _rms(data.reshape(data.shape[0], -1), axis=1)
    scalar = float(_rms(free_data, axis=None))
    acf = _spectral_acf(data)

    rep = _build_representation(data, sigma_xy, sigma_v, scalar, acf)
    return rep, NoiseQuality.MEASURED


def _rms(data: np.ndarray, *, axis) -> np.ndarray | float:
    """NaN-ignoring root-mean-square along *axis* (``numpy`` only — astrolyze stays thin)."""
    return np.sqrt(np.nanmean(np.square(data), axis=axis))


def _signal_free_channels(data: np.ndarray) -> np.ndarray:
    """A boolean per-channel mask: ``True`` where the channel carries no line signal.

    Each channel's robust level is its σ-clipped mean over the plane. A channel is line-bearing
    when that level sits more than ``RMS_SIGNAL_FREE_SIGMA`` robust-σ above the median channel
    level; the rest are signal-free. Channels with no finite data are not signal-free."""
    nz = data.shape[0]
    flat = data.reshape(nz, -1)
    levels = np.full(nz, np.nan)
    for k in range(nz):
        column = flat[k]
        if not np.isfinite(column).any():
            continue
        mean, _, _ = sigma_clipped_stats(
            column, sigma=SIGMA_CLIP_SIGMA, maxiters=SIGMA_CLIP_ITERS, stdfunc="mad_std"
        )
        levels[k] = mean
    finite = np.isfinite(levels)
    if not finite.any():
        return np.zeros(nz, dtype=bool)
    median_level = np.nanmedian(levels)
    spread = mad_std(levels, ignore_nan=True)
    if not np.isfinite(spread) or spread == 0.0:
        # No channel-to-channel spread (e.g. flat noise, no line): every finite channel is free.
        return finite
    threshold = median_level + RMS_SIGNAL_FREE_SIGMA * spread
    return finite & (levels <= threshold)


def _block_sigma_maps(data: np.ndarray, nblocks: int) -> np.ndarray:
    """A spatial σ map per spectral block: the matrix ``M[block, pixel]`` (NaN-ignoring).

    Splitting the band into a few blocks and measuring a robust σ map in each gives enough
    samples per pixel for the *spatial shape* of the noise to emerge above estimation scatter —
    which is what the separability test needs (does that shape change with frequency?)."""
    nz = data.shape[0]
    edges = np.linspace(0, nz, nblocks + 1).astype(int)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi - lo < 2:
            continue
        rows.append(mad_std(data[lo:hi], axis=0, ignore_nan=True).reshape(-1))
    return np.asarray(rows, dtype="float64")


def _spectral_acf(data: np.ndarray) -> np.ndarray:
    """The per-cube spectral noise autocorrelation (normalised to 1.0 at zero lag).

    Averaged over the spatial plane: for each pixel's spectrum we take the biased
    autocorrelation, then average across pixels and normalise. White noise -> a spike at lag 0."""
    nz = data.shape[0]
    flat = data.reshape(nz, -1)
    acc = np.zeros(nz)
    count = 0
    for j in range(flat.shape[1]):
        column = flat[:, j]
        if not np.all(np.isfinite(column)):
            continue
        column = column - np.nanmean(column)
        full = np.correlate(column, column, mode="full")
        acc += full[nz - 1 :]  # non-negative lags
        count += 1
    if count == 0:
        return np.full(nz, np.nan)
    acf = acc / count
    if acf[0] != 0.0:
        acf = acf / acf[0]
    return acf


def _build_representation(data, sigma_xy, sigma_v, scalar, acf):
    """Pick the separable form when the per-block spatial-σ matrix is effectively rank-1;
    otherwise keep a full 3D σ-cube (ADR-0003: no silent over-compression of structured noise).

    The separable σ-field is the rank-1 outer product ``σ_v ⊗ σ_xy / σ_0``. The full fallback
    builds a dense per-voxel field whose spatial shape is allowed to vary with frequency,
    interpolated from the per-block spatial σ maps and scaled to the per-channel σ_v."""
    separable = SeparableNoise(
        sigma_xy=np.asarray(sigma_xy, dtype="float64"),
        sigma_v=np.asarray(sigma_v, dtype="float64"),
        scalar=scalar,
        acf=np.asarray(acf, dtype="float64"),
    )
    block_maps = _block_sigma_maps(data, SEPARABILITY_NBLOCKS)
    if _is_separable(block_maps):
        return separable
    sigma_field = _dense_sigma_field(data, block_maps, separable.sigma_v)
    return FullNoise(
        sigma_field=sigma_field,
        sigma_xy=separable.sigma_xy,
        sigma_v=separable.sigma_v,
        scalar=scalar,
        acf=separable.acf,
    )


def _is_separable(block_maps: np.ndarray) -> bool:
    """Whether the per-block spatial-σ matrix is effectively rank-1 (one frequency-independent
    spatial shape). Measured by the ratio of the second to the first singular value: a small
    ratio means a single spatial pattern explains every block, so σ_xy·σ_v is sufficient."""
    finite_rows = block_maps[np.all(np.isfinite(block_maps), axis=1)]
    if finite_rows.shape[0] < 2 or finite_rows.shape[1] < 2:
        # Too few blocks/pixels to detect a frequency-dependent shape: treat as separable
        # (the compact form is the conservative default; it cannot misstate a 1-block cube).
        return True
    singular_values = np.linalg.svd(finite_rows, compute_uv=False)
    if singular_values[0] == 0.0:
        return True
    return bool(singular_values[1] / singular_values[0] <= SEPARABILITY_SV_RATIO)


def _dense_sigma_field(data, block_maps, sigma_v) -> np.ndarray:
    """A dense per-voxel σ field for the non-separable fallback (cube's intensity unit).

    The spatial shape is allowed to change with frequency: each channel takes its spectral
    block's spatial σ map, renormalised so the channel's spatial mean matches the per-channel
    σ_v — so the field carries both the spectral profile (σ_v) and the frequency-dependent
    spatial pattern the compact form cannot."""
    nz, ny, nx = data.shape
    maps = block_maps.reshape(block_maps.shape[0], ny, nx)
    edges = np.linspace(0, nz, maps.shape[0] + 1).astype(int)
    field = np.empty((nz, ny, nx), dtype="float64")
    for block_index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        block_map = maps[block_index]
        block_level = np.nanmean(block_map)
        for k in range(lo, hi):
            level = sigma_v[k]
            if block_level and np.isfinite(block_level):
                field[k] = block_map * (level / block_level)
            else:
                field[k] = block_map
    return field


def _representation_from_field(field: np.ndarray, like) -> SeparableNoise | FullNoise:
    """Rebuild a representation from an already-converted σ-*field* (used by :meth:`NoiseModel.to`).

    The unit conversion is a per-voxel rescale, so it preserves the rank of the σ field; we keep
    the same branch (separable vs full) as the source representation and just rescale the stored
    arrays by the field's overall scale, so products stay consistent in the new unit."""
    scalar = float(mad_std(field, ignore_nan=True))
    sigma_xy = mad_std(field, axis=0, ignore_nan=True)
    sigma_v = mad_std(field.reshape(field.shape[0], -1), axis=1, ignore_nan=True)
    if isinstance(like, SeparableNoise):
        return SeparableNoise(
            sigma_xy=np.asarray(sigma_xy, dtype="float64"),
            sigma_v=np.asarray(sigma_v, dtype="float64"),
            scalar=scalar,
            acf=np.asarray(like.acf, dtype="float64"),  # ACF is unitless: unchanged
        )
    return FullNoise(
        sigma_field=np.asarray(field, dtype="float64"),
        sigma_xy=np.asarray(sigma_xy, dtype="float64"),
        sigma_v=np.asarray(sigma_v, dtype="float64"),
        scalar=scalar,
        acf=np.asarray(like.acf, dtype="float64"),
    )


# -- map-ingest helpers (issue #28) ----------------------------------------------------
def _map_values(map_or_array, unit: u.UnitBase) -> np.ndarray:
    """The bare 2D values of *map_or_array* in *unit* (a :class:`Map`, ``Quantity``, or array).

    Accepts a first-class :class:`~astrolyze.core.Map` (``.data``), a bare ``Quantity``, or a
    plain ndarray (assumed already in *unit*). Converting to *unit* keeps the no-silent-physics
    contract: a survey map in mismatched units is converted, not silently misread."""
    data = getattr(map_or_array, "data", map_or_array)
    if isinstance(data, u.Quantity):
        return np.asarray(data.to_value(unit), dtype="float64")
    return np.asarray(data, dtype="float64")


def _representation_from_sigma_map(cube, sigma_xy: np.ndarray):
    """A separable representation from a 2D σ map ``sigma_xy`` (the survey-map ingest path).

    A published σ map is a 2D field: σ is constant along the spectral axis, so the model is
    separable by construction — σ_v is flat at the scalar level and the reconstructed σ-cube is
    ``σ_xy`` broadcast over every channel. With no finite pixel the map carries no σ and the
    quality is :data:`NoiseQuality.UNRELIABLE` (ADR-0003)."""
    sigma_xy = np.asarray(sigma_xy, dtype="float64")
    nz = cube.shape[0]
    if not np.isfinite(sigma_xy).any():
        return _unreliable_representation(cube.shape), NoiseQuality.UNRELIABLE

    scalar = float(np.nanmean(sigma_xy))
    sigma_v = np.full(
        nz, scalar
    )  # flat spectral profile: one σ per pixel, all channels
    # ACF of pure (white, channel-independent) noise: a unit spike at zero lag.
    acf = np.zeros(nz)
    acf[0] = 1.0
    rep = SeparableNoise(
        sigma_xy=sigma_xy,
        sigma_v=sigma_v,
        scalar=scalar,
        acf=acf,
    )
    return rep, NoiseQuality.MEASURED


# -- built-in estimator registration ---------------------------------------------------
# Registering rather than hard-coding an if/elif keeps the routing open: a caller adds an
# estimator via register_estimator() and it is selectable by estimate_noise(method=...) with no
# core edit (the pluggability seam, issue #28).
register_estimator(
    DEFAULT_METHOD, _estimate_mad_std
)  # "mad_std" — robust default (#27)
register_estimator("rms", _estimate_rms)  # signal-free-channel RMS
register_estimator("mad", _estimate_mad_sigma_clip)  # robust σ-clipped MAD
# ======================================================================================
# Noise propagation through matching + correlated-noise synthesis (issue #32)
# ======================================================================================
# These are astrolyze's value-add: the propagation LAWS (correlation-area scaling, M_eff from
# the ACF) and a synthesis utility. numpy/astropy/FFT do the maths; the laws are the physics.
# Everything below is ADDITIVE — it does not touch #27's estimator or representation code.
#
# Why analytic, not re-estimated, on the hot path (ADR-0003/0004): matching is the fast path,
# and a real science cube is *full of source* — re-measuring σ from a matched cube would fold
# bright extended emission into the estimate and overstate the noise. The matching op therefore
# propagates the stored, signal-free σ through the known correlation structure; re-estimation is
# kept as the validation oracle / non-stationary fallback, never called inside the op.


def _spatial_rms_factor(beam_in, beam_out) -> float:
    """The per-pixel RMS scaling for convolving from *beam_in* to a larger *beam_out*.

    The correlation area is the **beam** solid angle (instrument-limited), NOT the pixel area:
    noise correlated over a beam averages down by the ratio of correlation footprints when
    smoothed to a larger beam. The RMS scales by ``√(Ω_corr,in / Ω_out)`` — a larger output beam
    averages over more correlated area, so σ drops (the factor is < 1). This is the spatial half
    of issue #32's propagation law; ``beam.sr`` is radio_beam's beam solid angle (Ω)."""
    omega_in = float(beam_in.sr.to_value(u.sr))
    omega_out = float(beam_out.sr.to_value(u.sr))
    if omega_out <= 0.0 or omega_in <= 0.0:
        return float("nan")
    return float(np.sqrt(omega_in / omega_out))


def _m_eff(acf: np.ndarray, m: int) -> float:
    """The *effective* number of independent samples in an ``m``-channel average, ``M²/Σρ``.

    For white noise averaging M channels gives σ/√M; for **correlated** noise fewer samples are
    independent, so the reduction is weaker. The autocorrelation ``acf`` (normalised to 1 at lag
    0) sets the effective count via the standard variance-of-the-mean result:

        Var(mean of M) = σ²/M² · Σ_{i,j} ρ(|i-j|)
                       = σ²/M · [ ρ(0) + 2 Σ_{k=1}^{M-1} (1 - k/M) ρ(k) ]

    so the variance reduction is ``M / Σ``, i.e. ``M_eff = M² / Σ`` with ``Σ`` the (triangular-
    weighted) two-sided lag sum out to lag M-1. White noise (ρ = δ) -> Σ = 1 -> M_eff = M. A
    broad ACF -> Σ > 1 -> M_eff < M (less averaging-down). This is the spectral half of the law:
    σ/√M_eff, NOT the naive σ/√M (ADR-0003 — do not silently assume whiteness)."""
    acf = np.asarray(acf, dtype="float64")
    m = int(m)
    if m <= 1 or acf.size == 0 or not np.isfinite(acf[0]) or acf[0] == 0.0:
        return float(m)
    rho = acf / acf[0]
    lags = np.arange(1, min(m, rho.size))
    triangular = rho[lags] * (1.0 - lags / m)
    lag_sum = float(rho[0] + 2.0 * np.nansum(triangular))
    if lag_sum <= 0.0:
        return float(m)
    return float((m * m) / (m * lag_sum))


def _scaled_representation(rep, factor: float):
    """Return a copy of *rep* with every σ array scaled by *factor* (the ACF is unitless: kept).

    Propagation is a uniform rescale of the σ field — the spatial/spectral shape and the
    correlation structure are unchanged, only the level moves — so we keep the same branch
    (separable vs full) and rescale the stored arrays. Reusing the #27 representation classes
    keeps every downstream product (sigma_cube/map/spectrum/scalar) consistent for free."""
    from dataclasses import replace

    if isinstance(rep, SeparableNoise):
        return replace(
            rep,
            sigma_xy=np.asarray(rep.sigma_xy) * factor,
            sigma_v=np.asarray(rep.sigma_v) * factor,
            scalar=float(rep.scalar) * factor,
        )
    return replace(
        rep,
        sigma_field=np.asarray(rep.sigma_field) * factor,
        sigma_xy=np.asarray(rep.sigma_xy) * factor,
        sigma_v=np.asarray(rep.sigma_v) * factor,
        scalar=float(rep.scalar) * factor,
    )


# -- spectral rebin (issue #87): rebin the noise spectral axis onto the binned data grid -----
# spectral_factor is variance-only (smoothing keeps the channel grid); spectral_bin instead
# coarsens the spectral axis to match a k-channel averaging of the DATA so the σ companion stays
# paired voxel-for-voxel with the binned cube. The block structure mirrors spectral-cube's
# downsample_axis(k, axis=0) (truncate=False): contiguous k-channel blocks with a trailing
# PARTIAL block, so nz_out = ceil(N/k). The σ-level reduction follows the same σ/√M_eff law the
# spectral_factor path uses (the variance-of-the-mean result through the stored ACF); per-channel
# σ is rebinned in quadrature through the same normalised ACF ρ.


def _block_bounds(n: int, k: int):
    """The [start, stop) channel index pairs of an ``k``-channel binning of ``n`` channels.

    Contiguous blocks of ``k`` with a trailing partial block (truncate=False) — matches
    spectral-cube's ``downsample_axis`` exactly, including the remainder block for non-divisors."""
    starts = np.arange(0, n, k)
    return [(int(s), int(min(s + k, n))) for s in starts]


def _normalised_acf(acf: np.ndarray) -> np.ndarray:
    """The stored ACF normalised to ρ(0)=1; empty/degenerate -> a white spike ``[1.0]``."""
    acf = np.asarray(acf, dtype="float64")
    if acf.size == 0 or not np.isfinite(acf[0]) or acf[0] == 0.0:
        return np.array([1.0])
    return acf / acf[0]


def _rho_at(rho: np.ndarray, lag: int) -> float:
    """ρ at integer ``lag``; symmetric (ρ(-l)=ρ(l)), zero beyond the stored length."""
    lag = abs(int(lag))
    if lag < rho.size:
        v = float(rho[lag])
        return v if np.isfinite(v) else 0.0
    return 0.0


def _block_quadrature(sigma: np.ndarray, rho: np.ndarray, lo: int, hi: int) -> float:
    """The block-averaged σ over channels ``[lo:hi)``: ``sqrt((1/m²) ΣΣ σ_p σ_q ρ(|p-q|))``.

    The variance-of-the-block-mean through the noise correlation ρ. White noise (ρ=δ) collapses
    to rms/√m -> σ/√k for a full block of equal σ; correlation raises it (fewer independent
    samples), so this is the exact per-channel analogue of σ/√M_eff (issue #87)."""
    idx = np.arange(lo, hi)
    m = idx.size
    if m == 0:
        return 0.0
    s = sigma[idx]
    # Lag matrix |p-q| over the block; ρ(|p-q|) via the symmetric, zero-padded lookup.
    lag = np.abs(idx[:, None] - idx[None, :])
    rho_mat = np.array([[_rho_at(rho, int(lg)) for lg in row] for row in lag])
    var = float((s[:, None] * s[None, :] * rho_mat).sum()) / (m * m)
    return float(np.sqrt(var)) if var > 0.0 else 0.0


def _rebin_sigma_v(sigma_v: np.ndarray, rho: np.ndarray, k: int) -> np.ndarray:
    """Rebin the 1D per-channel σ (length N -> ceil(N/k)) by the block quadrature above."""
    sigma_v = np.asarray(sigma_v, dtype="float64")
    bounds = _block_bounds(sigma_v.size, k)
    return np.array([_block_quadrature(sigma_v, rho, lo, hi) for lo, hi in bounds])


def _rebin_sigma_field(sigma_field: np.ndarray, rho: np.ndarray, k: int) -> np.ndarray:
    """Rebin the full 3D σ field along the spectral axis (nz -> ceil(nz/k)) per spatial pixel."""
    sigma_field = np.asarray(sigma_field, dtype="float64")
    nz, ny, nx = sigma_field.shape
    bounds = _block_bounds(nz, k)
    out = np.empty((len(bounds), ny, nx))
    for j, (lo, hi) in enumerate(bounds):
        m = hi - lo
        block = sigma_field[lo:hi]  # (m, ny, nx)
        idx = np.arange(lo, hi)
        lag = np.abs(idx[:, None] - idx[None, :])
        rho_mat = np.array([[_rho_at(rho, int(lg)) for lg in row] for row in lag])
        # ΣΣ σ_p σ_q ρ(|p-q|) per pixel via the (m,m) ρ contraction over the two block axes.
        var = np.einsum("pyx,qyx,pq->yx", block, block, rho_mat) / (m * m)
        out[j] = np.sqrt(np.where(var > 0.0, var, 0.0))
    return out


def _coarse_acf(acf: np.ndarray, k: int, nz_out: int) -> np.ndarray:
    """The ACF of the binned noise (length ``nz_out``) from the fine ACF and bin factor ``k``.

    A k-channel block average is a triangular (Bartlett) filter on the channel grid; the binned
    autocovariance at coarse lag L is the fine autocovariance summed against the overlap weight
    ``W(δ)=k-|δ|`` over the fine lags spanned by that coarse lag. Normalising by the L=0 sum gives
    ρ_out(0)=1. White noise -> a spike at 0. Degenerate normalisation -> white fallback (#87)."""
    rho = _normalised_acf(acf)
    deltas = np.arange(-(k - 1), k)
    weights = k - np.abs(deltas)
    ndenom = float(
        sum(_rho_at(rho, int(d)) * float(w) for d, w in zip(deltas, weights))
    )
    out = np.zeros(int(nz_out))
    if not np.isfinite(ndenom) or ndenom <= 0.0:
        out[0] = 1.0
        return out
    for L in range(int(nz_out)):
        nnum = sum(
            float(w) * _rho_at(rho, int(L * k + d)) for d, w in zip(deltas, weights)
        )
        out[L] = nnum / ndenom
    return out


def _rebinned_representation(rep, k: int, r: float):
    """Rebin *rep* onto a ``k``-channel-binned spectral grid (issue #87).

    σ_v is rebinned in quadrature through the normalised ACF; σ_xy and the scalar take the same
    global √M_eff level reduction ``r`` the variance-only path applies (so the reconstruct ratio
    σ_xy/σ_0 — hence sigma_cube — is preserved); the ACF is re-derived on the coarse grid. For the
    full (non-separable) case the dense σ field is block-rebinned per pixel as well."""
    from dataclasses import replace

    rho = _normalised_acf(rep.acf)
    sigma_v_out = _rebin_sigma_v(rep.sigma_v, rho, k)
    acf_out = _coarse_acf(rep.acf, k, sigma_v_out.size)
    if isinstance(rep, SeparableNoise):
        return replace(
            rep,
            sigma_xy=np.asarray(rep.sigma_xy) * r,
            sigma_v=sigma_v_out,
            scalar=float(rep.scalar) * r,
            acf=acf_out,
        )
    return replace(
        rep,
        sigma_field=_rebin_sigma_field(rep.sigma_field, rho, k),
        sigma_xy=np.asarray(rep.sigma_xy) * r,
        sigma_v=sigma_v_out,
        scalar=float(rep.scalar) * r,
        acf=acf_out,
    )


def propagate(
    model: NoiseModel,
    *,
    beam_out=None,
    spectral_factor: int | None = None,
    spectral_bin: int | None = None,
    per_channel_beam: bool = False,
) -> NoiseModel:
    """Propagate *model* through a matching op analytically; return a new :class:`NoiseModel`.

    The fast path (ADR-0003/0004): given the matching *geometry* — a larger output beam and/or a
    spectral averaging factor — scale the stored, signal-free σ through the known correlation
    structure rather than re-measuring it:

    - **spatial** (``beam_out``) — σ × ``√(Ω_corr,in / Ω_out)`` with the correlation area the
      **beam** solid angle (not the pixel area); the new beam is recorded as context;
    - **spectral, variance-only** (``spectral_factor`` M) — σ ÷ ``√M_eff`` with ``M_eff = M²/Σρ``
      from the stored ACF (the naive σ/√M only when the noise is white). The channel grid is
      UNCHANGED — this is the smoothing path (:meth:`Cube.spectral_smooth_to`);
    - **spectral, rebinning** (``spectral_bin`` k) — REBINS the noise spectral axis onto a
      k-channel binning of the data (:meth:`Cube.spectral_bin`, issue #87). σ_v is block-averaged
      in quadrature through the ACF and the axis coarsens to ceil(N/k) channels; σ_xy / σ_0 take
      the same √M_eff level reduction; the ACF is re-derived on the coarse grid. Distinct from
      ``spectral_factor`` (which never changes the channel count) — the two must not be conflated.

    The result carries :data:`NoiseQuality.PROPAGATED`. It degrades to
    :data:`NoiseQuality.APPROXIMATE` (no silent physics) when the single-kernel / stationary
    assumption breaks: a *per-channel beam* (``per_channel_beam=True``) or a non-stationary
    source model (a full-3D, non-separable input — propagated as a uniform rescale, which is only
    approximate). Re-estimation is **never** invoked here — that is :func:`reestimate`."""
    if beam_out is None and spectral_factor is None and spectral_bin is None:
        raise ValueError(
            "propagate() needs a geometry to propagate through: pass beam_out (spatial), "
            "spectral_factor (spectral, variance-only), and/or spectral_bin (spectral rebin). "
            "With none there is nothing to propagate (astrolyze does not silently return the "
            "input model unchanged, ADR-0003)."
        )
    if spectral_bin is not None and spectral_factor is not None:
        raise ValueError(
            "propagate() got both spectral_bin and spectral_factor: they are different spectral "
            "ops — spectral_factor is a variance-only rescale (smoothing keeps the channel grid), "
            "spectral_bin rebins the spectral axis (#87). Pass one, not both."
        )
    if spectral_bin is not None and (
        not isinstance(spectral_bin, (int, np.integer)) or int(spectral_bin) < 2
    ):
        raise ValueError(
            f"propagate(spectral_bin=...) must be an integer >= 2, got {spectral_bin!r} "
            "(factor 1 is a no-op, < 1 invalid — mirrors Cube.spectral_bin's guard)."
        )

    rep = model._rep
    approximate = per_channel_beam or not model.is_separable

    if beam_out is not None:
        rep = _scaled_representation(rep, _spatial_rms_factor(model.beam, beam_out))
    if spectral_factor is not None:
        m_eff = _m_eff(np.asarray(model._rep.acf, dtype="float64"), spectral_factor)
        rep = _scaled_representation(rep, 1.0 / float(np.sqrt(m_eff)))
    if spectral_bin is not None:
        k = int(spectral_bin)
        # Same global √M_eff level reduction the variance-only path applies, so the scalar lands
        # at σ_0/√M_eff and the reconstruct ratio σ_xy/σ_0 (hence sigma_cube) is preserved (#87).
        r = 1.0 / float(np.sqrt(_m_eff(np.asarray(rep.acf, dtype="float64"), k)))
        rep = _rebinned_representation(rep, k, r)

    quality = NoiseQuality.APPROXIMATE if approximate else NoiseQuality.PROPAGATED
    new_cube = model._cube
    new_beam = beam_out if beam_out is not None else model.beam
    if spectral_bin is not None:
        # Downsample the parent cube through the SAME spectral-cube call the data binning uses, so
        # the returned model's coarse spectral axis / WCS / channel count match the binned data
        # EXACTLY (trailing partial block included) — the model's sigma_spectrum/sigma_cube build
        # off this cube's axis, so it must already be the coarse grid (#87). Lazy import keeps the
        # core import order clean (cube imports noise).
        from .cube import Cube

        binned_sc = model._cube._masked_sc().downsample_axis(int(spectral_bin), axis=0)
        new_cube = Cube(binned_sc, model._cube._metadata_with_beam(new_beam))
    elif beam_out is not None:
        # Carry the new beam onto the model's parent cube so the propagated products report the
        # post-match resolution. Lazy import keeps core import order clean (cube imports noise).
        from .cube import Cube

        new_cube = Cube(new_cube._sc, new_cube._metadata_with_beam(beam_out))
    _emit(
        "noise.propagate",
        params={
            "beam_out": str(beam_out) if beam_out is not None else None,
            "spectral_factor": spectral_factor,
            "spectral_bin": spectral_bin,
            "quality": quality.value,
        },
    )
    return NoiseModel(
        new_cube, rep, method=model.method, quality=quality, version=model.version
    )


def reestimate(cube) -> NoiseModel:
    """The validation oracle / non-stationary fallback: re-measure σ from *cube* via ``mad_std``.

    This is exactly the #27 :func:`estimate` (robust ``mad_std`` on signal-free voxels) surfaced
    under a propagation-validation name. It is the **oracle** that an analytic propagation can be
    checked against on a synthetic stationary-noise cube, and the fallback for the genuinely
    non-stationary case. It is **never** called inside a matching op — on a real science cube
    bright extended emission fills the field of view and would corrupt the estimate (ADR-0003),
    so the hot path stays analytic (:func:`propagate`)."""
    return estimate(cube)


def _beam_sigma_pixels(beam, pixel_scale: u.Quantity) -> tuple[float, float, float]:
    """The beam's Gaussian σ (major, minor in pixels) and position angle (rad), from *beam*.

    Synthesis convolves white noise by the beam, which is a Gaussian of these widths on the pixel
    grid; ``beam.major``/``minor`` are FWHM, converted to σ in pixels via the pixel scale."""
    pix = pixel_scale.to(u.deg)
    major_fwhm = (beam.major.to(u.deg) / pix).to_value(u.dimensionless_unscaled)
    minor_fwhm = (beam.minor.to(u.deg) / pix).to_value(u.dimensionless_unscaled)
    fwhm_to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    pa = float(beam.pa.to_value(u.rad))
    return major_fwhm * fwhm_to_sigma, minor_fwhm * fwhm_to_sigma, pa


def synthesize_correlated_noise(
    shape,
    *,
    beam,
    pixel_scale: u.Quantity,
    spectral_acf: np.ndarray | None = None,
    sigma: float = 1.0,
    rng=None,
) -> np.ndarray:
    """Synthesize realistic correlated noise at a target resolution; return a bare ``ndarray``.

    Draw white Gaussian noise of *shape* ``(nz, ny, nx)`` and impose the target correlation:

    - **spatial** — convolve each channel by the *beam* (an elliptical Gaussian of the beam's
      FWHM on the pixel grid set by *pixel_scale*), so the spatial autocorrelation matches the
      beam (instrument-limited correlation);
    - **spectral** — if a target ``spectral_acf`` is given, shape the spectral axis by it via the
      FFT: the noise power spectrum is the FFT of the ACF (Wiener-Khinchin), so multiplying the
      white spectrum's amplitude by ``√|FFT(acf)|`` yields noise with that autocorrelation.

    The realization is renormalised to the target per-voxel ``sigma`` (a robust scale), so it is
    correlated noise *at a stated level*. Stays thin: numpy/astropy convolution + numpy FFT do
    the maths; the recipe (beam = spatial kernel, ACF = spectral shaper) is the value-add. The
    inverse of :func:`Cube.estimate_noise`'s analysis — useful for testing the propagation laws
    and for injecting realistic noise into mock cubes."""
    from astropy.convolution import Gaussian2DKernel, convolve_fft
    from astropy.stats import mad_std

    rng = rng if rng is not None else np.random.default_rng()
    nz, ny, nx = shape
    white = rng.normal(0.0, 1.0, size=(nz, ny, nx))

    sigma_major, sigma_minor, pa = _beam_sigma_pixels(beam, pixel_scale)
    # Gaussian2DKernel's theta is measured from the +x axis; radio_beam PA is from +y (North),
    # so the kernel is rotated by pa + 90deg. The exact angle only matters for an elliptical
    # beam; the recovered correlation width is dominated by the (major, minor) widths.
    kernel = Gaussian2DKernel(
        x_stddev=sigma_minor,
        y_stddev=sigma_major,
        theta=pa,
    )
    spatial = np.empty_like(white)
    for k in range(nz):
        spatial[k] = convolve_fft(
            white[k], kernel, normalize_kernel=True, boundary="wrap"
        )

    if spectral_acf is not None:
        spatial = _shape_spectral_acf(
            spatial, np.asarray(spectral_acf, dtype="float64")
        )

    current = float(mad_std(spatial))
    if current > 0.0 and np.isfinite(current):
        spatial = spatial * (sigma / current)
    return spatial


def _shape_spectral_acf(field: np.ndarray, target_acf: np.ndarray) -> np.ndarray:
    """Shape *field*'s spectral axis to have autocorrelation *target_acf* via the FFT.

    Wiener-Khinchin: the power spectral density is the FFT of the autocorrelation. We build a
    symmetric (two-sided) ACF of the spectral length, take its real FFT as the target PSD, and
    multiply each pixel's white spectrum by ``√PSD`` — the result has the requested ACF. Negative
    numerical PSD values (from a non-positive-definite truncated ACF) are clipped to 0."""
    nz = field.shape[0]
    acf = np.zeros(nz)
    n = min(nz, target_acf.size)
    acf[:n] = target_acf[:n]
    if acf[0] != 0.0:
        acf = acf / acf[0]
    # Symmetric two-sided ACF -> its FFT is the (real, non-negative) power spectrum.
    symmetric = np.concatenate([acf, acf[1:][::-1]])
    psd = np.fft.rfft(symmetric).real
    psd = np.clip(psd, 0.0, None)
    amplitude = np.sqrt(psd)
    flat = field.reshape(nz, -1)
    out = np.empty_like(flat)
    nsym = symmetric.size
    for j in range(flat.shape[1]):
        spec = np.fft.rfft(flat[:, j], n=nsym)
        shaped = np.fft.irfft(spec * amplitude, n=nsym)
        out[:, j] = shaped[:nz]
    return out.reshape(field.shape)


def _measure_spectral_acf(field: np.ndarray) -> np.ndarray:
    """The mean per-pixel spectral autocorrelation of *field* (normalised to 1 at lag 0).

    A test/validation helper mirroring the #27 estimator's ACF: it is how a synthesis is checked
    against its target spectral correlation."""
    return _spectral_acf(np.asarray(field, dtype="float64"))


def _spatial_acf_fwhm_pixels(field: np.ndarray) -> float:
    """The FWHM (pixels) of *field*'s mean spatial autocorrelation — a synthesis check helper.

    Averages the per-channel 2D autocorrelation (via the FFT power spectrum), takes a central
    radial profile, and returns the full width at half maximum in pixels. Used only by tests to
    confirm a synthesized realization's spatial correlation tracks the target beam."""
    nz, ny, nx = field.shape
    acc = np.zeros((ny, nx))
    for k in range(nz):
        plane = field[k] - np.mean(field[k])
        power = np.abs(np.fft.fft2(plane)) ** 2
        acc += np.fft.ifft2(power).real
    acf2d = np.fft.fftshift(acc / nz)
    peak = acf2d.max()
    if peak <= 0.0:
        return float("nan")
    acf2d = acf2d / peak
    cy, cx = ny // 2, nx // 2
    row = acf2d[cy, cx:]  # radial profile along +x from the centre
    half = 0.5
    below = np.where(row < half)[0]
    if below.size == 0:
        return float("nan")
    # linear interpolation to the half-maximum crossing -> HWHM, doubled for FWHM.
    i = below[0]
    if i == 0:
        return 0.0
    x0, x1 = i - 1, i
    y0, y1 = row[x0], row[x1]
    hwhm = x0 + (y0 - half) / (y0 - y1) if y0 != y1 else float(x0)
    return 2.0 * hwhm


__all__ = [
    "NoiseModel",
    "NoiseQuality",
    "SeparableNoise",
    "FullNoise",
    "estimate",
    "register_estimator",
    "available_estimators",
    "propagate",
    "reestimate",
    "synthesize_correlated_noise",
    "NOISE_SCHEMA_VERSION",
    "DEFAULT_METHOD",
]
