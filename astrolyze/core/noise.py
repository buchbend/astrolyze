"""The context-carrying ``NoiseModel`` companion + its one default estimator (issue #27).

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

import numpy as np
import astropy.units as u
from astropy.stats import mad_std

from astrolyze.io import Metadata

from ._base import _emit

#: Schema version of the on-disk companion-group layout (provenance, bumped on layout change).
NOISE_SCHEMA_VERSION = 1

#: The default (and, for issue #27, only) estimator. The pluggable suite is issue #28.
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

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        kind = "separable" if self.is_separable else "full"
        return f"<NoiseModel {obj} {kind} [{self.unit}] {self.quality.value}>"


# -- the default estimator -------------------------------------------------------------
def estimate(cube, *, method: str = DEFAULT_METHOD) -> NoiseModel:
    """Estimate a :class:`NoiseModel` for *cube* with the given *method* (issue #27).

    One estimator is built — ``"mad_std"`` (the robust median-absolute-deviation σ); the
    pluggable suite is issue #28. An unknown *method* is refused rather than silently
    substituted (ADR-0003)."""
    if method != DEFAULT_METHOD:
        raise ValueError(
            f"unknown noise estimator method {method!r}: issue #27 ships only "
            f"{DEFAULT_METHOD!r} (the pluggable suite is issue #28); astrolyze never silently "
            "substitutes an estimator (ADR-0003)"
        )
    data = np.asarray(cube._data_quantity.value, dtype="float64")
    rep, quality = _estimate_mad_std(data)
    _emit("estimate_noise", params={"method": method, "quality": quality.value})
    return NoiseModel(cube, rep, method=method, quality=quality)


def _estimate_mad_std(data: np.ndarray):
    """The robust σ estimator: per-voxel scatter via :func:`astropy.stats.mad_std`.

    Returns ``(representation, quality)``. With no usable (finite) data the quality is
    :data:`NoiseQuality.UNRELIABLE` and σ is left ``NaN`` — never fabricated (ADR-0003)."""
    finite = np.isfinite(data)
    if not finite.any():
        # No usable signal-free data: refuse to invent a σ (ADR-0003).
        shape_xy = data.shape[1:]
        nz = data.shape[0]
        rep = SeparableNoise(
            sigma_xy=np.full(shape_xy, np.nan),
            sigma_v=np.full(nz, np.nan),
            scalar=np.nan,
            acf=np.array([np.nan]),
        )
        return rep, NoiseQuality.UNRELIABLE

    # Robust scatter on each axis-collapse (NaN-ignoring): spatial σ per pixel (over v),
    # spectral σ per channel (over the plane), and a single scalar over the whole cube.
    sigma_xy = mad_std(data, axis=0, ignore_nan=True)
    sigma_v = mad_std(data.reshape(data.shape[0], -1), axis=1, ignore_nan=True)
    scalar = float(mad_std(data, ignore_nan=True))
    acf = _spectral_acf(data)

    rep = _build_representation(data, sigma_xy, sigma_v, scalar, acf)
    return rep, NoiseQuality.MEASURED


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


__all__ = [
    "NoiseModel",
    "NoiseQuality",
    "SeparableNoise",
    "FullNoise",
    "estimate",
    "NOISE_SCHEMA_VERSION",
    "DEFAULT_METHOD",
]
