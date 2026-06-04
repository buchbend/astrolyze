"""Per-axis physical coordinates + a validity descriptor, surfaced off the core objects.

This is a **read-only surfacing** of what the WCS / spectral axis astrolyze already parsed
(no header reparse, no new geometry — ADR-0004): it turns the upstream-parsed axes into
first-class, context-carrying value objects (never a bare array). Two value objects:

- :class:`AxisCoordinates` — the per-axis physical coordinate arrays as ``Quantity`` (each
  carries its own unit): the **authoritative** absolute frequency, the per-line Δv, the sky
  longitude/latitude maps, and the celestial pixel scale. ``Cube`` exposes all of them;
  ``Map`` exposes only the sky/pixel subset; ``Spectrum`` only the spectral subset (the
  absent fields are simply ``None``).
- :class:`Validity` — *where the data is real*: the values with blanked / edge / outside-
  coverage voxels as ``NaN`` plus a boolean finite-data ``mask``. The mask is derived from
  the data itself, so it travels with the data through a subcube slice for free
  (mask-of-a-slice == slice-of-the-mask).

No-silent-physics (ADR-0003): the absolute frequency is *authoritative*, read from the
spectral axis by the object's stated velocity convention + rest frequency. When the spectral
axis is a velocity and that context is absent, deriving an absolute frequency would require
guessing the convention — so we **raise** :class:`~astrolyze.units.MissingContextError`
rather than assume one. The per-line Δv is the matching inverse, derived from that absolute
frequency relative to the metadata's ``rest_frequency`` under its ``velocity_convention``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import astropy.units as u
from astropy.coordinates import SpectralCoord
from astropy.wcs.utils import proj_plane_pixel_scales

from astrolyze.units import MissingContextError, doppler


@dataclass(frozen=True)
class AxisCoordinates:
    """Per-axis physical coordinate arrays, each a ``Quantity`` carrying its own unit.

    A field is ``None`` when the axis is absent for the wrapper that produced it (a ``Map``
    has no spectral axis; a ``Spectrum`` has no sky axis)."""

    #: Absolute frequency per channel (authoritative), or ``None`` for a non-spectral object.
    frequency: u.Quantity | None = None
    #: Per-line velocity offset relative to the rest frequency, or ``None``.
    delta_v: u.Quantity | None = None
    #: Sky longitude map (e.g. RA), or ``None`` for a non-spatial object.
    longitude: u.Quantity | None = None
    #: Sky latitude map (e.g. Dec), or ``None`` for a non-spatial object.
    latitude: u.Quantity | None = None
    #: Celestial pixel scale (one value per spatial axis), or ``None``.
    pixel_scale: u.Quantity | None = None


#: Schema version of the on-disk validity companion-group layout (provenance, bumped on a
#: layout change — the sibling of :data:`~astrolyze.core.noise.NOISE_SCHEMA_VERSION`).
VALIDITY_SCHEMA_VERSION = 1

#: How the mask was derived: the finite-data locations of the array (blanked/edge/coverage
#: voxels are ``NaN``). The provenance ``method`` of a companion, like the noise estimator name.
VALIDITY_METHOD = "finite_mask"


@dataclass(frozen=True)
class Validity:
    """Where the data is real: blanked voxels as ``NaN`` + a boolean finite-data mask.

    ``mask`` is ``True`` exactly where ``data`` is finite (blanked / edge / outside-coverage
    voxels are ``False`` and appear as ``NaN`` in ``data``). Because both are read from the
    same array, slicing the wrapper slices both consistently — the mask travels with the data
    (mask-of-a-slice == slice-of-the-mask).

    Like the :class:`~astrolyze.core.noise.NoiseModel`, the descriptor persists as a *companion
    group* inside the cube's #23 Zarr store (:meth:`to_zarr_companion` / :meth:`from_zarr_companion`):
    the blanked voxels already sit as ``NaN`` in the data array, so the companion adds only the
    explicit ``uint8`` mask the loader's lean core reads (ADR-0008)."""

    #: The values with blanked voxels exposed as ``NaN`` (carries the data's unit).
    data: u.Quantity
    #: Boolean finite-data mask, the shape of the data (``True`` where the data is real).
    mask: np.ndarray

    def to_zarr_companion(self, store):
        """Write this descriptor as a ``validity`` companion group in the cube's Zarr *store*.

        Mirrors :meth:`NoiseModel.to_zarr_companion`: the boolean mask lands as a ``uint8``
        array in a ``validity`` subgroup carrying its derivation method + schema version as
        provenance. Delegates the byte I/O to :mod:`astrolyze.io.zarr_backend`."""
        from astrolyze.io.zarr_backend import _save_validity_companion

        return _save_validity_companion(self, store)

    @classmethod
    def from_zarr_companion(cls, store) -> "Validity":
        """Reconstruct a :class:`Validity` from the ``validity`` companion group in *store*.

        The ``mask`` is read from the group (back to the cube's array order); the NaN-exposing
        ``data`` is the cube's intensity reloaded from the store, so ``mask == isfinite(data)``
        holds on the round-trip. The inverse of :meth:`to_zarr_companion`."""
        from astrolyze.io.zarr_backend import _load_validity_companion

        return _load_validity_companion(store)


def validity_of(data: u.Quantity) -> Validity:
    """Build the :class:`Validity` descriptor for a NaN-filled ``Quantity``.

    The mask is the finite locations of the data, so a consumer knows where the data is real
    without re-deriving it; both come from the one array, so they slice together."""
    quantity = u.Quantity(data)
    return Validity(data=quantity, mask=np.isfinite(quantity.value))


# -- spectral axis --------------------------------------------------------------------
def absolute_frequency(spectral_axis, *, rest_frequency, convention) -> u.Quantity:
    """The AUTHORITATIVE absolute frequency for an already-parsed spectral axis.

    If the spectral axis is already frequency-like it is returned in Hz directly (no context
    needed). If it is a velocity, mapping it to an absolute frequency needs the Doppler
    convention + rest frequency the object carries; absent, we never guess (ADR-0003) — the
    Doppler builder raises :class:`~astrolyze.units.MissingContextError`."""
    axis = u.Quantity(spectral_axis)
    if axis.unit.is_equivalent(u.Hz, equivalencies=u.spectral()):
        return axis.to(u.Hz, equivalencies=u.spectral())
    equiv = doppler(convention, rest_frequency)
    return axis.to(u.Hz, equivalencies=equiv + u.spectral())


def delta_v(spectral_axis, *, rest_frequency, convention) -> u.Quantity:
    """Per-line velocity offset of each channel from ``rest_frequency`` under ``convention``.

    Derived from the authoritative absolute frequency (see :func:`absolute_frequency`) by the
    object's stated convention — not assumed. Raises if the context is absent (ADR-0003)."""
    nu = absolute_frequency(
        spectral_axis, rest_frequency=rest_frequency, convention=convention
    )
    equiv = doppler(convention, rest_frequency)
    return nu.to(u.km / u.s, equivalencies=equiv + u.spectral())


def delta_v_relative_to(
    absolute_frequency_axis, *, line_rest_frequency, convention
) -> u.Quantity:
    """Per-channel Δv of an already-absolute frequency axis relative to a *chosen* line.

    The companion of :func:`delta_v` for broadband work: the authoritative absolute frequency
    is fixed (derived once from the primary line), and this measures each channel's velocity
    offset from *any* line's rest frequency under the same ``convention`` — so a second species
    in the band gets its own Δv. Raises if the convention is absent (ADR-0003)."""
    equiv = doppler(convention, line_rest_frequency)
    return u.Quantity(absolute_frequency_axis).to(
        u.km / u.s, equivalencies=equiv + u.spectral()
    )


# -- spectral-frame transform ---------------------------------------------------------
# FITS SPECSYS codes -> the astropy frame an observer is held stationary in (the spectral
# reference frame). Barycentric/heliocentric are ICRS-stationary to astropy; topocentric is
# the observer's own frame (a no-op shift). Only the frames astrolyze claims to transform
# between are listed; an unknown target raises rather than silently picking one (ADR-0003).
_SPECSYS_TO_ASTROPY_FRAME = {
    "LSRK": "lsrk",
    "LSRD": "lsrd",
    "BARYCENT": "icrs",
    "HELIOCEN": "hcrs",
    "TOPOCENT": "itrs",
    "GALACTOC": "galactocentric",
}


def to_spectral_frame(
    absolute_frequency_axis,
    target_frame,
    *,
    source_frame,
    location,
    obstime,
    target,
) -> u.Quantity:
    """Transform an absolute-frequency axis to a different spectral reference frame.

    Thin door over astropy :class:`~astropy.coordinates.SpectralCoord`: it builds the spectral
    coordinate with the observer (the telescope ``location`` at ``obstime``) and the ``target``
    sky position, then asks for the axis as seen by an observer stationary in ``target_frame``.

    Context-or-raise (ADR-0003): the observer ``location`` + ``obstime`` and the ``target``
    are all required — without them a frame shift cannot be computed, so we raise rather than
    silently leaving the frame unchanged. ``source_frame`` (the dataset's current ``SPECSYS``)
    must be known too: transforming *from* an unknown frame would be guessing, so an absent
    ``source_frame`` raises naming ``spectral_frame``."""
    if source_frame is None:
        raise MissingContextError(
            "spectral_frame is unknown (SPECSYS absent): a frame transform needs the source "
            "frame to interpret the axis; astrolyze never guesses it (ADR-0003)",
        )
    missing = [
        name
        for name, value in (
            ("location", location),
            ("obstime", obstime),
            ("target", target),
        )
        if value is None
    ]
    if missing:
        raise MissingContextError(
            "a spectral-frame transform requires observer location + obstime + target "
            "(SpectralCoord cannot shift the frame without them); missing: "
            + ", ".join(missing)
        )
    astropy_frame = _astropy_frame_for(target_frame)
    observer = location.get_itrs(obstime=obstime)
    coord = SpectralCoord(
        u.Quantity(absolute_frequency_axis).to(u.Hz, equivalencies=u.spectral()),
        observer=observer,
        target=target,
    )
    return coord.with_observer_stationary_relative_to(astropy_frame).to(u.Hz)


def _astropy_frame_for(specsys) -> str:
    """Map a FITS ``SPECSYS`` code to the astropy frame name; raise on an unknown code."""
    code = str(specsys).strip().upper()
    try:
        return _SPECSYS_TO_ASTROPY_FRAME[code]
    except KeyError as exc:
        valid = ", ".join(_SPECSYS_TO_ASTROPY_FRAME)
        raise MissingContextError(
            f"unknown spectral frame {specsys!r}; astrolyze knows: {valid} "
            "(it never guesses a frame, ADR-0003)"
        ) from exc


# -- sky axes -------------------------------------------------------------------------
def sky_coordinate_maps(spatial_coordinate_map) -> tuple[u.Quantity, u.Quantity]:
    """Split a spectral-cube ``spatial_coordinate_map`` into (longitude, latitude).

    spectral-cube returns the world maps in pixel-axis order ``[latitude, longitude]``; we
    return them as the named (longitude, latitude) pair so callers never juggle the order."""
    latitude, longitude = spatial_coordinate_map
    return u.Quantity(longitude), u.Quantity(latitude)


def sky_maps_from_wcs(wcs, shape) -> tuple[u.Quantity, u.Quantity]:
    """Per-pixel (longitude, latitude) maps for a 2D image from its already-parsed WCS.

    The :class:`~astrolyze.core.map.Map` holds a plain ``Quantity`` (no spectral-cube handle),
    so its sky maps are read straight off the stored celestial WCS — no header reparse. Returns
    the world coordinate of every pixel as angular ``Quantity`` arrays the shape of the image."""
    celestial = wcs.celestial
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    lon, lat = celestial.wcs_pix2world(xx, yy, 0)
    lon_unit, lat_unit = (u.Unit(c) for c in celestial.wcs.cunit)
    return u.Quantity(lon, lon_unit), u.Quantity(lat, lat_unit)


def pixel_scale(wcs) -> u.Quantity:
    """Celestial pixel scale (one value per spatial axis) as an angular ``Quantity``.

    Reads the already-parsed WCS (``proj_plane_pixel_scales`` returns degrees for a celestial
    WCS); no header reparse."""
    celestial = wcs.celestial
    scales = proj_plane_pixel_scales(celestial)
    return u.Quantity(scales, u.deg)


__all__ = [
    "AxisCoordinates",
    "Validity",
    "VALIDITY_SCHEMA_VERSION",
    "VALIDITY_METHOD",
    "validity_of",
    "absolute_frequency",
    "delta_v",
    "delta_v_relative_to",
    "to_spectral_frame",
    "sky_coordinate_maps",
    "pixel_scale",
    "MissingContextError",
]
