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


@dataclass(frozen=True)
class Validity:
    """Where the data is real: blanked voxels as ``NaN`` + a boolean finite-data mask.

    ``mask`` is ``True`` exactly where ``data`` is finite (blanked / edge / outside-coverage
    voxels are ``False`` and appear as ``NaN`` in ``data``). Because both are read from the
    same array, slicing the wrapper slices both consistently — the mask travels with the data
    (mask-of-a-slice == slice-of-the-mask)."""

    #: The values with blanked voxels exposed as ``NaN`` (carries the data's unit).
    data: u.Quantity
    #: Boolean finite-data mask, the shape of the data (``True`` where the data is real).
    mask: np.ndarray


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
    "validity_of",
    "absolute_frequency",
    "delta_v",
    "sky_coordinate_maps",
    "pixel_scale",
    "MissingContextError",
]
