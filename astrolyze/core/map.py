"""The 2D wrapper: :class:`Map` (an image, a channel, or a moment).

Thin — it holds a 2D :class:`~astropy.units.Quantity` plus its WCS and the astrolyze
:class:`~astrolyze.io.Metadata` (beam, rest frequency, convention). It is the type a
``Cube`` moment or single-channel slice transitions into, carrying that context for free.
"""

from __future__ import annotations

import astropy.units as u

from astrolyze.io import Metadata

from . import _coords
from ._base import ContextCarrier


class Map(ContextCarrier):
    """A 2D map: a wrapped ``Quantity`` + WCS that carries its physical context."""

    _viz_function = "plot_map"

    def __init__(self, data, wcs, metadata: Metadata):
        #: the 2D values + unit (a moment0 Projection is itself a Quantity and wraps cleanly).
        self._data = u.Quantity(data)
        self._wcs = wcs
        self.metadata = metadata

    @classmethod
    def from_loaded(cls, loaded) -> "Map":
        """Wrap an :class:`~astrolyze.io.LoadedData` (data + WCS + parsed Metadata)."""
        bunit = loaded.metadata.bunit or u.dimensionless_unscaled
        return cls(u.Quantity(loaded.data, bunit), loaded.wcs, loaded.metadata)

    # -- ContextCarrier seam ------------------------------------------------------------
    @property
    def _data_quantity(self) -> u.Quantity:
        return self._data

    def _with_data(self, new_quantity: u.Quantity) -> "Map":
        return Map(new_quantity, self._wcs, self._metadata_with_unit(new_quantity.unit))

    # -- coordinate-array + validity emission (issue #26) -------------------------------
    @property
    def coordinates(self) -> _coords.AxisCoordinates:
        """The sky / pixel coordinate subset, read from the stored WCS (no reparse). A 2D map
        has no spectral axis, so ``frequency`` / ``delta_v`` are ``None``."""
        longitude, latitude = _coords.sky_maps_from_wcs(self._wcs, self.shape)
        return _coords.AxisCoordinates(
            longitude=longitude,
            latitude=latitude,
            pixel_scale=_coords.pixel_scale(self._wcs),
        )

    @property
    def validity(self) -> _coords.Validity:
        """Validity descriptor: blanked voxels as ``NaN`` + a boolean finite-data mask."""
        return _coords.validity_of(self._data)

    # -- thin passthroughs --------------------------------------------------------------
    @property
    def data(self) -> u.Quantity:
        return self._data

    @property
    def wcs(self):
        return self._wcs

    @property
    def unit(self) -> u.UnitBase:
        return self._data.unit

    @property
    def shape(self) -> tuple[int, ...]:
        return self._data.shape

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        return f"<Map {obj} {self.shape} [{self.unit}]>"
