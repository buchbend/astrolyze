"""The 2D wrapper: :class:`Map` (an image, a channel, or a moment).

Thin — it holds a 2D :class:`~astropy.units.Quantity` plus its WCS and the astrolyze
:class:`~astrolyze.io.Metadata` (beam, rest frequency, convention). It is the type a
``Cube`` moment or single-channel slice transitions into, carrying that context for free.
"""

from __future__ import annotations

import astropy.units as u

from astrolyze.io import Metadata

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
