"""The PPV wrapper: :class:`Cube`, a thin composition over :class:`spectral_cube.SpectralCube`.

The heavy lifting — moments, the spectral axis, slicing — is delegated to ``spectral-cube``;
``Cube`` adds the astrolyze context (beam + rest frequency + velocity convention, sourced from
:class:`~astrolyze.io.Metadata`) and the type transitions that carry it: ``moment0() -> Map``
and ``cube[:, y, x] -> Spectrum`` (ADR-0004). It is **not** a ``SpectralCube`` subclass — it
holds one (``._sc``) and delegates, so we stay loosely coupled to upstream internals (ADR-0003).
"""

from __future__ import annotations

import astropy.units as u
from spectral_cube import SpectralCube

from astrolyze.io import Metadata

from ._base import ContextCarrier, _emit
from .map import Map
from .spectrum import Spectrum


class Cube(ContextCarrier):
    """A PPV cube: a wrapped :class:`~spectral_cube.SpectralCube` that carries context."""

    _viz_function = "plot_cube"

    def __init__(self, spectral_cube: SpectralCube, metadata: Metadata):
        self._sc = spectral_cube
        self.metadata = metadata

    @classmethod
    def from_loaded(cls, loaded) -> "Cube":
        """Build a :class:`Cube` from an :class:`~astrolyze.io.LoadedData` (the io seam).

        The beam is attached at construction rather than via ``with_beam`` — spectral-cube
        refuses ``with_beam`` on a Jy/beam cube (it would change the spatial resolution), but
        we are only recording the beam the header already states, not smoothing.
        """
        bunit = loaded.metadata.bunit or u.dimensionless_unscaled
        data = u.Quantity(loaded.data, bunit)
        sc = _build_cube(data, loaded.wcs, loaded.metadata.beam)
        return cls(sc, loaded.metadata)

    # -- delegated PPV operations (type transitions carry context) ----------------------
    def moment0(self) -> Map:
        """The zeroth moment (velocity-integrated intensity) as a :class:`Map`."""
        return self.moment(order=0)

    def moment(self, order: int = 0, axis: int = 0) -> Map:
        """Delegate the moment to spectral-cube and wrap the 2D result, carrying context.

        The moment changes the physical unit (e.g. Jy/beam -> Jy/beam km/s); the beam, rest
        frequency and convention are preserved on the resulting :class:`Map`."""
        proj = self._sc.moment(order=order, axis=axis)
        _emit("moment", params={"order": order, "axis": axis})
        return Map(
            proj, getattr(proj, "wcs", None), self._metadata_with_unit(proj.unit)
        )

    def __getitem__(self, key):
        """Slice the underlying cube and wrap by resulting dimensionality, carrying context:
        3D -> :class:`Cube`, 2D channel -> :class:`Map`, 1D spectrum -> :class:`Spectrum`.
        """
        result = self._sc[key]
        if isinstance(result, SpectralCube):
            return Cube(result, self.metadata)
        ndim = getattr(result, "ndim", None)
        if ndim == 1:  # spectral_cube.OneDSpectrum
            return Spectrum.from_oned(result, self.metadata)
        if ndim == 2:  # a channel map (slicing does not change the unit)
            return Map(
                result,
                getattr(result, "wcs", None),
                self._metadata_with_unit(result.unit),
            )
        return result  # 0D / anything else: hand back the raw upstream value

    # -- ContextCarrier seam ------------------------------------------------------------
    @property
    def _data_quantity(self) -> u.Quantity:
        return self._sc.filled_data[:]

    def _with_data(self, new_quantity: u.Quantity) -> "Cube":
        sc = _build_cube(new_quantity, self._sc.wcs, self.metadata.beam)
        return Cube(sc, self._metadata_with_unit(new_quantity.unit))

    # -- thin passthroughs --------------------------------------------------------------
    @property
    def spectral_axis(self):
        return self._sc.spectral_axis

    @property
    def unit(self) -> u.UnitBase:
        return self._sc.unit

    @property
    def shape(self) -> tuple[int, ...]:
        return self._sc.shape

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        return f"<Cube {obj} {self.shape} [{self.unit}]>"


def _build_cube(data: u.Quantity, wcs, beam) -> SpectralCube:
    """Construct a SpectralCube, attaching the beam at build time when present."""
    if beam is not None:
        return SpectralCube(data=data, wcs=wcs, beam=beam)
    return SpectralCube(data=data, wcs=wcs)
