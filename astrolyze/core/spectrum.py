"""The 1D wrapper: :class:`Spectrum`, a thin adaptor over :mod:`specutils` (ADR-0004 rider 2).

We do not reinvent spectral representation/fitting — the wrapper holds a
:class:`specutils.Spectrum` (the flux + spectral axis) plus the astrolyze
:class:`~astrolyze.io.Metadata`, and is the type a ``Cube`` spatial index transitions into.
"""

from __future__ import annotations

import astropy.units as u

# specutils 2.x exposes ``Spectrum`` (``Spectrum1D`` is a deprecated alias); fall back so the
# wrapper works across the 1.x/2.x boundary.
try:
    from specutils import Spectrum as _SpecutilsSpectrum
except ImportError:  # pragma: no cover - older specutils
    from specutils import Spectrum1D as _SpecutilsSpectrum

from astrolyze.io import Metadata

from ._base import ContextCarrier


class Spectrum(ContextCarrier):
    """A 1D spectrum: a wrapped :class:`specutils.Spectrum` that carries physical context."""

    _viz_function = "plot_spectrum"

    def __init__(self, spectrum, metadata: Metadata):
        self._spec = spectrum  # a specutils.Spectrum
        self.metadata = metadata

    @classmethod
    def from_oned(cls, oned, metadata: Metadata) -> "Spectrum":
        """Adapt a spectral-cube ``OneDSpectrum`` (e.g. from ``cube[:, y, x]``) into a
        specutils-backed :class:`Spectrum`, attaching the parent's context unchanged."""
        flux = oned.quantity if hasattr(oned, "quantity") else u.Quantity(oned)
        spec = _SpecutilsSpectrum(flux=flux, spectral_axis=oned.spectral_axis)
        return cls(spec, metadata)

    # -- ContextCarrier seam ------------------------------------------------------------
    @property
    def _data_quantity(self) -> u.Quantity:
        return self._spec.flux

    def _with_data(self, new_quantity: u.Quantity) -> "Spectrum":
        spec = _SpecutilsSpectrum(
            flux=new_quantity, spectral_axis=self._spec.spectral_axis
        )
        return Spectrum(spec, self._metadata_with_unit(new_quantity.unit))

    # -- thin passthroughs --------------------------------------------------------------
    @property
    def flux(self) -> u.Quantity:
        return self._spec.flux

    @property
    def spectral_axis(self):
        return self._spec.spectral_axis

    @property
    def unit(self) -> u.UnitBase:
        return self._spec.flux.unit

    def __repr__(self) -> str:
        obj = self.metadata.object or "?"
        return f"<Spectrum {obj} ({self.flux.size},) [{self.unit}]>"
