"""Shared context-carrying behaviour for the Cube / Map / Spectrum trio (ADR-0004).

The trio are **thin wrappers that compose** spectral-cube / astropy / specutils: each holds a
private upstream handle plus the astrolyze :class:`~astrolyze.io.Metadata` that is the single
source of truth for the physical context (beam + rest frequency + velocity convention +
provenance). This mixin factors out the two pieces of behaviour every wrapper shares — the
unit hub ``.to()`` (ADR-0003c) and the ``.plot()`` seam (ADR-0005) — over a small subclass
seam (``_data_quantity`` / ``_with_data``), so there is one implementation, three ergonomics.
"""

from __future__ import annotations

from dataclasses import replace

import astropy.units as u

from astrolyze.io import Metadata
from astrolyze.units import convert


class ContextCarrier:
    """Mixin giving a wrapper its context shortcuts, the unit hub, and the plot seam.

    Subclasses set ``self.metadata`` (an :class:`~astrolyze.io.Metadata`) and ``_viz_function``
    (the name of the free function in :mod:`astrolyze.viz` that plots them), and implement the
    data seam :pyattr:`_data_quantity` / :pymeth:`_with_data`.
    """

    metadata: Metadata
    #: Name of the free plotting function this type delegates to (filled by issue #6's viz).
    _viz_function: str = ""

    # -- context shortcuts (read from the metadata, never cached/duplicated) ------------
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
    def is_complete(self) -> bool:
        return self.metadata.is_complete

    @property
    def missing(self) -> list[str]:
        return self.metadata.missing

    def require_complete(self) -> None:
        """Raise :class:`~astrolyze.units.MissingContextError` if mandatory context is absent.

        The explicit door for velocity / spectral-axis operations, which always need *both*
        the rest frequency and the convention. (Intensity ``.to()`` does not call this: it
        delegates enforcement to the unit layer, which demands exactly — and only — the
        context the specific conversion needs, e.g. a beam-only Jy/beam<->Jy/sr conversion
        must not be blocked for want of a velocity convention it never uses.)"""
        self.metadata.ensure_complete()

    # -- the unit hub (ADR-0003c): supply this object's context, never guess ------------
    def to(
        self,
        unit,
        *,
        temperature_scale=None,
        convention=None,
        rest_frequency=None,
        beam=None,
    ):
        """Convert to ``unit`` through :func:`astrolyze.units.convert`, supplying the
        object's beam / rest frequency / convention so the caller need not.

        Returns a new wrapper of the same type carrying the converted data and an updated
        ``bunit``. The genuinely-ambiguous Rayleigh-Jeans-vs-Planck *scale* is deliberately
        **not** sourced from the object (the object cannot know it) and is never defaulted —
        a brightness-temperature conversion without ``temperature_scale`` raises (ADR-0003).
        Per-call keyword arguments override the object's context when given.
        """
        converted = convert(
            self._data_quantity,
            unit,
            rest_frequency=(
                rest_frequency
                if rest_frequency is not None
                else self.metadata.rest_frequency
            ),
            convention=(
                convention
                if convention is not None
                else self.metadata.velocity_convention
            ),
            beam=beam if beam is not None else self.metadata.beam,
            temperature_scale=temperature_scale,
        )
        return self._with_data(converted)

    # -- the display seam (ADR-0005): thin sugar over the free viz engine ---------------
    def plot(self, **kwargs):
        """Plot via the :mod:`astrolyze.viz` engine (the free ``plot_*`` functions).

        Object context (units, beam, WCS) makes the result auto-correct. The engine itself
        arrives in issue #6; until then this raises a clear :class:`NotImplementedError`
        rather than an obscure import/attribute error."""
        from astrolyze import viz  # lazy: keeps importing core free of matplotlib

        plotter = getattr(viz, self._viz_function, None)
        if plotter is None:
            raise NotImplementedError(
                f"plotting needs the viz layer ({self._viz_function}); it arrives in issue #6"
            )
        return plotter(self, **kwargs)

    # -- subclass seam ------------------------------------------------------------------
    @property
    def _data_quantity(self) -> u.Quantity:  # pragma: no cover - abstract
        """The wrapped values as a single :class:`~astropy.units.Quantity` (what ``.to()``
        converts)."""
        raise NotImplementedError

    def _with_data(self, new_quantity: u.Quantity):  # pragma: no cover - abstract
        """Build a new wrapper of this type from converted ``new_quantity``, preserving the
        rest of the context."""
        raise NotImplementedError

    def _metadata_with_unit(self, unit) -> Metadata:
        """This object's metadata with ``bunit`` set to ``unit`` (used after a conversion or a
        moment changes the physical unit)."""
        return replace(self.metadata, bunit=u.Unit(unit))
