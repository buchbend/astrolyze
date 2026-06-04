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
from astrolyze.units import ContextGap, Insufficiency, convert

# The unit-classification helpers that route a conversion live in the converter module; the
# probe reuses them (DRY) rather than re-deriving "which equivalency does this need?". Imported
# by full dotted path because the package re-exports the ``convert`` *function* under that name.
from astrolyze.units.convert import (
    _is_frequency_like,
    _is_per_beam,
    _is_temperature,
    _is_velocity,
    _shared_integration_axis,
)


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

    # -- the non-raising probe (issue #29): "what would I need to do X?" ----------------
    def can_convert_to(self, unit) -> Insufficiency:
        """Return the :class:`~astrolyze.units.Insufficiency` describing what is missing to
        convert this object to ``unit`` — **without raising** (issue #29).

        Where :meth:`to` *does* the conversion and raises on missing context (ADR-0003), this
        *probes* it: it returns the gaps as a structured, inspectable descriptor so a caller
        can branch (``if cube.can_convert_to("K"): ...``) instead of catching an exception. An
        empty (satisfied) descriptor means the conversion can proceed.

        The gap list is decided from the requested conversion, not from a fixed schema, so it
        can name parameters that have no :class:`~astrolyze.io.Metadata` field yet
        (``calibration_scale`` for the genuinely-ambiguous RJ-vs-Planck scale, arriving in
        #24/#25): converting brightness to/from temperature names ``rest_frequency`` +
        ``calibration_scale``; anything touching Jy/beam names ``beam``."""
        return _probe_conversion(self._data_quantity.unit, unit, self.metadata)

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
        _emit("to", params={"unit": str(unit), "temperature_scale": temperature_scale})
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
        result = plotter(self, **kwargs)
        _emit("plot", params={"viz_function": self._viz_function})
        return result

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


def _probe_conversion(src, target, metadata: Metadata) -> Insufficiency:
    """Decide the :class:`~astrolyze.units.Insufficiency` for converting ``src`` -> ``target``
    given an object's ``metadata`` — the non-raising counterpart to :func:`convert`'s routing.

    It mirrors the routing in :mod:`astrolyze.units.convert` (which equivalency does *this*
    conversion need?) but, instead of raising on the first missing piece, collects *all* the
    gaps the object cannot fill and returns them. The interpretation parameters that have no
    schema field yet (the RJ-vs-Planck ``calibration_scale``) are named from the conversion
    itself, so the probe integrates cleanly with #24/#25 without depending on them."""
    gaps: list[ContextGap] = []

    # An unparseable target is itself a no-silent-physics case: name it, don't raise opaquely.
    try:
        target_unit = u.Unit(target)
    except (ValueError, TypeError):
        return Insufficiency([ContextGap.for_parameter("unit", detail=str(target))])

    spectral = _is_frequency_like(src) or _is_velocity(src)
    spectral_target = _is_frequency_like(target_unit) or _is_velocity(target_unit)

    if spectral and spectral_target:
        # Spectral-axis conversion: only crossing frequency<->velocity needs context.
        crosses = _is_velocity(src) != _is_velocity(target_unit)
        if crosses:
            if metadata.rest_frequency is None:
                gaps.append(ContextGap.for_parameter("rest_frequency"))
            if metadata.velocity_convention is None:
                gaps.append(ContextGap.for_parameter("velocity_convention"))
        return Insufficiency(gaps)

    # Intensity conversion: factor out any shared velocity-integration axis, then classify.
    axis = _shared_integration_axis(src, target_unit)
    base_src = src / axis if axis is not None else src
    base_tgt = target_unit / axis if axis is not None else target_unit

    involves_temperature = _is_temperature(base_src) != _is_temperature(base_tgt)
    involves_beam = _is_per_beam(base_src) != _is_per_beam(base_tgt)

    if involves_temperature:
        # Brightness<->temperature needs the rest frequency and the genuinely-ambiguous
        # calibration scale (RJ vs Planck) — the latter has no schema field; it is a decision
        # the object cannot make, so it is always a gap here.
        if metadata.rest_frequency is None:
            gaps.append(ContextGap.for_parameter("rest_frequency"))
        gaps.append(ContextGap.for_parameter("calibration_scale"))
    if involves_beam and metadata.beam is None:
        gaps.append(ContextGap.for_parameter("beam"))

    return Insufficiency(gaps)


def _emit(op, **fields) -> None:
    """Append a run-log record for *op* through the always-on run-log seam (ADR-0010).

    Deferred to call time so ``core`` stays free of the experiment package on import; the seam
    (:mod:`astrolyze.experiment.runlog`) is a no-op when no run is active, so ``.to``/``.plot``
    behave identically outside an experiment."""
    from astrolyze.experiment.runlog import emit

    emit(op, **fields)
