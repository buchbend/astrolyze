"""Typed-insufficiency descriptor — a returnable "what I'd need" value object (issue #29).

The unit layer already refuses to guess the physics that silently biases radio/sub-mm work
(ADR-0003): a conversion that lacks the rest frequency, velocity convention, beam, or
brightness scale *raises* rather than defaulting. This module makes that refusal a
**first-class, inspectable, returnable** result so a caller can also ask "what would I need
to do X?" and branch on the answer without catching an exception (ADR-0013: results are
traceable and correctable).

Two value objects:

- :class:`ContextGap` — one missing-or-ambiguous interpretation parameter, as ``name`` +
  human ``reason`` (+ optional ``how_to_supply``). A gap is *just* name+reason: it can name a
  parameter that has no schema field yet (``spectral_frame`` / ``calibration_scale`` /
  ``eta_mb`` arrive in #24/#25), because naming a gap is decoupled from the schema shape.
- :class:`Insufficiency` — the descriptor: an ordered list of gaps, falsy when empty (no
  gaps == satisfied), so ``if cube.can_convert_to("K"): ...`` reads naturally.

``MissingContextError`` carries an :class:`Insufficiency` so a raise exposes the same
structured detail a probe returns.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Canned human reasons + how-to-supply text for every interpretation parameter the
# no-silent-physics guard can name (issue #29). Keyed by gap name. The three not yet on the
# Metadata schema (spectral_frame / calibration_scale / eta_mb, arriving in #24/#25) are
# listed here too: a gap names a parameter regardless of whether a field exists for it.
# ``{detail}`` is filled per occurrence (e.g. the offending unit, or the target).
_PARAMETER_GAPS: dict[str, dict[str, str]] = {
    "rest_frequency": {
        "reason": (
            "rest frequency is required (it sets the frequency of the Planck/RJ "
            "brightness law and the line frequency for velocity conversions); "
            "astrolyze never assumes one"
        ),
        "how_to_supply": (
            "set RESTFRQ on the FITS header, or pass rest_frequency=<Quantity> to the "
            "conversion"
        ),
    },
    "velocity_convention": {
        "reason": (
            "a velocity convention (radio | optical | relativistic) is required for "
            "frequency<->velocity conversions; the three give different velocities for the "
            "same line, so astrolyze never assumes one"
        ),
        "how_to_supply": (
            "set the ASTROLYZE VCONV header card (or a VRAD/VOPT/VELO CTYPE), or pass "
            "convention="
        ),
    },
    "beam": {
        "reason": (
            "a beam solid angle is required for any conversion involving Jy/beam (the "
            "per-beam<->per-sr geometry); a bare FWHM is not enough because it assumes a "
            "beam shape"
        ),
        "how_to_supply": (
            "set BMAJ/BMIN/BPA on the header, or pass beam=<radio_beam.Beam>"
        ),
    },
    "spectral_frame": {
        "reason": (
            "the spectral reference frame (SPECSYS, e.g. LSRK/BARYCENT) is required to "
            "interpret the spectral axis; astrolyze never assumes one"
        ),
        "how_to_supply": "set the SPECSYS header card",
    },
    "calibration_scale": {
        "reason": (
            "the brightness-temperature calibration scale (T_mb | T_A* | T_R*, i.e. "
            "Rayleigh-Jeans vs Planck) is genuinely ambiguous and the top silent-error trap "
            "in radio/sub-mm work, so astrolyze never picks one"
        ),
        "how_to_supply": (
            "state the scale on the conversion (e.g. temperature_scale='rayleigh_jeans' or "
            "'planck')"
        ),
    },
    "eta_mb": {
        "reason": (
            "the main-beam efficiency eta_mb is required to move between antenna and "
            "main-beam temperature scales; astrolyze never assumes a value"
        ),
        "how_to_supply": "supply eta_mb for this telescope/frequency",
    },
    "unit": {
        "reason": "the unit {detail!r} was not recognised",
        "how_to_supply": (
            "use a unit astropy/astrolyze understands (e.g. 'K', 'Jy/beam', 'MJy/sr', "
            "'K km/s')"
        ),
    },
}


@dataclass(frozen=True)
class ContextGap:
    """One missing-or-ambiguous interpretation parameter the guard can name.

    ``name`` is the machine token (e.g. ``"rest_frequency"``), ``reason`` the human
    explanation, and ``how_to_supply`` (optional) how a caller would provide it. A gap is
    deliberately *not* tied to a schema field — it can name a parameter that has no
    :class:`~astrolyze.io.Metadata` field yet.
    """

    name: str
    reason: str
    how_to_supply: str | None = None

    @classmethod
    def for_parameter(cls, name: str, *, detail: str | None = None) -> "ContextGap":
        """Build the gap for a known interpretation parameter with its canned reason.

        ``detail`` is interpolated into reasons that quote a value (currently the
        unrecognised ``unit``), so the message names the exact offending text. Unknown
        names still produce a usable, honest gap rather than raising — the descriptor's
        job is to *describe*, not to gatekeep its own vocabulary.
        """
        spec = _PARAMETER_GAPS.get(name)
        if spec is None:
            return cls(name=name, reason=f"{name} is required but was not supplied")
        reason = spec["reason"].format(detail=detail)
        return cls(name=name, reason=reason, how_to_supply=spec.get("how_to_supply"))


@dataclass(frozen=True)
class Insufficiency:
    """The descriptor: the ordered set of gaps that block (or qualify) an operation.

    Falsy when empty so ``if insufficiency:`` means "there are gaps". An empty descriptor is
    the *satisfied* result a non-raising probe returns when the operation can proceed.
    """

    gaps: list[ContextGap] = field(default_factory=list)

    @property
    def names(self) -> list[str]:
        """The gap names, in order (handy for membership checks and reporting)."""
        return [g.name for g in self.gaps]

    @property
    def is_satisfied(self) -> bool:
        """True when there are no gaps (the operation can proceed)."""
        return not self.gaps

    def get(self, name: str) -> ContextGap | None:
        """The gap with this ``name``, or ``None`` if it is not among the gaps."""
        return next((g for g in self.gaps if g.name == name), None)

    def __bool__(self) -> bool:
        return bool(self.gaps)

    def __iter__(self):
        return iter(self.gaps)

    def __len__(self) -> int:
        return len(self.gaps)

    def message(self) -> str:
        """A legible one-line-per-gap rendering (ADR-0013: surface the derivation)."""
        if not self.gaps:
            return "no missing context"
        return "; ".join(f"{g.name}: {g.reason}" for g in self.gaps)

    def __str__(self) -> str:
        return self.message()


__all__ = ["ContextGap", "Insufficiency"]
