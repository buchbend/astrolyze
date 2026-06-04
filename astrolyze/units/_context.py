"""Explicit-context types for the unit layer (ADR-0003).

The unit layer refuses to guess the physics that silently biases radio/sub-mm work:
the velocity/Doppler **convention**, the **rest frequency**, and the brightness-temperature
**scale** (Rayleigh-Jeans vs Planck). These types name those choices so call sites must
state them, and ``MissingContextError`` is what we raise when one is absent.
"""

from __future__ import annotations

from enum import Enum

from ._insufficiency import Insufficiency


class MissingContextError(ValueError):
    """A conversion needs physical context (convention / rest frequency / beam / scale)
    that was not supplied. Raised instead of silently defaulting — see ADR-0003.

    Carries an optional :class:`~astrolyze.units.Insufficiency` descriptor (``insufficiency``)
    so a raise exposes the same structured "what I'd need" detail a non-raising probe returns
    (issue #29; ADR-0013) — not only a message string. Backward compatible: existing raises
    that pass only a message keep ``insufficiency is None``."""

    def __init__(self, *args, insufficiency: Insufficiency | None = None):
        super().__init__(*args)
        self.insufficiency = insufficiency


class VelocityConvention(str, Enum):
    """Doppler convention relating frequency to velocity. There is no default: the radio,
    optical, and relativistic conventions give different velocities for the same line, and
    a silently assumed one is a classic, hard-to-detect bias."""

    RADIO = "radio"
    OPTICAL = "optical"
    RELATIVISTIC = "relativistic"


class BrightnessTemperatureScale(str, Enum):
    """Brightness-temperature scale. Rayleigh-Jeans is the linear low-frequency
    approximation (and the conventional radio intensity scale); Planck inverts the full
    Planck law and is correct into the Wien regime. astropy's built-in equivalency is
    RJ-only, so the choice is made explicit here (ADR-0003)."""

    RAYLEIGH_JEANS = "rayleigh_jeans"
    PLANCK = "planck"


class CalibrationScale(str, Enum):
    """The single-dish calibration *temperature scale* a Kelvin-valued cube is on (issue #25).

    This is a DIFFERENT concept from :class:`BrightnessTemperatureScale` (the Rayleigh-Jeans
    vs Planck brightness LAW): it names which observed temperature scale the calibration put
    the data on, and so whether a beam-efficiency correction is owed before the data is a
    main-beam brightness temperature:

    - ``T_mb`` — main-beam brightness temperature: already the source scale; convert directly.
    - ``T_A*`` — antenna temperature corrected for the atmosphere/spillover/ohmic losses but
      *not* the main-beam efficiency; ``T_mb = T_A*/eta_mb`` must be applied first.
    - ``T_R*`` — the forward-beam-corrected scale; treated as directly comparable to ``T_mb``
      for the point-vs-extended distinction astrolyze does not yet model (no extra factor).

    Kelvin is genuinely ambiguous about which of these it is, so astrolyze never assumes one
    (ADR-0003): a temperature cube must *declare* its scale before it can be harmonised."""

    T_MB = "T_mb"
    T_A_STAR = "T_A*"
    T_R_STAR = "T_R*"


def coerce_calibration_scale(value) -> CalibrationScale:
    """Accept a ``CalibrationScale`` or a string (case-insensitive on the leading ``T``);
    raise on anything else. The ``*``/``mb`` body is matched verbatim so ``T_A*`` round-trips."""
    if isinstance(value, CalibrationScale):
        return value
    text = str(value)
    try:
        return CalibrationScale(text)
    except ValueError:
        # Allow a forgiving case on the prefix (t_mb -> T_mb) without inventing new spellings.
        for scale in CalibrationScale:
            if text.lower() == scale.value.lower():
                return scale
        valid = ", ".join(s.value for s in CalibrationScale)
        raise ValueError(
            f"unknown calibration scale {value!r}; expected one of: {valid}"
        ) from None


def coerce_velocity_convention(value) -> VelocityConvention:
    """Accept a ``VelocityConvention`` or a case-insensitive string; raise on anything else."""
    if isinstance(value, VelocityConvention):
        return value
    try:
        return VelocityConvention(str(value).lower())
    except ValueError as exc:
        valid = ", ".join(c.value for c in VelocityConvention)
        raise ValueError(
            f"unknown velocity convention {value!r}; expected one of: {valid}"
        ) from exc


def coerce_temperature_scale(value) -> BrightnessTemperatureScale:
    """Accept a ``BrightnessTemperatureScale`` or a case-insensitive string; raise otherwise."""
    if isinstance(value, BrightnessTemperatureScale):
        return value
    try:
        return BrightnessTemperatureScale(str(value).lower())
    except ValueError as exc:
        valid = ", ".join(s.value for s in BrightnessTemperatureScale)
        raise ValueError(
            f"unknown brightness-temperature scale {value!r}; expected one of: {valid}"
        ) from exc
