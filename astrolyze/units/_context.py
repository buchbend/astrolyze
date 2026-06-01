"""Explicit-context types for the unit layer (ADR-0003).

The unit layer refuses to guess the physics that silently biases radio/sub-mm work:
the velocity/Doppler **convention**, the **rest frequency**, and the brightness-temperature
**scale** (Rayleigh-Jeans vs Planck). These types name those choices so call sites must
state them, and ``MissingContextError`` is what we raise when one is absent.
"""

from __future__ import annotations

from enum import Enum


class MissingContextError(ValueError):
    """A conversion needs physical context (convention / rest frequency / beam / scale)
    that was not supplied. Raised instead of silently defaulting — see ADR-0003."""


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
