"""Unit aliases, equivalency bundles, and converters (the deep module — ADR-0003).

Pure-astropy, no I/O, no import side effects. The velocity convention and rest frequency
(and the brightness-temperature scale) are explicit and mandatory where a conversion needs
them — astrolyze raises rather than assuming, because a silently assumed convention, rest
frequency, or Rayleigh-Jeans-vs-Planck scale is a classic hard-to-detect bias in this
science.

Two entry points over one implementation: the standalone :func:`convert`, and (later) the
data object's ``.to()`` which supplies beam/frequency/convention from its own context.
"""

from __future__ import annotations

from . import aliases, equivalencies
from ._context import (
    BrightnessTemperatureScale,
    MissingContextError,
    VelocityConvention,
    coerce_velocity_convention,
)
from .aliases import Jy_beam, Jy_sr, K_kms, MJy_sr, Ta, Tmb
from .convert import convert
from .equivalencies import (
    beam_angular_area,
    beam_solid_angle,
    brightness_temperature,
    doppler,
    spectral,
)

__all__ = [
    # converter
    "convert",
    # explicit-context types
    "VelocityConvention",
    "BrightnessTemperatureScale",
    "MissingContextError",
    "coerce_velocity_convention",
    # equivalency builders
    "brightness_temperature",
    "beam_angular_area",
    "beam_solid_angle",
    "doppler",
    "spectral",
    # named unit aliases
    "Tmb",
    "Ta",
    "Jy_beam",
    "Jy_sr",
    "MJy_sr",
    "K_kms",
    # submodules
    "aliases",
    "equivalencies",
]
