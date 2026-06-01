"""The radio unit converter (ADR-0003).

``convert`` is the single entry point that supplies the right equivalencies for a given
pair of units and *refuses to guess* the physics that silently biases this science: the
velocity convention, rest frequency, and brightness-temperature scale are explicit and
mandatory where a conversion depends on them (a missing one raises ``MissingContextError``,
it never defaults).

It delegates the actual maths to ``astropy`` equivalencies; the value it adds is routing
(which equivalency does *this* conversion need?) and the no-silent-physics contract.
"""

from __future__ import annotations

import astropy.units as u

from . import equivalencies as eq
from ._context import (
    BrightnessTemperatureScale,
    MissingContextError,
    coerce_temperature_scale,
)

# Velocity units that may ride along a brightness conversion as an integration axis.
_INTEGRATION_AXES = (u.km / u.s, u.m / u.s)


# -- unit classification ---------------------------------------------------------------
def _is_temperature(unit) -> bool:
    return unit.is_equivalent(u.K)


def _is_per_beam(unit) -> bool:
    return unit.is_equivalent(u.Jy / u.beam)


def _is_surface_brightness(unit) -> bool:
    return unit.is_equivalent(u.Jy / u.sr)


def _is_base_intensity(unit) -> bool:
    return _is_temperature(unit) or _is_per_beam(unit) or _is_surface_brightness(unit)


def _is_frequency_like(unit) -> bool:
    return unit.is_equivalent(u.Hz, equivalencies=u.spectral())


def _is_velocity(unit) -> bool:
    return unit.is_equivalent(u.km / u.s)


def _shared_integration_axis(src, target):
    """Return the velocity unit factored out of *both* src and target (a velocity-integrated
    intensity such as K km/s), or ``None`` if neither side is integrated."""
    if _is_base_intensity(src) and _is_base_intensity(target):
        return None
    for axis in _INTEGRATION_AXES:
        if _is_base_intensity(src / axis) and _is_base_intensity(target / axis):
            return axis
    return None


# -- conversion branches ---------------------------------------------------------------
def _convert_spectral_axis(quantity, target, rest_frequency, convention):
    src = quantity.unit
    if _is_velocity(src) == _is_velocity(target):
        # both frequency-like (freq<->wavelength<->energy) or both velocity: no convention.
        return quantity.to(target, equivalencies=u.spectral())
    # crossing the frequency<->velocity boundary needs the Doppler convention + rest freq.
    doppler = eq.doppler(convention, rest_frequency)
    return quantity.to(target, equivalencies=doppler + u.spectral())


def _convert_intensity(quantity, target, *, rest_frequency, beam, temperature_scale):
    src = quantity.unit
    axis = _shared_integration_axis(src, target)
    base_src = src / axis if axis is not None else src
    base_tgt = target / axis if axis is not None else target

    src_temp, tgt_temp = _is_temperature(base_src), _is_temperature(base_tgt)
    src_pb, tgt_pb = _is_per_beam(base_src), _is_per_beam(base_tgt)

    involves_temperature = src_temp != tgt_temp  # temperature on exactly one side
    involves_beam = src_pb != tgt_pb  # crossing per-beam <-> sr/K

    equivalencies = []
    if involves_temperature:
        if rest_frequency is None:
            raise MissingContextError(
                "rest_frequency is required to convert between brightness temperature and "
                "flux/surface brightness (it sets the frequency of the Planck/RJ law)"
            )
        if temperature_scale is None:
            raise MissingContextError(
                "temperature_scale is required (rayleigh_jeans | planck): RJ-vs-Planck is "
                "the top silent-error trap in radio/sub-mm work, so astrolyze never assumes it"
            )
        scale = coerce_temperature_scale(temperature_scale)
        if axis is not None and scale is BrightnessTemperatureScale.PLANCK:
            raise ValueError(
                "Planck brightness temperature is nonlinear and cannot be applied to a "
                "velocity-integrated intensity (the integral does not commute with B_nu); "
                "integrated intensities are on the Rayleigh-Jeans scale by construction"
            )
        beam_area = eq.beam_solid_angle(beam) if involves_beam else None
        equivalencies = eq.brightness_temperature(
            rest_frequency, scale=scale, beam_area=beam_area
        )
    elif involves_beam:
        # Pure geometry (per-beam <-> sr); beam_angular_area raises if beam is None.
        equivalencies = eq.beam_angular_area(beam)

    if axis is None:
        return quantity.to(target, equivalencies=equivalencies)

    # Velocity-integrated: the brightness step is linear here (Planck is rejected above and
    # pure geometry is linear), so it commutes with the velocity integral — convert the
    # per-channel factor and carry the velocity width through.
    factor = (1.0 * base_src).to(base_tgt, equivalencies=equivalencies)
    return (quantity.value * factor.value) * target


def convert(
    quantity,
    target,
    *,
    rest_frequency=None,
    convention=None,
    beam=None,
    temperature_scale=None,
):
    """Convert a radio/sub-mm ``Quantity`` to ``target`` units, supplying the right
    equivalencies and demanding the physical context each conversion needs.

    Parameters
    ----------
    quantity : astropy.units.Quantity
        The value to convert.
    target : astropy.units.Unit or str
        Desired output units.
    rest_frequency : astropy.units.Quantity, optional
        Line/observing frequency. **Required** for frequency<->velocity conversions and for
        brightness-temperature <-> flux/surface-brightness conversions. Missing -> raises.
    convention : VelocityConvention or {"radio","optical","relativistic"}, optional
        Doppler convention. **Required** for frequency<->velocity conversions. Missing -> raises.
    beam : radio_beam.Beam or astropy.units.Quantity, optional
        Beam (object) or beam solid angle. **Required** for any conversion involving Jy/beam.
    temperature_scale : BrightnessTemperatureScale or {"rayleigh_jeans","planck"}, optional
        Brightness-temperature scale. **Required** for brightness-temperature conversions.
        Missing -> raises. Planck is rejected for velocity-integrated intensities (nonlinear).

    Returns
    -------
    astropy.units.Quantity
        ``quantity`` expressed in ``target`` units.

    Raises
    ------
    MissingContextError
        When a required convention / rest frequency / beam / temperature scale is absent.
    """
    quantity = u.Quantity(quantity)
    target = u.Unit(target)
    src = quantity.unit

    if (_is_frequency_like(src) or _is_velocity(src)) and (
        _is_frequency_like(target) or _is_velocity(target)
    ):
        return _convert_spectral_axis(quantity, target, rest_frequency, convention)

    return _convert_intensity(
        quantity,
        target,
        rest_frequency=rest_frequency,
        beam=beam,
        temperature_scale=temperature_scale,
    )


__all__ = ["convert"]
