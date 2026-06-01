"""Equivalency bundles for radio/sub-mm unit conversion (ADR-0003).

Thin builders over ``astropy.units`` equivalencies, plus the one piece astropy does not
provide: an **exact Planck** brightness-temperature equivalency (astropy's built-in
``brightness_temperature`` is the Rayleigh-Jeans approximation only). Everything returns an
``astropy.units.Equivalency`` usable with ``Quantity.to(...)``.
"""

from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy import constants as const

from ._context import (
    BrightnessTemperatureScale,
    MissingContextError,
    coerce_temperature_scale,
    coerce_velocity_convention,
)


def spectral():
    """Frequency <-> wavelength <-> wavenumber <-> photon energy. Needs no rest frequency."""
    return u.spectral()


def beam_solid_angle(beam) -> u.Quantity:
    """Coerce a beam specification to a solid angle in steradian.

    Accepts a ``radio_beam.Beam`` (uses its ``.sr``) or an ``astropy`` solid-angle
    ``Quantity`` (e.g. ``sr``, ``arcsec**2``). A bare angular size (a FWHM) is rejected:
    turning a FWHM into a solid angle assumes a beam shape, which is exactly the kind of
    silent step ADR-0003 forbids — pass a ``radio_beam.Beam`` to make the shape explicit.
    """
    if beam is None:
        raise MissingContextError(
            "a beam is required for conversions involving Jy/beam; pass beam=<radio_beam.Beam> "
            "or a solid-angle Quantity (e.g. beam=2.66e-9*u.sr)"
        )
    # A radio_beam.Beam (or anything exposing a solid angle) is used directly.
    if hasattr(beam, "sr"):
        return beam.sr
    return u.Quantity(beam).to(u.sr)


def beam_angular_area(beam):
    """``beam`` <-> ``sr`` equivalency for a given beam (delegates to astropy)."""
    return u.beam_angular_area(beam_solid_angle(beam))


def doppler(convention, rest_frequency):
    """Doppler equivalency (frequency <-> velocity) for the given convention and rest
    frequency. Both are mandatory: a missing one raises rather than defaulting (ADR-0003).
    """
    if convention is None:
        raise MissingContextError(
            "a velocity convention (radio | optical | relativistic) is required for "
            "frequency<->velocity conversions; astrolyze never assumes one"
        )
    if rest_frequency is None:
        raise MissingContextError(
            "rest_frequency is required for frequency<->velocity conversions"
        )
    conv = coerce_velocity_convention(convention)
    nu0 = u.Quantity(rest_frequency)
    builder = {
        "radio": u.doppler_radio,
        "optical": u.doppler_optical,
        "relativistic": u.doppler_relativistic,
    }[conv.value]
    return builder(nu0)


def _planck_brightness_temperature(frequency, beam_area=None):
    """Exact Planck brightness-temperature equivalency.

    Inverts the full Planck law  I_nu = (2 h nu^3 / c^2) / (exp(h nu / k T) - 1)  rather than
    the Rayleigh-Jeans I_nu = 2 k nu^2 T / c^2 approximation. Registered between Jy/sr and K
    (or Jy/beam and K when a beam area is given), so astropy supplies the MJy<->Jy scaling.
    """
    nu = u.Quantity(frequency).to(u.Hz, equivalencies=u.spectral())
    # 2 h nu^3 / c^2 decomposes to spectral flux density (Jy); the per-steradian of specific
    # intensity is conventional (sr is dimensionless-valued), so the Jy value is the Jy/sr
    # numeric scale.
    A = (2 * const.h * nu**3 / const.c**2).to_value(u.Jy)  # I_nu scale [Jy/sr]
    xpl = (const.h * nu / const.k_B).to_value(u.K)  # h nu / k  [K]

    if beam_area is None:

        def to_K(intensity):  # intensity in Jy/sr
            return xpl / np.log1p(A / intensity)

        def from_K(temperature):  # temperature in K
            return A / np.expm1(xpl / temperature)

        return u.Equivalency(
            [(u.Jy / u.sr, u.K, to_K, from_K)],
            "planck_brightness_temperature",
            {"frequency": frequency},
        )

    omega = u.Quantity(beam_area).to_value(u.sr)

    def to_K(flux):  # flux in Jy/beam
        return xpl / np.log1p(A / (flux / omega))

    def from_K(temperature):  # temperature in K
        return (A / np.expm1(xpl / temperature)) * omega

    return u.Equivalency(
        [(u.Jy / u.beam, u.K, to_K, from_K)],
        "planck_brightness_temperature",
        {"frequency": frequency, "beam_area": beam_area},
    )


def brightness_temperature(frequency, *, scale, beam_area=None):
    """Brightness-temperature equivalency between surface brightness / per-beam flux and K.

    ``scale`` selects Rayleigh-Jeans (linear, astropy's built-in) or Planck (exact, this
    module's). It is mandatory — RJ-vs-Planck is the #1 silent-error trap in this science
    (ADR-0003), so there is no default. ``beam_area`` (a solid angle) is required only for
    the Jy/beam form.
    """
    scale = coerce_temperature_scale(scale)
    if scale is BrightnessTemperatureScale.RAYLEIGH_JEANS:
        return u.brightness_temperature(frequency, beam_area=beam_area)
    return _planck_brightness_temperature(frequency, beam_area=beam_area)


__all__ = [
    "spectral",
    "beam_solid_angle",
    "beam_angular_area",
    "doppler",
    "brightness_temperature",
]
