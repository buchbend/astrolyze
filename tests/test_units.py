"""Tests for astrolyze.units — the deep module (ADR-0003).

Written first (red/green TDD). These tests *are* the correctness obligation for radio
unit handling: round-trips across the unit zoo, one hand-checked Planck (non-RJ) value
computed independently from astropy constants, and the no-silent-physics contract
(missing velocity convention / rest frequency / temperature scale must raise, never
default).

The rest frequency used throughout is CO J=2-1 (230.538 GHz); the beam is a 10" circular
Gaussian — both realistic for the tracer-bullet PHANGS cube.
"""

import numpy as np
import pytest
import astropy.units as u
from astropy import constants as const

from astrolyze import units
from astrolyze.units import (
    BrightnessTemperatureScale,
    MissingContextError,
    VelocityConvention,
    convert,
)

REST = 230.538 * u.GHz  # CO(2-1)
BEAM_SR = (1.1330900354567984 * (10 * u.arcsec) ** 2).to(
    u.sr
)  # 10" Gaussian, pi/(4 ln2)


# --------------------------------------------------------------------------------------
# Aliases — named radio unit objects (plain astropy units; not globally registered)
# --------------------------------------------------------------------------------------
def test_aliases_are_plain_astropy_units():
    assert units.Jy_beam == u.Jy / u.beam
    assert units.MJy_sr == u.MJy / u.sr
    assert units.Jy_sr == u.Jy / u.sr
    assert units.K_kms == u.K * u.km / u.s
    # Tmb / Ta are the K (Rayleigh-Jeans) brightness-temperature scale, labelled for intent.
    assert units.Tmb == u.K
    assert units.Ta == u.K


def test_import_does_not_register_global_unit_names():
    # ADR-0003 / house rule: no global side effects. We never inject "Tmb" etc. into
    # astropy's global unit registry, so parsing the bare string must still fail.
    with pytest.raises(ValueError):
        u.Unit("Tmb")


# --------------------------------------------------------------------------------------
# Round-trips across the unit zoo (Rayleigh-Jeans = the standard radio intensity scale)
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "scale",
    [BrightnessTemperatureScale.RAYLEIGH_JEANS, BrightnessTemperatureScale.PLANCK],
)
def test_roundtrip_K_to_MJy_sr(scale):
    start = 12.0 * u.K
    out = convert(start, u.MJy / u.sr, rest_frequency=REST, temperature_scale=scale)
    assert out.unit.is_equivalent(u.MJy / u.sr)
    back = convert(out, u.K, rest_frequency=REST, temperature_scale=scale)
    assert u.isclose(back, start, rtol=1e-10)


@pytest.mark.parametrize(
    "scale",
    [BrightnessTemperatureScale.RAYLEIGH_JEANS, BrightnessTemperatureScale.PLANCK],
)
def test_roundtrip_K_to_Jy_beam(scale):
    start = 12.0 * u.K
    out = convert(
        start, u.Jy / u.beam, rest_frequency=REST, beam=BEAM_SR, temperature_scale=scale
    )
    assert out.unit.is_equivalent(u.Jy / u.beam)
    back = convert(out, u.K, rest_frequency=REST, beam=BEAM_SR, temperature_scale=scale)
    assert u.isclose(back, start, rtol=1e-10)


def test_roundtrip_Jy_beam_to_MJy_sr_is_pure_geometry():
    # Jy/beam <-> MJy/sr is a beam-solid-angle conversion only: no frequency, no
    # temperature scale required (and supplying them would be meaningless).
    start = 0.5 * u.Jy / u.beam
    out = convert(start, u.MJy / u.sr, beam=BEAM_SR)
    assert out.unit.is_equivalent(u.MJy / u.sr)
    back = convert(out, u.Jy / u.beam, beam=BEAM_SR)
    assert u.isclose(back, start, rtol=1e-12)


def test_roundtrip_full_chain_K_Jybeam_MJysr_K():
    start = 8.0 * u.K
    scale = BrightnessTemperatureScale.RAYLEIGH_JEANS
    a = convert(
        start, u.Jy / u.beam, rest_frequency=REST, beam=BEAM_SR, temperature_scale=scale
    )
    b = convert(a, u.MJy / u.sr, beam=BEAM_SR)
    c = convert(b, u.K, rest_frequency=REST, temperature_scale=scale)
    assert u.isclose(c, start, rtol=1e-10)


def test_integrated_intensity_K_kms_path():
    # Velocity-integrated intensity rides the velocity axis through a *linear* (RJ)
    # brightness-temperature conversion: K km/s <-> Jy/beam km/s.
    start = 5.0 * u.K * u.km / u.s
    scale = BrightnessTemperatureScale.RAYLEIGH_JEANS
    out = convert(
        start,
        u.Jy / u.beam * u.km / u.s,
        rest_frequency=REST,
        beam=BEAM_SR,
        temperature_scale=scale,
    )
    assert out.unit.is_equivalent(u.Jy / u.beam * u.km / u.s)
    back = convert(
        out,
        u.K * u.km / u.s,
        rest_frequency=REST,
        beam=BEAM_SR,
        temperature_scale=scale,
    )
    assert u.isclose(back, start, rtol=1e-10)

    # The integrated conversion must equal the per-channel factor times the velocity width
    # (proof that it is linear / commutes with the velocity integral).
    factor = convert(
        1.0 * u.K,
        u.Jy / u.beam,
        rest_frequency=REST,
        beam=BEAM_SR,
        temperature_scale=scale,
    )
    assert u.isclose(out, 5.0 * factor.value * u.Jy / u.beam * u.km / u.s, rtol=1e-10)


def test_integrated_geometric_path_needs_no_frequency():
    start = 2.0 * u.Jy / u.beam * u.km / u.s
    out = convert(start, u.MJy / u.sr * u.km / u.s, beam=BEAM_SR)
    back = convert(out, u.Jy / u.beam * u.km / u.s, beam=BEAM_SR)
    assert u.isclose(back, start, rtol=1e-12)


def test_beam_accepts_radio_beam_object():
    radio_beam = pytest.importorskip("radio_beam")
    beam = radio_beam.Beam(10 * u.arcsec)
    out_obj = convert(0.5 * u.Jy / u.beam, u.MJy / u.sr, beam=beam)
    out_sr = convert(0.5 * u.Jy / u.beam, u.MJy / u.sr, beam=beam.sr)
    assert u.isclose(out_obj, out_sr, rtol=1e-12)


# --------------------------------------------------------------------------------------
# Hand-checked Planck (non-RJ) value vs an independent calculation
# --------------------------------------------------------------------------------------
def test_planck_brightness_temperature_matches_independent_calc():
    T = 30.0 * u.K
    nu = REST.to(u.Hz, equivalencies=u.spectral())

    # Independent forward Planck law B_nu(T) = (2 h nu^3 / c^2) / (exp(h nu / k T) - 1).
    # The per-steradian of specific intensity is conventional (sr is dimensionless-valued),
    # so the Jy value is the Jy/sr numeric value.
    x = (const.h * nu / (const.k_B * T)).to_value(u.dimensionless_unscaled)
    B = (2 * const.h * nu**3 / const.c**2) / np.expm1(x)
    expected_MJy_sr = B.to_value(u.Jy) / 1e6

    out = convert(
        T,
        u.MJy / u.sr,
        rest_frequency=REST,
        temperature_scale=BrightnessTemperatureScale.PLANCK,
    )
    assert u.isclose(out, expected_MJy_sr * u.MJy / u.sr, rtol=1e-9)

    # Sanity: the hand-computed value (pinned literal) is what we expect.
    assert out.to_value(u.MJy / u.sr) == pytest.approx(40507.477, rel=1e-5)


def test_planck_differs_from_rayleigh_jeans_in_the_wien_regime():
    # At 230 GHz a 30 K source is past the RJ limit (h nu / k = 11 K), so the two scales
    # must disagree well outside round-off. This is the silent-error trap ADR-0003 targets.
    T = 30.0 * u.K
    nu = REST.to(u.Hz, equivalencies=u.spectral())
    x = (const.h * nu / (const.k_B * T)).to_value(u.dimensionless_unscaled)

    B_rj = 2 * const.k_B * T * nu**2 / const.c**2  # Rayleigh-Jeans intensity
    expected_rj = B_rj.to_value(u.Jy) / 1e6

    out_rj = convert(
        T,
        u.MJy / u.sr,
        rest_frequency=REST,
        temperature_scale=BrightnessTemperatureScale.RAYLEIGH_JEANS,
    )
    out_planck = convert(
        T,
        u.MJy / u.sr,
        rest_frequency=REST,
        temperature_scale=BrightnessTemperatureScale.PLANCK,
    )

    assert u.isclose(out_rj, expected_rj * u.MJy / u.sr, rtol=1e-9)
    # For T -> intensity the RJ line overestimates the Planck intensity in the Wien tail:
    # B_planck / B_rj = x / (exp(x) - 1) < 1. Demand the implementation match that exactly
    # (not merely "close to RJ"), and that the gap is clearly non-trivial (~17% here).
    assert u.isclose(out_planck / out_rj, x / np.expm1(x), rtol=1e-9)
    assert out_planck < out_rj * 0.9


# --------------------------------------------------------------------------------------
# No silent physics — missing required context raises (per case)
# --------------------------------------------------------------------------------------
def test_temperature_conversion_requires_rest_frequency():
    with pytest.raises(MissingContextError, match="rest_frequency"):
        convert(
            10 * u.K, u.MJy / u.sr, temperature_scale=BrightnessTemperatureScale.PLANCK
        )


def test_temperature_conversion_requires_temperature_scale():
    with pytest.raises(MissingContextError, match="temperature_scale"):
        convert(10 * u.K, u.MJy / u.sr, rest_frequency=REST)


def test_per_beam_temperature_conversion_requires_beam():
    with pytest.raises(MissingContextError, match="beam"):
        convert(
            10 * u.K,
            u.Jy / u.beam,
            rest_frequency=REST,
            temperature_scale=BrightnessTemperatureScale.RAYLEIGH_JEANS,
        )


def test_geometric_per_beam_conversion_requires_beam():
    with pytest.raises(MissingContextError, match="beam"):
        convert(0.5 * u.Jy / u.beam, u.MJy / u.sr)


def test_velocity_axis_requires_convention_and_rest_frequency():
    with pytest.raises(MissingContextError, match="convention"):
        convert(230.0 * u.GHz, u.km / u.s, rest_frequency=REST)
    with pytest.raises(MissingContextError, match="rest_frequency"):
        convert(230.0 * u.GHz, u.km / u.s, convention=VelocityConvention.RADIO)
    with pytest.raises(MissingContextError):
        convert(230.0 * u.GHz, u.km / u.s)


def test_planck_on_integrated_intensity_raises():
    # Planck B_nu(T) is nonlinear, so it cannot be applied to a velocity-integrated
    # intensity (the integral does not commute). Integrated intensities are RJ by
    # construction; asking for Planck here is a physics error, not a default to paper over.
    with pytest.raises(ValueError, match="(?i)planck.*integrated|integrated.*planck"):
        convert(
            5 * u.K * u.km / u.s,
            u.Jy / u.beam * u.km / u.s,
            rest_frequency=REST,
            beam=BEAM_SR,
            temperature_scale=BrightnessTemperatureScale.PLANCK,
        )


# --------------------------------------------------------------------------------------
# Velocity convention actually matters (radio != optical != relativistic)
# --------------------------------------------------------------------------------------
def test_velocity_conventions_are_distinct_and_match_astropy():
    sky = 230.0 * u.GHz
    expected = {
        VelocityConvention.RADIO: sky.to(
            u.km / u.s, equivalencies=u.doppler_radio(REST)
        ),
        VelocityConvention.OPTICAL: sky.to(
            u.km / u.s, equivalencies=u.doppler_optical(REST)
        ),
        VelocityConvention.RELATIVISTIC: sky.to(
            u.km / u.s, equivalencies=u.doppler_relativistic(REST)
        ),
    }
    got = {
        c: convert(sky, u.km / u.s, rest_frequency=REST, convention=c) for c in expected
    }
    for c in expected:
        assert u.isclose(got[c], expected[c], rtol=1e-12)
    # The three conventions must give genuinely different velocities.
    vals = [got[c].to_value(u.km / u.s) for c in expected]
    assert len({round(v, 3) for v in vals}) == 3


def test_velocity_roundtrip():
    v = convert(
        220.0 * u.GHz,
        u.km / u.s,
        rest_frequency=REST,
        convention=VelocityConvention.RADIO,
    )
    back = convert(v, u.GHz, rest_frequency=REST, convention=VelocityConvention.RADIO)
    assert u.isclose(back, 220.0 * u.GHz, rtol=1e-12)


# --------------------------------------------------------------------------------------
# Spectral-axis conversions that need no context
# --------------------------------------------------------------------------------------
def test_spectral_frequency_wavelength_needs_no_context():
    out = convert(REST, u.mm)
    assert out.unit.is_equivalent(u.mm)
    assert u.isclose(out, REST.to(u.mm, equivalencies=u.spectral()), rtol=1e-12)


# --------------------------------------------------------------------------------------
# String coercion (agent-native ergonomics) and result typing
# --------------------------------------------------------------------------------------
def test_string_context_is_accepted():
    out_enum = convert(
        10 * u.K,
        u.MJy / u.sr,
        rest_frequency=REST,
        temperature_scale=BrightnessTemperatureScale.RAYLEIGH_JEANS,
    )
    out_str = convert(
        10 * u.K, "MJy / sr", rest_frequency=REST, temperature_scale="rayleigh_jeans"
    )
    assert u.isclose(out_enum, out_str, rtol=1e-12)
    v_enum = convert(
        230.0 * u.GHz,
        u.km / u.s,
        rest_frequency=REST,
        convention=VelocityConvention.RADIO,
    )
    v_str = convert(230.0 * u.GHz, u.km / u.s, rest_frequency=REST, convention="radio")
    assert u.isclose(v_enum, v_str, rtol=1e-12)


def test_unknown_convention_string_raises():
    with pytest.raises(ValueError):
        convert(230.0 * u.GHz, u.km / u.s, rest_frequency=REST, convention="lsr")
    with pytest.raises(ValueError):
        convert(10 * u.K, u.MJy / u.sr, rest_frequency=REST, temperature_scale="bogus")


def test_convert_returns_quantity():
    out = convert(
        10 * u.K, u.MJy / u.sr, rest_frequency=REST, temperature_scale="rayleigh_jeans"
    )
    assert isinstance(out, u.Quantity)
