"""Tests for astrolyze.core — the Cube / Map / Spectrum context-carrying wrappers (ADR-0004).

Written first (red/green TDD). These tests *are* the correctness obligation for the object
model contract:

- the wrappers **compose** spectral-cube / astropy / specutils (a private upstream handle)
  rather than subclass them, and add the toolkit's value: object-carried context
  (beam + rest frequency + velocity convention, ADR-0003) sourced from ``io`` Metadata;
- type transitions carry context for free: ``Cube.moment0() -> Map`` keeps the beam and rest
  frequency; ``cube[:, y, x] -> Spectrum`` keeps the same context;
- ``.to(unit)`` is the unit hub (ADR-0003c) — it supplies the object's beam / rest frequency /
  convention so the caller need not, while still **never** defaulting the genuinely-ambiguous
  Rayleigh-Jeans-vs-Planck scale (that stays explicit);
- an operation that needs context the object does not have **raises** (lazy enforcement,
  ADR-0006 ii) rather than silently guessing.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention — here in Jy/beam
so the unit-hub conversions exercise real beam + rest-frequency context.
"""

import re
import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, Map, Spectrum
from astrolyze.io import Metadata, load
from astrolyze.units import MissingContextError, VelocityConvention, convert

# spectral-cube emits cosmetic warnings (e.g. no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


# --------------------------------------------------------------------------------------
# Synthetic FITS fixtures (a tiny PPV cube whose header carries the schema)
# --------------------------------------------------------------------------------------
def _cube_header():
    """A 3D Jy/beam cube on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "Jy/beam"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


@pytest.fixture
def cube(tmp_path):
    """A complete :class:`Cube` built through the real io.load seam (ADR-0004/0006)."""
    path = tmp_path / "ngc0628_co21.fits"
    data = np.ones((4, 3, 3), dtype="float32")
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture
def integrated_map(cube):
    """A velocity-integrated map (Jy/beam km/s), the moment0 of the reference cube."""
    return cube.moment0()


# --------------------------------------------------------------------------------------
# Composition, not subclassing (ADR-0004 rider 1)
# --------------------------------------------------------------------------------------
def test_wrappers_compose_and_do_not_subclass_upstream(cube):
    from spectral_cube import SpectralCube

    # holds a private handle and delegates; is NOT a SpectralCube itself.
    assert isinstance(cube, Cube)
    assert not isinstance(cube, SpectralCube)
    assert isinstance(cube._sc, SpectralCube)
    # context is sourced from the io Metadata (the single source of truth).
    assert isinstance(cube.metadata, Metadata)
    assert u.isclose(cube.rest_frequency, REST, rtol=1e-12)
    assert cube.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(cube.beam.major, 12 * u.arcsec, rtol=1e-9)


# --------------------------------------------------------------------------------------
# MANDATED: cube.moment0() -> Map carrying beam + rest frequency
# --------------------------------------------------------------------------------------
def test_moment0_returns_map_carrying_beam_and_frequency(cube):
    m = cube.moment0()
    assert isinstance(m, Map)
    # the 2D moment carries the cube's physical context for free.
    assert u.isclose(m.beam.major, cube.beam.major, rtol=1e-12)
    assert u.isclose(m.beam.minor, cube.beam.minor, rtol=1e-12)
    assert u.isclose(m.beam.pa, cube.beam.pa, rtol=1e-12)
    assert u.isclose(m.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert m.velocity_convention is cube.velocity_convention
    # moment0 of a Jy/beam cube over a velocity axis is a velocity-integrated intensity.
    assert m.shape == (3, 3)
    assert m.unit.is_equivalent(u.Jy / u.beam * u.km / u.s)


# --------------------------------------------------------------------------------------
# MANDATED: cube[:, y, x] -> Spectrum (carrying context)
# --------------------------------------------------------------------------------------
def test_spatial_index_returns_spectrum_carrying_context(cube):
    sp = cube[:, 1, 1]
    assert isinstance(sp, Spectrum)
    # the spectral axis comes through from the cube WCS.
    assert sp.flux.shape == (4,)
    assert sp.flux.unit.is_equivalent(u.Jy / u.beam)
    assert sp.spectral_axis.unit.is_equivalent(u.m / u.s)
    # same physical context as the parent cube.
    assert u.isclose(sp.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert sp.velocity_convention is cube.velocity_convention
    assert u.isclose(sp.beam.major, cube.beam.major, rtol=1e-12)


def test_spectrum_wraps_specutils(cube):
    from specutils import Spectrum as SpecutilsSpectrum

    sp = cube[:, 0, 0]
    assert isinstance(
        sp._spec, SpecutilsSpectrum
    )  # thin specutils adaptor (ADR-0004 rider 2)


# --------------------------------------------------------------------------------------
# MANDATED: map.to("K km/s") uses the object's context (beam + rest frequency)
# --------------------------------------------------------------------------------------
def test_map_to_uses_object_context(integrated_map):
    # The object supplies beam + rest frequency + convention; the caller supplies only the
    # Rayleigh-Jeans-vs-Planck scale, which astrolyze NEVER defaults (ADR-0003). For a
    # velocity-integrated intensity RJ is the only valid scale (Planck is non-linear).
    out = integrated_map.to("K km/s", temperature_scale="rayleigh_jeans")
    assert isinstance(out, Map)
    assert out.unit.is_equivalent(u.K * u.km / u.s)
    # identical to calling the unit layer directly with the object's context spelled out:
    # this is exactly the context the wrapper supplied on the caller's behalf.
    expected = convert(
        integrated_map.data,
        "K km/s",
        rest_frequency=REST,
        beam=BEAM,
        temperature_scale="rayleigh_jeans",
    )
    assert u.allclose(out.data, expected, rtol=1e-10)
    # context is preserved across the conversion.
    assert u.isclose(out.rest_frequency, REST, rtol=1e-12)
    assert u.isclose(out.beam.major, BEAM.major, rtol=1e-12)


def test_to_returns_same_wrapper_type_with_updated_unit(integrated_map):
    out = integrated_map.to("K km/s", temperature_scale="rayleigh_jeans")
    assert type(out) is Map
    assert out.metadata.bunit == out.unit


# --------------------------------------------------------------------------------------
# MANDATED: an operation needing absent context raises (no silent guess, ADR-0003/0006)
# --------------------------------------------------------------------------------------
def test_to_on_incomplete_object_raises():
    # A real-archival map with a beam but no rest frequency (lazy: holding it is fine).
    incomplete = Map(
        np.ones((3, 3)) * (u.Jy / u.beam),
        wcs=None,
        metadata=Metadata(bunit=u.Jy / u.beam, beam=BEAM),
    )
    assert incomplete.is_complete is False
    assert "rest_frequency" in incomplete.missing
    # Jy/beam -> K needs the rest frequency the object does not carry: it must raise.
    with pytest.raises(MissingContextError, match="rest_frequency"):
        incomplete.to("K", temperature_scale="rayleigh_jeans")


def test_to_never_defaults_the_brightness_scale(cube):
    # A per-channel Jy/beam -> K conversion is genuinely ambiguous (RJ vs Planck); the
    # object cannot know which, so .to() must refuse rather than pick one silently.
    channel = cube[0]  # a 2D channel map, still Jy/beam
    assert isinstance(channel, Map)
    with pytest.raises(MissingContextError, match="temperature_scale"):
        channel.to("K")
    # with the scale stated explicitly it succeeds (and both scales are reachable).
    assert channel.to("K", temperature_scale="planck").unit.is_equivalent(u.K)
    assert channel.to("K", temperature_scale="rayleigh_jeans").unit.is_equivalent(u.K)


# --------------------------------------------------------------------------------------
# Type transitions on slicing (Cube -> Cube / Map / Spectrum), context carried through
# --------------------------------------------------------------------------------------
def test_subcube_slice_returns_cube(cube):
    sub = cube[:, 0:2, 0:2]
    assert isinstance(sub, Cube)
    assert sub.shape == (4, 2, 2)
    assert u.isclose(sub.rest_frequency, cube.rest_frequency, rtol=1e-12)


def test_channel_slice_returns_map(cube):
    channel = cube[0]
    assert isinstance(channel, Map)
    assert channel.shape == (3, 3)
    assert channel.unit.is_equivalent(u.Jy / u.beam)
    assert u.isclose(channel.beam.major, cube.beam.major, rtol=1e-12)


# --------------------------------------------------------------------------------------
# .plot() is a thin seam onto the viz engine (ADR-0005); it arrives with issue #6
# --------------------------------------------------------------------------------------
def test_plot_is_a_seam_onto_the_unbuilt_viz_layer(integrated_map):
    # The sugar delegates to astrolyze.viz; until #6 ships that engine, calling it is a
    # clear NotImplementedError, not an obscure AttributeError.
    with pytest.raises(NotImplementedError):
        integrated_map.plot()


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code (ADR-0005 / the legacy sin)
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_library_code():
    from pathlib import Path

    import astrolyze.core as core_pkg

    pkg_dir = Path(core_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
