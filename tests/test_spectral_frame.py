"""Tests for issue #24 — the frequency-authoritative spectral axis + schema extensions.

Written first (red/green TDD). These tests *are* the correctness obligation for:

- the three new :class:`~astrolyze.io.Metadata` fields — ``spectral_frame`` (read from FITS
  ``SPECSYS``, ``None`` when absent — never guessed), ``systemic_velocity`` (optional
  ``Quantity``), and ``lines`` (an extensible list of ``(species, transition,
  rest_frequency)`` entries: one for a single-line cube, many for broadband);
- all three round-tripping through **both** projections — ``to_header``/``from_header``
  (HIERARCH ASTROLYZE cards) **and** ``to_attrs``/``from_attrs`` (#22);
- ``rest_frequency`` staying the optional *primary* line, so a single-line cube agrees with
  its single ``lines`` entry and existing single-line headers parse exactly as before;
- the spectral axis exposed/stored as absolute frequency (authoritative) with a derived
  per-line Δv (integrating #26's ``AxisCoordinates``);
- a frame transform (LSRK<->barycentric<->topocentric) being **context-or-raise**: it
  delegates to astropy ``SpectralCoord`` and raises without observer location + obstime +
  target, never silently shifting the frame (ADR-0003).

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, on the radio velocity convention, LSRK frame.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time

from astrolyze.core import Cube
from astrolyze.io import Line, Metadata, load
from astrolyze.units import MissingContextError, VelocityConvention

warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
REST_13CO = 220.39868 * u.GHz  # 13CO(2-1), the second broadband line


# --------------------------------------------------------------------------------------
# Synthetic FITS fixtures (a tiny PPV cube whose header carries the schema)
# --------------------------------------------------------------------------------------
def _cube_header(specsys: str | None = "LSRK"):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    if specsys is not None:
        h["SPECSYS"] = specsys
    return h


@pytest.fixture
def cube(tmp_path):
    path = tmp_path / "ngc0628_co21.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header())
    return Cube.from_loaded(load(path))


# ======================================================================================
# AC: Metadata gains spectral_frame (from SPECSYS, None if absent)
# ======================================================================================
def test_spectral_frame_read_from_specsys(tmp_path):
    path = tmp_path / "lsrk.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header("LSRK"))
    assert load(path).metadata.spectral_frame == "LSRK"


def test_spectral_frame_is_none_when_specsys_absent(tmp_path):
    # SPECSYS is read, never guessed: absent -> None (ADR-0003), even on a velocity axis.
    path = tmp_path / "no_specsys.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header(specsys=None))
    assert load(path).metadata.spectral_frame is None


@pytest.mark.parametrize("specsys", ["LSRK", "BARYCENT", "TOPOCENT", "HELIOCEN"])
def test_spectral_frame_codes_pass_through_verbatim(tmp_path, specsys):
    path = tmp_path / f"{specsys}.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header(specsys))
    assert load(path).metadata.spectral_frame == specsys


# ======================================================================================
# AC: systemic_velocity — optional Quantity
# ======================================================================================
def test_systemic_velocity_defaults_to_none():
    assert Metadata().systemic_velocity is None


def test_systemic_velocity_holds_a_quantity():
    m = Metadata(systemic_velocity=657.0 * u.km / u.s)
    assert u.isclose(m.systemic_velocity, 657.0 * u.km / u.s, rtol=1e-12)


# ======================================================================================
# AC: lines — extensible list; rest_frequency stays the optional PRIMARY line
# ======================================================================================
def test_lines_defaults_to_empty_list():
    assert Metadata().lines == []


def test_single_line_cube_primary_agrees_with_its_lines_entry():
    # A single-line cube: rest_frequency is the primary line; one lines entry agrees.
    line = Line(species="CO", transition="2-1", rest_frequency=REST)
    m = Metadata(rest_frequency=REST, lines=[line])
    assert len(m.lines) == 1
    assert u.isclose(m.lines[0].rest_frequency, m.rest_frequency, rtol=1e-12)


def test_multiple_lines_broadband():
    lines = [
        Line(species="CO", transition="2-1", rest_frequency=REST),
        Line(species="13CO", transition="2-1", rest_frequency=REST_13CO),
    ]
    m = Metadata(rest_frequency=REST, lines=lines)
    assert [line.species for line in m.lines] == ["CO", "13CO"]
    assert u.isclose(m.lines[1].rest_frequency, REST_13CO, rtol=1e-12)


# ======================================================================================
# AC: all three round-trip through to_header/from_header AND to_attrs/from_attrs
# ======================================================================================
def _rich_metadata():
    return Metadata(
        object="NGC0628",
        telescope="ALMA",
        species="CO21",
        rest_frequency=REST,
        velocity_convention=VelocityConvention.RADIO,
        spectral_frame="LSRK",
        systemic_velocity=657.0 * u.km / u.s,
        lines=[
            Line(species="CO", transition="2-1", rest_frequency=REST),
            Line(species="13CO", transition="2-1", rest_frequency=REST_13CO),
        ],
    )


def _assert_new_fields_equal(m0: Metadata, m1: Metadata):
    assert m1.spectral_frame == m0.spectral_frame
    assert u.isclose(m1.systemic_velocity, m0.systemic_velocity, rtol=1e-12)
    assert [line.species for line in m1.lines] == [line.species for line in m0.lines]
    assert [line.transition for line in m1.lines] == [
        line.transition for line in m0.lines
    ]
    for la, lb in zip(m0.lines, m1.lines):
        assert u.isclose(la.rest_frequency, lb.rest_frequency, rtol=1e-12)


def test_new_fields_roundtrip_through_header():
    m0 = _rich_metadata()
    m1 = Metadata.from_header(m0.to_header())
    _assert_new_fields_equal(m0, m1)


def test_new_fields_roundtrip_through_attrs():
    import json

    m0 = _rich_metadata()
    attrs = m0.to_attrs()
    json.dumps(attrs)  # the attrs projection must be JSON-serializable
    _assert_new_fields_equal(m0, Metadata.from_attrs(attrs))


def test_spectral_frame_written_as_standard_specsys_card():
    # spectral_frame round-trips via the standard FITS SPECSYS keyword (not a HIERARCH card),
    # so a plain FITS reader sees it.
    header = _rich_metadata().to_header()
    assert header["SPECSYS"] == "LSRK"


def test_single_line_header_parses_exactly_as_before(tmp_path):
    # Backward compatibility: a single-line header with only rest_frequency (no lines /
    # systemic_velocity cards) parses, and the primary line is preserved.
    path = tmp_path / "single.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header())
    m = load(path).metadata
    assert u.isclose(m.rest_frequency, REST, rtol=1e-12)
    assert m.systemic_velocity is None
    # An absent lines list does not invent broadband lines.
    assert m.lines == []


# ======================================================================================
# AC: spectral axis exposed/stored as absolute frequency + a derived per-line Δv
# ======================================================================================
def test_cube_exposes_absolute_frequency_and_per_line_delta_v(cube):
    coords = cube.coordinates
    assert coords.frequency.unit.is_equivalent(u.Hz)
    assert u.isclose(coords.frequency[0], REST, rtol=1e-9)
    assert coords.delta_v.unit.is_equivalent(u.km / u.s)


def test_per_line_delta_v_for_a_named_line(cube):
    # The derived Δv of every channel relative to a *chosen* line (here 13CO) — broadband
    # work needs Δv per line, not only relative to the primary rest frequency.
    dv = cube.delta_v_for(REST_13CO)
    assert dv.unit.is_equivalent(u.km / u.s)
    nu = cube.coordinates.frequency
    expected = nu.to(u.km / u.s, equivalencies=u.doppler_radio(REST_13CO))
    assert u.allclose(dv, expected, rtol=1e-9, atol=1e-6 * u.km / u.s)


# ======================================================================================
# AC: a frame transform without observer location + obstime + target RAISES
#     (delegates to SpectralCoord), never silently shifts frame
# ======================================================================================
def test_frame_transform_without_context_raises(cube):
    # No observer location / obstime / target supplied -> the transform must refuse (it would
    # otherwise silently shift the frame, ADR-0003). Delegates to astropy SpectralCoord.
    with pytest.raises((MissingContextError, ValueError, u.UnitsError)):
        cube.to_spectral_frame("BARYCENT")


def test_frame_transform_with_full_context_succeeds(cube):
    # With observer location + obstime + target, SpectralCoord does the maths and the cube's
    # spectral_frame is updated to the new frame.
    location = EarthLocation(lat=-23.0 * u.deg, lon=-67.75 * u.deg, height=5000 * u.m)
    obstime = Time("2020-01-01T00:00:00")
    target = SkyCoord(ra=24.174 * u.deg, dec=15.784 * u.deg, frame="icrs")
    shifted = cube.to_spectral_frame(
        "BARYCENT", location=location, obstime=obstime, target=target
    )
    assert shifted.metadata.spectral_frame == "BARYCENT"
    # The absolute frequency really moved (a barycentric shift is non-zero here).
    assert not u.allclose(
        shifted.coordinates.frequency, cube.coordinates.frequency, rtol=1e-12
    )


def test_frame_transform_when_frame_absent_raises_naming_spectral_frame(tmp_path):
    # Absent SPECSYS -> spectral_frame is None, flagged ONLY when a transform needs it: a
    # transform from an unknown source frame raises naming spectral_frame (lazy enforcement).
    path = tmp_path / "no_frame.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header(specsys=None))
    cube = Cube.from_loaded(load(path))
    assert cube.metadata.spectral_frame is None
    location = EarthLocation(lat=-23.0 * u.deg, lon=-67.75 * u.deg, height=5000 * u.m)
    with pytest.raises(MissingContextError, match="spectral_frame"):
        cube.to_spectral_frame(
            "BARYCENT",
            location=location,
            obstime=Time("2020-01-01T00:00:00"),
            target=SkyCoord(ra=24.174 * u.deg, dec=15.784 * u.deg, frame="icrs"),
        )
