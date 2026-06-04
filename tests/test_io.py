"""Tests for astrolyze.io — the metadata schema, lazy loader, and filename projection
(ADR-0006).

Written first (red/green TDD). These tests *are* the correctness obligation for the I/O
contract: the FITS header is authoritative; loading is lazy (an incomplete header opens but
is flagged, and asking for the missing context raises — never a silent guess, ADR-0003);
load -> save -> load preserves every schema field; and the filename is a derived projection
of the header (legacy astrolyze-1 grammar), regenerated, never hand-maintained.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import os

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.io import (
    CalibrationScale,
    LoadedData,
    Metadata,
    MissingContextError,
    VelocityConvention,
    load,
    project,
    save,
)

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


# --------------------------------------------------------------------------------------
# Synthetic FITS fixtures (tiny cubes; the header carries the full / partial schema)
# --------------------------------------------------------------------------------------
def _spatial_wcs(header):
    """Add a minimal 2D celestial WCS so astropy.wcs.WCS builds cleanly."""
    header["CTYPE1"], header["CRVAL1"] = "RA---SIN", 24.174
    header["CDELT1"], header["CRPIX1"], header["CUNIT1"] = -2e-4, 1.0, "deg"
    header["CTYPE2"], header["CRVAL2"] = "DEC--SIN", 15.784
    header["CDELT2"], header["CRPIX2"], header["CUNIT2"] = 2e-4, 1.0, "deg"


def _full_header():
    """A FITS header carrying the complete astrolyze metadata schema (a 3D cube)."""
    h = fits.Header()
    _spatial_wcs(h)
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    h["HIERARCH ASTROLYZE DISTANCE"] = 9.84
    h["HIERARCH ASTROLYZE DISTUNIT"] = "Mpc"
    h["HIERARCH ASTROLYZE CALERR"] = 0.1
    h["HIERARCH ASTROLYZE NAMETAG"] = "mom0"
    h["HIERARCH ASTROLYZE CALSCALE"] = "T_mb"
    h["HIERARCH ASTROLYZE ETAMB"] = 0.8
    return h


def _incomplete_header():
    """A 2D map header with an object but *no* rest frequency or velocity convention —
    the real-archival case ADR-0006 (ii) demands we open anyway."""
    h = fits.Header()
    _spatial_wcs(h)
    h["OBJECT"] = "NGC0628"
    h["BUNIT"] = "K"
    return h


@pytest.fixture
def full_fits(tmp_path):
    path = tmp_path / "raw" / "ngc0628_co21.fits"
    path.parent.mkdir(parents=True)
    data = np.arange(2 * 3 * 3, dtype="float32").reshape(2, 3, 3)
    fits.writeto(path, data, _full_header())
    return path


@pytest.fixture
def incomplete_fits(tmp_path):
    path = tmp_path / "raw" / "archival_map.fits"
    path.parent.mkdir(parents=True)
    data = np.arange(3 * 3, dtype="float32").reshape(3, 3)
    fits.writeto(path, data, _incomplete_header())
    return path


# --------------------------------------------------------------------------------------
# Parse the schema from a FITS header (the header is authoritative)
# --------------------------------------------------------------------------------------
def test_parse_required_keywords_from_synthetic_fits(full_fits):
    loaded = load(full_fits)
    assert isinstance(loaded, LoadedData)
    assert loaded.data.shape == (2, 3, 3)

    m = loaded.metadata
    assert m.object == "NGC0628"
    assert m.telescope == "ALMA"
    assert m.species == "CO21"
    assert u.isclose(m.rest_frequency, REST, rtol=1e-12)
    assert m.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(m.beam.major, 12 * u.arcsec, rtol=1e-9)
    assert u.isclose(m.beam.minor, 10 * u.arcsec, rtol=1e-9)
    assert u.isclose(m.beam.pa, 30 * u.deg, rtol=1e-9)
    assert m.bunit == u.K
    assert u.isclose(m.distance, 9.84 * u.Mpc, rtol=1e-9)
    assert m.calibration_error == pytest.approx(0.1)
    assert m.name_tag == "mom0"
    # The calibration TEMPERATURE SCALE + main-beam efficiency (issue #25).
    assert m.calibration_scale is CalibrationScale.T_MB
    assert m.eta_mb == pytest.approx(0.8)

    assert m.is_complete is True
    assert m.missing == []


# --------------------------------------------------------------------------------------
# Lazy enforcement — incomplete header loads, is flagged, and raises only on use
# --------------------------------------------------------------------------------------
def test_incomplete_header_loads_but_is_flagged(incomplete_fits):
    loaded = load(incomplete_fits)
    # The file opens and the data is there (never refuse a real file).
    assert loaded.data.shape == (3, 3)
    assert loaded.metadata.object == "NGC0628"

    m = loaded.metadata
    assert m.is_complete is False
    assert set(m.missing) == {"rest_frequency", "velocity_convention"}


def test_ensure_complete_raises_naming_the_missing_context(incomplete_fits):
    m = load(incomplete_fits).metadata
    with pytest.raises(MissingContextError) as exc:
        m.ensure_complete()
    msg = str(exc.value)
    assert "rest_frequency" in msg
    assert "velocity_convention" in msg


def test_complete_metadata_ensure_complete_is_noop(full_fits):
    # No exception, returns nothing useful — just must not raise.
    assert load(full_fits).metadata.ensure_complete() is None


# --------------------------------------------------------------------------------------
# Round-trip: load -> save -> load preserves every schema field (a test obligation)
# --------------------------------------------------------------------------------------
def test_roundtrip_load_save_load_preserves_all_fields(full_fits, tmp_path):
    loaded = load(full_fits)
    out_dir = tmp_path / "derived"
    written = save(loaded.data, loaded.metadata, out_dir, base_header=loaded.header)
    assert written.exists()

    m0, m1 = loaded.metadata, load(written).metadata
    assert m1.object == m0.object
    assert m1.telescope == m0.telescope
    assert m1.species == m0.species
    assert u.isclose(m1.rest_frequency, m0.rest_frequency, rtol=1e-12)
    assert m1.velocity_convention is m0.velocity_convention
    assert u.isclose(m1.beam.major, m0.beam.major, rtol=1e-12)
    assert u.isclose(m1.beam.minor, m0.beam.minor, rtol=1e-12)
    assert u.isclose(m1.beam.pa, m0.beam.pa, rtol=1e-12)
    assert m1.bunit == m0.bunit
    assert u.isclose(m1.distance, m0.distance, rtol=1e-12)
    assert m1.calibration_error == pytest.approx(m0.calibration_error)
    assert m1.name_tag == m0.name_tag
    assert m1.calibration_scale is m0.calibration_scale
    assert m1.eta_mb == pytest.approx(m0.eta_mb)
    assert m1.is_complete is True


def test_calibration_scale_and_eta_mb_round_trip_through_the_header(
    full_fits, tmp_path
):
    # The new calibration TEMPERATURE SCALE + eta_mb survive load -> save -> load through the
    # FITS header projection (issue #25), as the enum (not a bare string) and the float.
    loaded = load(full_fits)
    written = save(
        loaded.data, loaded.metadata, tmp_path / "cal", base_header=loaded.header
    )
    m = load(written).metadata
    assert m.calibration_scale is CalibrationScale.T_MB
    assert m.eta_mb == pytest.approx(0.8)


def test_distance_roundtrips_with_its_unit(full_fits, tmp_path):
    loaded = load(full_fits)
    written = save(
        loaded.data, loaded.metadata, tmp_path / "d", base_header=loaded.header
    )
    dist = load(written).metadata.distance
    assert dist.unit.is_equivalent(u.Mpc)
    assert u.isclose(dist, 9.84 * u.Mpc, rtol=1e-12)


# --------------------------------------------------------------------------------------
# Filename projection == header-derived name (legacy astrolyze-1 grammar)
# --------------------------------------------------------------------------------------
def test_filename_projection_matches_header_derived_name(full_fits, tmp_path):
    loaded = load(full_fits)
    expected = "NGC0628_ALMA_CO21_K_12.00x10.00a30.0.fits"
    assert project(loaded.metadata) == expected

    written = save(
        loaded.data, loaded.metadata, tmp_path / "out", base_header=loaded.header
    )
    assert written.name == expected


@pytest.mark.parametrize(
    "major, minor, pa, token",
    [
        (12 * u.arcsec, 12 * u.arcsec, 0 * u.deg, "12.00"),
        (12 * u.arcsec, 10 * u.arcsec, 0 * u.deg, "12.00x10.00"),
        (12 * u.arcsec, 12 * u.arcsec, 30 * u.deg, "12.00a30.0"),
        (12 * u.arcsec, 10 * u.arcsec, 30 * u.deg, "12.00x10.00a30.0"),
    ],
)
def test_resolution_grammar(major, minor, pa, token):
    beam = radio_beam.Beam(major=major, minor=minor, pa=pa)
    m = Metadata(object="S", telescope="T", species="L", bunit=u.K, beam=beam)
    assert project(m) == f"S_T_L_K_{token}.fits"


@pytest.mark.parametrize(
    "unit, token",
    [
        (u.K, "K"),
        (u.Jy / u.beam, "Jybeam"),
        (u.MJy / u.sr, "MJysr"),
        (u.Jy / u.sr, "Jysr"),
    ],
)
def test_flux_unit_is_sanitised_in_filename(unit, token):
    m = Metadata(object="S", telescope="T", species="L", bunit=unit, beam=BEAM)
    assert project(m) == f"S_T_L_{token}_12.00x10.00a30.0.fits"


def test_incomplete_metadata_still_projects_a_filename():
    # No field set: every token degrades to "unknown" rather than crashing.
    assert project(Metadata()) == "unknown_unknown_unknown_unknown_unknown.fits"


def test_projection_extension_is_configurable():
    m = Metadata(object="S", telescope="T", species="L", bunit=u.K, beam=BEAM)
    assert project(m, extension="fits.gz").endswith("_12.00x10.00a30.0.fits.gz")


# --------------------------------------------------------------------------------------
# Save writes a *derived* file and never touches the raw input (ADR-0006)
# --------------------------------------------------------------------------------------
def test_save_does_not_touch_the_source_raw_file(full_fits, tmp_path):
    raw_bytes = full_fits.read_bytes()
    raw_mtime = os.stat(full_fits).st_mtime

    written = save(
        load(full_fits).data,
        load(full_fits).metadata,
        tmp_path / "derived",
        base_header=load(full_fits).header,
    )

    assert written != full_fits
    assert full_fits.exists()
    assert full_fits.read_bytes() == raw_bytes
    assert os.stat(full_fits).st_mtime == raw_mtime


def test_save_refuses_to_overwrite_unless_asked(full_fits, tmp_path):
    loaded = load(full_fits)
    out = tmp_path / "out"
    save(loaded.data, loaded.metadata, out, base_header=loaded.header)
    with pytest.raises(OSError):
        save(loaded.data, loaded.metadata, out, base_header=loaded.header)
    # With overwrite=True it succeeds.
    again = save(
        loaded.data, loaded.metadata, out, base_header=loaded.header, overwrite=True
    )
    assert again.exists()


# --------------------------------------------------------------------------------------
# Velocity convention is read from the header, not guessed (header-authoritative)
# --------------------------------------------------------------------------------------
def test_convention_inferred_from_ctype3_when_keyword_absent(tmp_path):
    h = _full_header()
    del h["HIERARCH ASTROLYZE VCONV"]  # remove our keyword; CTYPE3='VRAD' remains
    path = tmp_path / "vrad.fits"
    fits.writeto(path, np.zeros((2, 3, 3), dtype="float32"), h)
    m = load(path).metadata
    # Reading what CTYPE3 *states* is not guessing — VRAD is the radio convention.
    assert m.velocity_convention is VelocityConvention.RADIO
    assert m.is_complete is True


@pytest.mark.parametrize(
    "ctype3, expected",
    [
        ("VRAD", VelocityConvention.RADIO),
        ("VOPT", VelocityConvention.OPTICAL),
        ("VELO", VelocityConvention.RELATIVISTIC),
    ],
)
def test_ctype3_velocity_codes_map_to_conventions(tmp_path, ctype3, expected):
    h = _full_header()
    del h["HIERARCH ASTROLYZE VCONV"]
    h["CTYPE3"] = ctype3
    path = tmp_path / f"{ctype3}.fits"
    fits.writeto(path, np.zeros((2, 3, 3), dtype="float32"), h)
    assert load(path).metadata.velocity_convention is expected


def test_explicit_vconv_keyword_wins_over_ctype3(tmp_path):
    h = _full_header()
    h["CTYPE3"] = "VOPT"  # disagrees with VCONV='radio'
    path = tmp_path / "conflict.fits"
    fits.writeto(path, np.zeros((2, 3, 3), dtype="float32"), h)
    # The dedicated astrolyze keyword is the authoritative source.
    assert load(path).metadata.velocity_convention is VelocityConvention.RADIO


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code (the legacy sin this module avoids)
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_library_code():
    import re
    from pathlib import Path

    import astrolyze.io as io_pkg

    pkg_dir = Path(io_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
