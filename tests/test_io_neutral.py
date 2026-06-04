"""Tests for the format-neutral I/O seam (issue #22, ADR-0006).

This file is the correctness obligation for *unwelding* ``io`` from FITS so a second
backend can slot in without touching the schema. It pins two contracts:

1. **The schema has two projections, not one.** :class:`~astrolyze.io.Metadata` is the
   in-memory truth; ``to_header``/``from_header`` is its FITS projection and
   ``to_attrs``/``from_attrs`` is a JSON-able projection of the *same* dataclass. The attrs
   dict maps each field to a primitive (``Quantity -> {"value","unit"}``,
   ``radio_beam.Beam -> {"bmaj","bmin","bpa"}``, enums/units -> ``str``) so a non-FITS
   backend (e.g. Zarr ``.attrs``) can carry the context verbatim. Attrs are *lazy like the
   header* (ADR-0006 ii): a partial ``Metadata`` projects and round-trips while staying
   flagged incomplete — we never refuse it and never silently complete it (ADR-0003).
2. **``LoadedData`` is backend-neutral.** A live ``fits.Header`` becomes optional; the exact
   WCS travels as the *verbatim FITS-WCS header string*, from which any backend reconstructs
   the WCS with ``WCS(Header.fromstring(s))`` (astropy still owns WCS reconstruction).

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import json

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from astropy.wcs import WCS

from astrolyze.io import LoadedData, Metadata, load
from astrolyze.io.schema import REQUIRED_CONTEXT
from astrolyze.units import VelocityConvention

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


def _full_metadata() -> Metadata:
    """A fully-populated schema — every field the AC enumerates is set."""
    return Metadata(
        object="NGC0628",
        telescope="ALMA",
        species="CO21",
        rest_frequency=REST.to(u.Hz),
        velocity_convention=VelocityConvention.RADIO,
        beam=BEAM,
        bunit=u.K,
        distance=9.84 * u.Mpc,
        calibration_error=0.1,
        name_tag="mom0",
    )


# --------------------------------------------------------------------------------------
# Synthetic FITS fixture (reused for the LoadedData / WCS contract)
# --------------------------------------------------------------------------------------
def _full_header() -> fits.Header:
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
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    h["HIERARCH ASTROLYZE DISTANCE"] = 9.84
    h["HIERARCH ASTROLYZE DISTUNIT"] = "Mpc"
    h["HIERARCH ASTROLYZE CALERR"] = 0.1
    h["HIERARCH ASTROLYZE NAMETAG"] = "mom0"
    return h


@pytest.fixture
def full_fits(tmp_path):
    path = tmp_path / "raw" / "ngc0628_co21.fits"
    path.parent.mkdir(parents=True)
    data = np.arange(2 * 3 * 3, dtype="float32").reshape(2, 3, 3)
    fits.writeto(path, data, _full_header())
    return path


# --------------------------------------------------------------------------------------
# to_attrs is JSON-serializable
# --------------------------------------------------------------------------------------
def test_to_attrs_is_json_serializable():
    attrs = _full_metadata().to_attrs()
    # The whole point of the attrs projection: it survives json.dumps unchanged.
    dumped = json.dumps(attrs)
    assert isinstance(dumped, str)
    # And the structured shapes are primitives, not astropy objects.
    assert attrs["rest_frequency"] == {"value": REST.to_value(u.Hz), "unit": "Hz"}
    assert set(attrs["beam"]) == {"bmaj", "bmin", "bpa"}
    assert attrs["velocity_convention"] == "radio"
    assert attrs["bunit"] == "K"


# --------------------------------------------------------------------------------------
# Full round-trip: from_attrs(to_attrs(m)) == m
# --------------------------------------------------------------------------------------
def test_full_metadata_round_trips_through_attrs():
    m = _full_metadata()
    back = Metadata.from_attrs(m.to_attrs())
    assert back == m


def test_units_survive_attrs_round_trip_exactly():
    m = _full_metadata()
    back = Metadata.from_attrs(m.to_attrs())

    # rest frequency in Hz
    assert back.rest_frequency.unit == u.Hz
    assert u.isclose(back.rest_frequency, REST, rtol=1e-12)
    # distance carries its unit
    assert back.distance.unit.is_equivalent(u.Mpc)
    assert u.isclose(back.distance, 9.84 * u.Mpc, rtol=1e-12)
    # beam bmaj/bmin/bpa
    assert u.isclose(back.beam.major, 12 * u.arcsec, rtol=1e-12)
    assert u.isclose(back.beam.minor, 10 * u.arcsec, rtol=1e-12)
    assert u.isclose(back.beam.pa, 30 * u.deg, rtol=1e-12)
    # bunit string
    assert back.bunit == u.K
    # velocity convention enum (not a bare string)
    assert back.velocity_convention is VelocityConvention.RADIO
    # calibration error
    assert back.calibration_error == pytest.approx(0.1)


# --------------------------------------------------------------------------------------
# Partial metadata: projects, round-trips, and stays flagged incomplete (ADR-0006 ii)
# --------------------------------------------------------------------------------------
def test_partial_metadata_round_trips_and_stays_incomplete():
    # An object/telescope but NO rest frequency or velocity convention — the archival case.
    partial = Metadata(object="NGC0628", telescope="ALMA", bunit=u.K)
    assert partial.is_complete is False
    missing_before = set(partial.missing)
    assert missing_before == set(REQUIRED_CONTEXT)

    attrs = partial.to_attrs()
    # The projection itself is still JSON-able even when context is absent.
    json.dumps(attrs)

    back = Metadata.from_attrs(attrs)
    assert back == partial
    # Lazy enforcement is preserved through the attrs projection: still incomplete, same gaps.
    assert back.is_complete is False
    assert set(back.missing) == missing_before


# --------------------------------------------------------------------------------------
# LoadedData is backend-neutral: header optional, WCS travels as a verbatim string
# --------------------------------------------------------------------------------------
def test_loaded_data_constructs_without_a_fits_header():
    data = np.zeros((3, 3), dtype="float32")
    wcs_string = _full_header().tostring()
    # No live fits.Header — a non-FITS backend has none, only the schema + a WCS string.
    loaded = LoadedData(
        data=data,
        wcs=WCS(fits.Header.fromstring(wcs_string)),
        metadata=Metadata(object="NGC0628"),
        path="NGC0628.zarr",
        header_string=wcs_string,
    )
    assert loaded.header is None
    assert isinstance(loaded.header_string, str)


def test_loaded_data_wcs_string_reconstructs_the_loaded_wcs(full_fits):
    loaded = load(full_fits)
    # The verbatim FITS-WCS header string is carried; reconstructing from it yields the
    # SAME WCS the loader produced (astropy owns the reconstruction).
    assert isinstance(loaded.header_string, str)
    reconstructed = WCS(fits.Header.fromstring(loaded.header_string))
    assert reconstructed.wcs.compare(loaded.wcs.wcs)


# --------------------------------------------------------------------------------------
# load() is behaviourally unchanged: still lazy, still populated header + wcs + metadata
# --------------------------------------------------------------------------------------
def test_load_still_populates_header_wcs_metadata(full_fits):
    loaded = load(full_fits)
    assert isinstance(loaded.header, fits.Header)
    assert isinstance(loaded.wcs, WCS)
    assert isinstance(loaded.metadata, Metadata)
    assert loaded.metadata.is_complete is True
    assert loaded.data.shape == (2, 3, 3)
