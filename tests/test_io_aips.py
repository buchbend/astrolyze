"""AIPS-era FITS header tolerance (THINGS, ifm#12).

Surveys imaged with classic AIPS (THINGS HI, Walter et al. 2008, AJ 136, 2563) write headers
modern parsers stumble over, in two specific ways the loader must absorb:

1. **All-caps unit strings.** AIPS writes ``BUNIT = 'JY/BEAM'``; astropy unit parsing is
   case-sensitive (wants ``Jy/beam``), so the verbatim string raises ``ValueError``.
2. **The restoring beam lives in HISTORY.** AIPS records the CLEAN beam as a HISTORY card
   (``AIPS   CLEAN BMAJ=  3.9274E-03 BMIN=  3.0983E-03 BPA= -61.66``, degrees) instead of
   BMAJ/BMIN/BPA keywords. ``radio_beam.Beam.from_fits_header`` can parse that line, but only
   if it is actually called — gating it on ``"BMAJ" in header`` skips AIPS headers entirely.

Real-world anchor: every staged THINGS cube (e.g. NGC_3521_NA_CUBE_THINGS.FITS) shows both.
"""

import astropy.units as u
import pytest
from astropy.io import fits

from astrolyze.io import Metadata

# Verbatim from the staged NGC 3521 THINGS cube: 3.9274e-3 deg = 14.139", 3.0983e-3 deg = 11.154".
AIPS_BEAM_HISTORY = "AIPS   CLEAN BMAJ=  3.9274E-03 BMIN=  3.0983E-03 BPA= -61.66"


def _aips_header(*, history_beam: bool) -> fits.Header:
    header = fits.Header()
    header["BUNIT"] = "JY/BEAM"
    header["OBJECT"] = "NGC3521H"
    header["TELESCOP"] = "VLA"
    header["RESTFREQ"] = 1420405750.0
    if history_beam:
        header.add_history("AIPS   IMAGR appears earlier in a real header")
        header.add_history(AIPS_BEAM_HISTORY)
    return header


def test_all_caps_jy_beam_bunit_parses():
    metadata = Metadata.from_header(_aips_header(history_beam=False))
    assert metadata.bunit == u.Jy / u.beam


def test_mixed_case_bunit_still_parses_strictly():
    header = _aips_header(history_beam=False)
    header["BUNIT"] = "Jy/beam"
    assert Metadata.from_header(header).bunit == u.Jy / u.beam


def test_garbage_bunit_still_raises():
    # Tolerance is for case, not for nonsense — an unrecognizable unit stays a loud error
    # (ADR-0003: never silently drop physical context).
    header = _aips_header(history_beam=False)
    header["BUNIT"] = "FURLONGS/FORTNIGHT"
    with pytest.raises(ValueError):
        Metadata.from_header(header)


def test_beam_recovered_from_aips_history():
    metadata = Metadata.from_header(_aips_header(history_beam=True))
    assert metadata.beam is not None
    assert metadata.beam.major.to_value(u.arcsec) == pytest.approx(14.139, abs=0.01)
    assert metadata.beam.minor.to_value(u.arcsec) == pytest.approx(11.154, abs=0.01)
    assert metadata.beam.pa.to_value(u.deg) == pytest.approx(-61.66, abs=0.01)


def test_bmaj_keyword_still_wins_when_present():
    header = _aips_header(history_beam=True)
    header["BMAJ"] = 12.0 / 3600.0
    header["BMIN"] = 10.0 / 3600.0
    header["BPA"] = 30.0
    metadata = Metadata.from_header(header)
    assert metadata.beam.major.to_value(u.arcsec) == pytest.approx(12.0, abs=0.01)


def test_no_beam_anywhere_stays_none():
    metadata = Metadata.from_header(_aips_header(history_beam=False))
    assert metadata.beam is None
