"""Tests for astrolyze.experiment.ingest — the merciless ingest gate (issue #12, ADR-0009).

Written first (red/green TDD). These are behavioural tests over the public surface — what a
caller observes from the :class:`IngestReport` and the resulting manifest — never internals:

- ``ingest`` walks ``raw/`` and partitions every dataset into *accepted* (schema-complete) and
  *rejected* (with the per-file list of missing mandatory fields);
- a schema-complete cube is accepted **and registered** in the manifest, with full provenance;
- an incomplete cube is rejected and the report **names the exact missing fields**
  (``rest_frequency``, ``velocity_convention``) — the merciless gate tells you what to fix;
- merciless = rejected files are **not** registered;
- ``raw/`` is sacred: ingest never renames or modifies anything under it, and drops no files in;
- re-ingest after a header is completed now accepts that file, idempotently (no duplicate row);
- the mandatory-context set is overridable (explicit ``required=`` and via ``[ingest]`` config).

Synthetic FITS are built header-first (mirroring tests/test_io.py); the array is a tiny stub —
ingest reads only the header, never the (potentially huge) raw cube body.
"""

import os

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.experiment import Experiment, Manifest
from astrolyze.experiment.ingest import (
    AcceptedDataset,
    IngestReport,
    RejectedDataset,
    ingest,
)

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


# --------------------------------------------------------------------------------------
# Synthetic FITS headers (the header carries the full / partial schema)
# --------------------------------------------------------------------------------------
def _spatial_wcs(header):
    header["CTYPE1"], header["CRVAL1"] = "RA---SIN", 24.174
    header["CDELT1"], header["CRPIX1"], header["CUNIT1"] = -2e-4, 1.0, "deg"
    header["CTYPE2"], header["CRVAL2"] = "DEC--SIN", 15.784
    header["CDELT2"], header["CRPIX2"], header["CUNIT2"] = 2e-4, 1.0, "deg"


def _complete_header(object_="NGC0628"):
    """A cube header carrying the complete mandatory context (rest freq + velocity conv)."""
    h = fits.Header()
    _spatial_wcs(h)
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = object_
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


def _incomplete_header(object_="NGC0628"):
    """A 2D map header with an object but no rest frequency or velocity convention —
    missing both mandatory-context fields."""
    h = fits.Header()
    _spatial_wcs(h)
    h["OBJECT"] = object_
    h["BUNIT"] = "K"
    return h


def _write(path, header, shape):
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.writeto(path, np.zeros(shape, dtype="float32"), header)
    return path


def _experiment(tmp_path):
    return Experiment.init(tmp_path / "study")


# --------------------------------------------------------------------------------------
# ingest returns a report partitioning raw/ into accepted and rejected
# --------------------------------------------------------------------------------------
def test_ingest_partitions_raw_into_accepted_and_rejected(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    _write(exp.raw / "bad.fits", _incomplete_header(), (3, 3))

    report = ingest(exp)

    assert isinstance(report, IngestReport)
    assert [a.source_path for a in report.accepted] == ["data/raw/good.fits"]
    assert [r.source_path for r in report.rejected] == ["data/raw/bad.fits"]
    assert report.n_accepted == 1
    assert report.n_rejected == 1


# --------------------------------------------------------------------------------------
# A schema-complete cube is accepted and registered (with full provenance)
# --------------------------------------------------------------------------------------
def test_complete_cube_is_accepted_and_registered(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))

    report = ingest(exp)
    assert len(report.accepted) == 1
    accepted = report.accepted[0]
    assert isinstance(accepted, AcceptedDataset)

    # Registered in the manifest the experiment config points at.
    manifest = Manifest.for_experiment(exp)
    rows = manifest.query(object="NGC0628")
    assert [r.source_path for r in rows] == ["data/raw/good.fits"]

    # The accepted entry carries the registered record with its reconstructed provenance.
    m = accepted.record.metadata
    assert m.object == "NGC0628"
    assert m.species == "CO21"
    assert u.isclose(m.rest_frequency, REST, rtol=1e-12)
    assert u.isclose(m.beam.major, 12 * u.arcsec, rtol=1e-9)


def test_accepted_record_is_retrievable_by_id(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    record = ingest(exp).accepted[0].record
    assert Manifest.for_experiment(exp).get(record.id).source_path == "data/raw/good.fits"


# --------------------------------------------------------------------------------------
# An incomplete cube is rejected, naming the EXACT missing mandatory fields, and is NOT
# registered (merciless)
# --------------------------------------------------------------------------------------
def test_incomplete_cube_is_rejected_naming_missing_fields(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "bad.fits", _incomplete_header(), (3, 3))

    report = ingest(exp)
    assert len(report.rejected) == 1
    rejected = report.rejected[0]
    assert isinstance(rejected, RejectedDataset)
    assert set(rejected.missing) == {"rest_frequency", "velocity_convention"}


def test_rejected_file_is_not_registered(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "bad.fits", _incomplete_header(), (3, 3))
    ingest(exp)
    # Merciless: nothing incomplete reaches the manifest.
    assert Manifest.for_experiment(exp).all() == []


# --------------------------------------------------------------------------------------
# raw/ is sacred: ingest never renames or modifies it, and drops no files into it
# --------------------------------------------------------------------------------------
def test_ingest_does_not_modify_or_rename_raw(tmp_path):
    exp = _experiment(tmp_path)
    good = _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    bad = _write(exp.raw / "bad.fits", _incomplete_header(), (3, 3))

    before = {
        p: (p.read_bytes(), os.stat(p).st_mtime)
        for p in (good, bad)
    }
    names_before = sorted(p.name for p in exp.raw.rglob("*") if p.is_file())

    ingest(exp)

    for path, (raw_bytes, mtime) in before.items():
        assert path.exists()
        assert path.read_bytes() == raw_bytes
        assert os.stat(path).st_mtime == mtime
    # No renames and nothing new written under raw/ (the manifest db lives at the root).
    assert sorted(p.name for p in exp.raw.rglob("*") if p.is_file()) == names_before
    assert not (exp.raw / "manifest.db").exists()


# --------------------------------------------------------------------------------------
# Re-ingest after completing a header now accepts the file (idempotent, no duplicate)
# --------------------------------------------------------------------------------------
def test_reingest_after_completing_header_accepts(tmp_path):
    exp = _experiment(tmp_path)
    path = _write(exp.raw / "ngc0628.fits", _incomplete_header(), (3, 3))

    first = ingest(exp)
    assert [r.source_path for r in first.rejected] == ["data/raw/ngc0628.fits"]
    assert first.accepted == []

    # Fix the header in place (a normal, iterative step) and re-ingest.
    path.unlink()
    _write(path, _complete_header(), (2, 3, 3))

    second = ingest(exp)
    assert [a.source_path for a in second.accepted] == ["data/raw/ngc0628.fits"]
    assert second.rejected == []
    assert len(Manifest.for_experiment(exp).all()) == 1  # registered once, no duplicate


def test_reingest_of_complete_cube_does_not_duplicate(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    ingest(exp)
    ingest(exp)  # idempotent on source path
    assert len(Manifest.for_experiment(exp).all()) == 1


# --------------------------------------------------------------------------------------
# Coverage details: recursion, non-FITS ignored, empty raw/
# --------------------------------------------------------------------------------------
def test_ingest_recurses_into_subdirectories(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "phangs" / "ngc0628.fits", _complete_header(), (2, 3, 3))
    report = ingest(exp)
    assert [a.source_path for a in report.accepted] == ["data/raw/phangs/ngc0628.fits"]


def test_ingest_ignores_non_fits_files(tmp_path):
    exp = _experiment(tmp_path)
    (exp.raw / "README.md").write_text("notes about these data\n")
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    report = ingest(exp)
    # The README is neither accepted nor rejected — ingest is about datasets, not every file.
    assert [a.source_path for a in report.accepted] == ["data/raw/good.fits"]
    assert report.rejected == []


def test_ingest_empty_raw_returns_empty_report(tmp_path):
    report = ingest(_experiment(tmp_path))
    assert report.accepted == []
    assert report.rejected == []


# --------------------------------------------------------------------------------------
# The mandatory-context set is overridable (explicit arg, and via [ingest] config)
# --------------------------------------------------------------------------------------
def test_required_override_argument_rejects_missing_beam(tmp_path):
    exp = _experiment(tmp_path)
    # A cube complete on the default mandatory context but with no beam.
    header = _complete_header()
    for key in ("BMAJ", "BMIN", "BPA"):
        del header[key]
    _write(exp.raw / "nobeam.fits", header, (2, 3, 3))

    report = ingest(exp, required=("rest_frequency", "velocity_convention", "beam"))
    assert [r.source_path for r in report.rejected] == ["data/raw/nobeam.fits"]
    assert report.rejected[0].missing == ["beam"]


def test_required_override_from_config(tmp_path):
    exp = _experiment(tmp_path)
    exp.config.write_text(
        '[manifest]\ndb_url = "sqlite:///manifest.db"\n\n'
        '[ingest]\nrequired_context = ["rest_frequency", "velocity_convention", "beam"]\n'
    )
    header = _complete_header()
    for key in ("BMAJ", "BMIN", "BPA"):
        del header[key]
    _write(exp.raw / "nobeam.fits", header, (2, 3, 3))

    report = ingest(exp)  # required set comes from the config
    assert report.rejected[0].missing == ["beam"]


def test_unknown_required_field_raises(tmp_path):
    exp = _experiment(tmp_path)
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))
    with pytest.raises(ValueError):
        ingest(exp, required=("not_a_schema_field",))


# --------------------------------------------------------------------------------------
# A corrupt / non-FITS-but-named-.fits file is rejected, not crashing the merciless pass
# --------------------------------------------------------------------------------------
def test_corrupt_fits_is_rejected_with_a_reason(tmp_path):
    exp = _experiment(tmp_path)
    (exp.raw / "junk.fits").write_text("this is not a FITS file\n")
    _write(exp.raw / "good.fits", _complete_header(), (2, 3, 3))

    report = ingest(exp)
    assert [a.source_path for a in report.accepted] == ["data/raw/good.fits"]
    junk = [r for r in report.rejected if r.source_path == "data/raw/junk.fits"]
    assert len(junk) == 1
    assert junk[0].error is not None  # the unreadable file is flagged, the pass continues
