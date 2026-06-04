"""Tests for astrolyze.experiment.manifest — the DB-backed dataset registry (issue #11,
ADR-0009).

Written first (red/green TDD). These are behavioural tests over the public surface — what a
caller observes from ``register`` / ``get`` / ``query`` / ``all`` — never the ORM internals:

- ``Manifest(db_url)`` opens (or creates) a registry; here on **in-memory sqlite**;
- ``register`` stores identity (source path + DOI) and the *full* ``io`` Metadata provenance;
- ``get`` / ``query`` / ``all`` round-trip that identity and provenance — crucially the
  reconstructed beam, rest frequency, velocity convention and species (we round-trip the
  *value*, persisting primitives, not the astropy/radio_beam objects);
- registration is **idempotent on source path** — re-registering the same path updates the
  existing row in place (no duplicate, stable id), the latest values winning;
- the backend is swappable: a real file-backed sqlite persists across two ``Manifest``
  instances pointed at the same DB.

The manifest is generic storage; the *merciless* completeness gate is ingest's job (#12), so
the manifest itself happily stores partial provenance as nulls.
"""

from datetime import datetime

import pytest
import astropy.units as u
import radio_beam

from astrolyze.experiment import Experiment
from astrolyze.experiment.manifest import DatasetRecord, Manifest
from astrolyze.io.schema import SCHEMA_VERSION, Metadata
from astrolyze.units import VelocityConvention

# The reference dataset mirrors the tracer-bullet PHANGS cube (NGC 0628, CO(2-1)).
REST = 230.538 * u.GHz
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


def _complete_metadata() -> Metadata:
    """A schema-complete Metadata carrying every provenance field (a fresh object each call)."""
    return Metadata(
        object="NGC0628",
        telescope="ALMA",
        species="CO21",
        rest_frequency=REST,
        velocity_convention=VelocityConvention.RADIO,
        beam=BEAM,
        bunit=u.K,
        distance=9.84 * u.Mpc,
        calibration_error=0.1,
        name_tag="mom0",
    )


@pytest.fixture
def manifest() -> Manifest:
    """A registry backed by in-memory sqlite (one shared connection for the test)."""
    return Manifest("sqlite://")


# --------------------------------------------------------------------------------------
# register -> get round-trips identity AND the full provenance
# --------------------------------------------------------------------------------------
def test_register_returns_a_record_with_a_stable_id(manifest):
    record = manifest.register(_complete_metadata(), "data/raw/ngc0628_co21.fits")
    assert isinstance(record, DatasetRecord)
    assert record.id is not None
    assert record.source_path == "data/raw/ngc0628_co21.fits"


def test_register_then_get_roundtrips_full_provenance(manifest):
    record = manifest.register(
        _complete_metadata(), "data/raw/ngc0628_co21.fits", doi="10.1234/abcd"
    )
    got = manifest.get(record.id)
    assert got is not None

    # Identity
    assert got.id == record.id
    assert got.source_path == "data/raw/ngc0628_co21.fits"
    assert got.doi == "10.1234/abcd"
    assert got.schema_version == SCHEMA_VERSION

    # Provenance — reconstructed from stored primitives back into the io schema types.
    m = got.metadata
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


def test_register_without_doi_stores_none(manifest):
    record = manifest.register(_complete_metadata(), "data/raw/a.fits")
    assert manifest.get(record.id).doi is None


def test_ingested_at_is_recorded(manifest):
    record = manifest.register(_complete_metadata(), "data/raw/a.fits")
    assert isinstance(record.ingested_at, datetime)


def test_get_unknown_id_returns_none(manifest):
    assert manifest.get(404) is None


# --------------------------------------------------------------------------------------
# query: find datasets by provenance fields (e.g. object / species)
# --------------------------------------------------------------------------------------
def test_query_by_object_and_by_species(manifest):
    manifest.register(_complete_metadata(), "data/raw/a.fits")
    other = _complete_metadata()
    other.object = "NGC1068"
    other.species = "HCN10"
    manifest.register(other, "data/raw/b.fits")

    by_object = manifest.query(object="NGC0628")
    assert [r.source_path for r in by_object] == ["data/raw/a.fits"]

    by_species = manifest.query(species="HCN10")
    assert [r.source_path for r in by_species] == ["data/raw/b.fits"]


def test_query_no_match_returns_empty(manifest):
    manifest.register(_complete_metadata(), "data/raw/a.fits")
    assert manifest.query(object="UNKNOWN") == []


def test_query_unknown_filter_raises(manifest):
    with pytest.raises(ValueError):
        manifest.query(not_a_column="x")


def test_all_returns_every_registered_dataset(manifest):
    manifest.register(_complete_metadata(), "data/raw/a.fits")
    manifest.register(_complete_metadata(), "data/raw/b.fits")
    assert {r.source_path for r in manifest.all()} == {
        "data/raw/a.fits",
        "data/raw/b.fits",
    }


# --------------------------------------------------------------------------------------
# Idempotent on source path: re-register updates in place (no duplicate, stable id)
# --------------------------------------------------------------------------------------
def test_reregister_same_path_does_not_duplicate(manifest):
    first = manifest.register(_complete_metadata(), "data/raw/a.fits", doi="10.1/first")
    second = manifest.register(
        _complete_metadata(), "data/raw/a.fits", doi="10.2/second"
    )

    assert second.id == first.id  # stable id across re-ingest
    assert len(manifest.all()) == 1
    assert manifest.get(first.id).doi == "10.2/second"  # latest values win


def test_reregister_updates_provenance_in_place(manifest):
    record = manifest.register(_complete_metadata(), "data/raw/a.fits")
    updated = _complete_metadata()
    updated.species = "13CO21"
    manifest.register(updated, "data/raw/a.fits")

    assert len(manifest.all()) == 1
    assert manifest.get(record.id).metadata.species == "13CO21"


# --------------------------------------------------------------------------------------
# The manifest stores whatever provenance is present (the merciless gate is ingest's, #12)
# --------------------------------------------------------------------------------------
def test_partial_metadata_persists_optional_fields_as_none(manifest):
    partial = Metadata(
        object="NGC0628",
        rest_frequency=REST,
        velocity_convention=VelocityConvention.RADIO,
    )
    got = manifest.get(manifest.register(partial, "data/raw/partial.fits").id)
    assert got.metadata.object == "NGC0628"
    assert got.metadata.beam is None
    assert got.metadata.distance is None
    assert got.metadata.species is None


# --------------------------------------------------------------------------------------
# Backend swappable: a real file-backed sqlite persists across Manifest instances
# --------------------------------------------------------------------------------------
def test_file_backend_persists_across_instances(tmp_path):
    url = f"sqlite:///{tmp_path / 'manifest.db'}"
    record = Manifest(url).register(
        _complete_metadata(), "data/raw/a.fits", doi="10.3/x"
    )

    reopened = Manifest(url)  # fresh instance, same DB file on disk
    got = reopened.get(record.id)
    assert got is not None
    assert got.source_path == "data/raw/a.fits"
    assert got.doi == "10.3/x"
    assert got.metadata.object == "NGC0628"
    assert u.isclose(got.metadata.rest_frequency, REST, rtol=1e-12)


# --------------------------------------------------------------------------------------
# for_experiment: the relative config db_url resolves against the experiment root, not CWD
# --------------------------------------------------------------------------------------
def test_for_experiment_places_db_inside_the_experiment(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    record = Manifest.for_experiment(experiment).register(
        _complete_metadata(), "data/raw/a.fits"
    )
    # The default config db_url is relative (sqlite:///manifest.db); it must land inside the
    # experiment regardless of the process CWD, and round-trip on reopen.
    assert (experiment.root / "manifest.db").is_file()
    assert (
        Manifest.for_experiment(experiment).get(record.id).source_path
        == "data/raw/a.fits"
    )


def test_for_experiment_honours_an_edited_config_db_url(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    experiment.config.write_text('[manifest]\ndb_url = "sqlite:///registry.db"\n')
    Manifest.for_experiment(experiment).register(
        _complete_metadata(), "data/raw/a.fits"
    )
    assert (experiment.root / "registry.db").is_file()


def test_for_experiment_leaves_in_memory_url_unanchored(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    # An explicit db_url overrides the config; in-memory must not spill a file into the tree.
    Manifest.for_experiment(experiment, db_url="sqlite://").register(
        _complete_metadata(), "data/raw/a.fits"
    )
    assert not (experiment.root / "manifest.db").exists()
