"""Tests for the ``collection`` drill-down + filtering layer — describe + query (issue #60).

Written first (red/green TDD). On top of the #57 facade (open / object-first ``list`` / lazy
records), #60 adds two consumer operations over the *same* catalog rows:

- ``Collection.describe(object)`` expands one source into per-store physical detail (beam,
  bunit, rest frequency, transition, provenance summary). Per the PRD the answer comes from the
  **catalog row first** — describing a remote corpus must stay cheap — and only ``deep=True``
  opens the live stores to fill the fields the catalog does not carry (channel width, velocity
  coverage). An unknown object raises a descriptive error, not an empty list.

- ``Collection.query(**filters)`` slices the catalog along any metadata axis (species, survey,
  telescope, transition, …) and returns a *sub-Collection* (composable: every facade method,
  including ``describe`` / ``query`` again, works on the result). An unknown filter key raises a
  clear error rather than silently matching nothing.

The synthetic corpus mirrors ``tests/test_collection.py``: NGC 3521 in HERACLES (CO 2-1) and
PHANGS-ALMA (CO 2-1), plus NGC 0628 in THINGS (HI 21cm) — the matched multi-survey source.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz
REST_HI = 1.420405751e9  # HI 21cm Hz


# --------------------------------------------------------------------------------------
# Synthetic corpus fixtures (the #57 patterns, reused)
# --------------------------------------------------------------------------------------
def _header(*, obj, rest_hz, bmaj_arcsec=13.0):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = obj
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (rest_hz, "Hz")
    beam = radio_beam.Beam(
        major=bmaj_arcsec * u.arcsec, minor=(bmaj_arcsec - 2) * u.arcsec, pa=30 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    return h, beam


def _loaded(data, header):
    from astropy.wcs import WCS

    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    return LoadedData(
        data=data,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )


def _write_store(root, rel_path, *, obj, rest_hz, telescope, species, bmaj_arcsec):
    header, beam = _header(obj=obj, rest_hz=rest_hz, bmaj_arcsec=bmaj_arcsec)
    header["TELESCOP"] = telescope
    header["HIERARCH ASTROLYZE SPECIES"] = species
    data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
    store = root / rel_path
    store.parent.mkdir(parents=True, exist_ok=True)
    Cube.from_loaded(_loaded(data, header)).to_zarr(store.parent)
    written = next(store.parent.glob("*.zarr"))
    return written, beam


def _catalog_table(rows, *, schema_version="1.0"):
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    columns = {name: [row.get(name) for row in rows] for name in CATALOG_COLUMNS}
    table = pa.table(columns)
    return table.replace_schema_metadata({"catalog_schema_version": schema_version})


_SPECS = [
    dict(
        rel_path="l1/heraclesa",
        obj="NGC3521",
        rest_hz=REST_CO21,
        telescope="IRAM30M",
        survey="HERACLES",
        species="CO",
        transition="2-1",
        bmaj_arcsec=13.0,
    ),
    dict(
        rel_path="l1/phangsa",
        obj="NGC3521",
        rest_hz=REST_CO21,
        telescope="ALMA",
        survey="PHANGS-ALMA",
        species="CO",
        transition="2-1",
        bmaj_arcsec=1.5,
    ),
    dict(
        rel_path="l1/thingsa",
        obj="NGC0628",
        rest_hz=REST_HI,
        telescope="VLA",
        survey="THINGS",
        species="HI",
        transition="21cm",
        bmaj_arcsec=11.0,
    ),
]


@pytest.fixture
def corpus(tmp_path):
    """A tiny published corpus: 3 stores across 2 objects / 3 surveys, plus catalog.parquet."""
    root = tmp_path / "ism_corpus"
    root.mkdir()

    rows = []
    for spec in _SPECS:
        store, beam = _write_store(
            root,
            spec["rel_path"],
            obj=spec["obj"],
            rest_hz=spec["rest_hz"],
            telescope=spec["telescope"],
            species=spec["species"],
            bmaj_arcsec=spec["bmaj_arcsec"],
        )
        rel = store.relative_to(root).as_posix()
        rows.append(
            {
                "object": spec["obj"],
                "survey": spec["survey"],
                "telescope": spec["telescope"],
                "species": spec["species"],
                "transition": spec["transition"],
                "rest_frequency_hz": float(spec["rest_hz"]),
                "beam_major_arcsec": beam.major.to_value(u.arcsec),
                "beam_minor_arcsec": beam.minor.to_value(u.arcsec),
                "beam_pa_deg": beam.pa.to_value(u.deg),
                "bunit": "K",
                "store_path": rel,
                "content_checksum": f"sha256:{spec['rel_path']}",
                "ra_deg": 170.0,
                "dec_deg": -0.04,
                "radius_deg": 0.01,
                "catalog_schema_version": "1.0",
            }
        )

    pq.write_table(_catalog_table(rows), root / "catalog.parquet")
    return root


# --------------------------------------------------------------------------------------
# describe(object): per-store physical detail, catalog-row-first
# --------------------------------------------------------------------------------------
def test_describe_lists_every_store_of_a_source(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    details = collection.describe("NGC3521")
    # Both of NGC3521's stores (HERACLES + PHANGS-ALMA), one detail each.
    assert len(details) == 2
    by_survey = {d.survey: d for d in details}
    assert set(by_survey) == {"HERACLES", "PHANGS-ALMA"}

    heracles = by_survey["HERACLES"]
    assert heracles.object == "NGC3521"
    assert heracles.telescope == "IRAM30M"
    assert heracles.species == "CO"
    assert heracles.transition == "2-1"
    assert heracles.bunit == "K"
    assert heracles.beam_major_arcsec == pytest.approx(13.0)
    assert heracles.beam_minor_arcsec == pytest.approx(11.0)
    assert heracles.beam_pa_deg == pytest.approx(30.0)
    assert heracles.rest_frequency_hz == pytest.approx(REST_CO21)
    assert heracles.store_path.endswith(".zarr")


def test_describe_unknown_object_raises_descriptive_error(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    with pytest.raises(KeyError) as exc:
        collection.describe("NGC9999")
    message = str(exc.value)
    assert "NGC9999" in message
    # The known objects are named so the message is actionable.
    assert "NGC3521" in message or "NGC0628" in message


def test_describe_shallow_does_not_open_any_store(corpus, monkeypatch):
    """The PRD rule: catalog-only describe issues no store reads — verifiable by making any
    record.open() blow up and confirming a shallow describe never reaches it."""
    from astrolyze.collection import Collection
    from astrolyze.collection import _facade as facade

    def _boom(self):
        raise AssertionError("shallow describe must not open a live store")

    monkeypatch.setattr(facade.Record, "open", _boom)

    collection = Collection.open(str(corpus))
    details = collection.describe("NGC3521")  # deep defaults False
    assert len(details) == 2
    # The deep-only fields are absent (None) on a shallow describe.
    assert all(d.channel_width_kms is None for d in details)
    assert all(d.velocity_min_kms is None for d in details)
    assert all(d.velocity_max_kms is None for d in details)


def test_describe_deep_reads_live_store_attrs(corpus):
    """deep=True opens each store and fills the fields the catalog does not carry: native
    channel width and velocity coverage, read from the store's own spectral axis."""
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    details = collection.describe("NGC3521", deep=True)
    assert len(details) == 2
    for detail in details:
        # The synthetic stores have a 3-channel axis at 0, 2, 4 km/s (CDELT3 = 2000 m/s).
        assert detail.channel_width_kms == pytest.approx(2.0)
        assert detail.velocity_min_kms == pytest.approx(0.0)
        assert detail.velocity_max_kms == pytest.approx(4.0)
        # The catalog-row fields are still present (deep augments, never replaces).
        assert detail.species == "CO"


# --------------------------------------------------------------------------------------
# query(**filters): metadata-axis filtering -> sub-Collection
# --------------------------------------------------------------------------------------
def test_query_single_filter_hits(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    hi = collection.query(species="HI")
    assert isinstance(hi, Collection)
    assert len(hi.records) == 1
    assert hi.records[0].object == "NGC0628"


def test_query_returns_composable_subcollection(corpus):
    """The query result is a Collection: list/describe/query again all work on it."""
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    co = collection.query(species="CO")
    # list() on the sub-collection groups only the CO stores.
    summaries = co.list()
    assert {s.object for s in summaries} == {"NGC3521"}
    # describe() on the sub-collection still expands the source.
    assert len(co.describe("NGC3521")) == 2
    # query() chains.
    heracles = co.query(survey="HERACLES")
    assert len(heracles.records) == 1


def test_query_miss_returns_empty_collection(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    empty = collection.query(species="SiO")
    assert isinstance(empty, Collection)
    assert empty.records == []
    assert empty.list() == []


def test_query_multi_filter_is_conjunctive(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    # CO AND PHANGS-ALMA -> the single interferometer store.
    hit = collection.query(species="CO", survey="PHANGS-ALMA")
    assert len(hit.records) == 1
    assert hit.records[0].row.telescope == "ALMA"
    # CO AND THINGS -> no store matches both (THINGS is HI) -> empty, not an error.
    miss = collection.query(species="CO", survey="THINGS")
    assert miss.records == []


def test_query_filters_on_object_and_telescope_and_transition(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    assert len(collection.query(object="NGC3521").records) == 2
    assert len(collection.query(telescope="VLA").records) == 1
    assert len(collection.query(transition="21cm").records) == 1


def test_query_unknown_filter_key_raises(corpus):
    """An unknown filter axis is a user error, not a silent no-op (ADR-0003)."""
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    with pytest.raises((TypeError, ValueError)) as exc:
        collection.query(colour="blue")
    message = str(exc.value)
    assert "colour" in message
    # The valid filter axes are named so the message is actionable.
    assert "species" in message or "survey" in message


def test_query_preserves_root_and_catalog_version(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    sub = collection.query(species="CO")
    assert sub.root_uri == collection.root_uri
    assert sub.catalog_version == collection.catalog_version
    # A record from the sub-collection still opens (its store URI resolves against the same root).
    cube = sub.query(survey="HERACLES").records[0].open()
    assert cube.metadata.origin_store_uri.endswith(".zarr")
