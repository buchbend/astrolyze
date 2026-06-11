"""Tests for the ``collection`` capability — catalog reader + Collection facade (issue #57).

Written first (red/green TDD). The collection tracer opens a *published corpus* — a directory
of astrolyze-flavoured Zarr stores with a single ``catalog.parquet`` index at its root — and
lists it object-first, opening any record into a lazy, dask-backed :class:`~astrolyze.core.Cube`
that carries its corpus origin in :class:`~astrolyze.io.Metadata`.

The contract under test is the consumer side of the catalog schema spec (the normative document
lives in the ifm repo): one Parquet row per published cube, columns in a fixed order, a
``catalog_schema_version`` of ``"1.0"`` (``MAJOR.MINOR``; an unknown MAJOR must be rejected
explicitly). The seams pinned here are the ones #60 (describe/query), #61 (scan-builder) and
#62 (covering) build on, so they are exercised directly:

1. **Catalog read + schema-validate** — ``read_catalog`` parses ``catalog.parquet`` into typed
   :class:`CatalogRow` records; an unknown/incompatible ``catalog_schema_version`` raises a
   descriptive :class:`CatalogSchemaError` (never a bare ``KeyError``).
2. **Object-first ``list()``** — one summary per source, its surveys / species / beam range
   aggregated across that source's stores.
3. **fsspec from day one** — ``Collection.open`` routes path handling through fsspec, so a
   local dir today and an ``s3://`` URL later share one call shape.
4. **Lazy Cubes with origin** — ``record.open()`` returns a dask-backed Cube whose Metadata
   carries the source store URI + catalog version, round-tripping through ``save``.

The reference corpus mirrors the real one: NGC 3521 observed in HERACLES (CO 2-1) and
PHANGS-ALMA (CO 2-1), plus NGC 0628 in THINGS (HI 21cm) — the matched multi-survey source is
exactly the #62/#64 acceptance fixture.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.io import load, save

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")  # the collection reader needs parquet
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz
REST_HI = 1.420405751e9  # HI 21cm Hz


# --------------------------------------------------------------------------------------
# Synthetic corpus fixtures: a couple of tiny astrolyze-flavoured Zarr stores + a catalog
# --------------------------------------------------------------------------------------
def _header(*, obj, rest_hz, ctype3="VRAD", bmaj_arcsec=13.0):
    """A small schema-complete cube header with a real celestial+spectral WCS."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = ctype3, 0.0
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


def _write_store(root, rel_path, *, obj, rest_hz, telescope, species, bmaj_arcsec):
    """Write one synthetic Zarr store under *root* and return its relative POSIX path + beam."""
    header, beam = _header(obj=obj, rest_hz=rest_hz, bmaj_arcsec=bmaj_arcsec)
    header["TELESCOP"] = telescope
    header["HIERARCH ASTROLYZE SPECIES"] = species
    data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
    cube = Cube.from_loaded(_loaded(data, header))
    store = root / rel_path
    store.parent.mkdir(parents=True, exist_ok=True)
    cube.to_zarr(store.parent)  # writes a store named by the header projection
    # Find the store just written (the projection name is opaque to the test).
    written = next(store.parent.glob("*.zarr"))
    return written, beam


def _loaded(data, header):
    """A LoadedData built straight from an in-memory array + header (no file)."""
    from astropy.wcs import WCS

    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    wcs = WCS(header)
    return LoadedData(
        data=data,
        wcs=wcs,
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )


def _catalog_table(rows, *, schema_version="1.0"):
    """Build a pyarrow table for *rows* in the contract column order, stamping the version."""
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    columns = {name: [row.get(name) for row in rows] for name in CATALOG_COLUMNS}
    table = pa.table(columns)
    # The version also lives in the file key-value metadata so a reader can gate from the header.
    return table.replace_schema_metadata({"catalog_schema_version": schema_version})


@pytest.fixture
def corpus(tmp_path):
    """A tiny published corpus: 3 stores across 2 objects / 3 surveys, plus catalog.parquet."""
    root = tmp_path / "ism_corpus"
    root.mkdir()

    specs = [
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

    rows = []
    for spec in specs:
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
# catalog deep module: read + schema-validate
# --------------------------------------------------------------------------------------
def test_read_catalog_returns_typed_rows(corpus):
    from astrolyze.collection.catalog import CatalogRow, read_catalog

    catalog = read_catalog(str(corpus))
    assert catalog.schema_version == "1.0"
    assert len(catalog.rows) == 3
    row = next(r for r in catalog.rows if r.survey == "HERACLES")
    assert isinstance(row, CatalogRow)
    assert row.object == "NGC3521"
    assert row.species == "CO"
    assert row.transition == "2-1"
    assert row.store_path == "l1/heraclesa.zarr" or row.store_path.endswith(".zarr")
    assert row.beam_major_arcsec == pytest.approx(13.0)


def test_unknown_schema_version_raises_descriptive_error(tmp_path):
    from astrolyze.collection.catalog import CatalogSchemaError, read_catalog

    root = tmp_path / "future_corpus"
    root.mkdir()
    rows = [
        {
            "object": "NGC3521",
            "survey": "HERACLES",
            "telescope": "IRAM30M",
            "species": "CO",
            "transition": "2-1",
            "rest_frequency_hz": REST_CO21,
            "beam_major_arcsec": 13.0,
            "beam_minor_arcsec": 11.0,
            "beam_pa_deg": 30.0,
            "bunit": "K",
            "store_path": "l1/x.zarr",
            "content_checksum": "sha256:x",
            "ra_deg": 170.0,
            "dec_deg": -0.04,
            "radius_deg": 0.01,
            "catalog_schema_version": "2.0",
        }
    ]
    pq.write_table(_catalog_table(rows, schema_version="2.0"), root / "catalog.parquet")

    with pytest.raises(CatalogSchemaError) as exc:
        read_catalog(str(root))
    message = str(exc.value)
    assert "2.0" in message  # the offending version is named
    assert "1" in message  # what astrolyze supports is named


def test_minor_bump_is_accepted(tmp_path):
    """A future 1.x MINOR (extra trailing columns) stays readable — additive evolution."""
    from astrolyze.collection.catalog import read_catalog

    root = tmp_path / "minor_corpus"
    root.mkdir()
    rows = [
        {
            "object": "NGC3521",
            "survey": "HERACLES",
            "telescope": "IRAM30M",
            "species": "CO",
            "transition": "2-1",
            "rest_frequency_hz": REST_CO21,
            "beam_major_arcsec": 13.0,
            "beam_minor_arcsec": 11.0,
            "beam_pa_deg": 30.0,
            "bunit": "K",
            "store_path": "l1/x.zarr",
            "content_checksum": "sha256:x",
            "ra_deg": 170.0,
            "dec_deg": -0.04,
            "radius_deg": 0.01,
            "catalog_schema_version": "1.7",
        }
    ]
    table = _catalog_table(rows, schema_version="1.7")
    # An additive trailing column a 1.0 reader has never heard of.
    table = table.append_column("v_sys_kms", pa.array([805.0]))
    pq.write_table(table, root / "catalog.parquet")

    catalog = read_catalog(str(root))
    assert catalog.schema_version == "1.7"
    assert len(catalog.rows) == 1


def test_read_catalog_missing_file_raises_filenotfound(tmp_path):
    from astrolyze.collection.catalog import read_catalog

    empty = tmp_path / "no_catalog"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        read_catalog(str(empty))


def test_null_store_path_raises_at_read_time(tmp_path):
    """The locator is non-null by contract: a row missing store_path is a non-conforming catalog
    and must fail descriptively at read time, not produce a broken URI that blows up at open()."""
    from astrolyze.collection.catalog import CatalogSchemaError, read_catalog

    root = tmp_path / "broken_corpus"
    root.mkdir()
    rows = [
        {
            "object": "NGC3521",
            "survey": "HERACLES",
            "telescope": "IRAM30M",
            "species": "CO",
            "transition": "2-1",
            "rest_frequency_hz": REST_CO21,
            "beam_major_arcsec": 13.0,
            "beam_minor_arcsec": 11.0,
            "beam_pa_deg": 30.0,
            "bunit": "K",
            "store_path": None,  # the non-conforming case
            "content_checksum": "sha256:x",
            "ra_deg": 170.0,
            "dec_deg": -0.04,
            "radius_deg": 0.01,
            "catalog_schema_version": "1.0",
        }
    ]
    pq.write_table(_catalog_table(rows), root / "catalog.parquet")
    with pytest.raises(CatalogSchemaError) as exc:
        read_catalog(str(root))
    assert "store_path" in str(exc.value)


# --------------------------------------------------------------------------------------
# Collection facade: open + object-first list
# --------------------------------------------------------------------------------------
def test_open_local_corpus(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    assert collection.catalog_version == "1.0"
    assert len(collection.records) == 3


def test_list_groups_object_first(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    summaries = collection.list()
    by_object = {s.object: s for s in summaries}
    assert set(by_object) == {"NGC3521", "NGC0628"}

    ngc3521 = by_object["NGC3521"]
    # NGC3521 appears in two surveys / two stores, both CO(2-1).
    assert set(ngc3521.surveys) == {"HERACLES", "PHANGS-ALMA"}
    assert set(ngc3521.species) == {"CO"}
    assert ngc3521.n_stores == 2
    # Beam range spans the coarse single-dish (13") and the fine interferometer (1.5") beams.
    lo, hi = ngc3521.beam_range_arcsec
    assert lo == pytest.approx(1.5)
    assert hi == pytest.approx(13.0)

    ngc0628 = by_object["NGC0628"]
    assert set(ngc0628.surveys) == {"THINGS"}
    assert set(ngc0628.species) == {"HI"}
    assert ngc0628.n_stores == 1


def test_records_for_object(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    records = [r for r in collection.records if r.object == "NGC3521"]
    assert len(records) == 2
    assert {r.survey for r in records} == {"HERACLES", "PHANGS-ALMA"}


# --------------------------------------------------------------------------------------
# Lazy dask-backed Cubes with origin provenance
# --------------------------------------------------------------------------------------
def test_record_open_returns_lazy_dask_cube(corpus):
    import dask

    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    record = next(r for r in collection.records if r.survey == "HERACLES")
    cube = record.open()
    assert isinstance(cube, Cube)
    # Lazy: the Zarr backend keeps the underlying array a dask collection (issue #23) — opening a
    # record never materialises the cube. Probed on the spectral-cube backing array, the same way
    # the io suite asserts a Zarr load stays lazy.
    assert dask.is_dask_collection(cube._sc._data)


def test_open_stamps_origin_metadata(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    record = next(r for r in collection.records if r.survey == "HERACLES")
    cube = record.open()
    assert cube.metadata.origin_catalog_version == "1.0"
    assert cube.metadata.origin_store_uri is not None
    assert cube.metadata.origin_store_uri.endswith(".zarr")
    # The URI points at the record's store (its physical home in the corpus).
    assert record.store_path.split("/")[-1] in cube.metadata.origin_store_uri


def test_origin_roundtrips_through_save(corpus, tmp_path):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    record = next(r for r in collection.records if r.survey == "HERACLES")
    cube = record.open()

    # Save the collection-born cube into a fresh store (the experiment-layer bridge): the origin
    # lineage must survive the round-trip so a saved product traces back to the corpus snapshot.
    out = save(
        np.asarray(cube._data_quantity.value),
        cube.metadata,
        tmp_path / "experiment",
        format="zarr",
    )
    back = load(out)
    assert back.metadata.origin_store_uri == cube.metadata.origin_store_uri
    assert back.metadata.origin_catalog_version == "1.0"
