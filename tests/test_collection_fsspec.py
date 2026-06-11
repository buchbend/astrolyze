"""Collection over fsspec URLs end-to-end (issue #63).

The same ``Collection.open(...)`` / ``Record.open()`` call shape must serve a local directory, a
``file://`` URL, and a non-local fsspec backend (``memory://``) identically — that is the whole
point of the catalog-as-contract corpus (PRD #56, user story 1). This file is the collection-side
companion to ``test_io_fsspec`` (the io-backend side) and reuses the published-corpus fixture
shape from ``test_collection``, but builds the corpus *through a fsspec URL* so the open path is
exercised over a remote-style root.

Assertions pinned here:

- a :class:`~astrolyze.collection.Collection` opened via a ``file://`` / ``memory://`` URL reads
  its catalog and lists its records exactly as the local-path open does;
- a record opened from a URL-rooted corpus returns a lazy, dask-backed Cube whose origin
  Metadata reflects the **URL** (the store URI is the fsspec URL, not a bare local path);
- a missing catalog *at a URL* raises a clear :class:`FileNotFoundError` (the deferred-error
  pitfall the #57 reader note flagged for remote backends), not an obscure parquet error.

No real network or S3 — ``file://`` and ``memory://`` are the CI remote stand-ins.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

import dask
import fsspec

from astrolyze.core import Cube
from astrolyze.io import LoadedData, Metadata

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz


# --------------------------------------------------------------------------------------
# Fixtures: a tiny published corpus written through a fsspec root (file:// or memory://)
# --------------------------------------------------------------------------------------
def _header(*, obj, bmaj_arcsec=13.0):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = obj
    h["TELESCOP"] = "IRAM30M"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST_CO21, "Hz")
    beam = radio_beam.Beam(
        major=bmaj_arcsec * u.arcsec, minor=(bmaj_arcsec - 2) * u.arcsec, pa=30 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO"
    return h, beam


def _loaded(data, header):
    from astropy.wcs import WCS

    return LoadedData(
        data=data,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )


def _catalog_table(rows, *, schema_version="1.0"):
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    columns = {name: [row.get(name) for row in rows] for name in CATALOG_COLUMNS}
    table = pa.table(columns)
    return table.replace_schema_metadata({"catalog_schema_version": schema_version})


@pytest.fixture(autouse=True)
def _clean_memory_fs():
    """A fresh fsspec in-memory filesystem per test (it is a process-global singleton)."""
    from fsspec.implementations.memory import MemoryFileSystem

    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]
    yield
    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]


def _build_corpus(root_url: str) -> str:
    """Write one store + a ``catalog.parquet`` under the fsspec *root_url*; return that root.

    The store is saved via the io backend through the URL (the #63 save path); the catalog is
    written with pyarrow through a fsspec handle so the corpus root is genuinely the given
    backend (``file://`` or ``memory://``) rather than a local directory."""
    header, beam = _header(obj="NGC3521", bmaj_arcsec=13.0)
    data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
    cube = Cube.from_loaded(_loaded(data, header))
    # Save the store under the URL root; its leaf name is the header projection (opaque here).
    saved = cube.to_zarr(root_url + "/l1")
    leaf = str(saved).rstrip("/").split("/")[-1]
    store_path = f"l1/{leaf}"

    rows = [
        {
            "object": "NGC3521",
            "survey": "HERACLES",
            "telescope": "IRAM30M",
            "species": "CO",
            "transition": "2-1",
            "rest_frequency_hz": float(REST_CO21),
            "beam_major_arcsec": beam.major.to_value(u.arcsec),
            "beam_minor_arcsec": beam.minor.to_value(u.arcsec),
            "beam_pa_deg": beam.pa.to_value(u.deg),
            "bunit": "K",
            "store_path": store_path,
            "content_checksum": "sha256:heracles",
            "ra_deg": 170.0,
            "dec_deg": -0.04,
            "radius_deg": 0.01,
            "catalog_schema_version": "1.0",
        }
    ]
    with fsspec.open(root_url + "/catalog.parquet", "wb") as handle:
        pq.write_table(_catalog_table(rows), handle)
    return root_url


@pytest.fixture
def file_corpus(tmp_path):
    """A published corpus rooted at a ``file://`` URL (the CI remote stand-in)."""
    return _build_corpus(f"file://{tmp_path / 'ism_corpus'}")


@pytest.fixture
def memory_corpus():
    """A published corpus rooted at a fsspec ``memory://`` URL (a non-local backend)."""
    return _build_corpus("memory://ism_corpus")


# --------------------------------------------------------------------------------------
# Collection.open + records over a URL root (file:// and memory://)
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("corpus_fixture", ["file_corpus", "memory_corpus"])
def test_open_collection_through_url(corpus_fixture, request):
    from astrolyze.collection import Collection

    root = request.getfixturevalue(corpus_fixture)
    collection = Collection.open(root)
    assert collection.catalog_version == "1.0"
    assert len(collection.records) == 1
    assert collection.records[0].object == "NGC3521"


@pytest.mark.parametrize("corpus_fixture", ["file_corpus", "memory_corpus"])
def test_record_open_through_url_is_lazy(corpus_fixture, request):
    from astrolyze.collection import Collection

    root = request.getfixturevalue(corpus_fixture)
    collection = Collection.open(root)
    cube = collection.records[0].open()
    assert isinstance(cube, Cube)
    # Lazy over the remote-style store — the URL open path must not eagerly download.
    assert dask.is_dask_collection(cube._sc._data)


@pytest.mark.parametrize("corpus_fixture", ["file_corpus", "memory_corpus"])
def test_origin_uri_reflects_the_url(corpus_fixture, request):
    from astrolyze.collection import Collection

    root = request.getfixturevalue(corpus_fixture)
    collection = Collection.open(root)
    record = collection.records[0]
    cube = record.open()
    # The origin store URI is the fsspec URL (scheme preserved), not a bare local path: a product
    # saved into an experiment must trace back to the corpus *as addressed*, on disk or on S3.
    uri = cube.metadata.origin_store_uri
    assert uri is not None
    scheme = root.split("://", 1)[0]
    assert uri.startswith(scheme + "://")
    assert uri.endswith(".zarr")
    assert cube.metadata.origin_catalog_version == "1.0"


# --------------------------------------------------------------------------------------
# Missing catalog at a URL raises a clear error (the deferred-error pitfall, #57 note)
# --------------------------------------------------------------------------------------
def test_missing_catalog_at_file_url_raises_filenotfound(tmp_path):
    from astrolyze.collection.catalog import read_catalog

    empty = tmp_path / "no_catalog"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        read_catalog(f"file://{empty}")


def test_missing_catalog_at_memory_url_raises_filenotfound():
    """A non-local backend that might defer its open must still fail clearly, not obscurely."""
    from astrolyze.collection.catalog import read_catalog

    with pytest.raises(FileNotFoundError):
        read_catalog("memory://empty_corpus")


def test_collection_open_missing_catalog_at_url_raises_filenotfound():
    from astrolyze.collection import Collection

    with pytest.raises(FileNotFoundError):
        Collection.open("memory://empty_corpus")
