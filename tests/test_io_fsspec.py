"""fsspec URLs end-to-end for the Zarr backend and the collection (issue #63).

One corpus, one call shape — on disk or on S3. This file pins the contract that the io Zarr
backend and the :class:`~astrolyze.collection.Collection` open path accept fsspec URLs (and
fsspec mappers) identically to a local directory path, so moving a corpus from a local directory
to an object store changes *nothing* in analysis code (PRD #56, user story 1).

No real network or S3 in CI. The remote backend is stood in for two ways, both fsspec
filesystems that are not the plain local path:

- ``file://`` — the explicit-scheme local filesystem (the CI remote stand-in the PRD names);
- ``memory://`` — fsspec's in-memory filesystem, a genuinely non-``LocalFileSystem`` backend
  whose bytes never touch disk, so a green test here proves the path layer is storage-neutral
  rather than quietly falling back to a local open.

The assertions mirror the local-path suite (``test_io_zarr`` / ``test_collection``): a store
saved and reloaded through a URL round-trips its data and stays lazy (dask-backed); a record
opened through a URL carries its origin store URI; and a missing catalog at a URL raises a
clear :class:`FileNotFoundError`, never an obscure deferred parquet error.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

import dask

from astrolyze.core import Cube
from astrolyze.io import LoadedData, Metadata, load, save

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz


# --------------------------------------------------------------------------------------
# Fixtures: a tiny schema-complete cube + helpers shared by the URL round-trip tests
# --------------------------------------------------------------------------------------
def _full_header(*, obj="NGC0628") -> fits.Header:
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = obj
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST_CO21, "Hz")
    beam = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


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


@pytest.fixture(autouse=True)
def _clean_memory_fs():
    """Each test gets a fresh fsspec in-memory filesystem (it is a process-global singleton)."""
    from fsspec.implementations.memory import MemoryFileSystem

    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]
    yield
    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]


def _file_url(path) -> str:
    """A ``file://`` URL for a local *path* (the CI remote stand-in the PRD names)."""
    return f"file://{path}"


# --------------------------------------------------------------------------------------
# io Zarr backend: save + load through a fsspec URL (file://) and a non-local backend
# --------------------------------------------------------------------------------------
def test_save_and_load_through_file_url(tmp_path):
    """A store saved and reloaded through a ``file://`` URL round-trips data + metadata."""
    loaded = _loaded(
        np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5), _full_header()
    )
    store = save(
        loaded.data,
        loaded.metadata,
        _file_url(tmp_path / "z"),
        format="zarr",
        base_header=loaded.header,
    )
    back = load(_file_url(str(store).removeprefix("file://")))
    assert isinstance(back, LoadedData)
    assert back.metadata.object == "NGC0628"
    np.testing.assert_array_equal(np.asarray(back.data), loaded.data)


def test_load_through_file_url_stays_lazy(tmp_path):
    """A URL-loaded Zarr store is still dask-backed — the URL path must not eagerly download."""
    loaded = _loaded(
        np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5), _full_header()
    )
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(2, 3, 5),
        base_header=loaded.header,
    )
    back = load(_file_url(store))
    assert dask.is_dask_collection(back.data)
    assert back.data.chunksize[0] == 2


def test_save_and_load_through_memory_url():
    """A genuinely non-local backend (fsspec ``memory://``) round-trips and stays lazy.

    memory:// never touches disk, so a green round-trip proves the backend opens a store from a
    fsspec mapper rather than quietly resolving a local path."""
    loaded = _loaded(
        np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5), _full_header()
    )
    store = save(
        loaded.data,
        loaded.metadata,
        "memory://corpus/z",
        format="zarr",
        base_header=loaded.header,
    )
    back = load("memory://corpus/z/" + str(store).rstrip("/").split("/")[-1])
    assert dask.is_dask_collection(back.data)
    np.testing.assert_array_equal(np.asarray(back.data.compute()), loaded.data)


def test_cube_from_zarr_through_memory_url():
    """``Cube.from_zarr`` accepts a fsspec URL identically to a local path (lazy cube)."""
    loaded = _loaded(
        np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5), _full_header()
    )
    store = save(
        loaded.data,
        loaded.metadata,
        "memory://corpus/z",
        format="zarr",
        base_header=loaded.header,
    )
    leaf = str(store).rstrip("/").split("/")[-1]
    cube = Cube.from_zarr("memory://corpus/z/" + leaf)
    assert isinstance(cube, Cube)
    assert dask.is_dask_collection(cube._sc._data)
