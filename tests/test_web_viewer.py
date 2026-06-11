"""Contract tests for the cube-viewer panel feeds (issue #67).

The viewer's four backend endpoints are *thin, lazy* wrappers over the public
:class:`~astrolyze.core.Cube` API. These tests pin each endpoint's JSON shape and values against a
**synthetic corpus** (the same two-tiny-Zarr-stores-plus-catalog shape the #66 explorer tests use),
driven through FastAPI's ``TestClient`` (no browser, no real network):

- ``GET /api/stores/{id}/cube`` — axis metadata: shape, n_channels, velocity axis + unit, bunit,
  sky extent. Read from headers/axes, no data plane.
- ``GET /api/stores/{id}/moment0`` — the integrated map: dims match the cube's sky plane; unit is
  the integrated unit (K·km/s); vmin/vmax bracket the data.
- ``GET /api/stores/{id}/channel/{i}`` — one channel's 2-D slice + its velocity (matches the axis);
  an out-of-range channel is a graceful ``422``.
- ``GET /api/stores/{id}/spectrum?x=&y=`` — a pixel spectrum of length == n_channels; an
  out-of-range pixel is a ``422``.

Three cross-cutting properties are pinned too:

1. **Laziness** (PRD #56, story 8): on a dask-backed store, a channel / spectrum request must NOT
   materialise the whole cube. Verified two ways — the opened cube's data stays a dask graph, and a
   spy on ``dask.array.Array.compute`` shows a single-channel request never computes an array as
   large as the full cube.
2. **fsspec parity**: the same contracts hold on a ``memory://`` collection as on disk (PRD #56,
   story 1) — the viewer rides the same fsspec opener as the rest of the collection API.
3. **Unknown store → 404**: an id not in the corpus is a clean 404 on every panel feed.
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
fastapi = pytest.importorskip("fastapi")
fsspec = pytest.importorskip("fsspec")
from fastapi.testclient import TestClient  # noqa: E402

REST_CO21 = 230.538e9

# The synthetic cube dimensions (n_channels, ny, nx) — small but enough to tell axes apart.
NCHAN, NY, NX = 6, 4, 5


def _header(*, obj, telescope, species, rest_hz, bmaj):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = obj, telescope, "K"
    h["RESTFRQ"] = (rest_hz, "Hz")
    beam = radio_beam.Beam(
        major=bmaj * u.arcsec, minor=(bmaj - 2) * u.arcsec, pa=30 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = species
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


def _build_corpus(root_url: str) -> str:
    """Write one CO 2-1 cube + a catalog.parquet under the fsspec *root_url*; return that root.

    Mirrors the proven ``test_collection_fsspec`` corpus builder (the store via the io backend
    through the URL, the catalog via a fsspec handle) so the same path serves a local ``file://``
    and a ``memory://`` root unchanged. **Per-channel zarr chunking** (``chunks=(1, NY, NX)``) is
    this fixture's one addition and its whole point: it makes a single-channel slice touch exactly
    ONE on-disk chunk, so the laziness assertion ("a channel request never computes the whole cube")
    is meaningful rather than trivially true on a one-chunk store."""
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    h, beam = _header(
        obj="NGC3521",
        telescope="IRAM30M",
        species="CO",
        rest_hz=REST_CO21,
        bmaj=13.0,
    )
    # A ramp so every channel / pixel has a distinct, checkable value.
    data = np.arange(NCHAN * NY * NX, dtype="float32").reshape(NCHAN, NY, NX)
    saved = Cube.from_loaded(_loaded(data, h)).to_zarr(
        root_url + "/l1", chunks=(1, NY, NX)
    )
    leaf = str(saved).rstrip("/").split("/")[-1]
    store_path = f"l1/{leaf}"

    row = {
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
        "content_checksum": "sha256:x",
        "ra_deg": 170.0,
        "dec_deg": -0.04,
        "radius_deg": 0.01,
        "catalog_schema_version": "1.0",
    }
    table = pa.table({k: [row.get(k)] for k in CATALOG_COLUMNS})
    table = table.replace_schema_metadata({"catalog_schema_version": "1.0"})
    with fsspec.open(root_url + "/catalog.parquet", "wb") as handle:
        pq.write_table(table, handle)
    return root_url, store_path


@pytest.fixture(autouse=True)
def _clean_memory_fs():
    """A fresh fsspec in-memory filesystem per test (it is a process-global singleton)."""
    from fsspec.implementations.memory import MemoryFileSystem

    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]
    yield
    MemoryFileSystem.store.clear()
    MemoryFileSystem.pseudo_dirs[:] = [""]


@pytest.fixture
def corpus(tmp_path):
    """An on-disk synthetic corpus: one CO 2-1 cube + catalog, per-channel chunked.

    Returns ``(root_url, store_path)`` where ``root_url`` is a ``file://`` URL (the local stand-in),
    so the on-disk and memory:// fixtures share one shape."""
    root_url, store_path = _build_corpus(f"file://{tmp_path / 'ism_corpus'}")
    return root_url, store_path


@pytest.fixture
def client(corpus):
    """A TestClient over the explorer app built on the on-disk synthetic corpus."""
    from astrolyze.web.api import create_app

    root_url, _ = corpus
    return TestClient(create_app(root_url))


def _store_id(store_path):
    from astrolyze.web.api import _encode_store_id

    return _encode_store_id(store_path)


# --------------------------------------------------------------------------------------
# GET /api/stores/{id}/cube — axis metadata (no data plane)
# --------------------------------------------------------------------------------------
def test_cube_endpoint_returns_axis_metadata(client, corpus):
    _, store_path = corpus
    body = client.get(f"/api/stores/{_store_id(store_path)}/cube").json()
    assert body["shape"] == [NCHAN, NY, NX]
    assert body["n_channels"] == NCHAN
    assert body["bunit"] == "K"
    # The velocity axis is the cube's own (CDELT3 = 2000 m/s → 0, 2, 4, ... km/s).
    assert body["velocity_unit"] == "km / s"
    assert body["velocity"] == pytest.approx([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    # The sky extent brackets the WCS reference pixel (RA 170°, Dec -0.04°).
    ext = body["sky_extent"]
    assert ext["ra_min_deg"] <= 170.0 <= ext["ra_max_deg"]
    assert ext["dec_min_deg"] <= -0.04 <= ext["dec_max_deg"]


def test_cube_endpoint_unknown_store_is_404(client):
    response = client.get(f"/api/stores/{_store_id('l1/nope/missing.zarr')}/cube")
    assert response.status_code == 404


# --------------------------------------------------------------------------------------
# GET /api/stores/{id}/moment0 — the integrated map
# --------------------------------------------------------------------------------------
def test_moment0_endpoint_dims_and_unit(client, corpus):
    _, store_path = corpus
    body = client.get(f"/api/stores/{_store_id(store_path)}/moment0").json()
    # The moment-0 map is the cube's SKY plane: (ny, nx).
    assert body["shape"] == [NY, NX]
    assert len(body["data"]) == NY and len(body["data"][0]) == NX
    # Integrating K over km/s gives K·km/s (spectral-cube's unit, carried through faithfully).
    assert "K" in body["unit"]
    assert body["downsample"] == 1  # tiny map, no decimation
    # vmin/vmax bracket the data plane.
    flat = [v for row in body["data"] for v in row if v is not None]
    assert body["vmin"] == pytest.approx(min(flat))
    assert body["vmax"] == pytest.approx(max(flat))


def test_moment0_matches_public_cube_moment0(client, corpus):
    """The served map is byte-for-byte the public Cube.moment0() (dogfooding, not reimplemented)."""
    from astrolyze.collection import Collection

    root_url, store_path = corpus
    body = client.get(f"/api/stores/{_store_id(store_path)}/moment0").json()
    record = next(
        r for r in Collection.open(root_url).records if r.store_path == store_path
    )
    expected = np.asarray(record.open().moment0().data.value, dtype="float64")
    served = np.array(
        [[np.nan if v is None else v for v in row] for row in body["data"]]
    )
    np.testing.assert_allclose(served, expected, rtol=1e-6)


# --------------------------------------------------------------------------------------
# GET /api/stores/{id}/channel/{i} — one channel + its velocity
# --------------------------------------------------------------------------------------
def test_channel_endpoint_returns_slice_and_velocity(client, corpus):
    _, store_path = corpus
    body = client.get(f"/api/stores/{_store_id(store_path)}/channel/2").json()
    assert body["index"] == 2
    assert body["shape"] == [NY, NX]
    # Channel 2 velocity == velocity axis[2] == 4 km/s (CDELT3 2000 m/s).
    assert body["velocity"] == pytest.approx(4.0)
    assert body["velocity_unit"] == "km / s"
    # The slice is the cube's channel-2 plane (ramp: channel c starts at c*NY*NX).
    served = np.array(
        [[np.nan if v is None else v for v in row] for row in body["data"]]
    )
    assert served[0, 0] == pytest.approx(2 * NY * NX)


def test_channel_out_of_range_is_422(client, corpus):
    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/channel/{NCHAN}")
    assert response.status_code == 422
    assert str(NCHAN) in response.json()["detail"]


def test_channel_unknown_store_is_404(client):
    response = client.get(f"/api/stores/{_store_id('l1/nope/missing.zarr')}/channel/0")
    assert response.status_code == 404


# --------------------------------------------------------------------------------------
# GET /api/stores/{id}/spectrum?x=&y= — a pixel spectrum
# --------------------------------------------------------------------------------------
def test_spectrum_endpoint_length_and_values(client, corpus):
    _, store_path = corpus
    body = client.get(f"/api/stores/{_store_id(store_path)}/spectrum?x=2&y=1").json()
    assert body["x"] == 2 and body["y"] == 1
    # Length == n_channels; velocity axis matches the cube endpoint's.
    assert len(body["value"]) == NCHAN
    assert len(body["velocity"]) == NCHAN
    assert body["velocity"] == pytest.approx([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    # Pixel (x=2, y=1): ramp value at channel c is c*NY*NX + 1*NX + 2.
    expected = [c * NY * NX + 1 * NX + 2 for c in range(NCHAN)]
    assert body["value"] == pytest.approx(expected)


def test_spectrum_out_of_range_pixel_is_422(client, corpus):
    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/spectrum?x={NX}&y=0")
    assert response.status_code == 422


def test_spectrum_missing_query_is_422(client, corpus):
    """x/y are required query params; omitting them is FastAPI's validation 422."""
    _, store_path = corpus
    assert (
        client.get(f"/api/stores/{_store_id(store_path)}/spectrum").status_code == 422
    )


# --------------------------------------------------------------------------------------
# Laziness (PRD #56, story 8): a slice request must not materialise the whole cube.
# --------------------------------------------------------------------------------------
def test_opened_cube_stays_a_dask_graph(corpus):
    """record.open() yields a dask-backed cube — nothing is materialised at open time."""
    import dask.array as da

    from astrolyze.collection import Collection

    root_url, store_path = corpus
    record = next(
        r for r in Collection.open(root_url).records if r.store_path == store_path
    )
    cube = record.open()
    assert isinstance(cube._sc._data, da.Array)


def test_channel_request_does_not_compute_full_cube(client, corpus, monkeypatch):
    """A single-channel request never computes an array as large as the full cube.

    Spy on ``dask.array.Array.compute``: with per-channel chunking, the channel-2 slice computes
    one channel's worth of elements (NY*NX), never the full NCHAN*NY*NX cube. The assertion is the
    inequality (no computed array reaches full-cube size), which is robust to spectral-cube's
    internal extra computes (masks etc.) while still failing loudly on a full materialisation."""
    import dask.array as da

    full_cube_size = NCHAN * NY * NX
    computed_sizes = []
    real_compute = da.Array.compute

    def spy_compute(self, *args, **kwargs):
        computed_sizes.append(int(np.prod(self.shape)))
        return real_compute(self, *args, **kwargs)

    monkeypatch.setattr(da.Array, "compute", spy_compute)

    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/channel/2")
    assert response.status_code == 200
    assert computed_sizes, "expected at least one dask compute for the channel slice"
    assert max(computed_sizes) < full_cube_size, (
        f"a channel request computed an array of {max(computed_sizes)} elements — "
        f"the full cube is {full_cube_size}; the slice was not lazy"
    )


def test_spectrum_request_does_not_compute_full_cube(client, corpus, monkeypatch):
    """A single-pixel spectrum computes one column (NCHAN), never the full cube."""
    import dask.array as da

    full_cube_size = NCHAN * NY * NX
    computed_sizes = []
    real_compute = da.Array.compute

    def spy_compute(self, *args, **kwargs):
        computed_sizes.append(int(np.prod(self.shape)))
        return real_compute(self, *args, **kwargs)

    monkeypatch.setattr(da.Array, "compute", spy_compute)

    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/spectrum?x=2&y=1")
    assert response.status_code == 200
    assert computed_sizes, "expected at least one dask compute for the spectrum slice"
    assert max(computed_sizes) < full_cube_size


# --------------------------------------------------------------------------------------
# fsspec parity (PRD #56, story 1): the same contracts on a memory:// collection.
# --------------------------------------------------------------------------------------
def test_viewer_works_on_memory_fsspec_collection():
    """Every panel feed serves identically over a memory:// collection (fsspec parity).

    The same builder + app on a ``memory://`` root (the CI non-local stand-in). The autouse
    ``_clean_memory_fs`` fixture keeps the process-global in-memory store fresh per test."""
    from astrolyze.web.api import create_app

    root_uri, store_path = _build_corpus("memory://ism_corpus_mem")
    client = TestClient(create_app(root_uri))
    sid = _store_id(store_path)

    cube = client.get(f"/api/stores/{sid}/cube").json()
    assert cube["shape"] == [NCHAN, NY, NX]

    moment = client.get(f"/api/stores/{sid}/moment0").json()
    assert moment["shape"] == [NY, NX]

    channel = client.get(f"/api/stores/{sid}/channel/2").json()
    assert channel["velocity"] == pytest.approx(4.0)

    spectrum = client.get(f"/api/stores/{sid}/spectrum?x=2&y=1").json()
    assert len(spectrum["value"]) == NCHAN
