"""Contract tests for the cube-viewer INTERACTIONS (issue #68).

The two #68 server-side reductions extend the #67 panel feeds and stay *thin, lazy* wrappers over
the public :class:`~astrolyze.core.Cube` API. These tests pin each against the same synthetic corpus
the #66/#67 tests use (one small per-channel-chunked CO 2-1 Zarr store + a catalog.parquet), driven
through FastAPI's ``TestClient``:

- ``POST /api/stores/{id}/region-spectrum`` — body ``{"vertices": [[x, y], …]}`` (image-pixel
  polygon) → the region-averaged spectrum: ``velocity`` + mean ``value`` (length == n_channels) +
  ``n_pixels``. The means match a hand-computed ``nanmean`` over a known box; a degenerate region
  (< 3 vertices) is a graceful ``422``.
- ``GET /api/stores/{id}/moment0?vmin=&vmax=`` — moment-0 recomputed over *exactly* the channels in
  the velocity window. The map equals a hand-computed integral over those channels (and differs from
  the full-band moment-0); an empty/half window is a ``422``.

The cross-cutting properties carry over: **laziness** (a windowed moment / region average never
materialises the whole cube — verified with a spy on ``dask.array.Array.compute``), and **fsspec
parity** (the same contracts on a ``memory://`` collection).

The synthetic cube is a ramp ``data[c, y, x] = c*NY*NX + y*NX + x`` on a velocity axis
``[0, 2, 4, 6, 8, 10] km/s`` (CDELT3 = 2000 m/s), so every region mean and windowed integral has a
closed form to check against.
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

# The synthetic cube dimensions (n_channels, ny, nx) — small but enough to check reductions.
NCHAN, NY, NX = 6, 4, 5
CHANNEL_WIDTH_MS = (
    2000.0  # CDELT3; moment-0 integrates value × Δv (m/s here → unit K·m/s)
)


def _header(*, obj, telescope, species, rest_hz, bmaj):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = CHANNEL_WIDTH_MS, 1.0, "m/s"
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


def _build_corpus(root_url: str):
    """Write one CO 2-1 ramp cube + catalog.parquet under *root_url*; return ``(root, store_path)``.

    Identical shape to the #67 viewer fixture (per-channel zarr chunking so a single-channel slab
    touches one chunk — the laziness assertions stay meaningful), reused here so both test modules
    speak the same synthetic corpus."""
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    h, beam = _header(
        obj="NGC3521",
        telescope="IRAM30M",
        species="CO",
        rest_hz=REST_CO21,
        bmaj=13.0,
    )
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
    """An on-disk synthetic corpus: one CO 2-1 ramp cube + catalog, per-channel chunked."""
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


def _ramp(c, y, x):
    """The synthetic cube's value at (channel c, row y, col x)."""
    return c * NY * NX + y * NX + x


# --------------------------------------------------------------------------------------
# POST /api/stores/{id}/region-spectrum — the region-averaged spectrum
# --------------------------------------------------------------------------------------
def test_region_spectrum_matches_hand_computed_box(client, corpus):
    """A polygon enclosing a known box of pixel centres averages to the hand-computed nanmean.

    The polygon ``[(0.5, -0.5), (2.5, -0.5), (2.5, 1.5), (0.5, 1.5)]`` encloses pixel centres with
    x ∈ {1, 2} and y ∈ {0, 1} (4 pixels). For each channel the expected mean is the ramp averaged
    over those four pixels."""
    _, store_path = corpus
    vertices = [[0.5, -0.5], [2.5, -0.5], [2.5, 1.5], [0.5, 1.5]]
    body = client.post(
        f"/api/stores/{_store_id(store_path)}/region-spectrum",
        json={"vertices": vertices},
    ).json()

    assert len(body["value"]) == NCHAN
    assert len(body["velocity"]) == NCHAN
    assert body["velocity"] == pytest.approx([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    assert body["n_pixels"] == 4
    assert body["value_unit"] == "K"

    pixels = [(y, x) for y in (0, 1) for x in (1, 2)]
    expected = [
        float(np.mean([_ramp(c, y, x) for (y, x) in pixels])) for c in range(NCHAN)
    ]
    assert body["value"] == pytest.approx(expected)


def test_region_spectrum_degenerate_is_422(client, corpus):
    """A region with < 3 vertices cannot enclose an area → a graceful 422 (not a 500)."""
    _, store_path = corpus
    response = client.post(
        f"/api/stores/{_store_id(store_path)}/region-spectrum",
        json={"vertices": [[1.0, 1.0], [2.0, 2.0]]},
    )
    assert response.status_code == 422
    assert "3 vertices" in response.json()["detail"]


def test_region_spectrum_outside_map_is_empty_not_error(client, corpus):
    """A polygon entirely off the map is an honest empty spectrum (n_pixels 0), not an error."""
    _, store_path = corpus
    far = [[100.0, 100.0], [110.0, 100.0], [110.0, 110.0]]
    body = client.post(
        f"/api/stores/{_store_id(store_path)}/region-spectrum",
        json={"vertices": far},
    ).json()
    assert body["n_pixels"] == 0
    assert body["value"] == [None] * NCHAN


def test_region_spectrum_request_does_not_compute_full_cube(
    client, corpus, monkeypatch
):
    """The region average reads only the polygon's bounding-box columns, never the whole cube."""
    import dask.array as da

    full_cube_size = NCHAN * NY * NX
    computed_sizes = []
    real_compute = da.Array.compute

    def spy_compute(self, *args, **kwargs):
        computed_sizes.append(int(np.prod(self.shape)))
        return real_compute(self, *args, **kwargs)

    monkeypatch.setattr(da.Array, "compute", spy_compute)

    _, store_path = corpus
    vertices = [[0.5, -0.5], [2.5, -0.5], [2.5, 1.5], [0.5, 1.5]]
    response = client.post(
        f"/api/stores/{_store_id(store_path)}/region-spectrum",
        json={"vertices": vertices},
    )
    assert response.status_code == 200
    assert computed_sizes, "expected at least one dask compute for the region slab"
    assert max(computed_sizes) < full_cube_size, (
        f"the region average computed {max(computed_sizes)} elements — the full cube is "
        f"{full_cube_size}; the bounding-box slice was not lazy"
    )


def test_region_spectrum_unknown_store_is_404(client):
    response = client.post(
        f"/api/stores/{_store_id('l1/nope/missing.zarr')}/region-spectrum",
        json={"vertices": [[0, 0], [1, 0], [1, 1]]},
    )
    assert response.status_code == 404


# --------------------------------------------------------------------------------------
# GET /api/stores/{id}/moment0?vmin=&vmax= — the velocity-window moment
# --------------------------------------------------------------------------------------
def test_windowed_moment0_matches_hand_computed_integral(client, corpus):
    """Moment-0 over the window [2, 6] km/s integrates EXACTLY channels 1, 2, 3 (and only those).

    Velocity axis is [0, 2, 4, 6, 8, 10] km/s, so v ∈ [2, 6] selects channels 1, 2, 3. spectral-cube
    integrates value × Δv with Δv = CDELT3 (2000 m/s here), so the windowed moment at (y, x) is
    ``2000 · Σ_{c∈{1,2,3}} ramp(c, y, x)``."""
    _, store_path = corpus
    sid = _store_id(store_path)
    body = client.get(f"/api/stores/{sid}/moment0?vmin=2&vmax=6").json()

    assert body["shape"] == [NY, NX]
    assert body["channel_start"] == 1
    assert body["channel_stop"] == 4  # exclusive
    assert body["vmin_window"] == pytest.approx(2.0)
    assert body["vmax_window"] == pytest.approx(6.0)

    served = np.array(
        [[np.nan if v is None else v for v in row] for row in body["data"]]
    )
    expected = np.array(
        [
            [
                CHANNEL_WIDTH_MS * sum(_ramp(c, y, x) for c in (1, 2, 3))
                for x in range(NX)
            ]
            for y in range(NY)
        ],
        dtype="float64",
    )
    np.testing.assert_allclose(served, expected, rtol=1e-6)


def test_windowed_moment0_differs_from_full_band(client, corpus):
    """The windowed map is a strict subset integral — it must differ from the full-band moment-0."""
    _, store_path = corpus
    sid = _store_id(store_path)
    full = client.get(f"/api/stores/{sid}/moment0").json()
    windowed = client.get(f"/api/stores/{sid}/moment0?vmin=2&vmax=6").json()
    full_arr = np.array(full["data"])
    win_arr = np.array(windowed["data"])
    assert not np.allclose(full_arr, win_arr)
    # The full band integrates all 6 channels; the window integrates 3 → smaller everywhere here.
    assert (win_arr < full_arr).all()
    assert (
        "channel_start" not in full
    )  # the unwindowed payload is unchanged (#67 compatible)


def test_windowed_moment0_empty_window_is_422(client, corpus):
    """A window between channels (no velocity inside) selects nothing → graceful 422."""
    _, store_path = corpus
    response = client.get(
        f"/api/stores/{_store_id(store_path)}/moment0?vmin=2.5&vmax=3.5"
    )
    assert response.status_code == 422
    assert "no channel" in response.json()["detail"]


def test_windowed_moment0_half_window_is_422(client, corpus):
    """One bound without the other is not a window → 422 (a window needs both edges)."""
    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/moment0?vmin=2")
    assert response.status_code == 422


def test_windowed_moment0_reversed_window_is_same(client, corpus):
    """A window dragged right-to-left (vmin > vmax) is the same window (bounds are ordered)."""
    _, store_path = corpus
    sid = _store_id(store_path)
    forward = client.get(f"/api/stores/{sid}/moment0?vmin=2&vmax=6").json()
    reversed_ = client.get(f"/api/stores/{sid}/moment0?vmin=6&vmax=2").json()
    np.testing.assert_allclose(
        np.array(forward["data"]), np.array(reversed_["data"]), rtol=1e-6
    )


def test_windowed_moment0_request_does_not_compute_full_cube(
    client, corpus, monkeypatch
):
    """A windowed moment reads only the slab's channels, never the full cube."""
    import dask.array as da

    full_cube_size = NCHAN * NY * NX
    computed_sizes = []
    real_compute = da.Array.compute

    def spy_compute(self, *args, **kwargs):
        computed_sizes.append(int(np.prod(self.shape)))
        return real_compute(self, *args, **kwargs)

    monkeypatch.setattr(da.Array, "compute", spy_compute)

    _, store_path = corpus
    response = client.get(f"/api/stores/{_store_id(store_path)}/moment0?vmin=2&vmax=6")
    assert response.status_code == 200
    assert computed_sizes, "expected at least one dask compute for the windowed slab"
    assert max(computed_sizes) < full_cube_size, (
        f"the windowed moment computed {max(computed_sizes)} elements — the full cube is "
        f"{full_cube_size}; the channel slab was not lazy"
    )


# --------------------------------------------------------------------------------------
# fsspec parity (PRD #56, story 1): the #68 reductions on a memory:// collection.
# --------------------------------------------------------------------------------------
def test_interactions_work_on_memory_fsspec_collection():
    """Region spectrum + windowed moment serve identically over a memory:// collection."""
    from astrolyze.web.api import create_app

    root_uri, store_path = _build_corpus("memory://ism_corpus_mem")
    client = TestClient(create_app(root_uri))
    sid = _store_id(store_path)

    region = client.post(
        f"/api/stores/{sid}/region-spectrum",
        json={"vertices": [[0.5, -0.5], [2.5, -0.5], [2.5, 1.5], [0.5, 1.5]]},
    ).json()
    assert region["n_pixels"] == 4
    assert len(region["value"]) == NCHAN

    windowed = client.get(f"/api/stores/{sid}/moment0?vmin=2&vmax=6").json()
    assert windowed["shape"] == [NY, NX]
    assert windowed["channel_start"] == 1 and windowed["channel_stop"] == 4
