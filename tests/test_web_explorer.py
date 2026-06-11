"""Contract tests for the web corpus explorer backend (issue #66).

The explorer's FastAPI endpoints are *thin wrappers* over the public
:class:`~astrolyze.collection.Collection` API. These tests pin the JSON contract each endpoint
serves against a **synthetic on-disk corpus** (the same fixture shape the library/CLI suites use —
two tiny Zarr stores + a ``catalog.parquet``), exercised through FastAPI's ``TestClient``. NO
browser automation, NO real network: TestClient drives the ASGI app in-process.

What is contracted here:

1. ``GET /api/collection`` mirrors :meth:`Collection.list` — the object-first overview (root URI,
   catalog version, one entry per source with surveys/species/store-count/beam-range).
2. ``GET /api/objects/{object}`` mirrors :meth:`Collection.describe` — per-store detail; an
   unknown object is a ``404`` (not an empty 200).
3. ``GET /api/stores/{store_id}`` mirrors a single :class:`~astrolyze.collection.Record`; an
   unknown id is a ``404``.
4. ``GET /api/stores/{store_id}/viewer`` is the reserved #67/#68 cube-viewer seam: a documented
   ``501`` (not a 404), so the follow-up slices fill a named hook.
5. The ``astrolyze[web]`` extra is opt-in: ``astrolyze explore`` (and ``astrolyze.web``'s gate)
   without fastapi installed raises a helpful, actionable error — simulated by monkeypatching the
   import, never by uninstalling fastapi.

The synthetic corpus mirrors the real one: NGC3521 in HERACLES (CO 2-1) and NGC0628 in THINGS
(HI 21cm) — the same two-store fixture the CLI smoke tests use.
"""

import builtins
import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.collection import Collection
from astrolyze.core import Cube

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")
# fastapi is the web extra; without it the backend tests cannot run (but the extra-gate test below
# does not need it — it monkeypatches the import).
fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

REST_CO21 = 230.538e9
REST_HI = 1.420405751e9


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


@pytest.fixture
def corpus(tmp_path):
    """Two stores: NGC3521 in HERACLES (CO 2-1) and NGC0628 in THINGS (HI 21cm)."""
    root = tmp_path / "ism_corpus"
    root.mkdir()
    (root / "l1").mkdir()

    specs = [
        dict(
            obj="NGC3521",
            survey="HERACLES",
            telescope="IRAM30M",
            species="CO",
            transition="2-1",
            rest_hz=REST_CO21,
            bmaj=13.0,
        ),
        dict(
            obj="NGC0628",
            survey="THINGS",
            telescope="VLA",
            species="HI",
            transition="21cm",
            rest_hz=REST_HI,
            bmaj=11.0,
        ),
    ]

    from astrolyze.collection.catalog import CATALOG_COLUMNS

    rows = []
    for spec in specs:
        h, beam = _header(
            obj=spec["obj"],
            telescope=spec["telescope"],
            species=spec["species"],
            rest_hz=spec["rest_hz"],
            bmaj=spec["bmaj"],
        )
        data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
        sub = root / "l1" / spec["survey"].lower()
        sub.mkdir()
        Cube.from_loaded(_loaded(data, h)).to_zarr(sub)
        store = next(sub.glob("*.zarr"))
        rows.append(
            {
                "object": spec["obj"],
                "survey": spec["survey"],
                "telescope": spec["telescope"],
                "species": spec["species"],
                "transition": spec["transition"],
                "rest_frequency_hz": spec["rest_hz"],
                "beam_major_arcsec": beam.major.to_value(u.arcsec),
                "beam_minor_arcsec": beam.minor.to_value(u.arcsec),
                "beam_pa_deg": beam.pa.to_value(u.deg),
                "bunit": "K",
                "store_path": store.relative_to(root).as_posix(),
                "content_checksum": "sha256:x",
                "ra_deg": 170.0,
                "dec_deg": -0.04,
                "radius_deg": 0.01,
                "catalog_schema_version": "1.0",
            }
        )

    table = pa.table({k: [row.get(k) for row in rows] for k in CATALOG_COLUMNS})
    table = table.replace_schema_metadata({"catalog_schema_version": "1.0"})
    pq.write_table(table, root / "catalog.parquet")
    return root


@pytest.fixture
def client(corpus):
    """A TestClient over the explorer app built on the synthetic corpus (in-process ASGI)."""
    from astrolyze.web.api import create_app

    return TestClient(create_app(str(corpus)))


# --------------------------------------------------------------------------------------
# GET /api/collection — object-first list (mirrors Collection.list())
# --------------------------------------------------------------------------------------
def test_collection_endpoint_returns_object_first_overview(client, corpus):
    response = client.get("/api/collection")
    assert response.status_code == 200
    body = response.json()

    # Root URI + catalog version come straight off the public Collection handle.
    expected = Collection.open(str(corpus))
    assert body["root_uri"] == expected.root_uri
    assert body["catalog_version"] == expected.catalog_version

    # The list mirrors Collection.list(): same objects, same order, same aggregates.
    summaries = expected.list()
    assert [o["object"] for o in body["objects"]] == [s.object for s in summaries]

    by_name = {o["object"]: o for o in body["objects"]}
    ngc3521 = by_name["NGC3521"]
    assert ngc3521["surveys"] == ["HERACLES"]
    assert ngc3521["species"] == ["CO"]
    assert ngc3521["n_stores"] == 1
    # Single store → beam_min == beam_max (the catalog's 13" major axis).
    assert ngc3521["beam_min_arcsec"] == pytest.approx(13.0)
    assert ngc3521["beam_max_arcsec"] == pytest.approx(13.0)


def test_collection_endpoint_list_mirrors_collection_list(client, corpus):
    """The endpoint's object set + aggregates are byte-for-byte the library's list()."""
    body = client.get("/api/collection").json()
    summaries = Collection.open(str(corpus)).list()

    for summary, payload in zip(summaries, body["objects"]):
        assert payload["object"] == summary.object
        assert payload["surveys"] == list(summary.surveys)
        assert payload["species"] == list(summary.species)
        assert payload["n_stores"] == summary.n_stores
        lo, hi = summary.beam_range_arcsec
        assert payload["beam_min_arcsec"] == pytest.approx(lo)
        assert payload["beam_max_arcsec"] == pytest.approx(hi)


# --------------------------------------------------------------------------------------
# GET /api/objects/{object} — per-store detail (mirrors Collection.describe())
# --------------------------------------------------------------------------------------
def test_object_endpoint_mirrors_describe(client, corpus):
    response = client.get("/api/objects/NGC3521")
    assert response.status_code == 200
    body = response.json()
    assert body["object"] == "NGC3521"

    details = Collection.open(str(corpus)).describe("NGC3521")
    assert len(body["stores"]) == len(details) == 1
    store = body["stores"][0]
    detail = details[0]

    assert store["survey"] == detail.survey == "HERACLES"
    assert store["telescope"] == detail.telescope == "IRAM30M"
    assert store["species"] == detail.species == "CO"
    assert store["transition"] == detail.transition == "2-1"
    assert store["beam_major_arcsec"] == pytest.approx(detail.beam_major_arcsec)
    assert store["bunit"] == detail.bunit == "K"
    assert store["store_path"] == detail.store_path
    assert store["store_uri"] == detail.store_uri
    assert store["content_checksum"] == detail.content_checksum


def test_object_endpoint_is_shallow_velocity_fields_null(client):
    """The list+detail slice never opens a store: the deep-only velocity fields are null."""
    store = client.get("/api/objects/NGC3521").json()["stores"][0]
    assert store["channel_width_kms"] is None
    assert store["velocity_min_kms"] is None
    assert store["velocity_max_kms"] is None


def test_object_endpoint_unknown_object_is_404(client):
    response = client.get("/api/objects/NGC9999")
    assert response.status_code == 404
    # The library's KeyError message (names the known objects) is carried through.
    assert "NGC9999" in response.json()["detail"]
    assert "NGC3521" in response.json()["detail"]


# --------------------------------------------------------------------------------------
# GET /api/stores/{store_id} — single record (mirrors a Record)
# --------------------------------------------------------------------------------------
def test_store_endpoint_returns_record(client, corpus):
    from astrolyze.web.api import _encode_store_id

    record = Collection.open(str(corpus)).records[0]
    store_id = _encode_store_id(record.store_path)
    response = client.get(f"/api/stores/{store_id}")
    assert response.status_code == 200
    body = response.json()
    assert body["store_path"] == record.store_path
    assert body["store_uri"] == record.store_uri
    assert body["object"] == record.object
    # The footprint columns from the catalog row are carried through.
    assert body["ra_deg"] == pytest.approx(170.0)
    assert body["catalog_schema_version"] == "1.0"


def test_store_endpoint_unknown_id_is_404(client):
    from astrolyze.web.api import _encode_store_id

    response = client.get(f"/api/stores/{_encode_store_id('l1/nope/missing.zarr')}")
    assert response.status_code == 404


# --------------------------------------------------------------------------------------
# GET /api/stores/{store_id}/viewer — reserved #67/#68 seam (501, not 404)
# --------------------------------------------------------------------------------------
def test_viewer_seam_is_not_implemented_not_404(client, corpus):
    """The cube-viewer route exists as a documented 501 seam, so #67/#68 fill a named hook."""
    from astrolyze.web.api import _encode_store_id

    record = Collection.open(str(corpus)).records[0]
    store_id = _encode_store_id(record.store_path)
    response = client.get(f"/api/stores/{store_id}/viewer")
    assert response.status_code == 501
    body = response.json()
    assert "67" in body["detail"] and "68" in body["detail"]
    assert body["store_id"] == store_id


# --------------------------------------------------------------------------------------
# astrolyze[web] is opt-in: a bare install errors helpfully (simulated import failure)
# --------------------------------------------------------------------------------------
def test_require_web_extra_raises_helpful_error_without_fastapi(monkeypatch):
    """Without fastapi importable, the gate raises an actionable error naming the extra.

    Simulated by monkeypatching __import__ to fail for fastapi/uvicorn — never by uninstalling the
    package — so the test is hermetic and reversible."""
    import astrolyze.web as web

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name in ("fastapi", "uvicorn"):
            raise ImportError(f"No module named {name!r}")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(web.WebExtraNotInstalled) as excinfo:
        web.require_web_extra()
    message = str(excinfo.value)
    assert "astrolyze[web]" in message
    assert "pip install" in message


def test_cli_explore_without_extra_exits_with_helpful_error(monkeypatch, corpus):
    """`astrolyze explore` on a bare install exits non-zero with the install hint, not a traceback.

    The web extra IS installed in the test venv, so we simulate its absence by patching the gate to
    raise as it would without fastapi. The CLI must catch it, print the hint, and exit code 1."""
    from typer.testing import CliRunner

    import astrolyze.web as web
    from astrolyze.cli import app

    def boom():
        raise web.WebExtraNotInstalled(web.WEB_EXTRA_HINT)

    # serve() gates via require_web_extra first; patch it to simulate the missing extra.
    monkeypatch.setattr(web, "require_web_extra", boom)

    runner = CliRunner()
    result = runner.invoke(app, ["explore", str(corpus)])
    assert result.exit_code == 1
    out = (result.stdout or "") + (result.stderr if hasattr(result, "stderr") else "")
    assert "astrolyze[web]" in out
