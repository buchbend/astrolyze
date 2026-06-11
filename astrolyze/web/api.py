"""The FastAPI app + JSON contracts for the corpus explorer (issue #66).

Every endpoint is a **thin wrapper over the public Collection API**: it calls
:meth:`~astrolyze.collection.Collection.list` / :meth:`~astrolyze.collection.Collection.describe`
/ :attr:`~astrolyze.collection.Collection.records` and serializes the returned dataclasses
field-for-field. It never touches ``_facade`` internals, the parquet reader, or store files
directly — the same public surface a Python user or the CLI uses (dogfooding, PRD #56). The wire
JSON therefore mirrors the Python shape: an ``ObjectSummary`` field maps one-to-one to a JSON key.

The collection is opened **once** at app-construction time (``Collection.open`` reads the catalog
on a local corpus, scans a catalog-less directory) and held on the app, so each request reuses the
parsed in-memory catalog rather than re-reading parquet per call. ``describe`` stays *shallow* (no
``deep=True``): the list+detail slice never opens a store, so a remote corpus explorer is cheap —
the velocity-axis fields the catalog does not carry surface as ``null`` (the same honesty as the
CLI's shallow describe). Opening a store is the cube viewer's job (#67/#68).

JSON-shape discipline:

- ``store_id`` is the URL-safe identity a store is addressed by in the API: the catalog
  ``store_path`` (relative to the corpus root, POSIX, unique per store) percent-quoted so a path
  with slashes is one path segment. :func:`_encode_store_id` / :func:`_decode_store_id` are the
  single round-trip seam (DRY) — the frontend treats it as opaque.
- a NaN/inf float never goes on the wire (it is not valid JSON): :func:`_clean_floats` maps it to
  ``null``, the same "unknown, not guessed" the dataclasses already encode with ``None``.

The endpoints (all under ``/api``):

- ``GET /api/collection`` — corpus root + catalog version + the object-first list (mirrors
  ``Collection.list()``); the list view's data source.
- ``GET /api/objects/{object}`` — one object expanded into its per-store details (mirrors
  ``Collection.describe(object)``); the detail view's data source. Unknown object → 404.
- ``GET /api/stores/{store_id}`` — one store's full catalog row + resolved URI (mirrors a single
  :class:`~astrolyze.collection.Record`); a deep-link target. Unknown id → 404.
- ``GET /api/stores/{store_id}/viewer`` — the reserved #67/#68 cube-viewer seam: a documented
  ``501 Not Implemented`` (not a 404) so the follow-up slices fill a named hook.

The built Vue frontend is mounted at ``/`` (``astrolyze/web/static/``) so one process serves the
SPA and its API; a request for a non-existent static path falls through to ``index.html`` (SPA
client-side routing), while ``/api/*`` 404s stay JSON.
"""

from __future__ import annotations

import base64

from astrolyze.collection import Collection

from . import VIEWER_STUB_ROUTE

# Where the built frontend bundle lands (npm run build emits here; the wheel ships it as package
# data). Resolved at import so a missing build is a clear, early condition rather than a 404 maze.
from pathlib import Path

STATIC_DIR = Path(__file__).parent / "static"


# -- dataclass -> JSON serialization (faithful, field-for-field) ------------------------
def _clean_floats(value):
    """Map a non-finite float (NaN/inf) to ``None`` so the value is valid JSON.

    JSON has no NaN/inf; emitting one yields invalid JSON some clients reject. A non-finite value
    here means "unknown" exactly as ``None`` does in the dataclasses, so collapsing it to ``null``
    is faithful, not lossy (the same honesty rule — ADR-0003, never guess a value)."""
    import math

    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _object_summary_json(summary) -> dict:
    """Serialize an :class:`~astrolyze.collection.ObjectSummary` faithfully (one key per field).

    ``beam_range_arcsec`` is the dataclass's ``(min, max)`` tuple split into named ``beam_min`` /
    ``beam_max`` keys (a JSON object is clearer on the wire than a positional 2-tuple, and each
    half may be ``null`` when no beam is stated). ``surveys`` / ``species`` stay arrays (the
    dataclass holds tuples — JSON has no tuple)."""
    beam_min, beam_max = summary.beam_range_arcsec
    return {
        "object": summary.object,
        "surveys": list(summary.surveys),
        "species": list(summary.species),
        "n_stores": summary.n_stores,
        "beam_min_arcsec": _clean_floats(beam_min),
        "beam_max_arcsec": _clean_floats(beam_max),
    }


def _store_detail_json(detail) -> dict:
    """Serialize a :class:`~astrolyze.collection.StoreDetail` faithfully, plus its ``store_id``.

    Every field is carried straight through (including the deep-only velocity fields, which are
    ``null`` on a shallow describe — never invented). ``store_id`` is added so the detail view can
    deep-link each store to ``GET /api/stores/{store_id}`` and the future viewer seam."""
    return {
        "store_id": _encode_store_id(detail.store_path),
        "object": detail.object,
        "survey": detail.survey,
        "telescope": detail.telescope,
        "species": detail.species,
        "transition": detail.transition,
        "rest_frequency_hz": _clean_floats(detail.rest_frequency_hz),
        "beam_major_arcsec": _clean_floats(detail.beam_major_arcsec),
        "beam_minor_arcsec": _clean_floats(detail.beam_minor_arcsec),
        "beam_pa_deg": _clean_floats(detail.beam_pa_deg),
        "bunit": detail.bunit,
        "store_path": detail.store_path,
        "store_uri": detail.store_uri,
        "content_checksum": detail.content_checksum,
        "channel_width_kms": _clean_floats(detail.channel_width_kms),
        "velocity_min_kms": _clean_floats(detail.velocity_min_kms),
        "velocity_max_kms": _clean_floats(detail.velocity_max_kms),
    }


def _record_json(record) -> dict:
    """Serialize a single :class:`~astrolyze.collection.Record` (its catalog row + resolved URI).

    The store deep-link payload: the full typed :class:`~astrolyze.collection.catalog.CatalogRow`
    field-for-field, plus the record's identity helpers (``store_id`` / ``store_uri``). Read off
    the public ``record.row`` / ``record.store_uri`` — no parquet or facade internals."""
    row = record.row
    return {
        "store_id": _encode_store_id(record.store_path),
        "object": row.object,
        "survey": row.survey,
        "telescope": row.telescope,
        "species": row.species,
        "transition": row.transition,
        "rest_frequency_hz": _clean_floats(row.rest_frequency_hz),
        "beam_major_arcsec": _clean_floats(row.beam_major_arcsec),
        "beam_minor_arcsec": _clean_floats(row.beam_minor_arcsec),
        "beam_pa_deg": _clean_floats(row.beam_pa_deg),
        "bunit": row.bunit,
        "store_path": row.store_path,
        "store_uri": record.store_uri,
        "content_checksum": row.content_checksum,
        "ra_deg": _clean_floats(row.ra_deg),
        "dec_deg": _clean_floats(row.dec_deg),
        "radius_deg": _clean_floats(row.radius_deg),
        "catalog_schema_version": row.catalog_schema_version,
    }


# -- store_id round-trip (the single encode/decode seam) --------------------------------
def _encode_store_id(store_path: str) -> str:
    """The URL-safe id for *store_path* — base64-urlsafe so a slashed path is ONE path segment.

    The catalog ``store_path`` (relative, POSIX, unique per store) is the natural store identity,
    but it carries ``/`` separators. Percent-quoting them is *not* enough: ASGI servers decode
    ``%2F`` back to ``/`` before routing, so a quoted slash would split the path param and miss the
    route. base64-urlsafe encoding has NO ``/`` in its alphabet (it uses ``-``/``_``), so the id is
    a single, route-safe segment with no server-decoding surprises. The frontend treats it as
    opaque. ``=`` padding is kept (urlsafe-safe and the decoder needs it)."""
    return base64.urlsafe_b64encode(store_path.encode()).decode()


def _decode_store_id(store_id: str) -> str:
    """Inverse of :func:`_encode_store_id`: recover the catalog ``store_path`` from a URL id."""
    return base64.urlsafe_b64decode(store_id.encode()).decode()


def create_app(collection_root: str):
    """Build the corpus-explorer FastAPI app over the corpus at *collection_root*.

    *collection_root* is a path or fsspec URL (``s3://…`` works — it rides
    :meth:`Collection.open`'s fsspec layer, so a local and a remote corpus serve identically). The
    collection is opened once here and held on the app; each request reuses the parsed catalog. The
    built frontend is mounted last so ``/api/*`` routes win over the static catch-all.

    fastapi/starlette are imported **inside** this function (not at module top) so importing
    ``astrolyze.web.api`` does not require the web extra — only building the app does. The
    extra-gate has already run by the time the public :func:`astrolyze.web.create_app` reaches
    here, so a missing extra surfaced the helpful error upstream, not an ``ImportError`` here."""
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import FileResponse, JSONResponse
    from fastapi.staticfiles import StaticFiles

    collection = Collection.open(collection_root)

    app = FastAPI(
        title="astrolyze corpus explorer",
        description=(
            "Read-only web explorer for a published astrolyze corpus — object-first list "
            "and per-store detail over the public Collection API (issue #66)."
        ),
        version="0.1",
    )

    @app.get("/api/collection")
    def get_collection() -> dict:
        """The corpus overview + object-first list (mirrors ``Collection.list()``).

        The list view's data source: the corpus root URI, the catalog version, and one entry per
        source object (surveys / species / store count / beam range). Object-first grouping is the
        library's — this endpoint only serializes :meth:`Collection.list`."""
        return {
            "root_uri": collection.root_uri,
            "catalog_version": collection.catalog_version,
            "objects": [_object_summary_json(s) for s in collection.list()],
        }

    @app.get("/api/objects/{object}")
    def get_object(object: str) -> dict:
        """One source expanded into its per-store details (mirrors ``Collection.describe``).

        The detail view's data source. Shallow (no ``deep=True``) so no store is opened — the
        velocity-axis fields the catalog does not carry are ``null`` (the cube viewer #67/#68 will
        read them on demand). An unknown object → ``404`` carrying the library's KeyError message
        (which names the known objects), not an empty 200 (an empty result would hide a typo)."""
        try:
            details = collection.describe(object)
        except KeyError as exc:
            # describe() raises KeyError(message); .args[0] is the clean message (no repr quotes).
            raise HTTPException(status_code=404, detail=exc.args[0]) from exc
        return {
            "object": object,
            "stores": [_store_detail_json(d) for d in details],
        }

    @app.get("/api/stores/{store_id}")
    def get_store(store_id: str) -> dict:
        """One store's full catalog row + resolved URI (mirrors a single ``Record``).

        A deep-link target keyed by the opaque ``store_id`` (the percent-quoted ``store_path``).
        Looked up against the public :attr:`Collection.records` by ``store_path`` — no facade
        internals. Unknown id → ``404`` (the store_path is not in this corpus)."""
        store_path = _decode_store_id(store_id)
        for record in collection.records:
            if record.store_path == store_path:
                return _record_json(record)
        raise HTTPException(
            status_code=404,
            detail=f"no store {store_path!r} in this corpus",
        )

    @app.get(VIEWER_STUB_ROUTE)
    def get_store_viewer(store_id: str):
        """The reserved cube-viewer seam (#67/#68): a documented ``501``, not a ``404``.

        The list+detail slice ships NO viewer logic. This stub gives the follow-up viewer slices a
        named hook (the route exists, returns ``501 Not Implemented`` with a pointer) rather than a
        route they must invent — a clean seam, in the spirit of the facade's documented #60/#61/#62
        extension points."""
        return JSONResponse(
            status_code=501,
            content={
                "detail": (
                    "the cube viewer (integrated map / channel maps / spectra) is not part of "
                    "this slice — it arrives in issues #67/#68. This route is the reserved seam."
                ),
                "store_id": store_id,
            },
        )

    # The built SPA, mounted last so /api/* always wins. html=True makes StaticFiles serve
    # index.html for a directory request; the explicit index route + catch-all below give the
    # SPA client-side routing (an unknown non-/api path returns index.html, not a 404), but only
    # when a build is present — a dev checkout without `npm run build` still serves the API.
    if STATIC_DIR.is_dir() and (STATIC_DIR / "index.html").exists():
        app.mount(
            "/assets",
            StaticFiles(directory=STATIC_DIR / "assets"),
            name="assets",
        ) if (STATIC_DIR / "assets").is_dir() else None

        @app.get("/")
        def index():
            """Serve the built SPA entry point."""
            return FileResponse(STATIC_DIR / "index.html")

        @app.get("/{full_path:path}")
        def spa_fallback(full_path: str):
            """SPA client-side routing: any non-/api path serves index.html.

            A deep link the browser opens directly (e.g. ``/objects/NGC3521``) has no static file,
            so it falls through to index.html and the Vue router takes over. ``/api/*`` is matched
            by the API routes above (registered first), so this never shadows them. A real static
            asset (``/favicon.ico``, ``/vite.svg``) is served if it exists, else index.html."""
            candidate = STATIC_DIR / full_path
            if candidate.is_file():
                return FileResponse(candidate)
            return FileResponse(STATIC_DIR / "index.html")

    # Expose the held collection on the app for tests/introspection (read-only handle).
    app.state.collection = collection
    return app


__all__ = ["create_app", "STATIC_DIR"]
