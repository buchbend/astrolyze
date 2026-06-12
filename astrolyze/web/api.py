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

The cube-viewer slice (#67) — the four basic panels. Each opens the store **lazily** into a
dask-backed :class:`~astrolyze.core.Cube` (via the public :meth:`Record.open`) and serves a
**bounded JSON slice** for client-side D3 (no full-cube load per request, no server-side
matplotlib; the serializers live in :mod:`astrolyze.web._viewer`):

- ``GET /api/stores/{store_id}/cube`` — axis metadata to drive the UI (shape, n_channels, velocity
  axis + unit, sky extent, bunit, value-range hint). Reads headers/axes, not the data plane.
- ``GET /api/stores/{store_id}/moment0`` — the integrated (moment-0) map as a 2-D JSON array +
  extent/unit/vmin/vmax (public :meth:`Cube.moment0`).
- ``GET /api/stores/{store_id}/channel/{index}`` — one channel's 2-D slice + its velocity (a single
  lazy ``cube[index]``). Out-of-range index → ``422``.
- ``GET /api/stores/{store_id}/spectrum?x=&y=`` — the spectrum at a pixel (velocity + value arrays,
  length == n_channels; a single lazy ``cube[:, y, x]``). Out-of-range pixel → ``422``.

The cube-viewer interactions (#68) extend that surface with two server-side reductions, still thin
over the public Cube API (the serializers live in :mod:`astrolyze.web._viewer`):

- ``GET /api/stores/{store_id}/moment0?vmin=&vmax=`` — the **same** moment-0 endpoint, now accepting
  an optional velocity window. With ``vmin``/``vmax`` it recomputes moment-0 over *exactly* the
  channels whose velocity falls in ``[vmin, vmax]`` (a public ``cube[i0:i1]`` slab → ``moment0``);
  without them it is the full-band map unchanged (backward compatible). An empty window (no channel
  inside) → ``422``.
- ``POST /api/stores/{store_id}/region-spectrum`` — body ``{"vertices": [[x, y], …]}`` (image-pixel
  polygon, the coordinates the maps report) → the region-averaged spectrum (``velocity`` + mean
  ``value`` arrays of length ``n_channels`` + ``n_pixels`` averaged). A degenerate region
  (< 3 vertices) → ``422``. POST (not a GET with an encoded polygon) because the vertex list is the
  request payload — many vertices, naturally a JSON body, not a URL.

The reserved ``GET /api/stores/{store_id}/viewer`` seam from #66 indexes the panel feeds (the #68
routes are listed there too so the frontend can discover them without hard-coding paths).

The built Vue frontend is mounted at ``/`` (``astrolyze/web/static/``) so one process serves the
SPA and its API; a request for a non-existent static path falls through to ``index.html`` (SPA
client-side routing), while ``/api/*`` 404s stay JSON.
"""

from __future__ import annotations

import base64
import threading
from collections import OrderedDict

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


def _find_record(collection, store_id: str):
    """The :class:`~astrolyze.collection.Record` for *store_id*, or ``None`` if absent.

    The single store-id → record lookup shared by the record and the viewer endpoints (DRY): decode
    the id to a ``store_path`` and match it against the public :attr:`Collection.records`. Returns
    ``None`` (the caller raises the 404) rather than raising here, so the message stays one place."""
    store_path = _decode_store_id(store_id)
    for record in collection.records:
        if record.store_path == store_path:
            return record
    return None


def create_app(collection_root: str):
    """Build the corpus-explorer FastAPI app over the corpus at *collection_root*.

    *collection_root* is a path or fsspec URL (``s3://…`` works — it rides
    :meth:`Collection.open`'s fsspec layer, so a local and a remote corpus serve identically). The
    collection is opened once here and held on the app; each request reuses the parsed catalog. The
    built frontend is mounted last so ``/api/*`` routes win over the static catch-all.

    The viewer panels share a small **app-scoped** cache (created here as a closure, also exposed on
    ``app.state``) of opened cubes plus the per-store ``cube_axes``/``moment0_map`` results, so a
    viewer session does not re-open the store and re-run the full-cube moment-0 integral on every
    request (each open is several round-trips over ``s3://``). It is app-scoped on purpose, NOT a
    module-level ``lru_cache`` keyed by ``store_id``: two apps built over different corpora can reuse
    the SAME ``store_id`` (a relative ``store_path`` like ``l1/<leaf>`` is identical across a
    ``file://`` and a ``memory://`` corpus), so a process-global key would hand one app the other's
    cube. A fresh cache per ``create_app`` keeps each app's cubes its own.

    fastapi/starlette are imported **inside** this function (not at module top) so importing
    ``astrolyze.web.api`` does not require the web extra — only building the app does. The
    extra-gate has already run by the time the public :func:`astrolyze.web.create_app` reaches
    here, so a missing extra surfaced the helpful error upstream, not an ``ImportError`` here."""
    from fastapi import Body, FastAPI, HTTPException, Query
    from fastapi.responses import FileResponse
    from fastapi.staticfiles import StaticFiles

    from . import _viewer

    collection = Collection.open(collection_root)

    # App-scoped viewer cache: store_id -> entry {"cube", "axes", "moment0"}. Bounded LRU so a
    # long-lived process serving many stores cannot grow without limit; an OrderedDict gives the
    # eviction (move-to-end on hit, popitem(last=False) when full) without a dependency. A single
    # lock guards both the read-then-insert of an entry and the lazy fill of its memoized fields so
    # concurrent threadpool requests (FastAPI runs sync endpoints in a threadpool) cannot double-open
    # or race the dict; the heavy work (open / moment0) is cheap to repeat at worst, so holding the
    # lock across it is acceptable here and keeps the logic simple and race-free.
    _CACHE_MAXSIZE = 16
    _cube_cache: OrderedDict[str, dict] = OrderedDict()
    _cube_cache_lock = threading.Lock()

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
        record = _find_record(collection, store_id)
        if record is None:
            raise HTTPException(
                status_code=404,
                detail=f"no store {_decode_store_id(store_id)!r} in this corpus",
            )
        return _record_json(record)

    # -- cube viewer (#67): open the store lazily, serve bounded JSON slices for D3 --------
    def _cache_entry(store_id: str) -> dict:
        """The cache entry for *store_id* (opening + inserting it on a miss), or 404.

        The single seam every viewer endpoint shares. On a hit the cube is reused (no re-open); on a
        miss the store is opened once via the public :meth:`Record.open` (``Cube.from_zarr`` under
        the hood, dask-backed — nothing materialises until a slice is read) and an entry holding the
        cube plus empty ``axes``/``moment0`` slots is inserted, evicting the LRU tail if full.

        Laziness is preserved (PRD #56, story 8): caching the *opened* cube caches a dask graph, not
        data — slice endpoints still compute exactly their one slice. The lock makes the
        get-or-create atomic so concurrent threadpool requests share one cube rather than racing."""
        with _cube_cache_lock:
            entry = _cube_cache.get(store_id)
            if entry is not None:
                _cube_cache.move_to_end(store_id)  # mark most-recently-used
                return entry
            record = _find_record(collection, store_id)
            if record is None:
                raise HTTPException(
                    status_code=404,
                    detail=f"no store {_decode_store_id(store_id)!r} in this corpus",
                )
            entry = {"cube": record.open(), "axes": None, "moment0": None}
            _cube_cache[store_id] = entry
            if len(_cube_cache) > _CACHE_MAXSIZE:
                _cube_cache.popitem(last=False)  # evict the oldest store
            return entry

    def _open_cube(store_id: str):
        """The cached lazy, dask-backed :class:`Cube` for *store_id* (or 404).

        Slice endpoints (channel/spectrum/region/windowed-moment) take the cube straight: each
        slices it lazily, so only their bounded slice computes — the cube must stay an unrealised
        dask graph, which the cache holds verbatim."""
        return _cache_entry(store_id)["cube"]

    def _cached_axes(store_id: str) -> dict:
        """The memoized :func:`_viewer.cube_axes` result for *store_id* (computed once per store).

        ``cube_axes`` reads a whole mid-channel for a colour-range hint (a single lazy 2-D slice);
        memoizing its result spares that read on repeat ``/cube`` requests within a viewer session.
        Filled under the lock (which also guards the entry's existence) so two concurrent
        first-requests do not both compute it."""
        entry = _cache_entry(
            store_id
        )  # ensures the entry exists + is most-recently-used
        with _cube_cache_lock:
            if entry["axes"] is None:
                entry["axes"] = _viewer.cube_axes(entry["cube"])
            return entry["axes"]

    def _cached_moment0(store_id: str) -> dict:
        """The memoized full-band :func:`_viewer.moment0_map` result for *store_id*.

        moment-0 is the one full-cube integral on the hot path; memoizing its RESULT (the whole
        point of this cache) keeps it one-per-store instead of one-per-request. Only the *full-band*
        map is memoized — windowed moment-0 (``vmin``/``vmax``) varies per request, so it is not
        cached and recomputes its slab each time."""
        entry = _cache_entry(store_id)
        with _cube_cache_lock:
            if entry["moment0"] is None:
                entry["moment0"] = _viewer.moment0_map(entry["cube"])
            return entry["moment0"]

    @app.get(VIEWER_STUB_ROUTE)
    def get_store_viewer(store_id: str) -> dict:
        """The cube-viewer index (#67): the four panel-feed routes for this store.

        The #66 reserved seam, now repointed from a ``501`` to a small route index so the frontend
        (or a script) can discover the panel feeds without hard-coding paths. The #68 panels
        (region-averaged spectrum, velocity-window moment, linked-zoom) will extend this list."""
        base = f"/api/stores/{store_id}"
        return {
            "store_id": store_id,
            "panels": {
                "cube": f"{base}/cube",
                "moment0": f"{base}/moment0",
                "channel": f"{base}/channel/{{index}}",
                "spectrum": f"{base}/spectrum?x={{x}}&y={{y}}",
                # #68 interactions
                "windowed_moment0": f"{base}/moment0?vmin={{vmin}}&vmax={{vmax}}",
                "region_spectrum": f"{base}/region-spectrum",
            },
        }

    @app.get("/api/stores/{store_id}/cube")
    def get_store_cube(store_id: str) -> dict:
        """Axis metadata to drive the viewer UI (shape, velocity axis, sky extent, units, range).

        Reads only headers/axes via :func:`_viewer.cube_axes` — no full data plane — so it is cheap
        even on a huge corpus cube. Memoized per store (the colour-range hint reads one mid-channel
        once, not per request). Unknown store id → ``404``."""
        return _cached_axes(store_id)

    @app.get("/api/stores/{store_id}/moment0")
    def get_store_moment0(
        store_id: str,
        vmin: float | None = Query(
            None,
            description="velocity-window lower bound (cube velocity unit); pair with vmax",
        ),
        vmax: float | None = Query(
            None,
            description="velocity-window upper bound (cube velocity unit); pair with vmin",
        ),
    ) -> dict:
        """The integrated (moment-0) map as a 2-D JSON array + extent/unit/vmin/vmax.

        Dogfoods the public :meth:`Cube.moment0`; the single 2-D result is decimated to a bounded
        payload if large (the response carries the ``downsample`` factor).

        With **both** ``vmin`` and ``vmax`` the map is recomputed over *exactly* the channels whose
        velocity falls in that window (the velocity-window panel, #68): a public ``cube[i0:i1]`` slab
        handed to the same ``moment0``, so the windowed map is the full map's integral restricted to
        the slab — not a reimplemented sum. An empty/invalid window (no channel inside) → ``422``;
        one bound without the other → ``422`` (a window needs both edges). Without either it is the
        full-band map, unchanged (backward compatible with #67)."""
        if vmin is None and vmax is None:
            # The full-band map is the heavy full-cube integral; memoize its result per store.
            return _cached_moment0(store_id)
        if vmin is None or vmax is None:
            raise HTTPException(
                status_code=422,
                detail="a velocity window needs both vmin and vmax (got only one bound)",
            )
        try:
            # Windowed maps vary per request — not memoized; the cached cube keeps them lazy.
            return _viewer.windowed_moment0_map(_open_cube(store_id), vmin, vmax)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.post("/api/stores/{store_id}/region-spectrum")
    def post_store_region_spectrum(
        store_id: str,
        vertices: list[tuple[float, float]] = Body(
            ...,
            embed=True,
            description="polygon vertices as [x, y] image-pixel pairs",
        ),
    ) -> dict:
        """The cube averaged over a polygon region → its region-averaged spectrum (#68).

        Body ``{"vertices": [[x, y], …]}`` — an image-pixel polygon (the coordinates the maps
        report; ``embed=True`` keeps the ``vertices`` key so the body is a JSON object, not a bare
        array). Returns ``velocity`` + mean ``value`` arrays (length == ``n_channels``) + the
        ``n_pixels`` averaged, via :func:`_viewer.region_spectrum` (a thin composition of public
        spatial slicing + nanmean). A degenerate region (< 3 vertices) → ``422``; a polygon enclosing
        no in-map pixel returns an honest empty spectrum (``n_pixels == 0``), not an error."""
        try:
            return _viewer.region_spectrum(_open_cube(store_id), vertices)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/stores/{store_id}/channel/{index}")
    def get_store_channel(store_id: str, index: int) -> dict:
        """One channel's 2-D slice + its velocity — a single lazy ``cube[index]``.

        Out-of-range *index* → ``422`` (a graceful client error naming the valid range), not a 500
        or a wrapped traceback."""
        try:
            return _viewer.channel_slice(_open_cube(store_id), index)
        except IndexError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

    @app.get("/api/stores/{store_id}/spectrum")
    def get_store_spectrum(
        store_id: str,
        x: int = Query(..., description="pixel column (image x)"),
        y: int = Query(..., description="pixel row (image y)"),
    ) -> dict:
        """The spectrum at pixel (x, y) — a single lazy ``cube[:, y, x]`` (length == n_channels).

        Out-of-range pixel → ``422`` (naming the valid x/y ranges)."""
        try:
            return _viewer.pixel_spectrum(_open_cube(store_id), x, y)
        except IndexError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc

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

    # Expose the held collection + the app-scoped viewer cache for tests/introspection (read-only
    # handles; the cache stays per-app — see the cube-cache note above for why it must not be global).
    app.state.collection = collection
    app.state.cube_cache = _cube_cache
    return app


__all__ = ["create_app", "STATIC_DIR"]
