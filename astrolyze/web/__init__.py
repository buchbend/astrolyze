"""The ``web`` capability: serve a published corpus as a browsable web app (PRD #56, issue #66).

The **corpus explorer**: ``astrolyze explore <collection>`` starts a local FastAPI server whose
endpoints are *thin wrappers* over the public :class:`~astrolyze.collection.Collection` API
(:meth:`~astrolyze.collection.Collection.list` / :meth:`~astrolyze.collection.Collection.describe`
/ records), and serves a single-page Vue frontend that renders the object-first **list** view and
the per-store **detail** view. It is the GUI sibling of the ``astrolyze collection …`` CLI group —
same read-only public API underneath, a browser on top instead of a rich table.

This is an **optional extra** (``astrolyze[web]`` → fastapi + uvicorn). The dependency is opt-in,
not a hard requirement: a bare install never imports fastapi/uvicorn (the heavy server stack is
irrelevant to the analysis library), and ``astrolyze explore`` without the extra fails with a
single actionable line — :func:`require_web_extra` — rather than an opaque ``ImportError``. The
imports of fastapi/uvicorn stay **lazy** (inside :func:`create_app` / :func:`serve`) so importing
``astrolyze.web`` itself, and listing CLI ``--help``, costs nothing.

Scope (issues #66 + #67): the explorer ships the object-first **list** + per-store **detail** views
(#66) and the cube viewer's **three basic panels** (#67) — an integrated (moment-0) map, a channel
map with a velocity slider + keyboard stepping, and a pixel spectrum at the clicked position. The
viewer's backend (:mod:`astrolyze.web._viewer`) opens each store *lazily* (dask-backed
:class:`~astrolyze.core.Cube`) and serves bounded JSON slices for client-side D3 — no full-cube load
per request. The remaining viewer panels (region-averaged spectrum, velocity-window moment,
linked-zoom) are #68; clean seams are left for them (``VIEWER_STUB_ROUTE`` is now the panel index,
and the Vuex viewer state reserves the slots) but no #68 logic lives here.

Layering, deliberately thin (the same "astrolyze stays thin" rule as the rest of the package):

- :mod:`astrolyze.web.api` builds the FastAPI app and the JSON contracts. Every endpoint
  *dogfoods* the public Collection API — it never reaches into ``_facade`` internals or the
  parquet reader. The dataclasses (:class:`~astrolyze.collection.ObjectSummary`,
  :class:`~astrolyze.collection.StoreDetail`, :class:`~astrolyze.collection.catalog.CatalogRow`)
  are serialized faithfully (field-for-field) so the wire shape mirrors the Python shape.
- :mod:`astrolyze.web.server` is the uvicorn launcher the CLI ``explore`` command calls.
- ``astrolyze/web/static/`` is where the **built** Vue frontend lands (``frontend/`` is the
  source; ``npm run build`` emits the bundle here, and the wheel ships it as package data). The
  app mounts this directory at ``/`` so one process serves both the API and the UI.
"""

from __future__ import annotations

# The single actionable install line surfaced everywhere the extra is missing — defined once so
# the CLI and any future entry point share the exact wording (DRY).
WEB_EXTRA_HINT = (
    "the web explorer needs the optional 'web' extra (fastapi + uvicorn). "
    "Install it with:  pip install 'astrolyze[web]'"
)

# The cube-viewer route base. #66 reserved it as a documented 501 seam; #67 repoints it to a small
# panel-index endpoint (the four panel-feed routes for a store) so the frontend can discover the
# feeds without hard-coding paths. The #68 panels extend that index. The name is kept stable so the
# api module and the tests share one constant (DRY).
VIEWER_STUB_ROUTE = "/api/stores/{store_id}/viewer"


class WebExtraNotInstalled(RuntimeError):
    """Raised when the corpus explorer is invoked without the ``astrolyze[web]`` extra.

    Carries the single actionable install line (:data:`WEB_EXTRA_HINT`) — never a bare
    ``ImportError`` traceback. The CLI catches it and prints the hint, so a user who typed
    ``astrolyze explore`` on a bare install learns exactly what to install in one line (the same
    opt-in-extra contract as ``astrolyze[s3]`` / ``astrolyze[coverage]``)."""


def require_web_extra() -> None:
    """Confirm fastapi + uvicorn are importable; raise :class:`WebExtraNotInstalled` if not.

    The single gate the ``explore`` command (and any future web entry point) calls before doing
    anything that needs the server stack. The check is an actual import attempt — not a guess from
    a version table — so it is true to what the environment can run. Lazy by design: it imports
    fastapi/uvicorn only when called, never at ``import astrolyze.web`` time, so a bare install is
    unaffected and ``astrolyze --help`` stays fast."""
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401
    except ImportError as exc:
        raise WebExtraNotInstalled(WEB_EXTRA_HINT) from exc


def create_app(collection_root: str):
    """Build the FastAPI corpus-explorer app for the corpus at *collection_root*.

    A thin re-export of :func:`astrolyze.web.api.create_app` with the extra-gate applied first, so
    a caller (the CLI, or a test) gets the helpful error rather than an ``ImportError`` when the
    web extra is absent. The fastapi import lives in :mod:`astrolyze.web.api`, reached only past the
    gate — keeping ``astrolyze.web`` import side-effect free."""
    require_web_extra()
    from .api import create_app as _create_app

    return _create_app(collection_root)


def serve(collection_root: str, *, host: str = "127.0.0.1", port: int = 8000) -> None:
    """Serve the corpus explorer for *collection_root* on *host*:*port* (blocking, via uvicorn).

    The entry point the CLI ``explore`` command calls. Gates on the web extra first (helpful error
    if absent), then hands the built app to uvicorn. uvicorn is imported lazily inside
    :mod:`astrolyze.web.server` so a bare install never pulls it."""
    require_web_extra()
    from .server import serve as _serve

    _serve(collection_root, host=host, port=port)


__all__ = [
    "WEB_EXTRA_HINT",
    "VIEWER_STUB_ROUTE",
    "WebExtraNotInstalled",
    "require_web_extra",
    "create_app",
    "serve",
]
