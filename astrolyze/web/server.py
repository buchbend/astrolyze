"""The uvicorn launcher for the corpus explorer (issue #66).

The blocking ``serve(...)`` the CLI ``explore`` command calls: build the FastAPI app over the
corpus and hand it to uvicorn. uvicorn is imported **inside** :func:`serve` (lazy) so importing
this module never requires the web extra — only actually serving does. The extra-gate has run
upstream (``astrolyze.web.serve`` calls :func:`~astrolyze.web.require_web_extra` first), so a
missing extra surfaced the helpful error before this module's uvicorn import is reached.
"""

from __future__ import annotations

from .api import create_app


def serve(collection_root: str, *, host: str = "127.0.0.1", port: int = 8000) -> None:
    """Serve the corpus explorer for *collection_root* on *host*:*port* (blocking).

    Builds the app (which opens the collection once) and runs it under uvicorn until interrupted.
    Defaults to loopback (``127.0.0.1``) — a local analysis tool, not a public server; bind
    ``0.0.0.0`` explicitly via the CLI ``--host`` to expose it. The app is constructed eagerly here
    (not passed as an import string) so a bad corpus root fails fast with the library's clear error
    before uvicorn binds the socket."""
    import uvicorn

    app = create_app(collection_root)
    uvicorn.run(app, host=host, port=port)


__all__ = ["serve"]
