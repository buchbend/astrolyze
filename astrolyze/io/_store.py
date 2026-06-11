"""fsspec URL resolution for the Zarr backend (issue #63, ADR-0006).

One corpus, one call shape — on disk or on S3. The Zarr backend opens and writes stores given a
plain local path, a ``file://`` URL, an in-memory ``memory://`` URL, or an object-store
``s3://`` URL, with no API change between them (PRD #56). This module is the single seam that
turns any of those into the concrete object xarray/zarr opens:

- a **local** target (a bare path or a ``file://`` URL) resolves to a :class:`~pathlib.Path`, so
  the on-disk behaviour — and the ``store / "noise"`` companion-group arithmetic the rest of the
  backend relies on — is byte-for-byte unchanged;
- a **remote** target (``memory://`` / ``s3://`` / any non-local fsspec scheme) resolves to a
  fsspec **mapper** (:meth:`AbstractFileSystem.get_mapper`). xarray opens a Zarr store straight
  from a mapper and the read stays **lazy** (dask chunks pull only the bytes touched), so a 100 GB
  S3 corpus is browsed a chunk at a time rather than downloaded (PRD #56, user story 8).

astrolyze stays thin: fsspec owns the filesystems and the credentials, zarr owns the bytes. The
only thing added here is the local-vs-remote routing and URL-space path joining. Remote caching
is **fsspec passthrough** — a ``simplecache::s3://…`` URL or ``storage_options`` flow through
:func:`fsspec.core.url_to_fs` untouched; astrolyze builds no cache layer of its own.

Imports stay lazy at the call boundary and no object-store driver (``s3fs``) is imported here, so
a bare install keeps working locally; ``s3fs`` is the opt-in :pep:`508` extra ``astrolyze[s3]``
that lights up ``s3://`` only when present (fsspec raises a clear ImportError otherwise).
"""

from __future__ import annotations

from pathlib import Path

from fsspec.core import url_to_fs

# fsspec protocol names that mean "the local filesystem" — these resolve to a Path so the on-disk
# code path (and companion-group ``store / name`` arithmetic) is unchanged. Everything else is a
# remote-style backend handed to xarray as a mapper.
_LOCAL_PROTOCOLS = ("file", "local")


def _is_local_protocol(fs) -> bool:
    """Whether a fsspec filesystem is the local one (its protocol may be a str or a tuple)."""
    protocol = fs.protocol
    names = (protocol,) if isinstance(protocol, str) else tuple(protocol)
    return any(name in _LOCAL_PROTOCOLS for name in names)


def resolve_store(target, **storage_options):
    """Resolve *target* (a path or fsspec URL) to what xarray opens: a Path or a fsspec mapper.

    A local target (bare path or ``file://``) returns a :class:`~pathlib.Path`; any remote scheme
    (``memory://`` / ``s3://`` / …) returns a fsspec ``FSMap`` whose reads stay lazy. *storage_
    options* pass straight through to fsspec (credentials, anon, caching) — astrolyze never
    inspects them (PRD #56: storage options are fsspec's, untouched)."""
    fs, path = url_to_fs(str(target), **storage_options)
    if _is_local_protocol(fs):
        return Path(path)
    return fs.get_mapper(path)


def store_uri(target) -> str:
    """The canonical fsspec URI for *target* (scheme preserved): ``file://…`` / ``s3://…``.

    A bare local path normalises to a ``file://`` URI; a URL keeps its scheme. Used to stamp the
    provenance path a store was opened from, so a saved store names its real home (on disk or on
    S3) rather than a scheme-stripped local guess."""
    fs, path = url_to_fs(str(target))
    return fs.unstrip_protocol(path)


def join_uri(root, name: str) -> str:
    """Join leaf *name* onto a path-or-URL *root* in URL space (POSIX forward-slash).

    Keeps a URL's scheme (``memory://corpus`` -> ``memory://corpus/store.zarr``) and normalises a
    bare local path the same way, so the save path can build a child store URI without collapsing
    a remote scheme to a local Path. The catalog module joins the same way (#57)."""
    text = str(root).rstrip("/")
    return f"{text}/{name}"


def is_local(target) -> bool:
    """Whether *target* names the local filesystem (a bare path or a ``file://`` URL)."""
    fs, _ = url_to_fs(str(target))
    return _is_local_protocol(fs)
