"""The catalog deep module: all ``catalog.parquet`` format knowledge behind one seam (#57).

A *published corpus* exposes one ``catalog.parquet`` at its root — a flat, per-store index over
the Zarr cubes that live under it. This module is the **only** place in astrolyze that knows that
file's shape: its column set, their order, and the ``catalog_schema_version`` discipline. Read it
through the single entry point :func:`read_catalog`; everything downstream (the :class:`~
astrolyze.collection.Collection` facade, and later the #61 scan-builder and #62 ``covering()``
reader) consumes typed :class:`CatalogRow` records and never touches parquet directly.

The format is the consumer contract normatively specified in the ifm repo
(``docs/catalog-schema.md``); astrolyze reads it, ifm writes it, and the two codebases never
import each other — the versioned spec is the whole interface (PRD #56). The version is
``MAJOR.MINOR``:

- a row order and column set fixed per MAJOR (§4 of the spec): order is part of the contract so a
  MINOR can append columns without breaking an older reader;
- an unknown **MAJOR** is rejected explicitly (:class:`CatalogSchemaError`, never a bare
  ``KeyError``) — a consumer must refuse rather than silently mis-read (ADR-0003, no silent
  guessing);
- a newer **MINOR** of a known MAJOR is accepted: the reader finds every column it knows at its
  position and ignores trailing additions (e.g. the future ``moc`` / ``v_sys_kms`` of a ``1.1``).

All path handling routes through fsspec from day one (:func:`fsspec.open`), so the same call
shape serves a local directory today and an ``s3://`` URL later with no API change (#63). The
**scan-builder** (#61) and **covering()/MOC reader** (#62) extend this module: a builder will add
a ``build_catalog(store_dir)`` producing the same :class:`Catalog`, and covering will read the
footprint columns already parsed here — both reuse :class:`CatalogRow` and the single read path
rather than re-deriving the format.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import fsspec

# The catalog format version this reader is written against. ``MAJOR.MINOR`` — see the module
# docstring for the gating discipline. Bump MINOR for an additive column; MAJOR for a break.
CATALOG_SCHEMA_VERSION = "1.0"

# The filename the index always carries at the corpus root (the spec's fixed container name).
CATALOG_FILENAME = "catalog.parquet"

# The v1.0 column set, in the contract order (the order is part of the contract — a MINOR
# appends, never reorders). A reader maps a row by these names, so a trailing column it has
# never heard of is simply not requested.
CATALOG_COLUMNS = (
    "object",
    "survey",
    "telescope",
    "species",
    "transition",
    "rest_frequency_hz",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "bunit",
    "store_path",
    "content_checksum",
    "ra_deg",
    "dec_deg",
    "radius_deg",
    "catalog_schema_version",
)


class CatalogSchemaError(Exception):
    """Raised when ``catalog.parquet`` declares a ``catalog_schema_version`` astrolyze cannot
    read — an unknown or incompatible MAJOR. Descriptive by contract (the offending version and
    the supported one are named), never a bare ``KeyError`` (ADR-0003: no silent mis-read)."""


@dataclass(frozen=True)
class CatalogRow:
    """One published cube's catalog entry — the typed view of a ``catalog.parquet`` row.

    The stable record other collection code passes around (the #60 query layer filters these,
    #62 covering reads the footprint fields, the :class:`~astrolyze.collection.Collection` facade
    groups them object-first). Every field except ``store_path`` and ``catalog_schema_version``
    may be ``None`` — a partial row is valid and the catalog never invents context (the same
    honesty rule as the store schema). Fields mirror :data:`CATALOG_COLUMNS` one-to-one.
    """

    object: str | None
    survey: str | None
    telescope: str | None
    species: str | None
    transition: str | None
    rest_frequency_hz: float | None
    beam_major_arcsec: float | None
    beam_minor_arcsec: float | None
    beam_pa_deg: float | None
    bunit: str | None
    store_path: str  # relative to the corpus root, POSIX; never null (the locator)
    content_checksum: str | None
    ra_deg: float | None
    dec_deg: float | None
    radius_deg: float | None
    catalog_schema_version: str

    @classmethod
    def from_record(cls, record: dict, *, version: str) -> "CatalogRow":
        """Build a row from a column->value mapping, taking only the fields this reader knows.

        Trailing columns a newer MINOR added (and this reader does not model) are ignored here —
        that is the additive-evolution contract in code: an unknown column never reaches the
        record, so a ``1.1`` file opens cleanly in a ``1.0`` reader.

        The non-null locator (``store_path``) is enforced *here*, at read time: a null is a
        non-conforming catalog and must fail with a descriptive :class:`CatalogSchemaError`
        rather than producing a broken store URI that only blows up opaquely at ``open()``
        (ADR-0003: refuse, never silently mis-read). ``catalog_schema_version`` is stamped from
        the file's gated version, so a per-row value never diverges from the header astrolyze
        validated against."""
        known = {f.name for f in fields(cls)}
        values = {name: record.get(name) for name in known}
        if values.get("store_path") is None:
            raise CatalogSchemaError(
                "catalog row is missing its store_path — the non-null locator the catalog "
                "format requires (a published corpus row always names its store)."
            )
        # The version is the file's gated value, not a per-row column: it cannot disagree with
        # the header astrolyze just version-checked (the spec stamps every row identically).
        values["catalog_schema_version"] = version
        return cls(**values)


@dataclass(frozen=True)
class Catalog:
    """A read ``catalog.parquet``: its rows plus the format version they declared.

    The return of :func:`read_catalog`. ``schema_version`` is the file's
    ``catalog_schema_version`` (already gated as a readable MAJOR); ``rows`` are the typed
    per-store records in file order. The :class:`~astrolyze.collection.Collection` facade is the
    consumer-facing object built on top of this — this stays a thin, format-pure value."""

    rows: tuple[CatalogRow, ...]
    schema_version: str


def read_catalog(root) -> Catalog:
    """Read + schema-validate the ``catalog.parquet`` at corpus *root*; return a :class:`Catalog`.

    *root* is a path or fsspec URL of the corpus root (the directory the Zarr stores live under).
    All access goes through fsspec (:func:`fsspec.open`), so a local directory and an ``s3://``
    prefix share this one path — the S3 case (#63) needs no change here. The version is read from
    the parquet file's key-value metadata when present (so a reader can gate without scanning a
    row) and otherwise from the ``catalog_schema_version`` column.

    Raises :class:`FileNotFoundError` when no catalog exists at *root*, and
    :class:`CatalogSchemaError` (descriptive) when the declared MAJOR is one astrolyze cannot
    read. A newer MINOR of a known MAJOR is accepted (additive evolution)."""
    import pyarrow.parquet as pq

    catalog_url = _join(root, CATALOG_FILENAME)
    # fsspec.open raises FileNotFoundError for a missing local/remote target — surface it as-is
    # (a clear, standard signal the corpus has no index), don't wrap it.
    with fsspec.open(catalog_url, "rb") as handle:
        table = pq.read_table(handle)

    version = _declared_version(table)
    _check_readable(version)

    records = table.to_pylist()
    rows = tuple(CatalogRow.from_record(record, version=version) for record in records)
    return Catalog(rows=rows, schema_version=version)


# -- version gating --------------------------------------------------------------------
def _declared_version(table) -> str:
    """The catalog's ``catalog_schema_version``: from the parquet key-value metadata if stamped
    there (the header path the spec mandates), else from the column. Defaults to the reader's own
    version only when neither is present (a hand-built table) — never invents a *different* one."""
    meta = table.schema.metadata or {}
    raw = meta.get(b"catalog_schema_version")
    if raw is not None:
        return raw.decode() if isinstance(raw, bytes) else str(raw)
    if "catalog_schema_version" in table.column_names:
        column = table.column("catalog_schema_version")
        if len(column) and column[0].is_valid:
            return str(column[0].as_py())
    return CATALOG_SCHEMA_VERSION


def _check_readable(version: str) -> None:
    """Gate *version* against the reader's MAJOR; raise :class:`CatalogSchemaError` on a mismatch.

    A newer MINOR of the supported MAJOR is fine (additive). An unknown MAJOR — or an
    unparseable version — is refused explicitly, naming both the offending and the supported
    version so the message is actionable (ADR-0003)."""
    supported_major = _major(CATALOG_SCHEMA_VERSION)
    try:
        major = _major(version)
    except (ValueError, AttributeError) as exc:
        raise CatalogSchemaError(
            f"catalog declares an unparseable catalog_schema_version {version!r}; "
            f"astrolyze reads catalog format {supported_major}.x (the version must be MAJOR.MINOR)"
        ) from exc
    if major != supported_major:
        raise CatalogSchemaError(
            f"catalog_schema_version {version!r} has an unsupported MAJOR ({major}); "
            f"astrolyze reads catalog format {supported_major}.x. A MAJOR bump is a breaking "
            "change to the catalog format — upgrade astrolyze (or read an older snapshot)."
        )


def _major(version: str) -> int:
    """The MAJOR of a ``MAJOR.MINOR`` version string."""
    return int(str(version).split(".", 1)[0])


# -- fsspec path joining ---------------------------------------------------------------
def _join(root, name: str) -> str:
    """Join *name* onto corpus *root* as an fsspec-resolvable URL.

    A URL keeps its scheme (``s3://bucket/corpus`` -> ``s3://bucket/corpus/catalog.parquet``); a
    bare local path is normalised the same way. POSIX forward-slash throughout — relative store
    paths in the catalog are POSIX by contract, so corpus relocation needs no rewrite."""
    text = str(root).rstrip("/")
    return f"{text}/{name}"


__all__ = [
    "Catalog",
    "CatalogRow",
    "CatalogSchemaError",
    "read_catalog",
    "CATALOG_SCHEMA_VERSION",
    "CATALOG_COLUMNS",
    "CATALOG_FILENAME",
]
