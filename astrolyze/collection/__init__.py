"""The ``collection`` capability: open a published corpus and browse it (PRD #56, issue #57).

A *collection* is a published corpus — a directory of astrolyze-flavoured Zarr stores indexed by
one ``catalog.parquet`` at its root — opened **read only** for analysis. astrolyze consumes the
catalog as a versioned data-format contract (normatively specified in the ifm repo); the two
codebases never import each other (PRD #56). This package is the consumer side:

- :mod:`astrolyze.collection.catalog` — the deep module owning *all* catalog-format knowledge:
  read + schema-validate ``catalog.parquet``, reject an unknown ``catalog_schema_version``
  explicitly. The :class:`~astrolyze.collection.catalog.CatalogRow` is the typed record everything
  downstream passes around.
- :class:`Collection` — the facade: :meth:`Collection.open` (fsspec, local today / ``s3://``
  later with no API change), :meth:`Collection.list` (object-first overview), and records that
  open into lazy dask-backed Cubes carrying their corpus origin in
  :class:`~astrolyze.io.Metadata`.

This is the keystone tracer the rest of the corpus story builds on: #60 (describe/query), #61
(scan-builder), #62 (covering), #63 (S3 end-to-end). The seams here — the single catalog read
entry point, the typed row, the ``Record.open`` opener, the origin Metadata bridge — are designed
so those extend it without rewrites.
"""

from __future__ import annotations

from . import catalog, scan
from ._facade import Collection, ObjectSummary, Record, StoreDetail
from .catalog import (
    Catalog,
    CatalogRow,
    CatalogSchemaError,
    read_catalog,
)
from .scan import (
    ScanResult,
    ScanWarning,
    build_catalog,
    scan_directory,
    write_catalog,
)
from .stack import HomogeneityReport, Selection, Stack, StackMember

__all__ = [
    "Collection",
    "Record",
    "ObjectSummary",
    "StoreDetail",
    "Catalog",
    "CatalogRow",
    "CatalogSchemaError",
    "read_catalog",
    "catalog",
    "scan",
    "build_catalog",
    "scan_directory",
    "write_catalog",
    "ScanResult",
    "ScanWarning",
    "Stack",
    "StackMember",
    "Selection",
    "HomogeneityReport",
]
