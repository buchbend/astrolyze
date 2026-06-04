"""The experiment layer (ADR-0009/0010).

An *experiment* is the fixed home of one analysis: a predictable skeleton on disk, a merciless
ingest gate, a DB-backed manifest, and an always-on run log. The surface grows as the slices
land — :mod:`layout` (the skeleton, #10), :mod:`manifest` (the dataset registry, #11) and
:mod:`ingest` (the merciless gate, #12) are here now.

Public surface::

    Experiment(root)            # value object resolving the skeleton paths
    Experiment.init(root)       # create the skeleton + default config.toml (idempotent)
    role_of(experiment, path)   # classify a path as raw / interim / processed / output
    Role                        # the four roles
    Manifest(db_url)            # DB-backed dataset registry (register / get / query / all)
    DatasetRecord               # a registered dataset (identity + reconstructed Metadata)
    ingest(experiment)          # validate raw/ and register what is complete (-> IngestReport)
    IngestReport                # accepted / rejected partition of one ingest pass
    AcceptedDataset             # a registered raw file (+ its DatasetRecord)
    RejectedDataset             # a refused raw file (+ the missing fields / read error)

Importing this subpackage pulls dynaconf (for the config seam) and SQLAlchemy (for the
manifest) but not matplotlib or spectral-cube; the CLI still defers the import into the command
body to keep ``--help`` fast.
"""

from __future__ import annotations

from . import ingest as ingest_module
from . import layout, manifest
from .ingest import AcceptedDataset, IngestReport, RejectedDataset, ingest
from .layout import Experiment, Role, role_of
from .manifest import DatasetRecord, Manifest

__all__ = [
    "Experiment",
    "Role",
    "role_of",
    "Manifest",
    "DatasetRecord",
    "ingest",
    "IngestReport",
    "AcceptedDataset",
    "RejectedDataset",
    "layout",
    "manifest",
    "ingest_module",
]
