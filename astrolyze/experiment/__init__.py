"""The experiment layer (ADR-0009/0010).

An *experiment* is the fixed home of one analysis: a predictable skeleton on disk, a merciless
ingest gate, a DB-backed manifest, and an always-on run log. The surface grows as the slices
land — :mod:`layout` (the skeleton, #10), :mod:`manifest` (the dataset registry, #11),
:mod:`ingest` (the merciless gate, #12) and :mod:`runlog` (the always-on run log, #13).

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
    RunLog.open(experiment)     # start a run; per-run append-only JSONL in logs/
    emit(op, …)                 # the seam the operations call (no-op when no run is active)
    narrate(experiment)         # scaffold/open the markdown narrative offer (never enforced)

Importing this subpackage pulls dynaconf (for the config seam) and SQLAlchemy (for the
manifest) but not matplotlib or spectral-cube; the CLI still defers the import into the command
body to keep ``--help`` fast. The run-log seam (:mod:`runlog`) is light on its own (stdlib +
the io schema version) — but the operations reach it through ``from astrolyze.experiment.runlog
import emit``, which runs this ``__init__`` and so pays the manifest/ingest import once, on the
first logged operation, not at ``import astrolyze.io`` time (the deferred-import pattern keeps
import-time side-effect-free; see ``io.access._emit`` / ``core/_base._emit``).
"""

from __future__ import annotations

from . import ingest as ingest_module
from . import layout, manifest, narrative, runlog
from .ingest import AcceptedDataset, IngestReport, RejectedDataset, ingest
from .layout import Experiment, Role, role_of
from .manifest import DatasetRecord, Manifest
from .narrative import narrate
from .runlog import RunLog, active_run, emit

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
    "RunLog",
    "emit",
    "active_run",
    "narrate",
    "layout",
    "manifest",
    "ingest_module",
    "runlog",
    "narrative",
]
