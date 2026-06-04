"""The always-on run log (ADR-0010).

"What ran" is captured automatically — never depending on a human or an agent remembering to
write it. A *run* is one session of work against an experiment: :meth:`RunLog.open` starts it,
creating a per-run append-only **JSONL** file in the experiment's ``logs/`` and installing the
run as the process-level *active run*. While a run is active, the existing astrolyze operations
(``io.load`` / ``io.save`` / ``Cube.moment*`` / ``.to`` / ``.plot``) append one record each
through the single :func:`emit` seam; when **no** run is active :func:`emit` is a silent no-op,
so the library stays fully usable outside an experiment.

Each record is one JSON object — schema version, run id, timestamp, operation, params, inputs
and outputs (paths/ids), the software versions (astrolyze + key deps) and the data version (io
schema version + any manifest ids). JSONL is machine-readable, diff-friendly, and append-only
by construction (we only ever open the file in append mode), so a run is fully reconstructable
and never silently overwritten. The format is kept deliberately flat and UI-readable: a future
frontend (ADR-0009/0010) can render and annotate it without rework.

Two design points worth stating:

- **The emit seam is decoupled from the operations.** Each op calls one tiny function —
  ``runlog.emit(op, …)`` — and knows nothing about JSONL, contextvars, or files. Emission is
  driven entirely by whether a run is active (a :class:`~contextvars.ContextVar`), so wiring a
  new op into the log is a one-line call and removing the log touches only this module.
- **This module is light on purpose.** It pulls only the standard library plus the io schema
  version — never SQLAlchemy or matplotlib — so ``io.load`` importing the seam stays cheap and
  the library's no-side-effect import contract holds (the heavier manifest/ingest siblings are
  imported lazily by :mod:`astrolyze.experiment`).
"""

from __future__ import annotations

import json
import uuid
from contextvars import ContextVar
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from astrolyze.io.schema import SCHEMA_VERSION

if TYPE_CHECKING:
    from .layout import Experiment

#: Run-log *record* format version (independent of the io metadata schema version, which rides
#: in each record's ``data`` block). Bump when the record shape changes.
RUNLOG_SCHEMA_VERSION = 1

#: Distribution names whose versions are stamped on every record (astrolyze + the key deps that
#: actually shape a result). Missing distributions are simply skipped — the set is informative,
#: not a hard dependency contract.
_SOFTWARE_DISTRIBUTIONS = (
    "astrolyze",
    "numpy",
    "astropy",
    "spectral-cube",
    "specutils",
    "radio-beam",
    "matplotlib",
)

# The active run for this process/context, or None when none is open. A ContextVar (not a plain
# global) so the active-run state is well-defined under threading/async and resets cleanly.
_ACTIVE_RUN: ContextVar["RunLog | None"] = ContextVar(
    "astrolyze_active_run", default=None
)

# Software versions are constant for the life of the process; resolve them once, lazily.
_SOFTWARE_VERSIONS_CACHE: dict[str, str] | None = None


class RunLog:
    """One run's append-only JSONL log in an experiment's ``logs/``.

    Open it with :meth:`open` (which also installs it as the active run); it doubles as a
    context manager so ``with RunLog.open(exp) as run:`` deactivates the run on exit. Records
    are written by :func:`emit` (the seam the operations call) or directly via :meth:`record`.
    """

    def __init__(self, path: Path, run_id: str) -> None:
        self.path = path
        self.run_id = run_id
        self._token = None  # ContextVar reset token, set while this run is active

    # -- lifecycle ---------------------------------------------------------------------
    @classmethod
    def open(cls, experiment: "Experiment") -> "RunLog":
        """Start a run for *experiment* and make it the active run.

        Creates ``logs/`` if needed and an empty per-run ``run-<id>.jsonl`` file (so the file
        exists even before the first operation), then installs the run as active so subsequent
        operations emit into it. Each call gets a fresh id and file — separate sessions never
        tangle (ADR-0010)."""
        logs_dir = experiment.logs
        logs_dir.mkdir(parents=True, exist_ok=True)
        run_id = _new_run_id()
        path = logs_dir / f"run-{run_id}.jsonl"
        path.touch()  # the run-log file exists from the moment the run opens, even with no ops
        run = cls(path, run_id)
        run._token = _ACTIVE_RUN.set(run)
        return run

    def close(self) -> None:
        """Deactivate this run (idempotent). Subsequent :func:`emit` calls become no-ops again
        unless another run is active."""
        if self._token is not None:
            _ACTIVE_RUN.reset(self._token)
            self._token = None

    def __enter__(self) -> "RunLog":
        return self

    def __exit__(self, *exc_info) -> bool:
        self.close()
        return False  # never suppress exceptions

    # -- writing / reading -------------------------------------------------------------
    def record(
        self,
        op: str,
        *,
        params: dict | None = None,
        inputs: Any = None,
        outputs: Any = None,
        manifest_ids: Any = None,
    ) -> dict:
        """Append one record for *op* to this run's file and return it.

        Append-only: the file is opened in append mode, so earlier records are never rewritten.
        Inputs/outputs are normalised to a list of reference strings/ids; *params* is JSON-
        encoded with a string fallback so a stray astropy unit/quantity can never break a run.
        """
        record = _build_record(self.run_id, op, params, inputs, outputs, manifest_ids)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, default=str) + "\n")
        return record

    def entries(self) -> list[dict]:
        """Parse this run's file back into a list of records (UI-/inspection-friendly)."""
        return read(self.path)


def read(path: Path | str) -> list[dict]:
    """Parse a run-log JSONL file into its list of records (empty if the file is absent).

    The reusable reader behind :meth:`RunLog.entries` and the narrative offer (#14): a run log
    is plain JSONL, so reading it back is just line-by-line :func:`json.loads`."""
    path = Path(path)
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8")
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def emit(
    op: str,
    *,
    params: dict | None = None,
    inputs: Any = None,
    outputs: Any = None,
    manifest_ids: Any = None,
) -> dict | None:
    """The seam operations call: append a record to the active run, or do nothing.

    Returns the written record while a run is active, or ``None`` when none is — that no-op is
    what keeps astrolyze fully usable outside an experiment (the tracer spine runs unchanged).
    Operations call this and know nothing about JSONL/files/contextvars.
    """
    run = _ACTIVE_RUN.get()
    if run is None:
        return None
    return run.record(
        op, params=params, inputs=inputs, outputs=outputs, manifest_ids=manifest_ids
    )


def active_run() -> "RunLog | None":
    """The run currently active in this context, or ``None`` (the seam :func:`emit` reads)."""
    return _ACTIVE_RUN.get()


# -- record construction ---------------------------------------------------------------
def _build_record(
    run_id: str,
    op: str,
    params: dict | None,
    inputs: Any,
    outputs: Any,
    manifest_ids: Any,
) -> dict:
    """Assemble one run-log record. Flat and UI-readable; values JSON-coercible (paths -> str,
    ids -> int), with the schema/run/timestamp envelope plus software and data versions."""
    return {
        "schema": RUNLOG_SCHEMA_VERSION,
        "run_id": run_id,
        "timestamp": _now_iso(),
        "op": op,
        "params": params or {},
        "inputs": _as_refs(inputs),
        "outputs": _as_refs(outputs),
        "software": _software_versions(),
        "data": {
            "schema_version": SCHEMA_VERSION,
            "manifest_ids": _as_refs(manifest_ids),
        },
    }


def _as_refs(value: Any) -> list:
    """Normalise an input/output reference (or list of them) to a JSON-friendly list: a
    :class:`~pathlib.Path` becomes its string form, an int (a manifest id) stays an int, and a
    bare scalar is wrapped into a one-element list. ``None`` is the empty list."""
    if value is None:
        return []
    if isinstance(value, (str, bytes, Path)) or not isinstance(value, Sequence):
        value = [value]
    return [_as_ref(item) for item in value]


def _as_ref(item: Any):
    if isinstance(item, Path):
        return str(item)
    if isinstance(item, int):  # a manifest id
        return item
    return str(item)


def _software_versions() -> dict[str, str]:
    """The versions of astrolyze and its key deps, resolved once and cached. A distribution that
    is not installed is omitted rather than raising — the stamp is informative, not a contract."""
    global _SOFTWARE_VERSIONS_CACHE
    if _SOFTWARE_VERSIONS_CACHE is None:
        versions: dict[str, str] = {}
        for distribution in _SOFTWARE_DISTRIBUTIONS:
            try:
                versions[distribution] = importlib_metadata.version(distribution)
            except importlib_metadata.PackageNotFoundError:
                continue
        _SOFTWARE_VERSIONS_CACHE = versions
    return dict(_SOFTWARE_VERSIONS_CACHE)


def _new_run_id() -> str:
    """A sortable, collision-resistant run id: a microsecond UTC timestamp plus a short random
    suffix. The microseconds make the id (and so the ``run-<id>.jsonl`` filename) sort
    chronologically — the narrative offer and any future UI can take the lexicographic max as
    the latest run — while the random suffix keeps ids distinct across processes."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%f")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


__all__ = ["RunLog", "emit", "active_run", "read", "RUNLOG_SCHEMA_VERSION"]
