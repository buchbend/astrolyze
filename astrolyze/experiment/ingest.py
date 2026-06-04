"""The merciless ingest gate (ADR-0009).

Ingest is the *strict* counterpart to the lazy :func:`astrolyze.io.load`. It walks an
experiment's ``raw/`` tree, parses each dataset's header through the same ``io``
:class:`~astrolyze.io.schema.Metadata` schema, and partitions every file into:

- **accepted** — the mandatory physical context is present; the file is **registered** into the
  manifest (one row, idempotent on its source path);
- **rejected** — at least one mandatory field is missing; the file is **not** registered, and
  the report names exactly which fields are absent so the user knows what to fix.

That is what *merciless* means: nothing incomplete is ever registered, so a dataset in the
manifest is always safe to compute on. Fixing a header and re-running ingest is the normal,
iterative path — re-ingest updates the same row in place (#11), never duplicating.

Two design points worth stating:

- **``raw/`` is sacred, and reads are header-only.** Ingest never renames or writes anything
  under ``raw/`` (the manifest db lives at the experiment root). It also reads *only* each
  file's header, never the (potentially gigabyte) cube body — scanning a whole ``raw/`` tree
  must stay cheap, and the mandatory context lives entirely in the header.
- **Validation is reused, not reimplemented.** Completeness is decided by the same
  ``REQUIRED_CONTEXT`` / field semantics the schema and lazy loader already use; ingest only
  *applies* that contract as a gate. The required set is overridable — explicitly via
  ``required=`` or per-experiment via ``[ingest] required_context`` in ``config.toml`` — for
  the studies that demand more than the schema baseline (e.g. a beam).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

from astropy.io import fits

from astrolyze.io.schema import REQUIRED_CONTEXT, SCHEMA_FIELDS, Metadata

from .manifest import DatasetRecord, Manifest

if TYPE_CHECKING:
    from .layout import Experiment

# Filename suffixes treated as FITS datasets. Anything else under ``raw/`` (READMEs,
# provenance notes, sidecars) is simply ignored — ingest is about datasets, not every file.
FITS_SUFFIXES = (".fits", ".fits.gz", ".fit", ".fit.gz", ".fts", ".fts.gz")


# -- report value objects --------------------------------------------------------------
@dataclass(frozen=True)
class AcceptedDataset:
    """A ``raw/`` file that passed the gate and was registered.

    :attr:`path` is the absolute location on disk; :attr:`source_path` is its key in the
    manifest (relative to the experiment root, posix); :attr:`record` is the registered
    :class:`~astrolyze.experiment.manifest.DatasetRecord`."""

    path: Path
    source_path: str
    record: DatasetRecord


@dataclass(frozen=True)
class RejectedDataset:
    """A ``raw/`` file the gate refused — and *why*.

    :attr:`missing` lists the mandatory-context fields that are absent (e.g.
    ``["rest_frequency", "velocity_convention"]``). :attr:`error` is set instead when the file
    could not be read as FITS at all (corrupt / not really FITS); such a file is reported, not
    allowed to abort the pass. A rejected dataset is never registered."""

    path: Path
    source_path: str
    missing: list[str] = field(default_factory=list)
    error: str | None = None


@dataclass(frozen=True)
class IngestReport:
    """The outcome of one ingest pass: the accepted and rejected datasets."""

    accepted: list[AcceptedDataset]
    rejected: list[RejectedDataset]

    @property
    def n_accepted(self) -> int:
        return len(self.accepted)

    @property
    def n_rejected(self) -> int:
        return len(self.rejected)


# -- the gate --------------------------------------------------------------------------
def ingest(
    experiment: "Experiment",
    *,
    required: Sequence[str] | None = None,
    manifest: Manifest | None = None,
) -> IngestReport:
    """Validate and register every dataset under *experiment*'s ``raw/``.

    Parameters
    ----------
    experiment:
        The experiment whose ``raw/`` tree is scanned.
    required:
        Mandatory-context field names a file must carry to be accepted. ``None`` (the default)
        resolves the set from ``[ingest] required_context`` in the experiment config, falling
        back to the schema baseline :data:`~astrolyze.io.schema.REQUIRED_CONTEXT`. An unknown
        field name raises :class:`ValueError`.
    manifest:
        The manifest to register accepted datasets into. Defaults to the experiment's own
        (resolved from its config); injectable for tests.

    Returns the :class:`IngestReport`. Accepted files are registered idempotently on their
    source path; rejected files are left untouched and unregistered.
    """
    required = _resolve_required(experiment, required)
    if manifest is None:
        manifest = Manifest.for_experiment(experiment)

    accepted: list[AcceptedDataset] = []
    rejected: list[RejectedDataset] = []
    for path in _raw_datasets(experiment):
        source_path = path.relative_to(experiment.root).as_posix()
        try:
            header = _read_header(path)
        except OSError as exc:
            # A corrupt / non-FITS file must not abort the merciless pass: flag it and move on.
            rejected.append(
                RejectedDataset(path=path, source_path=source_path, error=str(exc))
            )
            continue

        metadata = Metadata.from_header(header)
        missing = _missing(metadata, required)
        if missing:
            rejected.append(
                RejectedDataset(path=path, source_path=source_path, missing=missing)
            )
        else:
            record = manifest.register(metadata, source_path)
            accepted.append(
                AcceptedDataset(path=path, source_path=source_path, record=record)
            )
    return IngestReport(accepted=accepted, rejected=rejected)


# -- helpers ---------------------------------------------------------------------------
def _resolve_required(
    experiment: "Experiment", required: Sequence[str] | None
) -> tuple[str, ...]:
    """Decide the mandatory-context set: explicit *required* wins, then the ``[ingest]`` config
    override, then the schema baseline. Validate every name is a real schema field so a typo is
    a clear error rather than a silent always-reject."""
    if required is None:
        configured = experiment.settings.get("ingest.required_context", None)
        required = configured if configured else REQUIRED_CONTEXT
    required = tuple(required)
    unknown = [name for name in required if name not in SCHEMA_FIELDS]
    if unknown:
        raise ValueError(
            f"unknown required-context field(s): {unknown} "
            f"(choose from {list(SCHEMA_FIELDS)})"
        )
    return required


def _missing(metadata: Metadata, required: Sequence[str]) -> list[str]:
    """Mandatory-context fields absent from *metadata*, in the *required* order. With the
    default set this is exactly :attr:`Metadata.missing`; the function generalises it to an
    overridden required set without reimplementing the schema's notion of 'present'."""
    return [name for name in required if getattr(metadata, name) is None]


def _raw_datasets(experiment: "Experiment") -> list[Path]:
    """Every FITS file under ``raw/`` (recursively), sorted for a deterministic report. An
    absent ``raw/`` yields nothing rather than raising — an empty experiment ingests cleanly."""
    raw = experiment.raw
    if not raw.is_dir():
        return []
    return sorted(p for p in raw.rglob("*") if p.is_file() and _is_fits(p))


def _is_fits(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FITS_SUFFIXES)


def _read_header(path: Path) -> fits.Header:
    """Return the header of the first data-bearing HDU (else the primary), reading **only**
    headers — never the data array. Mirrors :func:`astrolyze.io.load`'s HDU choice via
    ``NAXIS`` so we don't pull a multi-gigabyte raw cube into memory just to validate it."""
    with fits.open(path) as hdul:
        for hdu in hdul:
            if int(hdu.header.get("NAXIS", 0)) > 0:
                return hdu.header.copy()
        return hdul[0].header.copy()


__all__ = [
    "ingest",
    "IngestReport",
    "AcceptedDataset",
    "RejectedDataset",
    "FITS_SUFFIXES",
]
