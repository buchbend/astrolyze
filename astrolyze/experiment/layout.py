"""The fixed experiment skeleton (ADR-0009).

Every astrolyze analysis lives in the same predictable shape on disk, so that work is
portable and legible — and so a coding agent and a human assume the identical layout:

    <experiment>/
      data/{raw,interim,processed}/
      outputs/{figures,tables}/
      logs/
      config.toml

:class:`Experiment` is the value object that resolves those paths; :meth:`Experiment.init`
creates the skeleton and writes a default ``config.toml``; :func:`role_of` classifies where a
path sits (raw / interim / processed / output).

Two design points worth stating:

- **``raw/`` is sacred.** Nothing in this subpackage ever writes to or renames anything under
  ``raw/`` — that would break provenance back to the upstream source. Derivation only ever
  flows raw -> interim -> processed, and the filename projection (ADR-0006) applies to derived
  files only; raw files keep their upstream names.
- **The layout is fixed, the behaviour is configured.** The skeleton path *names* are fixed by
  ADR-0009, so the structure cannot silently drift; they are nonetheless recorded in
  ``config.toml`` for legibility (a human or a future UI can read the layout without knowing
  this module). ``config.toml`` also holds the genuinely-configurable knobs that later slices
  read through dynaconf — the manifest DB URL and an optional required-context override.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

CONFIG_NAME = "config.toml"

# The skeleton directories, as POSIX-relative paths under the experiment root. ``data`` and
# ``outputs`` are created implicitly as parents (mkdir(parents=True)).
SKELETON_DIRS = (
    "data/raw",
    "data/interim",
    "data/processed",
    "outputs/figures",
    "outputs/tables",
    "logs",
)

# Default config.toml written by ``init``. The [paths] block records the fixed layout for
# legibility (it is not a relocation switch in v1 — the skeleton is fixed by ADR-0009);
# [manifest]/[ingest] are the configurable knobs later slices read through dynaconf. The
# manifest db_url is relative and resolves against the experiment root (the manifest slice
# owns that resolution).
DEFAULT_CONFIG = """\
# astrolyze experiment configuration (dynaconf). Written by `astrolyze init`.
# The layout below is the fixed ADR-0009 skeleton, recorded here for legibility.

[paths]
raw = "data/raw"
interim = "data/interim"
processed = "data/processed"
figures = "outputs/figures"
tables = "outputs/tables"
logs = "logs"

[manifest]
# Dataset-manifest database. sqlite to start (swappable); a relative path resolves against
# the experiment root.
db_url = "sqlite:///manifest.db"

[ingest]
# Optional override of the mandatory-context fields the merciless gate requires. Unset = use
# the schema default (io REQUIRED_CONTEXT: rest_frequency, velocity_convention).
# required_context = ["rest_frequency", "velocity_convention", "beam"]
"""


class Role(str, Enum):
    """Where a path sits in the experiment. ``raw`` is sacred (never written by astrolyze);
    derived products live in ``interim``/``processed``; figures and tables are ``output``."""

    RAW = "raw"
    INTERIM = "interim"
    PROCESSED = "processed"
    OUTPUT = "output"


@dataclass(frozen=True)
class Experiment:
    """The fixed experiment skeleton rooted at :attr:`root`.

    A value object: constructing it resolves the skeleton paths but touches no disk. Use
    :meth:`init` to create the skeleton and write a default ``config.toml``.
    """

    root: Path

    def __post_init__(self) -> None:
        # Accept a str or Path; store a Path (frozen, so set via object.__setattr__).
        if not isinstance(self.root, Path):
            object.__setattr__(self, "root", Path(self.root))

    # -- resolved skeleton paths -------------------------------------------------------
    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def raw(self) -> Path:
        return self.data / "raw"

    @property
    def interim(self) -> Path:
        return self.data / "interim"

    @property
    def processed(self) -> Path:
        return self.data / "processed"

    @property
    def outputs(self) -> Path:
        return self.root / "outputs"

    @property
    def figures(self) -> Path:
        return self.outputs / "figures"

    @property
    def tables(self) -> Path:
        return self.outputs / "tables"

    @property
    def logs(self) -> Path:
        return self.root / "logs"

    @property
    def config(self) -> Path:
        return self.root / CONFIG_NAME

    # -- scaffolding -------------------------------------------------------------------
    @classmethod
    def init(cls, root: Path | str) -> "Experiment":
        """Create the skeleton at *root* and write a default ``config.toml``.

        Idempotent: directories are created with ``exist_ok=True`` and an existing
        ``config.toml`` is left untouched (so a re-init never clobbers the user's edits).
        Returns the :class:`Experiment` for the created (or already-present) tree."""
        experiment = cls(root)
        for relative in SKELETON_DIRS:
            (experiment.root / relative).mkdir(parents=True, exist_ok=True)
        if not experiment.config.exists():
            experiment.config.write_text(DEFAULT_CONFIG)
        return experiment

    # -- classification ----------------------------------------------------------------
    def role_of(self, path: Path | str) -> Role | None:
        """Classify *path* as raw / interim / processed / output, or ``None`` when it is
        outside those four trees (e.g. ``logs/``, ``config.toml``, or anything not under this
        experiment). Paths need not exist; classification is purely structural."""
        target = Path(path).resolve()
        for role, base in (
            (Role.RAW, self.raw),
            (Role.INTERIM, self.interim),
            (Role.PROCESSED, self.processed),
            (Role.OUTPUT, self.outputs),
        ):
            if target.is_relative_to(base.resolve()):
                return role
        return None

    # -- configuration -----------------------------------------------------------------
    @property
    def settings(self):
        """The dynaconf settings for this experiment, read from ``config.toml``.

        A fresh object each access (cheap, and avoids caching a stale view after an edit).
        When no ``config.toml`` exists yet, dynaconf returns an empty settings object rather
        than raising — callers use ``.get(key, default)``."""
        from dynaconf import Dynaconf

        return Dynaconf(settings_files=[str(self.config)])


def role_of(experiment: Experiment, path: Path | str) -> Role | None:
    """Module-level mirror of :meth:`Experiment.role_of` — classify *path* within
    *experiment* as raw / interim / processed / output (or ``None``)."""
    return experiment.role_of(path)


__all__ = ["Experiment", "Role", "role_of", "SKELETON_DIRS", "CONFIG_NAME"]
