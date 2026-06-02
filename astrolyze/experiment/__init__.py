"""The experiment layer (ADR-0009/0010).

An *experiment* is the fixed home of one analysis: a predictable skeleton on disk, a merciless
ingest gate, a DB-backed manifest, and an always-on run log. This slice (issue #10) ships the
first piece — the :mod:`layout` module — and the surface grows as the later slices land.

Public surface::

    Experiment(root)            # value object resolving the skeleton paths
    Experiment.init(root)       # create the skeleton + default config.toml (idempotent)
    role_of(experiment, path)   # classify a path as raw / interim / processed / output
    Role                        # the four roles

Importing this subpackage pulls dynaconf (for the config seam) but not matplotlib or
spectral-cube; the CLI still defers the import into the command body to keep ``--help`` fast.
"""

from __future__ import annotations

from . import layout
from .layout import Experiment, Role, role_of

__all__ = ["Experiment", "Role", "role_of", "layout"]
