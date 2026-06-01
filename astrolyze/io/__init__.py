"""I/O and the metadata schema (ADR-0006).

The FITS header is authoritative; the filename is a derived, browsable projection of it.
Loading is lazy: an incomplete header still opens, but is flagged, and operations needing the
missing context raise (``MissingContextError``) rather than silently guessing (ADR-0003).

Public surface::

    load(path) -> LoadedData          # data + header + WCS + parsed Metadata (lazy)
    save(data, metadata, directory)   # write a derived file, named from the header
    project(metadata) -> str          # the header-derived filename projection
    Metadata                          # the typed header schema (is_complete / ensure_complete)

Byte I/O delegates to astropy; astrolyze adds the schema, not a FITS reader.
"""

from __future__ import annotations

# Re-exported from the unit layer: the I/O contract reuses the same no-silent-physics types
# (the velocity convention and the missing-context error) rather than defining parallel ones.
from astrolyze.units import MissingContextError, VelocityConvention

from . import naming, schema
from .access import LoadedData, load, save
from .naming import project
from .schema import Metadata

__all__ = [
    # loading / saving
    "load",
    "save",
    "LoadedData",
    # schema + projection
    "Metadata",
    "project",
    # reused context types
    "VelocityConvention",
    "MissingContextError",
    # submodules
    "schema",
    "naming",
]
