"""I/O and the metadata schema (ADR-0006).

The FITS header is authoritative; the filename is a derived, browsable projection of it.
Loading is lazy: an incomplete header still opens, but is flagged, and operations needing the
missing context raise (``MissingContextError``) rather than silently guessing (ADR-0003).

Two on-disk backends sit behind the seam (issue #23): **FITS** (eager) and an xarray-native
**Zarr v3** store (lazy, dask-backed). ``load`` dispatches on the target and ``save(...,
format=...)`` selects the writer; both meet the same backend-neutral ``LoadedData`` contract.

Public surface::

    load(path) -> LoadedData          # FITS file or Zarr store -> data + WCS + Metadata (lazy)
    save(data, metadata, directory, format="fits"|"zarr")   # write a derived dataset
    project(metadata) -> str          # the header-derived filename projection
    Metadata                          # the typed header schema (is_complete / ensure_complete)

Byte I/O delegates to astropy (FITS) and xarray/zarr/dask (Zarr); astrolyze adds the schema
and the dispatch, not a storage layer.
"""

from __future__ import annotations

# Re-exported from the unit layer: the I/O contract reuses the same no-silent-physics types
# (the velocity convention and the missing-context error) rather than defining parallel ones.
from astrolyze.units import CalibrationScale, MissingContextError, VelocityConvention

from . import naming, schema
from .access import LoadedData, load, save
from .naming import project
from .schema import Line, Metadata

__all__ = [
    # loading / saving
    "load",
    "save",
    "LoadedData",
    # schema + projection
    "Metadata",
    "Line",
    "project",
    # reused context types
    "VelocityConvention",
    "CalibrationScale",
    "MissingContextError",
    # submodules
    "schema",
    "naming",
]
