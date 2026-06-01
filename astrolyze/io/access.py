"""Lazy FITS loading and derived-file writing (ADR-0006).

``load`` reads a FITS file and parses the metadata schema; it is *lazy* — a header missing
mandatory context still opens, with the gap recorded on :attr:`LoadedData.metadata`.
``save`` writes a new file whose name is the header-derived projection; it writes **derived
files only and never touches the raw input** (ADR-0006).

Byte-level I/O and WCS parsing are delegated to ``astropy`` — astrolyze adds the schema, not
a FITS reader.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .naming import project
from .schema import Metadata


@dataclass
class LoadedData:
    """The result of :func:`load`: the array plus its header, WCS, and parsed schema.

    This is the seam the ``core`` layer (Cube/Map/Spectrum, issue #5) consumes; ``io`` stops
    at parsed bytes + metadata and does not build the data objects itself.
    """

    data: np.ndarray
    header: fits.Header
    wcs: WCS
    metadata: Metadata
    path: Path


def load(path) -> LoadedData:
    """Read a FITS file into a :class:`LoadedData`.

    Lazy by contract (ADR-0006 ii): an incomplete header opens and is flagged via
    ``result.metadata.is_complete`` / ``.missing`` — it never raises here. The first HDU that
    carries data is used.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        hdu = next((h for h in hdul if h.data is not None), hdul[0])
        header = hdu.header.copy()
        data = None if hdu.data is None else np.array(hdu.data)
    return LoadedData(
        data=data,
        header=header,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path=path,
    )


def save(
    data: np.ndarray,
    metadata: Metadata,
    directory,
    *,
    base_header: fits.Header | None = None,
    overwrite: bool = False,
    extension: str = "fits",
) -> Path:
    """Write *data* + *metadata* to a new FITS file under *directory*, named by the
    header-derived projection (:func:`~astrolyze.io.naming.project`).

    The schema is written onto a copy of *base_header* when given (preserving its WCS).
    Returns the path written. The source file of a prior :func:`load` is never modified —
    this only ever creates a derived file (ADR-0006).
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    out_path = directory / project(metadata, extension=extension)
    header = metadata.to_header(base_header)
    fits.writeto(out_path, data, header, overwrite=overwrite)
    return out_path


__all__ = ["LoadedData", "load", "save"]
