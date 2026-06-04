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
    """The result of :func:`load`: the array plus its WCS, parsed schema, and an optional
    live header.

    This is the seam the ``core`` layer (Cube/Map/Spectrum, issue #5) consumes; ``io`` stops
    at parsed bytes + metadata and does not build the data objects itself.

    Backend-neutral by design (issue #22, ADR-0006). A live :class:`~astropy.io.fits.Header`
    is *optional* — a non-FITS backend (Zarr, issue #23) has none. The exact WCS instead
    travels as :attr:`header_string`, the **verbatim FITS-WCS header string**, from which any
    backend reconstructs the WCS with ``WCS(fits.Header.fromstring(header_string))`` (astropy
    still owns the reconstruction). The FITS loader populates both ``header`` and
    ``header_string``; a non-FITS loader populates only ``header_string``.
    """

    data: np.ndarray
    wcs: WCS
    metadata: Metadata
    path: Path
    header_string: str | None = None
    header: fits.Header | None = None


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
    _emit("load", inputs=[path])
    return LoadedData(
        data=data,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path=path,
        # Carry the verbatim FITS-WCS header string so a non-FITS backend can reconstruct the
        # exact WCS without a live fits.Header (ADR-0006); the FITS path also keeps the header.
        header_string=header.tostring(),
        header=header,
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
    _emit("save", params={"extension": extension}, outputs=[out_path])
    return out_path


def _emit(op, **fields) -> None:
    """Append a run-log record for *op* through the always-on run-log seam (ADR-0010).

    The import is deferred to call time so ``io`` stays light — the seam lives in
    :mod:`astrolyze.experiment.runlog`, which pulls only the stdlib + the schema version (never
    SQLAlchemy/matplotlib). When no run is active the seam is a no-op, so ``load``/``save``
    behave identically outside an experiment."""
    from astrolyze.experiment.runlog import emit

    emit(op, **fields)


__all__ = ["LoadedData", "load", "save"]
