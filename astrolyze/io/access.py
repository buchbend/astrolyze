"""Backend dispatch + lazy FITS loading and derived-file writing (ADR-0006).

``load`` reads a dataset and parses the metadata schema; it is *lazy* — a header missing
mandatory context still opens, with the gap recorded on :attr:`LoadedData.metadata`.
``save`` writes a new derived dataset and **never touches the raw input** (ADR-0006).

Two on-disk backends sit behind this one seam (issue #23): **FITS** (eager, the default) and
an xarray-native **Zarr v3** store (lazy, dask-backed). ``load`` dispatches on the target — a
``.fits`` / ``.fits.gz`` path goes to the FITS reader, a Zarr store to the Zarr reader — and
``save(..., format="fits"|"zarr")`` selects the writer. Both return / consume the same
backend-neutral :class:`LoadedData`, whose contract (the verbatim FITS-WCS header string as the
WCS vehicle, an optional live ``fits.Header``) is what lets a non-FITS backend slot in (#22).

Byte-level I/O and WCS parsing are delegated to ``astropy`` / ``xarray`` / ``zarr`` — astrolyze
adds the schema and the dispatch, not a storage layer.
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


# FITS path suffixes the FITS backend claims; anything else (a directory store) is Zarr.
_FITS_SUFFIXES = (".fits", ".fits.gz", ".fit", ".fts")


def load(path) -> LoadedData:
    """Read a dataset into a :class:`LoadedData`, dispatching FITS vs Zarr on the target.

    A ``.fits`` / ``.fits.gz`` path routes to the eager FITS reader; anything else (a Zarr
    store directory) routes to the lazy, dask-backed Zarr reader (issue #23). Both return a
    :class:`LoadedData`.

    Lazy by contract (ADR-0006 ii): an incomplete header opens and is flagged via
    ``result.metadata.is_complete`` / ``.missing`` — it never raises here.
    """
    path = Path(path)
    if _is_fits_path(path):
        return _load_fits(path)
    from .zarr_backend import _load_zarr

    return _load_zarr(path)


def _is_fits_path(path: Path) -> bool:
    """Whether *path* names a FITS file (by suffix). A two-part ``.fits.gz`` counts."""
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in _FITS_SUFFIXES)


def _load_fits(path: Path) -> LoadedData:
    """Read a FITS file into a :class:`LoadedData` (the eager backend).

    The first HDU that carries data is used."""
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
    format: str = "fits",
    base_header: fits.Header | None = None,
    overwrite: bool = False,
    extension: str = "fits",
    chunks=None,
    shards=None,
    compressors=None,
) -> Path:
    """Write *data* + *metadata* to a new derived dataset under *directory*; return its path.

    *format* selects the backend (issue #23): ``"fits"`` (default, eager) writes a FITS file
    named by the header-derived projection (:func:`~astrolyze.io.naming.project`); ``"zarr"``
    writes an xarray-native Zarr v3 store under the same name (a directory). An unknown *format*
    raises rather than guessing (ADR-0003). The source file of a prior :func:`load` is never
    modified — this only ever creates a derived dataset (ADR-0006).

    The schema is written onto a copy of *base_header* when given (FITS) or carried as attrs +
    the verbatim WCS header string (Zarr), preserving the WCS either way. *chunks* / *shards* /
    *compressors* are the caller's Zarr layout choices, passed straight through to zarr — there
    is no hard-coded chunking policy (they are ignored by the FITS path).
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    if format == "fits":
        out_path = directory / project(metadata, extension=extension)
        header = metadata.to_header(base_header)
        fits.writeto(out_path, data, header, overwrite=overwrite)
        _emit(
            "save",
            params={"format": "fits", "extension": extension},
            outputs=[out_path],
        )
        return out_path
    if format == "zarr":
        from .zarr_backend import _save_zarr

        store = directory / project(metadata, extension="zarr")
        return _save_zarr(
            data,
            metadata,
            store,
            base_header=base_header,
            chunks=chunks,
            shards=shards,
            compressors=compressors,
            overwrite=overwrite,
        )
    raise ValueError(
        f"unknown save format {format!r}: astrolyze writes 'fits' or 'zarr' "
        "(it never guesses a storage backend — ADR-0003)"
    )


def _emit(op, **fields) -> None:
    """Append a run-log record for *op* through the always-on run-log seam (ADR-0010).

    The import is deferred to call time so ``io`` stays light — the seam lives in
    :mod:`astrolyze.experiment.runlog`, which pulls only the stdlib + the schema version (never
    SQLAlchemy/matplotlib). When no run is active the seam is a no-op, so ``load``/``save``
    behave identically outside an experiment."""
    from astrolyze.experiment.runlog import emit

    emit(op, **fields)


__all__ = ["LoadedData", "load", "save"]
