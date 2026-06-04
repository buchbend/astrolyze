"""The xarray-native Zarr backend (issue #23, ADR-0006).

A **second on-disk backend alongside FITS** — not a replacement. astrolyze stays thin: the
store *is* an xarray :class:`~xarray.Dataset` written to Zarr v3, so storage, coordinates,
chunking, sharding and compression are xarray/zarr/dask's job. astrolyze adds only

- the ``Metadata`` <-> ``attrs`` projection (issue #22): the physical context travels in the
  group's ``attrs`` as JSON primitives via :meth:`Metadata.to_attrs` / :meth:`from_attrs`;
- the **verbatim FITS-WCS header string** as the WCS round-trip vehicle (issue #22): any
  backend reconstructs the exact WCS with ``WCS(fits.Header.fromstring(s))`` (astropy owns the
  reconstruction). The Zarr store has no live ``fits.Header`` — only this string;
- a small **provenance** marker recording the backend.

The data variable is ``intensity`` with dims ``(sky_y, sky_x, freq)`` for a cube (a 2D map
drops ``freq``). The on-disk array axis order differs from the FITS array order
``(spectral, y, x)``; we transpose at the boundary so ``LoadedData.data`` keeps the FITS axis
order the WCS expects — the store layout is xarray-natural, the in-memory contract is unchanged.

Zarr cubes are **dask-backed (lazy)**: :func:`_load_zarr` opens with dask chunks so constructing
a cube and slicing a subcube never materialise the full array (FITS stays eager). Chunk / shard
/ compression are **caller parameters** on save — astrolyze hard-codes no chunking policy
(ADR-0003: no silent decisions about layout either).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from astropy.io import fits
from astropy.wcs import WCS

from .schema import Metadata

if (
    TYPE_CHECKING
):  # for type hints only; the real class lives in io.access (circular at runtime)
    from .access import LoadedData

# The single data variable. Dim order is xarray-natural (sky first, spectral last); the FITS
# numpy order is (spectral, y, x), so we transpose at the save/load boundary.
DATA_VAR = "intensity"
SPECTRAL_DIM = "freq"
SKY_DIMS = ("sky_y", "sky_x")
CUBE_DIMS = (*SKY_DIMS, SPECTRAL_DIM)

# attrs keys astrolyze owns on the group (beside the Metadata.to_attrs() projection).
ATTR_WCS_HEADER = "fits_wcs_header"
ATTR_PROVENANCE = "provenance"
PROVENANCE_BACKEND = "zarr"


def _save_zarr(
    data: np.ndarray,
    metadata: Metadata,
    store,
    *,
    header_string: str | None = None,
    base_header: fits.Header | None = None,
    chunks=None,
    shards=None,
    compressors=None,
    overwrite: bool = False,
) -> Path:
    """Write *data* + *metadata* to a Zarr v3 store at *store*; return its path.

    The store is an xarray ``Dataset`` whose group ``attrs`` carry the schema projection
    (:meth:`Metadata.to_attrs`), the verbatim FITS-WCS header string, and a provenance marker.
    *chunks* / *shards* / *compressors* are passed straight through to zarr (per data var) — the
    caller owns the layout. A spectral WCS, when present, fills the authoritative ``freq`` coord.
    """
    store = Path(store)
    # Prefer the explicit verbatim string; fall back to a base header (the FITS save path passes
    # one). Either is the authoritative WCS vehicle — astrolyze never re-derives the WCS.
    if header_string is None and base_header is not None:
        header_string = base_header.tostring()

    dataset = _to_dataset(data, metadata, header_string)
    # Caller layout (chunks/shards) is given in the *data* axis order (spectral, y, x), to match
    # LoadedData.data; the store axis order is (sky_y, sky_x, freq), so translate before zarr.
    ndim = np.asarray(data).ndim
    encoding = _encoding(
        chunks=_to_store_axes(chunks, ndim),
        shards=_to_store_axes(shards, ndim),
        compressors=compressors,
    )
    dataset.to_zarr(
        store,
        mode="w" if overwrite else "w-",
        zarr_format=3,
        encoding=encoding,
        consolidated=False,
    )
    _emit("save", params={"format": "zarr"}, outputs=[store])
    return store


def _load_zarr(store) -> "LoadedData":
    """Open a Zarr v3 store into a **dask-backed (lazy)** :class:`LoadedData`.

    ``chunks={}`` makes xarray hand back dask arrays sized to the on-disk chunks, so neither
    constructing the cube nor slicing a subcube reads the whole array. The schema is rebuilt
    from the group ``attrs`` (:meth:`Metadata.from_attrs`) and the exact WCS from the verbatim
    header string (astropy owns the reconstruction). No live ``fits.Header`` exists here.
    """
    from .access import (
        LoadedData,
    )  # local import: avoid an io.access <-> zarr_backend cycle

    store = Path(store)
    dataset = xr.open_zarr(store, zarr_format=3, consolidated=False, chunks={})
    array = dataset[DATA_VAR]
    # Restore the FITS array order (spectral, y, x) the WCS expects; stays lazy (dask transpose).
    data = _to_fits_order(array).data

    attrs = dict(dataset.attrs)
    header_string = attrs.get(ATTR_WCS_HEADER)
    wcs = (
        WCS(fits.Header.fromstring(header_string))
        if header_string is not None
        else WCS(naxis=array.ndim)
    )
    _emit("load", inputs=[store])
    return LoadedData(
        data=data,
        wcs=wcs,
        metadata=Metadata.from_attrs(attrs),
        path=store,
        header_string=header_string,
        header=None,  # a non-FITS backend has no live fits.Header (issue #22)
    )


# -- dataset <-> array projection ------------------------------------------------------
def _to_dataset(
    data: np.ndarray, metadata: Metadata, header_string: str | None
) -> xr.Dataset:
    """Project *data* + *metadata* onto an xarray ``Dataset`` (the on-disk shape)."""
    arr = np.asarray(data)
    if arr.ndim == 3:
        # FITS order (spectral, y, x) -> store order (sky_y, sky_x, freq).
        values = np.transpose(arr, (1, 2, 0))
        dims = CUBE_DIMS
    elif arr.ndim == 2:
        values = arr
        dims = SKY_DIMS
    else:
        # 1D / other: store as-is under the spectral dim so the round-trip stays lossless.
        values = arr
        dims = (SPECTRAL_DIM,) * arr.ndim

    coords = _coords(dims, values.shape, header_string)
    attrs = metadata.to_attrs()
    if header_string is not None:
        attrs[ATTR_WCS_HEADER] = header_string
    attrs[ATTR_PROVENANCE] = {"backend": PROVENANCE_BACKEND}
    return xr.Dataset({DATA_VAR: (dims, values)}, coords=coords, attrs=attrs)


def _coords(dims, shape, header_string: str | None) -> dict:
    """Build 1D coordinate arrays for *dims*.

    Sky axes carry pixel-index coords — astrolyze does not fabricate separable sky world
    coordinates for a non-separable projection (that would be silent physics, ADR-0003); the
    verbatim WCS remains the authoritative 2D sky mapping. The spectral axis carries its
    *authoritative* world values, read straight from the spectral sub-WCS when one is present.
    """
    coords: dict = {}
    for dim, size in zip(dims, shape):
        if dim == SPECTRAL_DIM:
            coords[dim] = _spectral_world(size, header_string)
        else:
            coords[dim] = np.arange(size)
    return coords


def _spectral_world(size: int, header_string: str | None) -> np.ndarray:
    """The spectral axis' world values (frequency or velocity, whatever the WCS states).

    Read from the WCS' spectral axis — this is authoritative, not guessed: it is what the
    header declares. Falls back to a pixel index if no spectral WCS axis is present."""
    if header_string is not None:
        wcs = WCS(fits.Header.fromstring(header_string))
        spectral_axes = [
            i for i, t in enumerate(wcs.wcs.ctype) if _is_spectral_ctype(t)
        ]
        if spectral_axes:
            # Reverse-engineer the numpy<->WCS axis flip: WCS axis index -> world coord.
            wcs_axis = spectral_axes[0]
            pix = np.zeros((size, wcs.naxis))
            pix[:, wcs_axis] = np.arange(size)
            world = wcs.all_pix2world(pix, 0)
            return world[:, wcs_axis]
    return np.arange(size, dtype=float)


def _is_spectral_ctype(ctype: str) -> bool:
    code = str(ctype).strip().upper()
    return code.startswith(("FREQ", "VRAD", "VOPT", "VELO", "WAVE", "ENER", "FELO"))


def _to_fits_order(array: xr.DataArray) -> xr.DataArray:
    """Transpose a stored DataArray back to FITS numpy order (spectral, y, x)."""
    if set(array.dims) == set(CUBE_DIMS):
        return array.transpose(SPECTRAL_DIM, *SKY_DIMS)
    return array  # 2D map / other: stored in FITS order already


# -- zarr encoding (caller-chosen layout) ----------------------------------------------
def _to_store_axes(layout, ndim: int):
    """Translate a per-axis layout tuple from data order (spectral, y, x) to store order
    (sky_y, sky_x, freq). Only the 3D cube case reorders; 2D/other passes through unchanged."""
    if layout is None or ndim != 3:
        return layout
    spectral, sky_y, sky_x = layout
    return (sky_y, sky_x, spectral)


def _encoding(*, chunks, shards, compressors) -> dict:
    """Build the per-variable zarr encoding from caller layout params (none hard-coded).

    Only keys the caller actually set are passed through, so xarray/zarr apply their own
    defaults for the rest — astrolyze decides no chunking policy of its own."""
    var_encoding: dict = {}
    if chunks is not None:
        var_encoding["chunks"] = tuple(chunks)
    if shards is not None:
        var_encoding["shards"] = tuple(shards)
    if compressors is not None:
        var_encoding["compressors"] = compressors
    return {DATA_VAR: var_encoding} if var_encoding else {}


def _emit(op, **fields) -> None:
    """Append a run-log record through the always-on run-log seam (ADR-0010); a no-op outside
    an experiment. Deferred to call time so ``io`` stays light on import."""
    from astrolyze.experiment.runlog import emit

    emit(op, **fields)


__all__ = ["_save_zarr", "_load_zarr"]
