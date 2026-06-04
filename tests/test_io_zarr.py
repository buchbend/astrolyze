"""Tests for the xarray-native Zarr backend behind the ``io`` seam (issue #23, ADR-0006).

This file is the correctness obligation for adding a **second on-disk backend alongside FITS**
(never replacing it). astrolyze stays thin: storage, coordinates, chunking and compression are
xarray/zarr/dask's job; astrolyze supplies only the ``Metadata`` <-> attrs projection (from
issue #22) and the FITS-vs-Zarr dispatch. The contracts pinned here:

1. **Round-trip through attrs.** ``save(format="zarr")`` writes a Zarr v3 store; ``load(store)``
   returns a :class:`~astrolyze.io.LoadedData` whose ``metadata`` equals the original — the
   physical context travels verbatim in the store's ``attrs`` via :meth:`Metadata.to_attrs` /
   :meth:`Metadata.from_attrs`.
2. **WCS survives FITS -> Zarr -> FITS.** The verbatim FITS-WCS header string
   (:attr:`LoadedData.header_string`, from #22) is carried in the store attrs, so the
   reconstructed WCS is *equal* to the original and a re-saved FITS preserves the full schema.
3. **Dispatch on the target.** ``load`` routes a ``.fits`` / ``.fits.gz`` path to the FITS
   backend and a Zarr store to the Zarr backend; both yield a ``LoadedData``.
4. **Lazy / dask-backed.** A Zarr-sourced cube is dask-backed — constructing it and slicing a
   subcube never materialises the full array. FITS stays eager.
5. **Caller-chosen layout.** Chunk / shard / compression are *parameters*; astrolyze hard-codes
   no chunking policy. A non-default chunking set on save is read back from the store.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import json

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from astropy.wcs import WCS

from astrolyze.io import LoadedData, Metadata, load, save

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


# --------------------------------------------------------------------------------------
# Synthetic FITS fixture (a small 3D PPV cube carrying the full schema + a real WCS)
# --------------------------------------------------------------------------------------
def _full_header() -> fits.Header:
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    h["HIERARCH ASTROLYZE DISTANCE"] = 9.84
    h["HIERARCH ASTROLYZE DISTUNIT"] = "Mpc"
    h["HIERARCH ASTROLYZE CALERR"] = 0.1
    h["HIERARCH ASTROLYZE NAMETAG"] = "mom0"
    return h


@pytest.fixture
def full_fits(tmp_path):
    """A tiny 3D cube on disk in FITS (spectral, y, x) order — the eager backend."""
    path = tmp_path / "raw" / "ngc0628_co21.fits"
    path.parent.mkdir(parents=True)
    # A distinctive ramp so a transposed/garbled axis order is detectable on round-trip.
    data = np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5)
    fits.writeto(path, data, _full_header())
    return path


def _assert_metadata_equal(m1: Metadata, m0: Metadata) -> None:
    """Every schema field round-tripped (Quantities compared by value+unit)."""
    assert m1.object == m0.object
    assert m1.telescope == m0.telescope
    assert m1.species == m0.species
    assert u.isclose(m1.rest_frequency, m0.rest_frequency, rtol=1e-12)
    assert m1.velocity_convention is m0.velocity_convention
    assert u.isclose(m1.beam.major, m0.beam.major, rtol=1e-12)
    assert u.isclose(m1.beam.minor, m0.beam.minor, rtol=1e-12)
    assert u.isclose(m1.beam.pa, m0.beam.pa, rtol=1e-12)
    assert m1.bunit == m0.bunit
    assert u.isclose(m1.distance, m0.distance, rtol=1e-12)
    assert m1.calibration_error == pytest.approx(m0.calibration_error)
    assert m1.name_tag == m0.name_tag


# --------------------------------------------------------------------------------------
# AC: save(format="zarr") writes a Zarr v3 store; load(store) round-trips the metadata
# --------------------------------------------------------------------------------------
def test_save_zarr_writes_a_v3_store(full_fits, tmp_path):
    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    assert store.exists()
    # zarr_format 3: the v3 group sentinel is zarr.json (not the v2 .zgroup).
    assert (store / "zarr.json").exists()


def test_zarr_roundtrip_metadata_equals_original(full_fits, tmp_path):
    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    back = load(store)
    assert isinstance(back, LoadedData)
    _assert_metadata_equal(back.metadata, loaded.metadata)
    assert back.metadata.is_complete is True
    # The whole reason for the attrs projection: it equals the in-memory truth (#22).
    assert back.metadata == Metadata.from_attrs(loaded.metadata.to_attrs())


def test_zarr_roundtrip_preserves_the_data(full_fits, tmp_path):
    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    back = load(store)
    # Axis order is preserved end-to-end (the distinctive ramp would expose a transpose bug).
    np.testing.assert_array_equal(np.asarray(back.data), loaded.data)


# --------------------------------------------------------------------------------------
# AC: FITS -> Zarr -> FITS round-trips the full schema AND the WCS (verbatim header string)
# --------------------------------------------------------------------------------------
def test_fits_zarr_fits_roundtrips_schema_and_wcs(full_fits, tmp_path):
    loaded = load(full_fits)

    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    from_zarr = load(store)

    # The verbatim FITS-WCS string is carried through the Zarr store (the #22 vehicle).
    assert isinstance(from_zarr.header_string, str)
    reconstructed = WCS(fits.Header.fromstring(from_zarr.header_string))
    assert reconstructed.wcs.compare(loaded.wcs.wcs)

    # ...and re-saving as FITS preserves the schema and the WCS.
    refits = save(
        from_zarr.data,
        from_zarr.metadata,
        tmp_path / "back",
        format="fits",
        base_header=fits.Header.fromstring(from_zarr.header_string),
    )
    reloaded = load(refits)
    _assert_metadata_equal(reloaded.metadata, loaded.metadata)
    assert reloaded.wcs.wcs.compare(loaded.wcs.wcs)


# --------------------------------------------------------------------------------------
# AC: load() dispatches on the target — .fits/.fits.gz -> FITS, Zarr store -> Zarr
# --------------------------------------------------------------------------------------
def test_load_dispatches_fits_path_to_fits_backend(full_fits):
    loaded = load(full_fits)
    # The FITS path stays eager (a live header is populated, data is a real ndarray).
    assert isinstance(loaded, LoadedData)
    assert isinstance(loaded.header, fits.Header)
    assert isinstance(loaded.data, np.ndarray)


def test_load_dispatches_zarr_store_to_zarr_backend(full_fits, tmp_path):
    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    back = load(store)
    assert isinstance(back, LoadedData)
    # A non-FITS backend has no live fits.Header — only the WCS string (#22).
    assert back.header is None
    assert isinstance(back.header_string, str)


@pytest.mark.parametrize("extension", ["fits", "fits.gz"])
def test_load_dispatches_both_fits_extensions(tmp_path, extension):
    path = tmp_path / f"cube.{extension}"
    data = np.arange(4 * 3 * 5, dtype="float32").reshape(4, 3, 5)
    fits.writeto(path, data, _full_header())
    loaded = load(path)
    assert isinstance(loaded, LoadedData)
    assert isinstance(loaded.header, fits.Header)


# --------------------------------------------------------------------------------------
# AC: a Zarr-backed cube is lazy / dask-backed — no full-array materialisation
# --------------------------------------------------------------------------------------
def test_zarr_load_is_dask_backed_and_lazy(full_fits, tmp_path):
    import dask

    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(2, 3, 5),
        base_header=loaded.header,
    )
    back = load(store)
    # Lazy: the loaded data is a dask collection, not a materialised ndarray.
    assert dask.is_dask_collection(back.data)
    # The chunking the caller chose is what the lazy array carries (storage is zarr's job).
    assert back.data.chunksize[0] == 2


def test_zarr_subcube_slice_does_not_materialise(full_fits, tmp_path):
    import dask

    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(2, 3, 5),
        base_header=loaded.header,
    )
    back = load(store)
    # Slicing a subcube stays lazy — it is still a dask graph, not a NumPy array.
    sub = back.data[:2]
    assert dask.is_dask_collection(sub)
    # Only when explicitly computed does it become concrete (and match the eager slice).
    np.testing.assert_array_equal(np.asarray(sub.compute()), loaded.data[:2])


# --------------------------------------------------------------------------------------
# AC: chunk / shard / compression are caller parameters (no hard-coded policy)
# --------------------------------------------------------------------------------------
def test_caller_chunking_is_read_back_from_the_store(full_fits, tmp_path):
    loaded = load(full_fits)
    # Caller chunks are in data axis order (spectral, y, x) — the same order as loaded.data.
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(1, 3, 5),
        base_header=loaded.header,
    )
    # The on-disk grid is in the store's axis order (sky_y, sky_x, freq); reading it back from
    # the array's v3 metadata (zarr owns the layout) shows the caller's choice, reordered.
    import zarr

    arr = zarr.open_array(store=str(store / "intensity"), mode="r")
    assert tuple(arr.chunks) == (3, 5, 1)
    # ...and round-tripping through load() returns the chunk grid in data order again.
    back = load(store)
    assert back.data.chunksize == (1, 3, 5)


def test_caller_shard_and_compression_are_honoured(full_fits, tmp_path):
    loaded = load(full_fits)
    from zarr.codecs import BloscCodec

    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(1, 3, 5),
        shards=(2, 3, 5),
        compressors=[BloscCodec(cname="zstd", clevel=3)],
        base_header=loaded.header,
    )
    back = load(store)
    # Shard/compression are storage details; the read-back value must still be exact.
    np.testing.assert_array_equal(np.asarray(back.data.compute()), loaded.data)


# --------------------------------------------------------------------------------------
# AC: the store carries the attrs projection + verbatim WCS string + provenance
# --------------------------------------------------------------------------------------
def test_store_attrs_carry_metadata_wcs_and_provenance(full_fits, tmp_path):
    loaded = load(full_fits)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    group_attrs = json.loads((store / "zarr.json").read_text())["attributes"]
    # The schema attrs travel verbatim (JSON primitives, #22) ...
    assert group_attrs["object"] == "NGC0628"
    assert group_attrs["velocity_convention"] == "radio"
    # ... alongside the verbatim FITS-WCS header string and a provenance marker.
    assert isinstance(group_attrs["fits_wcs_header"], str)
    assert group_attrs["provenance"]["backend"] == "zarr"


# --------------------------------------------------------------------------------------
# House rule: a partial schema still round-trips through Zarr and stays flagged incomplete
# --------------------------------------------------------------------------------------
def test_partial_metadata_roundtrips_through_zarr_and_stays_incomplete(tmp_path):
    # An archival-style header with no rest frequency / velocity convention.
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["OBJECT"] = "NGC0628"
    fits_path = tmp_path / "archival.fits"
    fits.writeto(fits_path, np.zeros((3, 5), dtype="float32"), h)

    loaded = load(fits_path)
    assert loaded.metadata.is_complete is False

    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        base_header=loaded.header,
    )
    back = load(store)
    # Lazy enforcement survives the Zarr projection: still incomplete, same gaps (ADR-0006 ii).
    assert back.metadata.is_complete is False
    assert set(back.metadata.missing) == set(loaded.metadata.missing)


# --------------------------------------------------------------------------------------
# House rule: save() rejects an unknown format rather than guessing (ADR-0003)
# --------------------------------------------------------------------------------------
def test_save_rejects_unknown_format(full_fits, tmp_path):
    loaded = load(full_fits)
    with pytest.raises(ValueError, match="format"):
        save(
            loaded.data,
            loaded.metadata,
            tmp_path / "out",
            format="hdf5",
            base_header=loaded.header,
        )
