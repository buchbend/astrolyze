"""Tests for ``Cube.tiles`` — the spatial-cutout / mosaic-tiling primitive (issue #30).

Written first (red/green TDD). ``Cube.tiles(size, stride, overlap)`` is a cutout *iterator*
yielding spatial sub-``Cube``s that keep the **full** velocity/frequency extent — pure spatial
slicing (ADR-0004): astrolyze adds no resampling/interpolation, it delegates the slice to
spectral-cube and carries the context (beam + rest frequency + velocity convention) and a
correct sub-image WCS onto each cutout. The correctness obligations pinned here:

1. **Geometry** — on a known cube shape the tile *count and positions* are exactly the strided
   grid the caller asked for (``size`` / ``stride`` / ``overlap`` are caller parameters).
2. **Full spectral extent** — every cutout keeps all channels; only the sky plane is cut.
3. **Pure slicing** — each cutout's voxels are byte-for-byte the parent's voxels at that window
   (no resampling/interpolation/padding invented).
4. **Context carry + sub-WCS** — each cutout is a ``Cube`` carrying the same context, and its
   WCS is the parent WCS shifted to the window origin (so world coordinates are still correct).
5. **Explicit edge rule** — a partial edge window (smaller than ``size``) is **dropped by
   default** and never silently padded; ``partial="keep"`` yields the truncated cutout instead.
6. **Lazy on dask** — on a dask-backed (Zarr) cube, iterating tiles and taking a cutout stays a
   dask graph: no full materialisation (assert the cutout is still a dask collection + chunked).

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.io import load, save
from astrolyze.units import VelocityConvention

# spectral-cube emits cosmetic warnings (no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)

# A known cube shape: 4 channels x 8 (sky_y) x 10 (sky_x).
N_CHAN, NY, NX = 4, 8, 10


def _cube_header() -> fits.Header:
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
    return h


def _ramp() -> np.ndarray:
    """A distinctive per-voxel ramp so a wrong window / transpose is detectable."""
    return np.arange(N_CHAN * NY * NX, dtype="float32").reshape(N_CHAN, NY, NX)


@pytest.fixture
def fits_cube(tmp_path):
    """An eager :class:`Cube` built through the FITS io seam (the eager reference path)."""
    path = tmp_path / "ngc0628_co21.fits"
    fits.writeto(path, _ramp(), _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture
def zarr_cube(tmp_path):
    """A dask-backed (lazy) :class:`Cube` via the Zarr backend (#23) — the laziness fixture."""
    path = tmp_path / "ngc0628_co21.fits"
    fits.writeto(path, _ramp(), _cube_header())
    loaded = load(path)
    store = save(
        loaded.data,
        loaded.metadata,
        tmp_path / "z",
        format="zarr",
        chunks=(N_CHAN, NY, NX),  # one chunk over the whole array
        base_header=loaded.header,
    )
    return Cube.from_zarr(store)


# --------------------------------------------------------------------------------------
# AC: tiles() yields Cubes; the COUNT and POSITIONS on a known shape are exactly the
# strided grid the caller asked for (size / stride / overlap are caller parameters).
# --------------------------------------------------------------------------------------
def test_non_overlapping_tile_count_and_positions(fits_cube):
    # size 4, default stride == size (overlap 0): a non-overlapping grid that tiles 8x10
    # exactly in y (2 rows of 4) and leaves a partial 2-wide column in x (10 = 2*4 + 2),
    # which is DROPPED by default — so only the full 4x4 windows survive.
    tiles = list(fits_cube.tiles(size=4))
    assert all(isinstance(t, Cube) for t in tiles)
    # y origins: 0, 4 (both full); x origins: 0, 4 (8 is partial -> dropped).
    assert all(t.shape == (N_CHAN, 4, 4) for t in tiles)
    assert len(tiles) == 4  # 2 (y) x 2 (x)


def test_overlap_changes_the_stride_and_the_count(fits_cube):
    # overlap=2 with size=4 -> stride 2. On y=8: origins 0,2,4,6; window 6:10->6:8 is partial
    # (height 2) -> dropped, so y in {0,2,4}. On x=10: origins 0,2,4,6,8; 6:10 is FULL (width 4),
    # 8:12->8:10 is partial -> dropped, so x in {0,2,4,6}. The full-window grid is 3x4 = 12.
    tiles = list(fits_cube.tiles(size=4, overlap=2))
    # Positions are decided by stride; recover them from the data ramp's first voxel.
    full = fits_cube._data_quantity.value
    seen = []
    for t in tiles:
        assert t.shape == (N_CHAN, 4, 4)
        block = t._data_quantity.value
        # find the (y0, x0) window whose top-left voxel matches this tile's top-left.
        match = np.argwhere(full[0] == block[0, 0, 0])
        assert match.size  # the value is unique in channel 0 (a ramp)
        seen.append(tuple(int(v) for v in match[0]))
    assert sorted(seen) == [(y, x) for y in (0, 2, 4) for x in (0, 2, 4, 6)]
    assert len(tiles) == 12


def test_explicit_stride_overrides_overlap(fits_cube):
    # When stride is given it wins; stride 5 over size 4 in a 10-wide axis -> x origins 0, 5.
    tiles = list(fits_cube.tiles(size=(8, 4), stride=(8, 5)))
    assert len(tiles) == 2
    assert all(t.shape == (N_CHAN, 8, 4) for t in tiles)


# --------------------------------------------------------------------------------------
# AC: every cutout keeps the FULL spectral axis (only the sky plane is cut).
# --------------------------------------------------------------------------------------
def test_every_tile_keeps_the_full_spectral_axis(fits_cube):
    for tile in fits_cube.tiles(size=4):
        assert tile.shape[0] == N_CHAN
        assert tile.spectral_axis.shape == fits_cube.spectral_axis.shape
        assert u.allclose(tile.spectral_axis, fits_cube.spectral_axis, rtol=1e-12)


# --------------------------------------------------------------------------------------
# AC: tiling is pure slicing — voxels are byte-for-byte the parent's at that window
# (no resampling / interpolation / padding invented).
# --------------------------------------------------------------------------------------
def test_tiles_are_pure_slices_of_the_parent(fits_cube):
    full = fits_cube._data_quantity.value
    first = next(iter(fits_cube.tiles(size=4)))
    # The first tile is the y[0:4], x[0:4] window of every channel, unmodified.
    np.testing.assert_array_equal(first._data_quantity.value, full[:, 0:4, 0:4])


# --------------------------------------------------------------------------------------
# AC: each tile carries context (beam, rest_frequency, velocity_convention) and a correct
# sub-image WCS (the parent WCS shifted to the window origin).
# --------------------------------------------------------------------------------------
def test_each_tile_carries_context(fits_cube):
    for tile in fits_cube.tiles(size=4):
        assert u.isclose(tile.beam.major, BEAM.major, rtol=1e-12)
        assert u.isclose(tile.beam.minor, BEAM.minor, rtol=1e-12)
        assert u.isclose(tile.beam.pa, BEAM.pa, rtol=1e-12)
        assert u.isclose(tile.rest_frequency, REST, rtol=1e-12)
        assert tile.velocity_convention is VelocityConvention.RADIO


def test_tile_has_a_correct_sub_image_wcs(fits_cube):
    # The tile at the (y=4, x=4) window: its pixel (0,0) must map to the SAME sky position
    # as the parent's pixel (x=4, y=4) — i.e. the WCS reference moved with the window, so
    # world coordinates of a cutout voxel are still correct (a sub-image WCS, not the parent's).
    tiles = list(fits_cube.tiles(size=4))
    # locate the (4,4) tile by its ramp value.
    full = fits_cube._data_quantity.value
    target = full[0, 4, 4]
    tile = next(t for t in tiles if t._data_quantity.value[0, 0, 0] == target)

    parent_world = fits_cube._sc.wcs.celestial.wcs_pix2world([[4, 4]], 0)
    tile_world = tile._sc.wcs.celestial.wcs_pix2world([[0, 0]], 0)
    np.testing.assert_allclose(tile_world, parent_world, rtol=1e-12)


# --------------------------------------------------------------------------------------
# AC: edge handling is EXPLICIT — a partial edge window is DROPPED by default (no silent
# padding); partial="keep" yields the truncated cutout instead.
# --------------------------------------------------------------------------------------
def test_partial_edge_tiles_are_dropped_by_default(fits_cube):
    # x is 10 wide; size 4 stride 4 -> windows at 0,4 full and 8:12 partial (only 2 wide).
    tiles = list(fits_cube.tiles(size=4))
    assert all(t.shape[2] == 4 for t in tiles)  # no 2-wide partial survived
    # the dropped column is never padded to 4: there is simply no tile covering x>=8.
    assert len(tiles) == 4


def test_partial_edge_tiles_kept_when_requested(fits_cube):
    # partial="keep": the truncated 2-wide windows are yielded as smaller cutouts (still pure
    # slices — never padded). Now y has 0,4 and x has 0,4,8(partial,width 2).
    tiles = list(fits_cube.tiles(size=4, partial="keep"))
    widths = sorted({t.shape[2] for t in tiles})
    assert widths == [2, 4]  # both full and the truncated edge appear, no padding to 4
    assert len(tiles) == 6  # 2 (y) x 3 (x, incl. the partial edge)
    # the kept edge tile is exactly the parent's 2-wide window, unpadded.
    full = fits_cube._data_quantity.value
    edge = next(t for t in tiles if t.shape[2] == 2)
    assert edge.shape[1] == 4
    np.testing.assert_array_equal(
        edge._data_quantity.value[:, :, :],
        full[:, edge_y(edge, full) : edge_y(edge, full) + 4, 8:10],
    )


def edge_y(tile, full):
    """Recover the y-origin of a width-2 edge tile from the ramp (test helper)."""
    match = np.argwhere(full[0] == tile._data_quantity.value[0, 0, 0])
    return int(match[0][0])


def test_partial_mode_rejects_unknown_value(fits_cube):
    # No silent guess on an unknown edge rule (ADR-0003): an explicit error.
    with pytest.raises(ValueError, match="partial"):
        list(fits_cube.tiles(size=4, partial="pad"))


# --------------------------------------------------------------------------------------
# AC: lazy on a dask-backed Zarr cube — iterating tiles and taking a cutout never
# materialises the full array (assert chunk / laziness).
# --------------------------------------------------------------------------------------
def test_tiles_on_dask_cube_stay_lazy(zarr_cube):
    import dask

    # The source cube must actually be dask-backed (else the laziness claim is vacuous).
    assert dask.is_dask_collection(zarr_cube._sc._data)

    tile = next(iter(zarr_cube.tiles(size=4)))
    assert isinstance(tile, Cube)
    # The cutout is STILL a dask graph — no full materialisation to slice a window.
    assert dask.is_dask_collection(tile._sc._data)
    # …and it is chunked (carries a chunk grid), confirming laziness end-to-end.
    assert tile._sc._data.chunks is not None
    # Only an explicit compute makes it concrete — and it equals the eager slice.
    np.testing.assert_array_equal(
        np.asarray(tile._sc._data.compute())[:, :4, :4],
        _ramp()[:, 0:4, 0:4],
    )


def test_dask_tile_keeps_full_spectral_extent_and_context(zarr_cube):
    tile = next(iter(zarr_cube.tiles(size=4)))
    assert tile.shape[0] == N_CHAN
    assert u.isclose(tile.rest_frequency, REST, rtol=1e-12)
    assert tile.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(tile.beam.major, BEAM.major, rtol=1e-12)


# --------------------------------------------------------------------------------------
# Guard rails: bad geometry parameters are explicit errors, not silent guesses.
# --------------------------------------------------------------------------------------
def test_nonpositive_size_and_stride_raise(fits_cube):
    with pytest.raises(ValueError):
        list(fits_cube.tiles(size=0))
    with pytest.raises(ValueError):
        list(fits_cube.tiles(size=4, stride=0))
    with pytest.raises(ValueError):
        list(
            fits_cube.tiles(size=4, overlap=4)
        )  # overlap >= size -> non-positive stride
