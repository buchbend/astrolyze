"""Tests for ``Cube.cutout(SkyCoord, size)`` — the sky-coordinate postage stamp (issue #58).

Written first (red/green TDD). ``Cube.cutout`` is a *spatial* postage-stamp primitive: a
window centred on a :class:`~astropy.coordinates.SkyCoord` with a requested angular ``size``,
keeping the **full** spectral axis (spectra are never cut) and a correct sub-image WCS, while
carrying the cube's context (beam + rest frequency + velocity convention + Metadata). It stays
thin (ADR-0004): the center+size -> pixel-window maths is astropy's :class:`~astropy.nddata.Cutout2D`,
the slice itself is the existing ``self[:, y0:y1, x0:x1]`` path, so spectral-cube does the
sub-WCS and ``__getitem__`` carries the context. The correctness obligations pinned here:

1. **Sky-coordinate window** — on a known WCS the data window equals the hand-computed pixel
   bracket around the SkyCoord for the requested angular size, byte-for-byte the parent's voxels.
2. **Sub-WCS** — the cutout's pixel (0,0) maps to the SAME sky position as the parent's pixel at
   the window origin (the WCS moved with the window; world coordinates stay correct).
3. **Full spectral extent** — every cutout keeps all channels; only the sky plane is cut.
4. **Context carry** — the cutout is a ``Cube`` carrying the same beam / rest frequency /
   convention / Metadata.
5. **Off-footprint** — a center SkyCoord outside the cube's sky footprint raises a descriptive
   error (no empty/garbage stamp).
6. **Edge overlap is explicit** — a window straddling the edge raises by default (no silent
   clipping); ``partial="trim"`` opts in to the trimmed (smaller) stamp instead.
7. **Lazy on dask** — on a dask-backed (Zarr) cube the cutout stays a dask graph (no full
   materialisation).

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at 230.538 GHz,
a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings
from pathlib import Path

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs.utils import pixel_to_skycoord

from astrolyze.core import Cube
from astrolyze.io import load, save
from astrolyze.units import VelocityConvention

# spectral-cube emits cosmetic warnings (no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)

# A known cube shape: 4 channels x 8 (sky_y) x 10 (sky_x). The pixel scale is |CDELT| = 2e-4
# deg = 0.72 arcsec/pixel, so an angular size of N*0.72 arcsec brackets N pixels.
N_CHAN, NY, NX = 4, 8, 10
PIX_DEG = 2e-4
PIX_ARCSEC = PIX_DEG * 3600.0  # 0.72 arcsec/pixel
CUTOUT_FIXTURE = Path(__file__).parent / "data" / "ngc0628_co21_cutout.fits.gz"


def _cube_header() -> fits.Header:
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -PIX_DEG, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = PIX_DEG, 1.0, "deg"
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


def _sky_at_pixel(cube: Cube, x: int, y: int) -> SkyCoord:
    """The sky position of integer pixel ``(x, y)`` on the cube's celestial WCS."""
    return pixel_to_skycoord(x, y, cube._sc.wcs.celestial, origin=0)


# --------------------------------------------------------------------------------------
# AC 1: the data window equals the hand-computed pixel bracket around the SkyCoord for the
# requested angular size, byte-for-byte the parent's voxels (pure slice, no resampling).
# --------------------------------------------------------------------------------------
def test_cutout_window_matches_hand_computed_pixels(fits_cube):
    # Centre on pixel (x=5, y=4); request an odd 3x3 (px) stamp => exactly pixels
    # x in [4,6], y in [3,5] (centre +/- 1). Hand-computed slice: [3:6, 4:7].
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    assert isinstance(stamp, Cube)
    assert stamp.shape == (N_CHAN, 3, 3)
    parent = fits_cube._data_quantity.value
    np.testing.assert_array_equal(stamp._data_quantity.value, parent[:, 3:6, 4:7])


def test_cutout_accepts_a_rectangular_size_pair(fits_cube):
    # size=(height, width) as a Quantity pair: a 3 (y) x 5 (x) px stamp centred on (x=5,y=4)
    # => y in [3,5], x in [3,7] -> slice [3:6, 3:8].
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(
        centre, size=(3 * PIX_ARCSEC * u.arcsec, 5 * PIX_ARCSEC * u.arcsec)
    )
    assert stamp.shape == (N_CHAN, 3, 5)
    parent = fits_cube._data_quantity.value
    np.testing.assert_array_equal(stamp._data_quantity.value, parent[:, 3:6, 3:8])


def test_cutout_is_a_pure_slice_no_resampling(fits_cube):
    # The stamp's voxels are exactly the parent's at the window — never interpolated/padded.
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    parent = fits_cube._data_quantity.value
    np.testing.assert_array_equal(stamp._data_quantity.value, parent[:, 3:6, 4:7])


# --------------------------------------------------------------------------------------
# AC 2: the cutout's sub-WCS is the parent WCS shifted to the window origin (world
# coordinates of a cutout voxel are still correct).
# --------------------------------------------------------------------------------------
def test_cutout_has_a_correct_sub_image_wcs(fits_cube):
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    # The window origin is parent pixel (x=4, y=3); stamp pixel (0,0) must map to its sky.
    parent_world = fits_cube._sc.wcs.celestial.wcs_pix2world([[4, 3]], 0)
    stamp_world = stamp._sc.wcs.celestial.wcs_pix2world([[0, 0]], 0)
    np.testing.assert_allclose(stamp_world, parent_world, rtol=1e-12)
    # And the requested SkyCoord lands back at the stamp's centre pixel (x=1, y=1).
    centre_px = stamp._sc.wcs.celestial.world_to_pixel(centre)
    np.testing.assert_allclose(
        [float(centre_px[0]), float(centre_px[1])], [1.0, 1.0], atol=1e-9
    )


# --------------------------------------------------------------------------------------
# AC 3: full spectral axis preserved (only the sky plane is cut).
# --------------------------------------------------------------------------------------
def test_cutout_preserves_the_full_spectral_axis(fits_cube):
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    assert stamp.shape[0] == N_CHAN
    assert stamp.spectral_axis.shape == fits_cube.spectral_axis.shape
    assert u.allclose(stamp.spectral_axis, fits_cube.spectral_axis, rtol=1e-12)


# --------------------------------------------------------------------------------------
# AC 4: the cutout carries context (beam, rest frequency, convention, Metadata).
# --------------------------------------------------------------------------------------
def test_cutout_carries_context(fits_cube):
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    stamp = fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    assert u.isclose(stamp.beam.major, BEAM.major, rtol=1e-12)
    assert u.isclose(stamp.beam.minor, BEAM.minor, rtol=1e-12)
    assert u.isclose(stamp.beam.pa, BEAM.pa, rtol=1e-12)
    assert u.isclose(stamp.rest_frequency, REST, rtol=1e-12)
    assert stamp.velocity_convention is VelocityConvention.RADIO
    assert stamp.metadata.object == "NGC0628"


# --------------------------------------------------------------------------------------
# AC 5: a center coordinate outside the cube footprint raises a descriptive error.
# --------------------------------------------------------------------------------------
def test_cutout_off_footprint_centre_raises(fits_cube):
    # A SkyCoord ~50 pixels off the corner: well outside the 8x10 footprint.
    off = pixel_to_skycoord(-50, -50, fits_cube._sc.wcs.celestial, origin=0)
    with pytest.raises(ValueError, match="footprint"):
        fits_cube.cutout(off, size=3 * PIX_ARCSEC * u.arcsec)


# --------------------------------------------------------------------------------------
# AC 6: edge-overlap behaviour is EXPLICIT — raise by default (no silent clip), trim on opt-in.
# --------------------------------------------------------------------------------------
def test_cutout_partial_overlap_raises_by_default(fits_cube):
    # Centre on the (x=0, y=0) corner with a 3px stamp -> the window straddles the edge.
    corner = _sky_at_pixel(fits_cube, x=0, y=0)
    with pytest.raises(ValueError, match="partial|edge|overlap"):
        fits_cube.cutout(corner, size=3 * PIX_ARCSEC * u.arcsec)


def test_cutout_partial_trim_opt_in_yields_smaller_stamp(fits_cube):
    # partial="trim": the same edge stamp is returned trimmed to the in-bounds region (2x2),
    # never padded out to 3x3. Centre (0,0), 3px -> overlap is pixels x in [0,1], y in [0,1].
    corner = _sky_at_pixel(fits_cube, x=0, y=0)
    stamp = fits_cube.cutout(corner, size=3 * PIX_ARCSEC * u.arcsec, partial="trim")
    assert stamp.shape == (N_CHAN, 2, 2)
    parent = fits_cube._data_quantity.value
    np.testing.assert_array_equal(stamp._data_quantity.value, parent[:, 0:2, 0:2])


def test_cutout_rejects_unknown_partial_value(fits_cube):
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    with pytest.raises(ValueError, match="partial"):
        fits_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec, partial="pad")


def test_cutout_requires_angular_size(fits_cube):
    # No silent guess of pixels-vs-angle: a bare number / wrong-dimension size raises.
    centre = _sky_at_pixel(fits_cube, x=5, y=4)
    with pytest.raises((u.UnitsError, ValueError, TypeError)):
        fits_cube.cutout(centre, size=3)
    with pytest.raises((u.UnitsError, ValueError)):
        fits_cube.cutout(centre, size=3 * u.km / u.s)


# --------------------------------------------------------------------------------------
# AC 7: lazy on a dask-backed Zarr cube — the cutout never materialises the full array.
# --------------------------------------------------------------------------------------
def test_cutout_on_dask_cube_stays_lazy(zarr_cube):
    import dask

    assert dask.is_dask_collection(
        zarr_cube._sc._data
    )  # the laziness claim is not vacuous
    centre = _sky_at_pixel(zarr_cube, x=5, y=4)
    stamp = zarr_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    assert isinstance(stamp, Cube)
    assert dask.is_dask_collection(stamp._sc._data)  # still a dask graph
    assert stamp._sc._data.chunks is not None  # carries a chunk grid
    # Only an explicit compute makes it concrete — and it equals the eager slice.
    np.testing.assert_array_equal(
        np.asarray(stamp._sc._data.compute()),
        _ramp()[:, 3:6, 4:7],
    )


def test_cutout_dask_keeps_full_spectral_extent_and_context(zarr_cube):
    centre = _sky_at_pixel(zarr_cube, x=5, y=4)
    stamp = zarr_cube.cutout(centre, size=3 * PIX_ARCSEC * u.arcsec)
    assert stamp.shape[0] == N_CHAN
    assert u.isclose(stamp.rest_frequency, REST, rtol=1e-12)
    assert stamp.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(stamp.beam.major, BEAM.major, rtol=1e-12)


# --------------------------------------------------------------------------------------
# Realism smoke: the committed real-data PHANGS NGC 628 cutout (always-on, mirrors the
# tracer's real-data pattern). Cutting a stamp at the cube's own centre keeps the spectral
# axis and yields real, finite signal.
# --------------------------------------------------------------------------------------
@pytest.mark.skipif(
    not CUTOUT_FIXTURE.exists(), reason="committed cutout fixture missing"
)
def test_cutout_on_real_phangs_cube():
    cube = Cube.from_loaded(load(CUTOUT_FIXTURE))
    ny, nx = cube.shape[1], cube.shape[2]
    centre = pixel_to_skycoord(nx // 2, ny // 2, cube._sc.wcs.celestial, origin=0)
    # A 20 arcsec stamp comfortably inside the ~92 arcsec field (128 px * 0.72 arcsec).
    stamp = cube.cutout(centre, size=20 * u.arcsec)
    assert isinstance(stamp, Cube)
    # full spectral axis preserved, smaller sky plane.
    assert stamp.shape[0] == cube.shape[0]
    assert stamp.shape[1] < ny and stamp.shape[2] < nx
    # context carried through (the fixture parses as complete: beam + rest freq + convention).
    assert stamp.rest_frequency is not None
    assert stamp.beam is not None
    # real CO emission: the stamp's moment-0 has genuine, finite signal.
    mom0 = stamp.moment0().to("K km/s")
    values = u.Quantity(mom0.data).value
    assert np.isfinite(values).all()
    assert np.nanmax(values) > 1.0  # K km/s
