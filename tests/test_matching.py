"""Tests for guarded beam/channel matching on :class:`~astrolyze.core.Cube` (issue #31).

Written first (red/green TDD). These tests *are* the correctness obligation for the
matching geometry + its lossy-direction guard:

- the data maths is delegated wholesale to the ecosystem — spectral-cube
  (``convolve_to`` / ``downsample_axis`` / ``spectral_smooth``) does the convolution and
  binning, radio_beam (``Beam.deconvolve`` / solid angle) decides which direction is
  resolution-*losing*, and reproject puts two cubes on a common grid;
- astrolyze adds exactly one thing on top: **no silent physics** (ADR-0003). Every op only
  ever *degrades* resolution. The inverse — deconvolving to a smaller beam, up-sampling to
  finer channels — would invent information that is not in the data, so it **raises** a
  clear, named :class:`~astrolyze.core.LossyDirectionError` rather than super-resolving;
- every op is a first-class transition: it returns a :class:`~astrolyze.core.Cube` carrying
  updated context — the new, *larger* beam is recorded on the metadata (ADR-0004);
- ``match_to`` brings two cubes to a common beam and *optionally* reprojects them onto a
  common spatial grid. The reproject is **explicit** (an opt-in flag), never a silent side
  effect — line-ratio work must not silently regrid one map onto another.

Noise propagation through these ops is a separate concern (issue #32) and is deliberately
**not** exercised here.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention, in Jy/beam.
"""

import re
import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, LossyDirectionError

# spectral-cube emits cosmetic warnings (no-beam, stokes, big-cube) we do not care about.
warnings.filterwarnings("ignore", module="spectral_cube")
warnings.filterwarnings("ignore", module="astropy")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)
LARGER_BEAM = radio_beam.Beam(major=20 * u.arcsec, minor=18 * u.arcsec, pa=30 * u.deg)
SMALLER_BEAM = radio_beam.Beam(major=6 * u.arcsec, minor=5 * u.arcsec, pa=30 * u.deg)
CHANNEL_WIDTH = 2000.0 * u.m / u.s  # CDELT3 of the reference header


# --------------------------------------------------------------------------------------
# Synthetic FITS fixtures (a small PPV cube whose header carries the full schema)
# --------------------------------------------------------------------------------------
def _cube_header(crval1=24.174, crval2=15.784):
    """A 3D Jy/beam cube on the radio velocity convention with the full astrolyze schema.

    ``crval1`` / ``crval2`` are exposed so a second cube can be placed on a *shifted* grid,
    which is what makes the reproject step observable (a no-op reproject would prove nothing).
    """
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", crval1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", crval2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = (
        CHANNEL_WIDTH.to_value(u.m / u.s),
        1.0,
        "m/s",
    )
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "Jy/beam"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


def _make_cube(tmp_path, name="ngc0628_co21.fits", nspec=8, nx=5, ny=5, **header_kw):
    from astrolyze.io import load

    path = tmp_path / name
    data = np.ones((nspec, ny, nx), dtype="float32")
    fits.writeto(path, data, _cube_header(**header_kw))
    return Cube.from_loaded(load(path))


@pytest.fixture
def cube(tmp_path):
    """A complete :class:`Cube` built through the real io.load seam."""
    return _make_cube(tmp_path)


# --------------------------------------------------------------------------------------
# convolve_to_beam: the ALLOWED direction (smooth to a LARGER beam) delegates + carries
# --------------------------------------------------------------------------------------
def test_convolve_to_larger_beam_smooths_and_records_new_beam(cube):
    out = cube.convolve_to_beam(LARGER_BEAM)
    # first-class transition: a Cube back, with the *new, larger* beam recorded as context.
    assert isinstance(out, Cube)
    assert u.isclose(out.beam.major, LARGER_BEAM.major, rtol=1e-9)
    assert u.isclose(out.beam.minor, LARGER_BEAM.minor, rtol=1e-9)
    # the underlying SpectralCube was actually convolved to the new beam (delegation).
    assert u.isclose(out._sc.beam.major, LARGER_BEAM.major, rtol=1e-9)
    # the rest of the physical context flows through unchanged.
    assert u.isclose(out.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert out.velocity_convention is cube.velocity_convention
    # spatial-only smoothing: shape is unchanged.
    assert out.shape == cube.shape
    # and the input cube is not mutated (first-class return, no side effect).
    assert u.isclose(cube.beam.major, BEAM.major, rtol=1e-9)


def test_convolve_to_beam_delegates_to_spectral_cube(cube, monkeypatch):
    # astrolyze stays thin: the convolution maths is spectral-cube's ``convolve_to``. We do
    # not reimplement it — we guard the direction and then hand off.
    from spectral_cube import SpectralCube

    captured = {}
    real = SpectralCube.convolve_to

    def spy(self, beam, *args, **kwargs):
        captured["beam"] = beam
        return real(self, beam, *args, **kwargs)

    monkeypatch.setattr(SpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM)
    assert captured["beam"] is LARGER_BEAM


# --------------------------------------------------------------------------------------
# convolve_to_beam: the FFT (not direct) convolution path on dask-backed cubes (issue #51)
#
# spectral-cube's per-class default differs — SpectralCube.convolve_to defaults to convolve_fft,
# DaskSpectralCube.convolve_to to the ~65x slower direct convolve. A Zarr-backed (dask) Cube hits
# the dask class, so without an explicit convolve= it silently ran direct. These guard the fix.
# --------------------------------------------------------------------------------------
def _make_dask_cube(tmp_path, name="dask_ngc0628.fits"):
    """A dask-backed :class:`Cube` (the lazy-Zarr path) with non-trivial, structured data.

    Built through the real FITS->Zarr io seam so it exercises the same backend production hits.
    Random (not uniform) data makes the convolution non-degenerate, so FFT vs direct is a real
    numerical comparison rather than a constant passing through. float64 (not the fixtures' float32)
    so FFT and direct agree at machine precision — the point is that FFT is exact here, not an
    approximation, and float32 round-off (~1e-8) would mask that."""
    from astrolyze.io import load

    rng = np.random.default_rng(51)
    data = rng.normal(size=(8, 16, 16)).astype("float64")
    path = tmp_path / name
    fits.writeto(path, data, _cube_header())
    eager = Cube.from_loaded(load(path))
    return Cube.from_zarr(eager.to_zarr(tmp_path / "z"))


def test_dask_convolve_to_beam_matches_direct(tmp_path):
    # The FFT path is bit-identical to direct here (a pure performance bug, not an approximation).
    # Convolve the dask cube the production way (FFT) and compare to an explicit direct reference.
    from astropy.convolution import convolve as direct_convolve

    cube = _make_dask_cube(tmp_path)
    out = cube.convolve_to_beam(LARGER_BEAM)
    reference = cube._masked_sc().convolve_to(LARGER_BEAM, convolve=direct_convolve)
    np.testing.assert_allclose(
        np.asarray(out._sc.unmasked_data[:]),
        np.asarray(reference.unmasked_data[:]),
        atol=1e-10,
        rtol=0,
    )


def test_dask_convolve_to_beam_uses_fft_not_direct(tmp_path, monkeypatch):
    # The no-silent-regression guard: assert the production code passes convolve_fft on the dask
    # path, so the backend's ~65x slower direct default can never silently creep back in.
    from spectral_cube import DaskSpectralCube
    from astropy.convolution import convolve_fft

    cube = _make_dask_cube(tmp_path)
    captured = {}
    real = DaskSpectralCube.convolve_to

    def spy(self, beam, *args, convolve=None, **kwargs):
        captured["convolve"] = convolve
        return real(self, beam, *args, convolve=convolve, **kwargs)

    monkeypatch.setattr(DaskSpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM)
    assert captured["convolve"] is convolve_fft


def test_convolve_to_beam_forwards_save_to_tmp_dir_on_dask(tmp_path, monkeypatch):
    # save_to_tmp_dir is a DaskSpectralCube-only eager-materialisation control; it must reach the
    # dask convolve_to when requested (the caller's lever against OOM on big cubes).
    from spectral_cube import DaskSpectralCube

    cube = _make_dask_cube(tmp_path)
    captured = {}
    real = DaskSpectralCube.convolve_to

    def spy(self, beam, *args, save_to_tmp_dir=False, **kwargs):
        captured["save_to_tmp_dir"] = save_to_tmp_dir
        return real(self, beam, *args, save_to_tmp_dir=save_to_tmp_dir, **kwargs)

    monkeypatch.setattr(DaskSpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM, save_to_tmp_dir=True)
    assert captured["save_to_tmp_dir"] is True


def test_convolve_to_beam_ignores_save_to_tmp_dir_on_eager(cube):
    # The eager base-class convolve_to has no save_to_tmp_dir (it would error if forwarded into
    # convolve_fft). The kwarg must be silently dropped on the in-memory path, not crash.
    out = cube.convolve_to_beam(LARGER_BEAM, save_to_tmp_dir=True)
    assert isinstance(out, Cube)
    assert u.isclose(out.beam.major, LARGER_BEAM.major, rtol=1e-9)


# --------------------------------------------------------------------------------------
# convolve_to_beam: allow_huge threads through to convolve_fft's size guard (issue #88)
#
# convolve_fft refuses to allocate kernels above a hard "Arrays will be N.NG" size. allow_huge=True
# is the explicit caller lever that overrides that guard; it must thread down through convolve_to
# into convolve_fft. Default False must never inject it — today's behaviour, backward compatible.
# --------------------------------------------------------------------------------------
def test_convolve_to_beam_threads_allow_huge(cube, monkeypatch):
    # allow_huge is the lever that overrides convolve_fft's hard size guard (issue #88): when the
    # caller asks for it, astrolyze must forward it through convolve_to into convolve_fft.
    from spectral_cube import SpectralCube

    captured = {}
    real = SpectralCube.convolve_to

    def spy(self, beam, *args, **kwargs):
        captured["kwargs"] = kwargs
        return real(self, beam, *args, **kwargs)

    monkeypatch.setattr(SpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM, allow_huge=True)
    assert captured["kwargs"].get("allow_huge") is True


def test_convolve_to_beam_threads_allow_huge_on_dask(tmp_path, monkeypatch):
    # The big-cube case lives on the dask path; allow_huge must reach the dask convolve_to too
    # (overriding convolve_fft's size guard is exactly the lever a huge Zarr cube needs, issue #88).
    from spectral_cube import DaskSpectralCube

    cube = _make_dask_cube(tmp_path)
    captured = {}
    real = DaskSpectralCube.convolve_to

    def spy(self, beam, *args, **kwargs):
        captured["kwargs"] = kwargs
        return real(self, beam, *args, **kwargs)

    monkeypatch.setattr(DaskSpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM, allow_huge=True)
    assert captured["kwargs"].get("allow_huge") is True


def test_convolve_to_beam_defaults_no_allow_huge(cube, monkeypatch):
    # Regression guard: without allow_huge, astrolyze must NOT inject it into convolve_to's kwargs.
    # Default False == today's behaviour, so existing callers keep convolve_fft's size guard intact.
    from spectral_cube import SpectralCube

    captured = {}
    real = SpectralCube.convolve_to

    def spy(self, beam, *args, **kwargs):
        captured["kwargs"] = kwargs
        return real(self, beam, *args, **kwargs)

    monkeypatch.setattr(SpectralCube, "convolve_to", spy)
    cube.convolve_to_beam(LARGER_BEAM)
    assert "allow_huge" not in captured["kwargs"]


# --------------------------------------------------------------------------------------
# convolve_to_beam: the REFUSED direction (a smaller / deconvolving target) raises
# --------------------------------------------------------------------------------------
def test_convolve_to_smaller_beam_raises_clear_named_error(cube):
    # Smoothing to a *smaller* beam would be super-resolution: it would invent structure the
    # data never contained. astrolyze refuses rather than letting spectral-cube produce it.
    with pytest.raises(LossyDirectionError) as exc:
        cube.convolve_to_beam(SMALLER_BEAM)
    msg = str(exc.value)
    assert "beam" in msg.lower()
    # the error is informative about *why* (never super-resolve / cannot deconvolve / smaller).
    assert re.search(r"smaller|deconvolv|super-resolv|larger", msg.lower())


def test_convolve_to_same_beam_raises(cube):
    # The identity is not a degradation either; there is nothing to gain and a same-size beam
    # is on the refused side of the guard (it cannot be deconvolved from itself meaningfully).
    with pytest.raises(LossyDirectionError):
        cube.convolve_to_beam(BEAM)


# --------------------------------------------------------------------------------------
# spectral_bin: the ALLOWED direction (coarser channels) degrades + carries context
# --------------------------------------------------------------------------------------
def test_spectral_bin_degrades_and_returns_cube(cube):
    out = cube.spectral_bin(2)
    assert isinstance(out, Cube)
    # binning by 2 halves the number of channels (coarser spectral resolution).
    assert out.shape[0] == cube.shape[0] // 2
    assert out.shape[1:] == cube.shape[1:]
    # the new channel width is coarser (degraded) by the bin factor.
    new_width = abs(out.spectral_axis[1] - out.spectral_axis[0])
    assert u.isclose(new_width, 2 * CHANNEL_WIDTH, rtol=1e-6)
    # the beam (spatial context) is untouched and carried through.
    assert u.isclose(out.beam.major, cube.beam.major, rtol=1e-12)
    assert u.isclose(out.rest_frequency, cube.rest_frequency, rtol=1e-12)


def test_spectral_bin_factor_one_or_less_raises(cube):
    # factor 1 is a no-op and < 1 would *finen* the grid (up-sampling): both are refused —
    # spectral_bin only ever coarsens.
    with pytest.raises(LossyDirectionError):
        cube.spectral_bin(1)
    with pytest.raises(LossyDirectionError):
        cube.spectral_bin(0)


# --------------------------------------------------------------------------------------
# spectral_smooth_to: the ALLOWED direction (broader target width) degrades + carries
# --------------------------------------------------------------------------------------
def test_spectral_smooth_to_broader_width_degrades_and_returns_cube(cube):
    target = 6 * u.km / u.s  # broader than the 2 km/s native channel
    out = cube.spectral_smooth_to(target)
    assert isinstance(out, Cube)
    # smoothing in the spectral direction does not change the channel grid/shape.
    assert out.shape == cube.shape
    # context carried; beam untouched.
    assert u.isclose(out.beam.major, cube.beam.major, rtol=1e-12)
    assert u.isclose(out.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert out.velocity_convention is cube.velocity_convention


def test_spectral_smooth_to_finer_width_raises(cube):
    # A target resolution *finer* than the native channel width is up-sampling: refused.
    with pytest.raises(LossyDirectionError) as exc:
        cube.spectral_smooth_to(0.5 * u.km / u.s)
    assert re.search(r"finer|width|up-?sampl|coarser", str(exc.value).lower())


# --------------------------------------------------------------------------------------
# match_to: common beam, and an EXPLICIT (never silent) reproject to a common grid
# --------------------------------------------------------------------------------------
def test_match_to_brings_both_cubes_to_a_common_beam(tmp_path):
    fine = _make_cube(tmp_path, name="fine.fits")  # the 12x10 BEAM
    # a second cube with a coarser beam — the common beam must be at least the larger of the two.
    coarse = _make_cube(tmp_path, name="coarse.fits")
    coarse = coarse.convolve_to_beam(LARGER_BEAM)

    matched_self, matched_other = fine.match_to(coarse)
    assert isinstance(matched_self, Cube)
    assert isinstance(matched_other, Cube)
    # both now share a common beam (the larger of the two — never super-resolving either).
    assert u.isclose(matched_self.beam.major, matched_other.beam.major, rtol=1e-6)
    assert u.isclose(matched_self.beam.minor, matched_other.beam.minor, rtol=1e-6)
    assert matched_self.beam.major >= LARGER_BEAM.major - 1e-6 * u.arcsec


def test_match_to_does_not_reproject_by_default(tmp_path):
    # Line-ratio work must not silently regrid: without an explicit opt-in, match_to leaves
    # each cube on its own spatial grid (only the beam is brought to common).
    a = _make_cube(tmp_path, name="a.fits", nx=5, ny=5)
    # a partially overlapping but distinct grid (shifted ~2 px) and a different shape.
    b = _make_cube(tmp_path, name="b.fits", nx=7, ny=7, crval1=24.1744, crval2=15.7844)

    matched_a, matched_b = a.match_to(b)
    # grids untouched: the reproject did NOT happen as a side effect.
    assert matched_a.shape[1:] == a.shape[1:]
    assert matched_b.shape[1:] == b.shape[1:]


def test_match_to_reprojects_to_common_grid_when_asked(tmp_path):
    # With the explicit flag, the first cube is reprojected onto the second's spatial grid so
    # the two share one pixel grid (and one beam) — ready for a line ratio.
    a = _make_cube(tmp_path, name="a.fits", nx=5, ny=5)
    # an overlapping but shifted, differently-shaped grid so the reproject is observable.
    b = _make_cube(tmp_path, name="b.fits", nx=7, ny=7, crval1=24.1744, crval2=15.7844)

    matched_a, matched_b = a.match_to(b, reproject=True)
    assert isinstance(matched_a, Cube)
    assert isinstance(matched_b, Cube)
    # a now lives on b's spatial grid (common grid) while keeping its own spectral axis.
    assert matched_a.shape[1:] == matched_b.shape[1:] == b.shape[1:]
    # and still a common beam.
    assert u.isclose(matched_a.beam.major, matched_b.beam.major, rtol=1e-6)


# --------------------------------------------------------------------------------------
# Every op records the new beam as context (ADR-0004): the metadata is updated, not stale.
# --------------------------------------------------------------------------------------
def test_ops_record_updated_beam_on_metadata(cube):
    out = cube.convolve_to_beam(LARGER_BEAM)
    # the context the wrappers read from is the io Metadata — it must carry the new beam.
    assert u.isclose(out.metadata.beam.major, LARGER_BEAM.major, rtol=1e-9)
    # spectral ops leave the (spatial) beam recorded unchanged.
    binned = cube.spectral_bin(2)
    assert u.isclose(binned.metadata.beam.major, BEAM.major, rtol=1e-12)


# --------------------------------------------------------------------------------------
# The named error is exported and is a sensible Python error subclass.
# --------------------------------------------------------------------------------------
def test_lossy_direction_error_is_exported_value_error_subclass():
    assert issubclass(LossyDirectionError, ValueError)
