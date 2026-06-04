"""Tests for astrolyze.core — the Cube / Map / Spectrum context-carrying wrappers (ADR-0004).

Written first (red/green TDD). These tests *are* the correctness obligation for the object
model contract:

- the wrappers **compose** spectral-cube / astropy / specutils (a private upstream handle)
  rather than subclass them, and add the toolkit's value: object-carried context
  (beam + rest frequency + velocity convention, ADR-0003) sourced from ``io`` Metadata;
- type transitions carry context for free: ``Cube.moment0() -> Map`` keeps the beam and rest
  frequency; ``cube[:, y, x] -> Spectrum`` keeps the same context;
- ``.to(unit)`` is the unit hub (ADR-0003c) — it supplies the object's beam / rest frequency /
  convention so the caller need not, while still **never** defaulting the genuinely-ambiguous
  Rayleigh-Jeans-vs-Planck scale (that stays explicit);
- an operation that needs context the object does not have **raises** (lazy enforcement,
  ADR-0006 ii) rather than silently guessing.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention — here in Jy/beam
so the unit-hub conversions exercise real beam + rest-frequency context.
"""

import re
import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube, Map, Spectrum
from astrolyze.io import Metadata, load
from astrolyze.units import MissingContextError, VelocityConvention, convert

# spectral-cube emits cosmetic warnings (e.g. no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


# --------------------------------------------------------------------------------------
# Synthetic FITS fixtures (a tiny PPV cube whose header carries the schema)
# --------------------------------------------------------------------------------------
def _cube_header():
    """A 3D Jy/beam cube on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = "Jy/beam"
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


@pytest.fixture
def cube(tmp_path):
    """A complete :class:`Cube` built through the real io.load seam (ADR-0004/0006)."""
    path = tmp_path / "ngc0628_co21.fits"
    data = np.ones((4, 3, 3), dtype="float32")
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture
def integrated_map(cube):
    """A velocity-integrated map (Jy/beam km/s), the moment0 of the reference cube."""
    return cube.moment0()


# --------------------------------------------------------------------------------------
# Composition, not subclassing (ADR-0004 rider 1)
# --------------------------------------------------------------------------------------
def test_wrappers_compose_and_do_not_subclass_upstream(cube):
    from spectral_cube import SpectralCube

    # holds a private handle and delegates; is NOT a SpectralCube itself.
    assert isinstance(cube, Cube)
    assert not isinstance(cube, SpectralCube)
    assert isinstance(cube._sc, SpectralCube)
    # context is sourced from the io Metadata (the single source of truth).
    assert isinstance(cube.metadata, Metadata)
    assert u.isclose(cube.rest_frequency, REST, rtol=1e-12)
    assert cube.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(cube.beam.major, 12 * u.arcsec, rtol=1e-9)


# --------------------------------------------------------------------------------------
# MANDATED: cube.moment0() -> Map carrying beam + rest frequency
# --------------------------------------------------------------------------------------
def test_moment0_returns_map_carrying_beam_and_frequency(cube):
    m = cube.moment0()
    assert isinstance(m, Map)
    # the 2D moment carries the cube's physical context for free.
    assert u.isclose(m.beam.major, cube.beam.major, rtol=1e-12)
    assert u.isclose(m.beam.minor, cube.beam.minor, rtol=1e-12)
    assert u.isclose(m.beam.pa, cube.beam.pa, rtol=1e-12)
    assert u.isclose(m.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert m.velocity_convention is cube.velocity_convention
    # moment0 of a Jy/beam cube over a velocity axis is a velocity-integrated intensity.
    assert m.shape == (3, 3)
    assert m.unit.is_equivalent(u.Jy / u.beam * u.km / u.s)


# --------------------------------------------------------------------------------------
# MANDATED: cube[:, y, x] -> Spectrum (carrying context)
# --------------------------------------------------------------------------------------
def test_spatial_index_returns_spectrum_carrying_context(cube):
    sp = cube[:, 1, 1]
    assert isinstance(sp, Spectrum)
    # the spectral axis comes through from the cube WCS.
    assert sp.flux.shape == (4,)
    assert sp.flux.unit.is_equivalent(u.Jy / u.beam)
    assert sp.spectral_axis.unit.is_equivalent(u.m / u.s)
    # same physical context as the parent cube.
    assert u.isclose(sp.rest_frequency, cube.rest_frequency, rtol=1e-12)
    assert sp.velocity_convention is cube.velocity_convention
    assert u.isclose(sp.beam.major, cube.beam.major, rtol=1e-12)


def test_spectrum_wraps_specutils(cube):
    from specutils import Spectrum as SpecutilsSpectrum

    sp = cube[:, 0, 0]
    assert isinstance(
        sp._spec, SpecutilsSpectrum
    )  # thin specutils adaptor (ADR-0004 rider 2)


# --------------------------------------------------------------------------------------
# MANDATED: map.to("K km/s") uses the object's context (beam + rest frequency)
# --------------------------------------------------------------------------------------
def test_map_to_uses_object_context(integrated_map):
    # The object supplies beam + rest frequency + convention; the caller supplies only the
    # Rayleigh-Jeans-vs-Planck scale, which astrolyze NEVER defaults (ADR-0003). For a
    # velocity-integrated intensity RJ is the only valid scale (Planck is non-linear).
    out = integrated_map.to("K km/s", temperature_scale="rayleigh_jeans")
    assert isinstance(out, Map)
    assert out.unit.is_equivalent(u.K * u.km / u.s)
    # identical to calling the unit layer directly with the object's context spelled out:
    # this is exactly the context the wrapper supplied on the caller's behalf.
    expected = convert(
        integrated_map.data,
        "K km/s",
        rest_frequency=REST,
        beam=BEAM,
        temperature_scale="rayleigh_jeans",
    )
    assert u.allclose(out.data, expected, rtol=1e-10)
    # context is preserved across the conversion.
    assert u.isclose(out.rest_frequency, REST, rtol=1e-12)
    assert u.isclose(out.beam.major, BEAM.major, rtol=1e-12)


def test_to_returns_same_wrapper_type_with_updated_unit(integrated_map):
    out = integrated_map.to("K km/s", temperature_scale="rayleigh_jeans")
    assert type(out) is Map
    assert out.metadata.bunit == out.unit


# --------------------------------------------------------------------------------------
# MANDATED: an operation needing absent context raises (no silent guess, ADR-0003/0006)
# --------------------------------------------------------------------------------------
def test_to_on_incomplete_object_raises():
    # A real-archival map with a beam but no rest frequency (lazy: holding it is fine).
    incomplete = Map(
        np.ones((3, 3)) * (u.Jy / u.beam),
        wcs=None,
        metadata=Metadata(bunit=u.Jy / u.beam, beam=BEAM),
    )
    assert incomplete.is_complete is False
    assert "rest_frequency" in incomplete.missing
    # Jy/beam -> K needs the rest frequency the object does not carry: it must raise.
    with pytest.raises(MissingContextError, match="rest_frequency"):
        incomplete.to("K", temperature_scale="rayleigh_jeans")


def test_to_never_defaults_the_brightness_scale(cube):
    # A per-channel Jy/beam -> K conversion is genuinely ambiguous (RJ vs Planck); the
    # object cannot know which, so .to() must refuse rather than pick one silently.
    channel = cube[0]  # a 2D channel map, still Jy/beam
    assert isinstance(channel, Map)
    with pytest.raises(MissingContextError, match="temperature_scale"):
        channel.to("K")
    # with the scale stated explicitly it succeeds (and both scales are reachable).
    assert channel.to("K", temperature_scale="planck").unit.is_equivalent(u.K)
    assert channel.to("K", temperature_scale="rayleigh_jeans").unit.is_equivalent(u.K)


# --------------------------------------------------------------------------------------
# Type transitions on slicing (Cube -> Cube / Map / Spectrum), context carried through
# --------------------------------------------------------------------------------------
def test_subcube_slice_returns_cube(cube):
    sub = cube[:, 0:2, 0:2]
    assert isinstance(sub, Cube)
    assert sub.shape == (4, 2, 2)
    assert u.isclose(sub.rest_frequency, cube.rest_frequency, rtol=1e-12)


def test_channel_slice_returns_map(cube):
    channel = cube[0]
    assert isinstance(channel, Map)
    assert channel.shape == (3, 3)
    assert channel.unit.is_equivalent(u.Jy / u.beam)
    assert u.isclose(channel.beam.major, cube.beam.major, rtol=1e-12)


# --------------------------------------------------------------------------------------
# .plot() is a thin seam onto the viz engine (ADR-0005): it delegates to the matching
# free function, passing the object + kwargs. (The engine itself is tested in test_viz.)
# --------------------------------------------------------------------------------------
def test_plot_delegates_to_the_matching_viz_function(integrated_map, monkeypatch):
    import astrolyze.viz as viz

    captured = {}

    def fake_plot_map(obj, **kwargs):
        captured["obj"], captured["kwargs"] = obj, kwargs
        return ("fig", "ax")

    monkeypatch.setattr(viz, "plot_map", fake_plot_map)
    result = integrated_map.plot(cmap="magma")

    assert result == ("fig", "ax")
    assert captured["obj"] is integrated_map  # the object is handed to the engine
    assert captured["kwargs"] == {
        "cmap": "magma"
    }  # caller kwargs flow straight through


def test_plot_raises_clearly_when_its_viz_function_is_absent(
    integrated_map, monkeypatch
):
    # Defensive branch: a wrapper pointing at a viz function that does not exist gets a
    # clear NotImplementedError, not an obscure AttributeError.
    monkeypatch.setattr(type(integrated_map), "_viz_function", "plot_does_not_exist")
    with pytest.raises(NotImplementedError):
        integrated_map.plot()


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code (ADR-0005 / the legacy sin)
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_library_code():
    from pathlib import Path

    import astrolyze.core as core_pkg

    pkg_dir = Path(core_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"


# ======================================================================================
# Issue #26 — per-axis physical coordinates + a validity descriptor, READ from the WCS /
# spectral axis astrolyze already parses (no reparse, no new geometry).
#
# Correctness obligation: a scientist (and any downstream consumer) must be able to read,
# off the core object,
#   - the per-axis physical coordinate arrays as Quantity (each carrying its own unit):
#     absolute frequency (authoritative), per-line Δv, sky coordinates, pixel scale;
#   - *where the data is real* — a validity descriptor exposing blanked / edge / outside-
#     coverage voxels as NaN plus a boolean finite-data mask;
# without re-deriving any of it from the header. The frequency is AUTHORITATIVE: it is read
# from the spectral axis (already velocity here) by the same convention + rest frequency the
# object carries, so the toolkit never silently guesses the convention (ADR-0003). The
# values are checked against the WCS directly, and the mask is checked on a cube containing
# NaNs (and that it travels through a subcube slice: mask-of-a-slice == slice-of-the-mask).
# ======================================================================================
@pytest.fixture
def cube_with_nans(tmp_path):
    """The reference cube with two blanked voxels — the validity / NaN-mask fixture."""
    path = tmp_path / "ngc0628_co21_blanked.fits"
    data = np.ones((4, 3, 3), dtype="float32")
    data[0, 0, 0] = np.nan  # a single blanked voxel
    data[3, 2, 2] = (
        np.nan
    )  # a second, in a different corner (so a 2x2 slice can exclude it)
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


def _wcs_spectral_velocity(cube):
    """The per-channel velocity, straight off the cube WCS (the VRAD spectral axis)."""
    return cube._sc.spectral_axis.to(u.km / u.s)


# --------------------------------------------------------------------------------------
# AC: Cube exposes per-axis physical coordinate arrays as Quantity (each with its unit),
# read from the existing WCS / spectral axis. Values checked against the WCS.
# --------------------------------------------------------------------------------------
def test_cube_coordinates_are_quantities_read_from_the_wcs(cube):
    coords = cube.coordinates

    # absolute frequency is AUTHORITATIVE and carries a frequency unit.
    assert isinstance(coords.frequency, u.Quantity)
    assert coords.frequency.unit.is_equivalent(u.Hz)
    assert coords.frequency.shape == (4,)
    # the channel at v=0 (CRVAL3) is the rest frequency on the radio convention.
    assert u.isclose(coords.frequency[0], REST, rtol=1e-9)

    # per-line Δv carries a velocity unit and equals the WCS velocity axis (radio + rest).
    assert isinstance(coords.delta_v, u.Quantity)
    assert coords.delta_v.unit.is_equivalent(u.km / u.s)
    assert u.allclose(
        coords.delta_v, _wcs_spectral_velocity(cube), rtol=1e-9, atol=1e-6 * u.km / u.s
    )

    # sky coordinates carry an angular unit and match the WCS world map.
    assert isinstance(coords.longitude, u.Quantity)
    assert isinstance(coords.latitude, u.Quantity)
    assert coords.longitude.unit.is_equivalent(u.deg)
    assert coords.latitude.unit.is_equivalent(u.deg)
    lat_wcs, lon_wcs = (
        cube._sc.spatial_coordinate_map
    )  # spectral-cube returns [lat, lon]
    assert u.allclose(coords.longitude, lon_wcs, rtol=1e-12)
    assert u.allclose(coords.latitude, lat_wcs, rtol=1e-12)

    # pixel scale carries an angular unit and matches CDELT (0.0002 deg per pixel).
    assert isinstance(coords.pixel_scale, u.Quantity)
    assert coords.pixel_scale.unit.is_equivalent(u.deg)
    assert u.allclose(coords.pixel_scale, [2e-4, 2e-4] * u.deg, rtol=1e-9)


def test_cube_coordinates_do_not_reparse_the_header(cube, monkeypatch):
    # The accessor must read the *already-parsed* WCS / spectral axis, never re-open or
    # re-parse a FITS header. We make any header parse explode and assert the read still works.
    import astrolyze.io.schema as schema

    def _boom(*a, **k):  # pragma: no cover - only fires on a (forbidden) reparse
        raise AssertionError("coordinates must not reparse the FITS header")

    monkeypatch.setattr(
        schema.Metadata, "from_header", classmethod(lambda cls, h: _boom())
    )
    coords = cube.coordinates
    assert coords.frequency.shape == (4,)
    assert coords.longitude.shape == (3, 3)


# --------------------------------------------------------------------------------------
# AC: per-line Δv is derived from the absolute frequency relative to the metadata's
# rest_frequency by its velocity_convention — and, when that context is absent, the
# toolkit does NOT silently guess (ADR-0003): the Δv / frequency accessors raise.
# --------------------------------------------------------------------------------------
def test_delta_v_uses_the_objects_convention_and_rest_frequency(cube):
    # The optical convention gives a *different* Δv for the same line than radio — proving
    # Δv is derived through the object's stated convention, not assumed.
    coords = cube.coordinates
    radio = coords.delta_v
    nu = coords.frequency
    expected_radio = nu.to(
        u.km / u.s,
        equivalencies=u.doppler_radio(REST),
    )
    assert u.allclose(radio, expected_radio, rtol=1e-9, atol=1e-6 * u.km / u.s)


def test_velocity_axis_cube_without_rest_or_convention_refuses_to_guess(tmp_path):
    # A velocity-axis cube whose header omits the rest frequency + convention: it loads
    # (lazy), but turning the velocity axis into an AUTHORITATIVE absolute frequency needs
    # the convention astrolyze never assumes (ADR-0003) — so frequency / Δv must raise.
    h = _cube_header()
    del h["RESTFRQ"]
    del h["HIERARCH ASTROLYZE VCONV"]
    path = tmp_path / "ngc0628_no_context.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), h)
    incomplete = Cube.from_loaded(load(path))
    assert incomplete.is_complete is False
    with pytest.raises(MissingContextError):
        _ = incomplete.coordinates.frequency
    with pytest.raises(MissingContextError):
        _ = incomplete.coordinates.delta_v


# --------------------------------------------------------------------------------------
# AC: a validity descriptor exposes blanked / edge / coverage voxels as NaN + a boolean
# finite-data mask. Checked on a cube containing NaNs.
# --------------------------------------------------------------------------------------
def test_validity_exposes_nan_and_a_boolean_finite_mask(cube_with_nans):
    validity = cube_with_nans.validity

    # the finite-data mask is a boolean array the shape of the cube.
    assert validity.mask.dtype == np.bool_
    assert validity.mask.shape == cube_with_nans.shape
    # exactly the two blanked voxels are flagged not-finite; everything else is real.
    assert not validity.mask[0, 0, 0]
    assert not validity.mask[3, 2, 2]
    assert validity.mask.sum() == cube_with_nans.shape[0] * 9 - 2

    # the data view exposes the blanked voxels as NaN (and carries its unit).
    assert isinstance(validity.data, u.Quantity)
    assert validity.data.unit.is_equivalent(u.Jy / u.beam)
    assert np.isnan(validity.data.value[0, 0, 0])
    assert np.isnan(validity.data.value[3, 2, 2])
    # mask and NaN agree: the mask is exactly the finite locations of the data.
    assert np.array_equal(validity.mask, np.isfinite(validity.data.value))


def test_validity_mask_travels_through_a_subcube_slice(cube_with_nans):
    # The mask of a slice must equal the slice of the mask (it travels WITH the data).
    full_mask = cube_with_nans.validity.mask
    sub = cube_with_nans[
        :, 0:2, 0:2
    ]  # excludes the [.,2,2] blank, keeps the [0,0,0] blank
    assert isinstance(sub, Cube)
    mask_of_slice = sub.validity.mask
    slice_of_mask = full_mask[:, 0:2, 0:2]
    assert np.array_equal(mask_of_slice, slice_of_mask)
    # and the surviving blank is still flagged in the sliced descriptor.
    assert not mask_of_slice[0, 0, 0]
    assert mask_of_slice.sum() == slice_of_mask.sum()


# --------------------------------------------------------------------------------------
# AC: the validity descriptor round-trips as a Zarr companion group alongside the #23 cube
# store (mirroring the NoiseModel companion, issue #27): a uint8 finite-data mask carrying
# its derivation method + a schema version, which the loader's lean core reads (ADR-0008).
# --------------------------------------------------------------------------------------
def test_validity_zarr_companion_roundtrip(cube_with_nans, tmp_path):
    from astrolyze.core._coords import Validity

    store = cube_with_nans.to_zarr(tmp_path / "z")
    validity = cube_with_nans.validity

    group = validity.to_zarr_companion(store)
    assert group.exists()

    back = Validity.from_zarr_companion(store)
    # The mask round-trips in the cube's array order...
    assert back.mask.shape == cube_with_nans.shape
    assert np.array_equal(back.mask, validity.mask)
    # ...and stays the finite-data mask: blanked voxels are False AND NaN in the data.
    assert np.array_equal(back.mask, np.isfinite(back.data.value))
    assert not back.mask[0, 0, 0]
    assert not back.mask[3, 2, 2]


def test_validity_companion_stored_low_cost_as_uint8(cube_with_nans, tmp_path):
    # The mask is a 1-byte/voxel uint8 array in the store axis order (not a float copy of the
    # data), carrying method + version provenance like the noise companion.
    import xarray as xr

    store = cube_with_nans.to_zarr(tmp_path / "z")
    cube_with_nans.validity.to_zarr_companion(store)

    ds = xr.open_zarr(store, group="validity", zarr_format=3, consolidated=False)
    assert ds["mask"].dtype == np.uint8
    assert set(ds["mask"].dims) == {"sky_y", "sky_x", "freq"}
    assert ds.attrs["method"] == "finite_mask"
    assert int(ds.attrs["version"]) == 1


# --------------------------------------------------------------------------------------
# AC: Map exposes the sky/pixel subset; Spectrum exposes the spectral subset.
# --------------------------------------------------------------------------------------
def test_map_exposes_the_sky_pixel_coordinate_subset(cube):
    channel = cube[0]  # a 2D channel Map
    assert isinstance(channel, Map)
    coords = channel.coordinates
    # sky + pixel scale present...
    assert coords.longitude.unit.is_equivalent(u.deg)
    assert coords.latitude.unit.is_equivalent(u.deg)
    assert u.allclose(coords.pixel_scale, [2e-4, 2e-4] * u.deg, rtol=1e-9)
    assert coords.longitude.shape == (3, 3)
    # ...and no spectral axis on a 2D map.
    assert not hasattr(coords, "frequency") or coords.frequency is None


def test_map_validity_exposes_nan_and_mask(cube_with_nans):
    channel = cube_with_nans[0]  # the channel holding the [0,0] blank
    validity = channel.validity
    assert validity.mask.shape == (3, 3)
    assert not validity.mask[0, 0]
    assert np.isnan(validity.data.value[0, 0])


def test_spectrum_exposes_the_spectral_coordinate_subset(cube):
    sp = cube[:, 1, 1]
    assert isinstance(sp, Spectrum)
    coords = sp.coordinates
    # absolute frequency (authoritative) + per-line Δv present...
    assert coords.frequency.unit.is_equivalent(u.Hz)
    assert coords.delta_v.unit.is_equivalent(u.km / u.s)
    assert coords.frequency.shape == (4,)
    assert u.isclose(coords.frequency[0], REST, rtol=1e-9)
    # ...and no sky/pixel subset on a 1D spectrum.
    assert not hasattr(coords, "longitude") or coords.longitude is None


# ======================================================================================
# Issue #25 — Cube.harmonize_to_surface_brightness() → I_ν in MJy/sr as the authoritative
# representation, with derived-view accessors for K and Jy/beam, and the new physics: the
# calibration TEMPERATURE SCALE (T_mb / T_A* / T_R*) astrolyze did not model.
#
# Correctness obligation:
#   - harmonize returns a Cube in MJy/sr (specific intensity I_ν), context preserved;
#   - K↔I_ν is the PER-CHANNEL Rayleigh-Jeans relation I_ν = 2kν²/c²·T_mb at the per-channel
#     ν (NOT band-center), and K → I_ν → K is bit-stable per channel;
#   - a Kelvin cube with NO declared calibration_scale RAISES (no silent K = T_mb);
#   - T_A* without eta_mb RAISES; with eta_mb it applies T_mb = T_A*/eta_mb then converts;
#   - derived accessors return K (needs the scale + RJ/Planck law) and Jy/beam (needs a beam)
#     WITHOUT mutating the authoritative I_ν cube;
#   - a flux↔temperature crossing without the required context raises MissingContextError.
# ======================================================================================
def _kelvin_header(calibration_scale=None, eta_mb=None):
    """The reference cube header, but in Kelvin (a temperature-valued cube)."""
    h = _cube_header()
    h["BUNIT"] = "K"
    if calibration_scale is not None:
        h["HIERARCH ASTROLYZE CALSCALE"] = calibration_scale
    if eta_mb is not None:
        h["HIERARCH ASTROLYZE ETAMB"] = eta_mb
    return h


def _kelvin_cube(tmp_path, *, calibration_scale=None, eta_mb=None, fill=12.0):
    path = tmp_path / "ngc0628_co21_K.fits"
    data = np.full((4, 3, 3), fill, dtype="float64")
    fits.writeto(path, data, _kelvin_header(calibration_scale, eta_mb), overwrite=True)
    return Cube.from_loaded(load(path))


def _per_channel_MJy_sr(temperatures, frequencies):
    """The expected per-channel RJ I_ν for a stack of per-channel T_mb values, computed
    independently through the unit layer at EACH channel's own absolute frequency."""
    out = np.empty_like(np.asarray(temperatures, dtype="float64"))
    for k, nu in enumerate(frequencies):
        out[k] = convert(
            temperatures[k] * u.K,
            u.MJy / u.sr,
            rest_frequency=nu,
            temperature_scale="rayleigh_jeans",
        ).to_value(u.MJy / u.sr)
    return out * (u.MJy / u.sr)


# --------------------------------------------------------------------------------------
# AC: harmonize returns a Cube in MJy/sr (I_ν), context preserved
# --------------------------------------------------------------------------------------
def test_harmonize_returns_cube_in_MJy_sr_with_context(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb")
    harmonized = cube.harmonize_to_surface_brightness()
    assert isinstance(harmonized, Cube)
    assert harmonized.unit.is_equivalent(u.MJy / u.sr)
    assert harmonized.unit == u.MJy / u.sr
    # physical context travels through the harmonisation.
    assert u.isclose(harmonized.rest_frequency, REST, rtol=1e-12)
    assert harmonized.velocity_convention is VelocityConvention.RADIO
    assert u.isclose(harmonized.beam.major, BEAM.major, rtol=1e-12)
    assert harmonized.shape == cube.shape


# --------------------------------------------------------------------------------------
# AC: K↔I_ν uses the PER-CHANNEL RJ relation at the per-channel ν (not band-center)
# --------------------------------------------------------------------------------------
def test_harmonize_uses_per_channel_frequency_not_band_center(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb", fill=12.0)
    harmonized = cube.harmonize_to_surface_brightness()

    freqs = cube.coordinates.frequency  # authoritative per-channel ν
    expected = _per_channel_MJy_sr([12.0] * 4, freqs)
    got = harmonized._data_quantity.to(u.MJy / u.sr)
    for k in range(4):
        assert u.allclose(got[k], expected[k], rtol=1e-12)

    # The channels are at DIFFERENT frequencies, so a per-channel (not band-center) law
    # gives DIFFERENT I_ν for the same 12 K — prove the per-channel ν is actually used.
    assert not np.isclose(
        expected[0].to_value(u.MJy / u.sr), expected[3].to_value(u.MJy / u.sr)
    )


def test_harmonize_K_to_Inu_to_K_round_trip_is_bit_stable_per_channel(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb", fill=8.0)
    harmonized = cube.harmonize_to_surface_brightness()
    back = harmonized.as_kelvin(temperature_scale="rayleigh_jeans")
    # bit-stable per channel: harmonisation multiplies K by the per-channel RJ factor and the
    # Kelvin view divides I_ν by the SAME per-channel factor, so the round-trip recovers the
    # original Kelvin values exactly (np.array_equal, not just close) — issue #25.
    original = cube._data_quantity.to_value(u.K)
    recovered = back._data_quantity.to_value(u.K)
    assert np.array_equal(recovered, original)


# --------------------------------------------------------------------------------------
# AC: a temperature cube with NO declared calibration_scale RAISES (no silent K = T_mb)
# --------------------------------------------------------------------------------------
def test_harmonize_temperature_cube_without_calibration_scale_raises(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale=None)
    assert cube.unit == u.K
    with pytest.raises(MissingContextError, match="calibration_scale"):
        cube.harmonize_to_surface_brightness()


# --------------------------------------------------------------------------------------
# AC: T_A* without eta_mb RAISES; with eta_mb it applies T_mb = T_A*/eta_mb then converts
# --------------------------------------------------------------------------------------
def test_harmonize_Tastar_without_eta_mb_raises(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale="T_A*", eta_mb=None)
    with pytest.raises(MissingContextError, match="eta_mb"):
        cube.harmonize_to_surface_brightness()


def test_harmonize_Tastar_with_eta_mb_applies_main_beam_correction(tmp_path):
    eta = 0.8
    ta_star = 12.0
    cube = _kelvin_cube(tmp_path, calibration_scale="T_A*", eta_mb=eta, fill=ta_star)
    harmonized = cube.harmonize_to_surface_brightness()

    # T_mb = T_A*/eta_mb is applied BEFORE the per-channel RJ conversion.
    freqs = cube.coordinates.frequency
    expected = _per_channel_MJy_sr([ta_star / eta] * 4, freqs)
    got = harmonized._data_quantity.to(u.MJy / u.sr)
    for k in range(4):
        assert u.allclose(got[k], expected[k], rtol=1e-12)

    # And it is genuinely the corrected (larger) value, not the raw T_A*.
    raw = _per_channel_MJy_sr([ta_star] * 4, freqs)
    assert np.all(got[0] > raw[0])


# --------------------------------------------------------------------------------------
# AC: derived accessors (K, Jy/beam) view the authoritative I_ν WITHOUT mutating it
# --------------------------------------------------------------------------------------
def test_derived_accessors_do_not_mutate_the_authoritative_Inu_cube(tmp_path):
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb", fill=12.0)
    harmonized = cube.harmonize_to_surface_brightness()
    before = harmonized._data_quantity.to_value(u.MJy / u.sr).copy()

    k_view = harmonized.as_kelvin(temperature_scale="rayleigh_jeans")
    jb_view = harmonized.as_jy_per_beam()

    assert isinstance(k_view, Cube)
    assert k_view.unit.is_equivalent(u.K)
    assert isinstance(jb_view, Cube)
    assert jb_view.unit.is_equivalent(u.Jy / u.beam)

    # The authoritative I_ν cube is unchanged by either derived view.
    assert harmonized.unit == u.MJy / u.sr
    assert np.array_equal(harmonized._data_quantity.to_value(u.MJy / u.sr), before)


def test_as_kelvin_requires_the_brightness_temperature_law(tmp_path):
    # The K view needs the RJ-vs-Planck law stated; the I_ν cube has no schema field for it
    # (it is the genuinely-ambiguous one), so omitting it RAISES rather than defaulting.
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb")
    harmonized = cube.harmonize_to_surface_brightness()
    with pytest.raises(MissingContextError, match="temperature_scale"):
        harmonized.as_kelvin()


def test_as_jy_per_beam_requires_a_beam(tmp_path):
    # The Jy/beam view needs a beam; an I_ν cube with no beam must refuse rather than guess.
    cube = _kelvin_cube(tmp_path, calibration_scale="T_mb")
    harmonized = cube.harmonize_to_surface_brightness()
    # Strip the beam from the harmonized cube's metadata to model a beam-less map.
    from dataclasses import replace

    no_beam = Cube(harmonized._sc, replace(harmonized.metadata, beam=None))
    with pytest.raises(MissingContextError, match="beam"):
        no_beam.as_jy_per_beam()


# --------------------------------------------------------------------------------------
# AC: a flux↔temperature crossing without the required context raises MissingContextError
# --------------------------------------------------------------------------------------
def test_harmonize_Jy_beam_cube_uses_beam_geometry(cube):
    # The reference fixture cube is Jy/beam: harmonising it to I_ν is pure beam geometry
    # (no temperature scale needed) and yields MJy/sr.
    harmonized = cube.harmonize_to_surface_brightness()
    assert harmonized.unit == u.MJy / u.sr
    expected = convert(0.5 * u.Jy / u.beam, u.MJy / u.sr, beam=BEAM)  # value irrelevant
    assert expected.unit.is_equivalent(u.MJy / u.sr)


def test_harmonize_temperature_cube_without_rest_frequency_raises(tmp_path):
    # Crossing temperature -> surface brightness needs the rest frequency for the law; an
    # archival Kelvin cube without it must raise MissingContextError, not guess.
    h = _kelvin_header(calibration_scale="T_mb")
    del h["RESTFRQ"]
    del h["HIERARCH ASTROLYZE VCONV"]
    path = tmp_path / "ngc0628_K_no_rest.fits"
    fits.writeto(path, np.full((4, 3, 3), 12.0), h)
    cube = Cube.from_loaded(load(path))
    with pytest.raises(MissingContextError):
        cube.harmonize_to_surface_brightness()
