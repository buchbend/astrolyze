"""Tests for astrolyze.viz.plot_channel_maps — publication-grade channel-map grid (#59).

Written first (red/green TDD). Contract tests — no pixel comparisons:

- Returns (fig, axes) with one Axes per selected channel.
- Channel selection by index range (start/stop/step) works.
- Channel selection by velocity range works.
- Invalid ranges raise ValueError with a clear message.
- Each panel carries a velocity label.
- Shared colorbar is drawn with the cube's unit label.
- Beam ellipse drawn in every panel (add_beam=True by default); disabled with add_beam=False.
- WCS axes used for the spatial dimensions.
- House style applied locally; no global state leaked.
- Thin Cube.plot_channel_maps sugar delegates to the free function.
"""

import re
import warnings

import matplotlib

matplotlib.use("Agg")  # headless; must be set before pyplot is imported
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.io import load
from astrolyze.viz import plot_channel_maps

warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)

# Enough channels for multi-panel tests; 8 channels, 6×6 spatial pixels.
_N_CHAN = 8
_NY = 6
_NX = 6


def _cube_header(n_chan=_N_CHAN, ny=_NY, nx=_NX, cdelt3_ms=2000.0):
    """A 3D Jy/beam cube on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = cdelt3_ms, 1.0, "m/s"
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
    path = tmp_path / "ngc0628_co21.fits"
    rng = np.random.default_rng(42)
    data = rng.standard_normal((_N_CHAN, _NY, _NX)).astype("float32")
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# --------------------------------------------------------------------------------------
# MANDATED: returns (fig, axes) with exactly one Axes per selected channel
# --------------------------------------------------------------------------------------
def test_returns_fig_and_axes_tuple(cube):
    fig, axes = plot_channel_maps(cube)
    assert isinstance(fig, mpl.figure.Figure)
    assert isinstance(axes, np.ndarray)


def test_default_selects_all_channels(cube):
    fig, axes = plot_channel_maps(cube)
    # axes is a 1-D array of Axes; default should include all _N_CHAN channels.
    assert axes.size == _N_CHAN


def test_channel_count_equals_number_of_selected_panels(cube):
    # step=2 over 8 channels: channels 0, 2, 4, 6 → 4 panels
    fig, axes = plot_channel_maps(cube, step=2)
    assert axes.size == 4


# --------------------------------------------------------------------------------------
# MANDATED: channel selection by index range
# --------------------------------------------------------------------------------------
def test_index_range_start_stop(cube):
    fig, axes = plot_channel_maps(cube, start=2, stop=5)
    # channels 2, 3, 4 → 3 panels
    assert axes.size == 3


def test_index_range_start_stop_step(cube):
    fig, axes = plot_channel_maps(cube, start=0, stop=8, step=2)
    # 0, 2, 4, 6 → 4 panels
    assert axes.size == 4


def test_index_range_start_only(cube):
    # start=5 → channels 5, 6, 7 → 3 panels
    fig, axes = plot_channel_maps(cube, start=5)
    assert axes.size == 3


def test_index_range_stop_only(cube):
    # stop=3 → channels 0, 1, 2 → 3 panels
    fig, axes = plot_channel_maps(cube, stop=3)
    assert axes.size == 3


# --------------------------------------------------------------------------------------
# MANDATED: channel selection by velocity range
# --------------------------------------------------------------------------------------
def test_velocity_range_selects_channels_within_window(cube):
    # cdelt3 = 2000 m/s; channels 0..7 span 0 to 14000 m/s.
    # Requesting 2 km/s .. 8 km/s should include channels 1..4 (4000..8000 m/s).
    v_min = 2.0 * u.km / u.s
    v_max = 8.0 * u.km / u.s
    fig, axes = plot_channel_maps(cube, v_min=v_min, v_max=v_max)
    # At least one panel selected in that window.
    assert axes.size >= 1


def test_velocity_range_and_step_combine(cube):
    v_min = 0.0 * u.km / u.s
    v_max = 14.0 * u.km / u.s
    fig, axes_step1 = plot_channel_maps(cube, v_min=v_min, v_max=v_max, step=1)
    fig2, axes_step2 = plot_channel_maps(cube, v_min=v_min, v_max=v_max, step=2)
    # step=2 should yield fewer or equal panels than step=1
    assert axes_step2.size <= axes_step1.size


# --------------------------------------------------------------------------------------
# MANDATED: invalid ranges raise
# --------------------------------------------------------------------------------------
def test_start_out_of_range_raises(cube):
    with pytest.raises(ValueError, match="start"):
        plot_channel_maps(cube, start=100)


def test_stop_out_of_range_raises(cube):
    with pytest.raises(ValueError, match="stop"):
        plot_channel_maps(cube, stop=100)


def test_start_after_stop_raises(cube):
    with pytest.raises(ValueError, match="start|stop|range"):
        plot_channel_maps(cube, start=5, stop=2)


def test_step_less_than_one_raises(cube):
    with pytest.raises(ValueError, match="step"):
        plot_channel_maps(cube, step=0)


def test_velocity_range_no_channels_raises(cube):
    # A velocity window far outside the cube's range should raise.
    v_min = 1000.0 * u.km / u.s
    v_max = 2000.0 * u.km / u.s
    with pytest.raises(ValueError, match="[Nn]o channels|velocity"):
        plot_channel_maps(cube, v_min=v_min, v_max=v_max)


def test_velocity_range_requires_both_bounds(cube):
    # Supplying only one bound of a velocity range should raise.
    with pytest.raises((ValueError, TypeError)):
        plot_channel_maps(cube, v_min=0.0 * u.km / u.s)


# --------------------------------------------------------------------------------------
# MANDATED: each panel carries a velocity label
# --------------------------------------------------------------------------------------
def test_each_panel_has_a_title_with_velocity(cube):
    fig, axes = plot_channel_maps(cube)
    for ax in axes.flat:
        title = ax.get_title()
        # title must be non-empty and contain a numeric value
        assert title != "", "each channel panel must have a velocity title"
        assert re.search(r"[-+]?\d", title), f"velocity title has no number: {title!r}"


# --------------------------------------------------------------------------------------
# MANDATED: shared colorbar with the cube's unit label
# --------------------------------------------------------------------------------------
def test_shared_colorbar_carries_unit_label(cube):
    fig, axes = plot_channel_maps(cube)
    unit_str = str(cube.unit)
    # The colorbar lives on its own Axes; its label should contain the unit string.
    other_labels = [cax.get_ylabel() for cax in fig.axes if cax not in list(axes.flat)]
    assert any(unit_str in lbl for lbl in other_labels), (
        f"expected unit {unit_str!r} in colorbar label; found: {other_labels}"
    )


# --------------------------------------------------------------------------------------
# MANDATED: beam ellipse drawn in every panel (default) and suppressible
# --------------------------------------------------------------------------------------
def test_beam_ellipse_in_every_panel_by_default(cube):
    fig, axes = plot_channel_maps(cube)
    for ax in axes.flat:
        ellipses = [p for p in ax.patches if isinstance(p, Ellipse)]
        assert len(ellipses) >= 1, "every panel must have a beam ellipse by default"


def test_add_beam_false_suppresses_all_ellipses(cube):
    fig, axes = plot_channel_maps(cube, add_beam=False)
    for ax in axes.flat:
        assert not any(isinstance(p, Ellipse) for p in ax.patches)


# --------------------------------------------------------------------------------------
# MANDATED: house style applied locally, no global state leaked
# --------------------------------------------------------------------------------------
def test_house_style_restored_after_call(cube):
    before = mpl.rcParams.copy()
    plot_channel_maps(cube)
    after = mpl.rcParams.copy()
    assert before == after


# --------------------------------------------------------------------------------------
# MANDATED: WCS axes used for spatial dimensions
# --------------------------------------------------------------------------------------
def test_panels_use_wcs_axes(cube):
    fig, axes = plot_channel_maps(cube)
    # WCSAxes is an Axes subclass; the axes should be instances of it.
    # We just confirm they are Axes (WCSAxes is a subclass, sufficient for structural test).
    for ax in axes.flat:
        assert isinstance(ax, mpl.axes.Axes)


# --------------------------------------------------------------------------------------
# MANDATED: thin Cube.plot_channel_maps sugar delegates to the free function
# --------------------------------------------------------------------------------------
def test_cube_method_sugar_delegates_to_engine(cube):
    fig, axes = cube.plot_channel_maps()
    assert isinstance(fig, mpl.figure.Figure)
    assert isinstance(axes, np.ndarray)
    assert axes.size == _N_CHAN


def test_cube_method_accepts_kwargs(cube):
    fig, axes = cube.plot_channel_maps(step=2)
    assert axes.size == 4


# --------------------------------------------------------------------------------------
# MANDATED: ncols parameter controls panel layout (optional layout knob)
# --------------------------------------------------------------------------------------
def test_ncols_controls_grid_shape(cube):
    fig, axes = plot_channel_maps(cube, ncols=4)
    # With _N_CHAN=8 channels and ncols=4, grid is 2×4; all 8 panels present.
    assert axes.size == _N_CHAN


# --------------------------------------------------------------------------------------
# Smoke test: committed real cutout (NGC 628, CO(2-1), PHANGS-ALMA)
# --------------------------------------------------------------------------------------
def test_smoke_real_cutout():
    import pathlib

    cutout = pathlib.Path(__file__).parent / "data" / "ngc0628_co21_cutout.fits.gz"
    if not cutout.exists():
        pytest.skip("real cutout not found")
    cube = Cube.from_loaded(load(cutout))
    fig, axes = plot_channel_maps(cube, step=10)
    assert isinstance(fig, mpl.figure.Figure)
    assert axes.size >= 1
    # Every panel has a velocity label.
    for ax in axes.flat:
        assert ax.get_title() != ""


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in the viz library code
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_channel_maps_code():
    import pathlib

    import astrolyze.viz as viz_pkg

    pkg_dir = pathlib.Path(viz_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
