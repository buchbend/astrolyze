"""Tests for astrolyze.viz — the plotting engine + house style (ADR-0005).

Written first (red/green TDD). These are **contract** tests, not pixel comparisons (issue #6):

- the engine is free functions (`plot_map(map, ax=...)`, `plot_spectrum`, `plot_cube`) that
  return ``(fig, ax)`` and compose into any figure the caller supplies;
- a map is drawn on its WCS, with the **beam ellipse** taken from the object's beam and a
  colorbar **labelled with the object's unit** (``str(unit)``);
- the house style is applied **locally**: importing astrolyze / astrolyze.viz never mutates
  matplotlib's global ``rcParams``, and a plot call leaves them exactly as it found them;
- the matplotlib backend is the only one built; the reserved backend seam refuses others
  (YAGNI — ADR-0005 rider 1).

The reference dataset mirrors the tracer-bullet PHANGS cube (NGC 0628, CO(2-1), 12"x10" beam,
radio convention, Jy/beam) so the beam + unit context the plots read is real.
"""

import re
import warnings

import matplotlib

matplotlib.use("Agg")  # headless; set before pyplot is imported via astrolyze.viz
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
import astrolyze.viz as viz

warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)


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
    path = tmp_path / "ngc0628_co21.fits"
    # a non-trivial gradient so the colorbar / line have real range to draw.
    data = np.arange(4 * 5 * 5, dtype="float32").reshape((4, 5, 5))
    fits.writeto(path, data, _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture
def integrated_map(cube):
    return cube.moment0()


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# --------------------------------------------------------------------------------------
# MANDATED: plot_map returns (fig, ax)
# --------------------------------------------------------------------------------------
def test_plot_map_returns_fig_and_ax(integrated_map):
    fig, ax = viz.plot_map(integrated_map)
    assert isinstance(fig, mpl.figure.Figure)
    assert isinstance(ax, mpl.axes.Axes)  # WCSAxes is an Axes subclass
    assert ax.get_figure() is fig


# --------------------------------------------------------------------------------------
# MANDATED: a beam artist is drawn from the object's beam
# --------------------------------------------------------------------------------------
def test_plot_map_draws_beam_ellipse_from_object_beam(integrated_map):
    fig, ax = viz.plot_map(integrated_map)
    ellipses = [p for p in ax.patches if isinstance(p, Ellipse)]
    assert len(ellipses) == 1  # exactly the beam, drawn once


def test_plot_map_beam_can_be_disabled(integrated_map):
    fig, ax = viz.plot_map(integrated_map, add_beam=False)
    assert not any(isinstance(p, Ellipse) for p in ax.patches)


# --------------------------------------------------------------------------------------
# MANDATED: colorbar label == str(object's unit)
# --------------------------------------------------------------------------------------
def test_plot_map_colorbar_labelled_with_object_unit(integrated_map):
    fig, ax = viz.plot_map(integrated_map)
    unit_str = str(integrated_map.unit)
    assert unit_str != ""  # a moment0 of Jy/beam has a real composite unit
    # the colorbar lives on its own Axes; its label is exactly the object's unit string.
    other_labels = [cax.get_ylabel() for cax in fig.axes if cax is not ax]
    assert unit_str in other_labels


# --------------------------------------------------------------------------------------
# MANDATED: default cmap is cividis, always overridable
# --------------------------------------------------------------------------------------
def test_plot_map_default_cmap_is_cividis_and_overridable(integrated_map):
    _, ax = viz.plot_map(integrated_map)
    assert ax.images[0].get_cmap().name == "cividis"
    _, ax2 = viz.plot_map(integrated_map, cmap="magma")
    assert ax2.images[0].get_cmap().name == "magma"


# --------------------------------------------------------------------------------------
# MANDATED: import never mutates global rcParams; style is applied locally + restored
# --------------------------------------------------------------------------------------
def test_import_does_not_mutate_global_rcparams():
    # importing astrolyze / astrolyze.viz must leave the house cmap OUT of the global state.
    assert mpl.rcParams["image.cmap"] == mpl.rcParamsDefault["image.cmap"]
    assert mpl.rcParams["image.cmap"] != "cividis"


def test_plot_restores_global_rcparams(integrated_map):
    before = mpl.rcParams.copy()
    viz.plot_map(integrated_map)
    after = mpl.rcParams.copy()
    assert before == after  # style.context restored every key on exit


def test_style_applies_house_sheet_only_inside_the_context():
    from astrolyze.viz import style

    assert mpl.rcParams["image.cmap"] != "cividis"
    with style():
        assert mpl.rcParams["image.cmap"] == "cividis"  # house default, locally
    assert mpl.rcParams["image.cmap"] != "cividis"  # restored on exit


def test_astrolyze_style_is_reachable_from_top_level():
    import astrolyze

    with astrolyze.style():
        assert mpl.rcParams["image.cmap"] == "cividis"


# --------------------------------------------------------------------------------------
# DONE WHEN: composing two plot_map calls into one figure works
# --------------------------------------------------------------------------------------
def test_two_plot_map_calls_compose_into_one_figure(integrated_map):
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={"projection": integrated_map.wcs})
    f1, a1 = viz.plot_map(integrated_map, ax=ax1)
    f2, a2 = viz.plot_map(integrated_map, ax=ax2)
    # both drew onto the caller's figure/axes — no new figure was created.
    assert f1 is fig and f2 is fig
    assert a1 is ax1 and a2 is ax2
    assert any(isinstance(p, Ellipse) for p in ax1.patches)
    assert any(isinstance(p, Ellipse) for p in ax2.patches)


# --------------------------------------------------------------------------------------
# plot_spectrum: free function over the 1D wrapper
# --------------------------------------------------------------------------------------
def test_plot_spectrum_returns_fig_ax_and_draws_a_line(cube):
    sp = cube[:, 2, 2]
    fig, ax = viz.plot_spectrum(sp)
    assert isinstance(fig, mpl.figure.Figure)
    assert isinstance(ax, mpl.axes.Axes)
    assert len(ax.lines) == 1
    # the y-axis carries the flux unit; the x-axis the spectral-axis unit.
    assert str(sp.unit) in ax.get_ylabel()
    assert str(sp.spectral_axis.unit) in ax.get_xlabel()


def test_plot_spectrum_composes_onto_supplied_axes(cube):
    sp = cube[:, 1, 1]
    fig, ax = plt.subplots()
    f, a = viz.plot_spectrum(sp, ax=ax)
    assert f is fig and a is ax


# --------------------------------------------------------------------------------------
# plot_cube: quick-look that reuses plot_map on the moment-0 map (DRY)
# --------------------------------------------------------------------------------------
def test_plot_cube_quicklooks_the_moment0_map(cube):
    fig, ax = viz.plot_cube(cube)
    assert isinstance(fig, mpl.figure.Figure)
    # it is the integrated map: one image + the beam ellipse, unit is the moment0 unit.
    assert len(ax.images) == 1
    assert any(isinstance(p, Ellipse) for p in ax.patches)
    unit_str = str(cube.moment0().unit)
    assert unit_str in [cax.get_ylabel() for cax in fig.axes if cax is not ax]


# --------------------------------------------------------------------------------------
# the object-method sugar now lights up via the engine (ADR-0005)
# --------------------------------------------------------------------------------------
def test_object_plot_sugar_delegates_to_engine(integrated_map):
    fig, ax = integrated_map.plot()
    assert isinstance(fig, mpl.figure.Figure)
    assert any(isinstance(p, Ellipse) for p in ax.patches)


# --------------------------------------------------------------------------------------
# the reserved backend seam refuses anything but matplotlib (YAGNI, ADR-0005 rider 1)
# --------------------------------------------------------------------------------------
def test_unbuilt_backend_is_refused(integrated_map):
    with pytest.raises(NotImplementedError, match="backend"):
        viz.plot_map(integrated_map, backend="plotly")


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in the viz library code (ADR-0005 / the legacy sin)
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_viz_library_code():
    from pathlib import Path

    import astrolyze.viz as viz_pkg

    pkg_dir = Path(viz_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
