"""The plotting engine: free functions that draw a wrapper onto an Axes (ADR-0005).

Decision (c) — engine + sugar: the real work lives in free functions (``plot_map(map,
ax=...)``, ``plot_spectrum``, ``plot_cube``) that take a wrapper plus an optional ``ax`` and
return ``(fig, ax)``, so they compose into any multi-panel figure the caller builds. The thin
``obj.plot()`` methods on Cube/Map/Spectrum just delegate here (see ``core/_base.py``). Object
context makes the result auto-correct: the WCS comes off the map, the colorbar is labelled
with the object's unit, and the beam ellipse is the object's beam — the caller spells out none
of it. The house style is applied locally via :func:`astrolyze.viz.style`.

A ``backend`` keyword is reserved on every signature for a future interactive backend, but
only matplotlib is built (YAGNI — ADR-0005 rider 1); any other value raises.
"""

from __future__ import annotations

import warnings

import numpy as np
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.pyplot as plt

from .style import style

#: Default colormap — perceptually uniform and legible printed on white (ADR-0005 rider 2).
DEFAULT_CMAP = "cividis"


def plot_map(
    map_obj,
    *,
    ax=None,
    cmap=DEFAULT_CMAP,
    add_beam=True,
    add_colorbar=True,
    backend="matplotlib",
    **imshow_kwargs,
):
    """Draw a :class:`~astrolyze.core.Map` on its WCS and return ``(fig, ax)``.

    The image is shown on a WCSAxes built from ``map_obj.wcs`` (a plain Axes if it has none),
    the colorbar is labelled with ``str(map_obj.unit)``, and — unless ``add_beam=False`` — the
    object's beam is drawn as an ellipse in the lower-left corner. Pass ``ax`` to compose into
    an existing figure; extra keywords flow through to :meth:`~matplotlib.axes.Axes.imshow`.
    """
    _require_backend(backend)
    with style():
        fig, ax = _resolve_axes(ax, map_obj.wcs)
        data = u.Quantity(map_obj.data)
        image = ax.imshow(np.asarray(data.value), cmap=cmap, **imshow_kwargs)
        if add_colorbar:
            cbar = fig.colorbar(image, ax=ax)
            cbar.set_label(str(map_obj.unit))
        if add_beam:
            _draw_beam(ax, map_obj)
    return fig, ax


def plot_spectrum(
    spec,
    *,
    ax=None,
    drawstyle="steps-mid",
    backend="matplotlib",
    **plot_kwargs,
):
    """Draw a :class:`~astrolyze.core.Spectrum` (flux vs spectral axis) and return ``(fig, ax)``.

    Axes are labelled from the object: the y-axis with the flux unit, the x-axis named for the
    spectral axis' physical type (Velocity / Frequency / Wavelength) and its unit. Pass ``ax``
    to overplot onto existing axes; extra keywords flow through to
    :meth:`~matplotlib.axes.Axes.plot`.
    """
    _require_backend(backend)
    with style():
        fig, ax = _resolve_axes(ax, None)
        x = spec.spectral_axis
        y = u.Quantity(spec.flux)
        ax.plot(
            np.asarray(x.value),
            np.asarray(y.value),
            drawstyle=drawstyle,
            **plot_kwargs,
        )
        ax.set_xlabel(_bracket_label(_spectral_axis_name(x.unit), x.unit))
        ax.set_ylabel(_bracket_label("Flux", spec.unit))
    return fig, ax


def plot_cube(cube, **kwargs):
    """Quick-look a :class:`~astrolyze.core.Cube` as its velocity-integrated (moment-0) map.

    A cube has no single 2D view, so ``cube.plot()`` shows the canonical summary image — the
    zeroth moment — by reusing :func:`plot_map` (DRY). For any other view, take the slice or
    moment you want and plot that (e.g. ``cube[k].plot()``, ``cube[:, y, x].plot()``).
    """
    return plot_map(cube.moment0(), **kwargs)


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def _resolve_axes(ax, wcs):
    """Return ``(fig, ax)``: reuse the caller's ``ax`` if given, else open a new figure whose
    sole Axes carries the map's WCS projection (a plain Axes when there is no WCS)."""
    if ax is not None:
        return ax.get_figure(), ax
    fig = plt.figure()
    if wcs is not None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = fig.add_subplot(111)
    return fig, ax


def _draw_beam(ax, map_obj):
    """Add the object's beam as an ellipse in the lower-left corner, sized via the WCS.

    The beam is a physical angle; turning it into pixels needs the map's pixel scale, so a map
    without a WCS cannot place it correctly — we warn and skip rather than draw a wrong beam.
    """
    beam = map_obj.beam
    if beam is None:
        return None
    wcs = map_obj.wcs
    if wcs is None:
        warnings.warn(
            "map has no WCS; cannot size the beam ellipse — omitting it",
            stacklevel=3,
        )
        return None
    pixscale = _pixel_scale(wcs)
    ny, nx = map_obj.data.shape[-2:]
    pad = 0.12  # corner inset as a fraction of the axes extent
    ellipse = beam.ellipse_to_plot(
        pad * nx,
        pad * ny,
        pixscale,
        facecolor="0.85",
        edgecolor="0.15",
        linewidth=1.25,
    )
    ax.add_patch(ellipse)
    return ellipse


def _pixel_scale(wcs) -> u.Quantity:
    """The mean celestial pixel scale as an angle Quantity (degrees per pixel), the form
    :meth:`radio_beam.Beam.ellipse_to_plot` expects."""
    celestial = wcs.celestial if hasattr(wcs, "celestial") else wcs
    scales = proj_plane_pixel_scales(celestial)  # in the WCS' own world units
    cunit = u.Unit(celestial.wcs.cunit[0] or u.deg)
    return (float(np.mean(scales)) * cunit).to(u.deg)


def _spectral_axis_name(unit) -> str:
    """A human label for the spectral axis from its physical type."""
    physical_type = str(u.Unit(unit).physical_type)
    return {
        "speed": "Velocity",
        "frequency": "Frequency",
        "length": "Wavelength",
    }.get(physical_type, "Spectral axis")


def _bracket_label(name: str, unit) -> str:
    """``"Flux [Jy / beam]"``; just ``name`` when the unit is dimensionless."""
    unit = u.Unit(unit)
    if unit == u.dimensionless_unscaled:
        return name
    return f"{name} [{unit}]"


def _require_backend(backend: str) -> None:
    """Guard the reserved backend seam: only matplotlib is built (ADR-0005 rider 1)."""
    if backend != "matplotlib":
        raise NotImplementedError(
            f"viz backend {backend!r} is not built; only 'matplotlib' is available. "
            "The backend seam is reserved for a future interactive backend (ADR-0005)."
        )
