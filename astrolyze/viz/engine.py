"""The plotting engine: free functions that draw a wrapper onto an Axes (ADR-0005).

Decision (c) — engine + sugar: the real work lives in free functions (``plot_map(map,
ax=...)``, ``plot_spectrum``, ``plot_cube``, ``plot_channel_maps``) that take a wrapper plus
an optional ``ax`` and return ``(fig, ax)`` / ``(fig, axes)``, so they compose into any
multi-panel figure the caller builds. The thin ``obj.plot()`` / ``obj.plot_channel_maps()``
methods on Cube/Map/Spectrum just delegate here (see ``core/_base.py``). Object context makes
the result auto-correct: the WCS comes off the map, the colorbar is labelled with the
object's unit, and the beam ellipse is the object's beam — the caller spells out none of it.
The house style is applied locally via :func:`astrolyze.viz.style`.

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


def plot_noise(model, **kwargs):
    """Quick-look a :class:`~astrolyze.core.NoiseModel` as its spatial σ map (issue #27).

    A noise model has several products; the canonical 2D summary is the spatial noise pattern
    σ_xy, so ``model.plot()`` shows that by reusing :func:`plot_map` (DRY). For the spectral
    σ_v or the autocorrelation, plot the corresponding product (``model.sigma_spectrum.plot()``,
    ``model.spectral_acf.plot()``).
    """
    return plot_map(model.sigma_map, **kwargs)


def plot_channel_maps(
    cube,
    *,
    start=None,
    stop=None,
    step=1,
    v_min=None,
    v_max=None,
    ncols=None,
    cmap=DEFAULT_CMAP,
    add_beam=True,
    vmin=None,
    vmax=None,
    backend="matplotlib",
):
    """Draw a publication-grade channel-map grid for a :class:`~astrolyze.core.Cube`.

    Returns ``(fig, axes)`` where ``axes`` is a 1-D :class:`numpy.ndarray` of Axes, one per
    selected channel. Every panel shows the channel image on WCS axes, labelled with the
    channel velocity (or frequency) as its title. A shared colorbar (on the right side of
    the figure) is labelled with the cube's unit. Unless ``add_beam=False``, the beam ellipse
    is drawn in the lower-left corner of every panel.

    Channel selection
    -----------------
    There are two mutually exclusive selection modes:

    **Index mode** — specify a Python-style integer range with any combination of ``start``,
    ``stop``, ``step``. The range is ``cube[start:stop:step]``. Defaults to all channels
    when no selection keyword is given.

    **Velocity mode** — specify ``v_min`` AND ``v_max`` (both required together) as
    :class:`~astropy.units.Quantity` with velocity units. All channels whose spectral-axis
    value falls within ``[v_min, v_max]`` (inclusive) are selected; ``step`` is then applied
    to thin that subset further.

    Parameters
    ----------
    cube :
        A :class:`~astrolyze.core.Cube` to display.
    start : int, optional
        First channel index (inclusive). Default: 0.
    stop : int, optional
        Last channel index (exclusive). Default: ``len(spectral_axis)``.
    step : int, optional
        Step between selected channels within the chosen range. Must be ≥ 1. Default: 1.
    v_min, v_max : :class:`~astropy.units.Quantity`, optional
        Velocity (or frequency) bounds for velocity-mode selection. Both must be supplied
        together; omitting only one raises :class:`ValueError`.
    ncols : int, optional
        Number of panel columns in the grid. Default: ``min(n_panels, 5)``.
    cmap : str, optional
        Colormap name. Default: ``"cividis"`` (the house default — ADR-0005 rider 2).
    add_beam : bool, optional
        Draw the beam ellipse in each panel's lower-left corner. Default: ``True``.
    vmin, vmax : float, optional
        Shared display range for all panels. When ``None`` (the default), the range is set
        from the full data of all selected channels so the colorbar is consistent across
        panels.
    backend : str, optional
        Reserved seam for a future interactive backend (ADR-0005 rider 1). Only
        ``"matplotlib"`` is built; any other value raises :class:`NotImplementedError`.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    axes : :class:`numpy.ndarray` of :class:`matplotlib.axes.Axes`
        One element per selected channel, in selection order.

    Raises
    ------
    ValueError
        ``start`` or ``stop`` out of range; ``start >= stop``; ``step < 1``; a velocity
        window that selects no channels; ``v_min`` supplied without ``v_max`` or vice versa.
    NotImplementedError
        ``backend`` is not ``"matplotlib"``.
    """
    _require_backend(backend)

    spectral_axis = cube.spectral_axis  # Quantity array, length n_chan
    n_chan = len(spectral_axis)

    # -- resolve channel indices ----------------------------------------------------------
    if v_min is not None or v_max is not None:
        # Velocity mode — both bounds must be given together.
        if v_min is None or v_max is None:
            raise ValueError(
                "v_min and v_max must both be supplied for velocity-range selection; "
                "provide both bounds or use index mode (start/stop/step)"
            )
        channel_indices = _channels_in_velocity_range(spectral_axis, v_min, v_max, step)
    else:
        # Index mode — validate and build the slice.
        channel_indices = _channels_from_index_range(n_chan, start, stop, step)

    n_panels = len(channel_indices)

    # -- layout ---------------------------------------------------------------------------
    if ncols is None:
        ncols = min(n_panels, 5)
    nrows = int(np.ceil(n_panels / ncols))

    # Determine shared display range from the selected channels when not caller-supplied.
    # This makes the shared colorbar physically meaningful — the same colour maps to the
    # same intensity in every panel.
    if vmin is None or vmax is None:
        selected_data = np.stack(
            [np.asarray(cube[k].data) for k in channel_indices], axis=0
        )
        finite_data = selected_data[np.isfinite(selected_data)]
        _vmin = float(finite_data.min()) if vmin is None else vmin
        _vmax = float(finite_data.max()) if vmax is None else vmax
    else:
        _vmin, _vmax = float(vmin), float(vmax)

    # Spatial WCS from a representative channel — all channels share the celestial axes.
    spatial_wcs = _spatial_wcs(cube)

    # One channel-map panel per subplot; the colorbar gets its own inset Axes via fig.colorbar.
    with style():
        # Each panel is square-ish; leave a little extra width on the right for the colorbar.
        panel_size = 2.2  # inches per panel
        fig_width = ncols * panel_size + 0.8  # 0.8 for the colorbar column
        fig_height = nrows * panel_size
        fig = plt.figure(figsize=(fig_width, fig_height))

        axes = []
        for i, k in enumerate(channel_indices):
            ax = fig.add_subplot(
                nrows,
                ncols,
                i + 1,
                projection=spatial_wcs,
            )
            channel_map = cube[k]  # a Map — carries unit + beam
            data = np.asarray(channel_map.data)
            img = ax.imshow(
                data,
                cmap=cmap,
                vmin=_vmin,
                vmax=_vmax,
                interpolation="nearest",
            )
            # Velocity label: round to 2 decimal places, state the unit.
            vel = spectral_axis[k]
            ax.set_title(_velocity_label(vel))
            # Suppress redundant axis labels on interior panels; keep them on the bottom-left
            # panel only. WCSAxes always shows the coordinate label; hide tick labels on
            # non-edge panels to avoid clutter in dense grids.
            if i % ncols != 0:
                ax.coords[1].set_ticklabel_visible(False)
            if i < (nrows - 1) * ncols:
                ax.coords[0].set_ticklabel_visible(False)
            if add_beam:
                _draw_beam(ax, channel_map)
            axes.append(ax)

        # Shared colorbar: attach to the rightmost column of panels. matplotlib's fig.colorbar
        # with ax=[...] builds one colorbar Axes spanning all supplied panels — that is the
        # "shared" colorbar pattern (single Axes, single scale for all panels).
        cbar = fig.colorbar(img, ax=axes, fraction=0.046, pad=0.04)
        cbar.set_label(str(cube.unit))

        # tight_layout is not compatible with WCSAxes (it warns and may misplace panels).
        # subplots_adjust gives the same breathing room without triggering the warning.
        fig.subplots_adjust(hspace=0.35, wspace=0.15)

    return fig, np.array(axes)


def plot_stack_grid(
    stack,
    *,
    ncols=None,
    cmap=DEFAULT_CMAP,
    add_beam=True,
    backend="matplotlib",
):
    """Draw one house-style panel per :class:`~astrolyze.collection.stack.Stack` member (#64).

    The browse-everything view (PRD #56 user story 12): a grid with one panel per stack member,
    each the member cube's canonical 2-D summary — its moment-0 map, the same view ``cube.plot()``
    shows for a cube — on its own WCS axes, titled with the member's identity (object / survey /
    species / transition). It reuses :func:`plot_map` (DRY) and the shipped house style, so a
    stack grid matches the channel-map grid's look.

    Crucially this works on a **heterogeneous** stack: members may carry different brightness
    units, so each panel gets its **own** colorbar labelled with that member's unit (a single
    shared scale across mixed units would be physically meaningless) — that is the deliberate
    difference from :func:`plot_channel_maps`, where all panels share one scale.

    Parameters
    ----------
    stack : ~astrolyze.collection.stack.Stack
        The stack to display. Must have at least one member.
    ncols : int, optional
        Number of panel columns. Default: ``min(n_members, 4)``.
    cmap : str, optional
        Colormap name. Default: ``"cividis"`` (the house default — ADR-0005 rider 2).
    add_beam : bool, optional
        Draw each member's beam ellipse in its panel's lower-left corner. Default: ``True``.
    backend : str, optional
        Reserved seam for a future interactive backend (ADR-0005 rider 1). Only ``"matplotlib"``.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    axes : :class:`numpy.ndarray` of :class:`matplotlib.axes.Axes`
        One element per member, in member order.

    Raises
    ------
    ValueError
        The stack is empty (nothing to plot is a caller mistake, not a silent empty figure).
    NotImplementedError
        ``backend`` is not ``"matplotlib"``.
    """
    _require_backend(backend)
    members = list(stack)
    n_panels = len(members)
    if n_panels == 0:
        raise ValueError(
            "cannot plot_grid an empty stack (no members to draw); a stack with no covering "
            "cutouts has nothing to show — filter or re-stack before plotting"
        )

    if ncols is None:
        ncols = min(n_panels, 4)
    nrows = int(np.ceil(n_panels / ncols))

    with style():
        panel_size = 2.6  # inches per panel — a touch larger than channel maps (own colorbar each)
        fig = plt.figure(figsize=(ncols * panel_size + 0.8, nrows * panel_size))
        axes = []
        for i, member in enumerate(members):
            # The member's canonical 2-D summary is its moment-0 map (same as cube.plot()); each
            # panel is an independent plot_map so it carries its own WCS, unit colorbar, and beam.
            summary = member.cube.moment0()
            ax = fig.add_subplot(nrows, ncols, i + 1, projection=summary.wcs)
            plot_map(
                summary,
                ax=ax,
                cmap=cmap,
                add_beam=add_beam,
                add_colorbar=True,
            )
            ax.set_title(member._label())
            axes.append(ax)
        # subplots_adjust (not tight_layout) — WCSAxes warns under tight_layout; generous spacing
        # because every panel carries its own colorbar.
        fig.subplots_adjust(hspace=0.4, wspace=0.45)

    return fig, np.array(axes)


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


def _channels_from_index_range(
    n_chan: int,
    start,
    stop,
    step: int,
) -> list[int]:
    """Validate index-mode channel selection and return a list of channel indices.

    Raises :class:`ValueError` on out-of-range or inverted bounds, or a non-positive step.
    """
    if not isinstance(step, (int, np.integer)) or step < 1:
        raise ValueError(
            f"step must be a positive integer, got {step!r}; "
            "use step >= 1 to select every Nth channel"
        )
    _start = 0 if start is None else int(start)
    _stop = n_chan if stop is None else int(stop)

    if _start < 0 or _start >= n_chan:
        raise ValueError(
            f"start={_start} is out of range for a cube with {n_chan} channels "
            f"(valid: 0 .. {n_chan - 1})"
        )
    if _stop < 1 or _stop > n_chan:
        raise ValueError(
            f"stop={_stop} is out of range for a cube with {n_chan} channels "
            f"(valid: 1 .. {n_chan})"
        )
    if _start >= _stop:
        raise ValueError(
            f"start={_start} must be less than stop={_stop}; the selected range is empty"
        )
    return list(range(_start, _stop, int(step)))


def _channels_in_velocity_range(
    spectral_axis,
    v_min,
    v_max,
    step: int,
) -> list[int]:
    """Return channel indices whose spectral-axis value falls within [v_min, v_max].

    Both bounds are converted to the spectral axis' own unit before comparison, so the
    caller can pass any velocity/frequency unit compatible with the axis. After filtering,
    ``step`` thins the matching subset. Raises :class:`ValueError` when no channels match or
    when the step is non-positive.
    """
    if not isinstance(step, (int, np.integer)) or step < 1:
        raise ValueError(f"step must be a positive integer, got {step!r}")
    axis_unit = spectral_axis.unit
    try:
        lo = v_min.to_value(axis_unit, equivalencies=u.spectral())
        hi = v_max.to_value(axis_unit, equivalencies=u.spectral())
    except u.UnitConversionError as exc:
        raise ValueError(
            f"v_min/v_max units {v_min.unit}/{v_max.unit} are not compatible with the "
            f"cube's spectral axis unit {axis_unit}: {exc}"
        ) from exc

    # Velocity axes may be in ascending or descending order; normalise the comparison.
    _lo, _hi = min(lo, hi), max(lo, hi)
    values = spectral_axis.to_value(axis_unit)
    matching = [k for k, v in enumerate(values) if _lo <= v <= _hi]
    if not matching:
        raise ValueError(
            f"No channels found in velocity range [{v_min}, {v_max}]; the cube's spectral "
            f"axis spans [{values.min():.4g} .. {values.max():.4g}] {axis_unit}"
        )
    return matching[:: int(step)]


def _velocity_label(vel) -> str:
    """A compact title string for a channel, e.g. ``"0.00 km/s"`` or ``"230.54 GHz"``.

    The spectral axis unit determines the display unit and number of significant figures.
    Velocities are shown in km/s; frequencies in GHz; anything else in its native unit.
    """
    unit = vel.unit
    physical_type = str(u.Unit(unit).physical_type)

    if "speed" in physical_type or "velocity" in physical_type:
        # Display in km/s for compact, human-readable velocity labels.
        display_val = vel.to_value(u.km / u.s)
        return f"{display_val:.2f} km/s"
    if "frequency" in physical_type:
        display_val = vel.to_value(u.GHz)
        return f"{display_val:.4f} GHz"
    # Fallback: use the native unit as-is.
    return f"{vel.value:.4g} {unit}"


def _spatial_wcs(cube):
    """Extract the 2D celestial WCS from the cube's 3D WCS.

    Uses the ``celestial`` attribute of the cube's internal SpectralCube WCS when available
    (WCS objects carry a ``.celestial`` property that drops the spectral axis). Falls back
    to ``None`` (plain Axes) when the celestial sub-WCS cannot be extracted.
    """
    wcs_3d = cube._sc.wcs
    if wcs_3d is None:
        return None
    celestial = getattr(wcs_3d, "celestial", None)
    if celestial is not None:
        return celestial
    # A SpectralCube's slice WCS is accessible via a representative channel map.
    try:
        chan = cube[0]
        return getattr(chan, "wcs", None)
    except Exception:
        return None
