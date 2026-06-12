"""Cube-viewer JSON serializers — lazy slices of a corpus cube for D3 (issues #67/#68).

The detail-view "open viewer" affordance opens a store's Zarr into a lazy, dask-backed
:class:`~astrolyze.core.Cube` (via the public :meth:`Record.open`) and the four viewer endpoints
serve **bounded JSON slices** of it for client-side D3 rendering — never a full-cube payload, never
server-side matplotlib. Every helper here is a *thin* wrapper over the public Cube API (the same
surface a Python user has), so astrolyze stays thin (AGENTS.md): the physics (moment-0, the
velocity axis, the WCS sky extents, single-pixel spectra) is the library's, this module only
projects the result onto a JSON shape D3 can draw.

Laziness is the contract (PRD #56, user story 8). A 100 GB corpus cube must never materialise per
request:

- :func:`cube_axes` reads only headers/axes (shape, velocity axis, WCS sky extents, units, value
  range hint) — no data plane is touched.
- :func:`channel_slice` computes **one** channel (``cube[index]`` — a 2-D dask slice), not the
  cube.
- :func:`pixel_spectrum` computes **one** pixel column (``cube[:, y, x]`` — a 1-D dask slice).
- :func:`moment0_map` is the one endpoint that reduces over the spectral axis; it goes through the
  public :meth:`Cube.moment0`, which is itself lazy over the dask graph until the (single, bounded)
  2-D result is read.

JSON-shape discipline (shared with :mod:`astrolyze.web.api`):

- numpy arrays serialize to plain nested lists; a NaN/inf becomes ``null`` (:func:`_array_to_json`)
  — invalid JSON never goes on the wire, and ``null`` reads as "blanked", the same honesty the rest
  of the package uses (ADR-0003: an unknown value is never guessed).
- a 2-D map that would be larger than :data:`MAX_MAP_DIM` on a side is **decimated by striding**
  (every k-th pixel) to keep the payload bounded; the response says so (``downsample`` factor) and
  the sky ``extent`` still spans the full map, so the client draws the true footprint. Striding (not
  averaging) keeps it a pure, cheap slice — no new resampling maths here.

The #68 interactions (region-averaged spectrum, velocity-window moment, linked-zoom) build *on top*
of the four basic reads without reshaping them. Two add server-side reductions, still thin over the
public Cube API:

- :func:`region_spectrum` averages the cube over a user-drawn polygon — composed from public
  spatial slicing (``cube[:, y0:y1, x0:x1]`` → a :class:`Cube`) plus a pure-numpy point-in-polygon
  mask and ``numpy.nanmean`` over the masked spatial pixels per channel. astrolyze has no
  first-class polygon-average method, so this is the documented composition; it stays lazy by
  slicing to the polygon's bounding box first (dask pulls only those columns).
- :func:`windowed_moment0_map` recomputes moment-0 over exactly the channels whose velocity falls in
  ``[vmin, vmax]`` — a public **spectral slice by channel index** (``cube[i0:i1]`` → a
  :class:`Cube`) handed to the same public :meth:`Cube.moment0`. Slicing by channel index (not
  spectral_slab) means it works on a cube without a rest-frequency velocity axis too (it reuses the
  viewer's own velocity-or-native axis), and it stays lazy (only the slab's channels are read).

The polygon coordinate convention is **image pixels** (full-resolution ``x`` = column, ``y`` = row),
the same coordinates the maps' click handler reports and the spectrum/channel endpoints take — so a
region drawn on a (possibly decimated) heatmap maps back to the same full-resolution grid the rest
of the viewer speaks, with no sky-WCS round-trip on the hot path.
"""

from __future__ import annotations

import math

import numpy as np

# A 2-D map larger than this on either side is decimated by striding before it goes on the wire, so
# a multi-thousand-pixel corpus map stays a bounded JSON payload the browser can hold and D3 can
# draw. 512 is a generous heatmap resolution on screen; the sky extent always spans the full map, so
# decimation costs detail, never the footprint. A square stride keeps pixels isotropic on screen.
MAX_MAP_DIM = 512


def _finite_or_none(value):
    """A finite float, else ``None`` — NaN/inf is not valid JSON and reads as "unknown"."""
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _array_to_json(array) -> list:
    """A numpy array → nested Python lists with every non-finite value mapped to ``None``.

    JSON has no NaN/inf; a blanked voxel (the cube's ``NaN``) must serialize as ``null`` (valid
    JSON, and the same "blanked, not zero" meaning the validity machinery uses).

    Fast path: when the map has NO non-finite pixel (the common case — most channel maps are fully
    finite) ``ndarray.tolist()`` is the C-level converter to plain Python floats (valid JSON, no
    ``None`` needed). Only when a blank exists do we pay for the ``object``-dtype ``where`` that
    substitutes ``None`` for the non-finite mask — a slow path because object dtype defeats the C
    converter, so it is reserved for the case that actually needs it. Output is identical either
    way; this only avoids the object-dtype detour when nothing is blanked."""
    array = np.asarray(array, dtype="float64")
    if np.isfinite(array).all():
        return array.tolist()
    return np.where(np.isfinite(array), array, None).tolist()


def _decimate_factor(ny: int, nx: int) -> int:
    """The integer stride that brings a (ny, nx) map under :data:`MAX_MAP_DIM` on each side.

    ``1`` when the map already fits (the common corpus-cutout case). Otherwise the smallest factor
    ``k`` such that ``ceil(dim / k) <= MAX_MAP_DIM`` for the larger side — a cheap, isotropic
    decimation that bounds the payload without inventing interpolated values."""
    largest = max(ny, nx)
    if largest <= MAX_MAP_DIM:
        return 1
    return math.ceil(largest / MAX_MAP_DIM)


def _sky_extent(wcs, ny: int, nx: int) -> dict:
    """The map's sky bounding box in degrees, from the *full*-resolution corner pixels.

    Two opposite corners (pixel ``(0, 0)`` and ``(nx-1, ny-1)``) projected through the celestial
    WCS give the RA/Dec span the heatmap covers — enough to label the D3 axes without shipping a
    per-pixel coordinate grid. Read off the public ``Map.wcs`` (the moment/channel WCS), so the
    extent is the library's WCS, not a reimplementation. ``None`` when no celestial WCS is present
    (a cube without sky calibration — honest, not guessed)."""
    if wcs is None:
        return None
    try:
        celestial = wcs.celestial
        corners = celestial.wcs_pix2world([[0, 0], [nx - 1, ny - 1]], 0)
    except Exception:
        # A non-celestial or malformed WCS yields no sky extent rather than a wrong one (ADR-0003).
        return None
    (ra0, dec0), (ra1, dec1) = corners
    return {
        "ra_min_deg": _finite_or_none(min(ra0, ra1)),
        "ra_max_deg": _finite_or_none(max(ra0, ra1)),
        "dec_min_deg": _finite_or_none(min(dec0, dec1)),
        "dec_max_deg": _finite_or_none(max(dec0, dec1)),
    }


def _velocity_values(cube) -> tuple[list, str]:
    """The cube's velocity axis as (values, unit) — km/s where the cube carries the context.

    Goes through the public :meth:`Cube.velocity_axis` (frequency-authoritative, raises if the rest
    frequency / convention is absent — astrolyze never guesses a velocity axis, ADR-0003). If that
    context is missing the spectral axis is reported in its own native unit instead, so a cube
    without a rest line still drives the slider (by channel, with a frequency/native label) rather
    than failing the whole viewer."""
    try:
        axis = cube.velocity_axis()
        return [_finite_or_none(v) for v in np.asarray(axis.value)], str(axis.unit)
    except Exception:
        axis = cube.spectral_axis
        return [_finite_or_none(v) for v in np.asarray(axis.value)], str(axis.unit)


# -- the four lazy slice serializers ---------------------------------------------------
def cube_axes(cube) -> dict:
    """Axis metadata to drive the viewer UI — **no data plane is read** (laziness, story 8).

    Shape, channel count, the velocity axis (values + unit), the sky extent for the map, the data
    unit (bunit), and a *coarse* value-range hint. Everything here comes from headers / WCS / the
    spectral axis, so opening a 100 GB corpus cube and asking for its axes touches no chunk. The
    value range is intentionally a coarse hint (sampled from the middle channel only, itself a
    single lazy 2-D slice) so the client can seed a colour scale without a full-cube reduction; the
    moment-0 endpoint returns the authoritative per-map vmin/vmax."""
    nchan, ny, nx = cube.shape
    velocity, vel_unit = _velocity_values(cube)
    # A cheap range hint: one mid channel (a single lazy 2-D slice), not the whole cube.
    mid = nchan // 2
    mid_map = np.asarray(cube[mid].data.value, dtype="float64")
    finite = mid_map[np.isfinite(mid_map)]
    vmin = _finite_or_none(finite.min()) if finite.size else None
    vmax = _finite_or_none(finite.max()) if finite.size else None
    return {
        "object": cube.metadata.object,
        "shape": [int(nchan), int(ny), int(nx)],
        "n_channels": int(nchan),
        "velocity": velocity,
        "velocity_unit": vel_unit,
        "bunit": str(cube.unit),
        "sky_extent": _sky_extent(cube[mid].wcs, ny, nx),
        "value_min_hint": vmin,
        "value_max_hint": vmax,
    }


def moment0_map(cube) -> dict:
    """The integrated (moment-0) map as a 2-D JSON array + extent/unit/vmin/vmax.

    Dogfoods the public :meth:`Cube.moment0` (spectral-cube does the integral; the unit becomes
    e.g. K·km/s). The single 2-D result is read once, decimated by striding if it would exceed
    :data:`MAX_MAP_DIM`, and serialized with NaN→null. The sky ``extent`` is computed at *full*
    resolution so the decimated array still spans the true footprint on screen."""
    moment = cube.moment0()
    full = np.asarray(moment.data.value, dtype="float64")
    ny, nx = full.shape
    factor = _decimate_factor(ny, nx)
    shown = full[::factor, ::factor]
    finite = shown[np.isfinite(shown)]
    return {
        "data": _array_to_json(shown),
        "shape": [int(shown.shape[0]), int(shown.shape[1])],
        "unit": str(moment.unit),
        "extent": _sky_extent(moment.wcs, ny, nx),
        "downsample": int(factor),
        "vmin": _finite_or_none(finite.min()) if finite.size else None,
        "vmax": _finite_or_none(finite.max()) if finite.size else None,
    }


def channel_slice(cube, index: int) -> dict:
    """One channel as a 2-D JSON array + its velocity — a **single lazy 2-D slice** (``cube[index]``).

    ``cube[index]`` is a 2-D dask slice through the public ``__getitem__``; only that channel's
    chunks are pulled. Decimated/serialized exactly like the moment-0 map so the two heatmaps share
    a pixel grid and a sky extent. The channel velocity is read off the (already-computed) velocity
    axis so the slider/keyboard panel can display "v = … km/s" without another request. An
    out-of-range *index* raises :class:`IndexError` (the API maps it to a graceful 4xx)."""
    nchan = cube.shape[0]
    if not 0 <= index < nchan:
        raise IndexError(
            f"channel {index} out of range for a cube with {nchan} channels (0..{nchan - 1})"
        )
    channel = cube[index]
    full = np.asarray(channel.data.value, dtype="float64")
    ny, nx = full.shape
    factor = _decimate_factor(ny, nx)
    shown = full[::factor, ::factor]
    velocity, vel_unit = _velocity_values(cube)
    return {
        "index": int(index),
        "data": _array_to_json(shown),
        "shape": [int(shown.shape[0]), int(shown.shape[1])],
        "unit": str(channel.unit),
        "extent": _sky_extent(channel.wcs, ny, nx),
        "downsample": int(factor),
        "velocity": velocity[index] if index < len(velocity) else None,
        "velocity_unit": vel_unit,
    }


def pixel_spectrum(cube, x: int, y: int) -> dict:
    """The spectrum at pixel (x, y) as (velocity, value) arrays — a **single lazy 1-D slice**.

    ``cube[:, y, x]`` is a 1-D dask slice (one pixel column) through the public ``__getitem__``;
    only that column's chunks are pulled, never the cube. Returned as parallel ``velocity`` /
    ``value`` arrays (length == ``n_channels``) plus units, the natural shape for a D3 line. Pixel
    coordinates are image-order ``x`` (column) / ``y`` (row), matching the click position the maps
    report. An out-of-range pixel raises :class:`IndexError` (→ graceful 4xx)."""
    nchan, ny, nx = cube.shape
    if not (0 <= x < nx and 0 <= y < ny):
        raise IndexError(
            f"pixel (x={x}, y={y}) out of range for a {nx}×{ny} map (x 0..{nx - 1}, y 0..{ny - 1})"
        )
    spectrum = cube[:, y, x]
    values = np.asarray(spectrum.flux.value, dtype="float64")
    velocity, vel_unit = _velocity_values(cube)
    return {
        "x": int(x),
        "y": int(y),
        "velocity": velocity,
        "velocity_unit": vel_unit,
        "value": [_finite_or_none(v) for v in values],
        "value_unit": str(spectrum.unit),
    }


# -- #68 interactions: region-averaged spectrum + velocity-window moment ----------------
def _polygon_pixel_mask(vertices, ny: int, nx: int):
    """A boolean (ny, nx) mask of pixels whose CENTRE lies inside the polygon *vertices*.

    *vertices* is a sequence of ``(x, y)`` image-pixel pairs (x = column, y = row); the polygon is
    implicitly closed (last vertex → first). Membership is the standard even-odd ray-casting test
    (a horizontal ray from each pixel centre crosses the polygon edges an odd number of times when
    inside) — implemented in pure numpy so the web path pulls in no extra geometry dependency and
    the rule is auditable here. A pixel centre exactly on an edge is treated consistently by the
    half-open ``(y0 <= y) != (y1 <= y)`` edge test; the result is intended as a fill mask, not a
    sub-pixel-exact aperture (region averaging over whole pixels, the honest map resolution)."""
    verts = np.asarray(vertices, dtype="float64")
    xs = verts[:, 0]
    ys = verts[:, 1]
    # Pixel-centre coordinate grids (centre of pixel (col=c, row=r) is at (c, r) in this convention).
    col = np.arange(nx, dtype="float64")[None, :]
    row = np.arange(ny, dtype="float64")[:, None]
    inside = np.zeros((ny, nx), dtype=bool)
    n = len(verts)
    j = n - 1
    for i in range(n):
        xi, yi = xs[i], ys[i]
        xj, yj = xs[j], ys[j]
        # Edge straddles the pixel row (half-open in y so a shared vertex is counted once).
        straddles = (yi <= row) != (yj <= row)
        # x of the edge at this row; guard a horizontal edge (yj == yi) where straddles is False.
        denom = yj - yi
        denom = np.where(denom == 0.0, 1.0, denom)
        x_cross = xi + (row - yi) / denom * (xj - xi)
        inside ^= straddles & (col < x_cross)
        j = i
    return inside


def region_spectrum(cube, vertices) -> dict:
    """The cube averaged over a polygon *region* → a (velocity, mean value, n_pixels) spectrum.

    *vertices* is a list of ``[x, y]`` image-pixel pairs (full-resolution column/row, the same
    coordinates the maps report). Composed from the public Cube API (astrolyze has no first-class
    polygon-average): slice the cube to the polygon's integer bounding box (``cube[:, y0:y1, x0:x1]``
    — a :class:`Cube`, so dask pulls only those columns and the average stays lazy), build the
    even-odd :func:`_polygon_pixel_mask` over that box, and take ``numpy.nanmean`` over the masked
    spatial pixels per channel (NaN voxels are excluded — a blanked pixel does not drag the mean to
    zero, ADR-0003). The velocity axis is the viewer's own (km/s where the cube carries the line
    context, else native), so the region spectrum shares the pixel spectrum's x-axis exactly.

    Returns ``velocity`` / ``value`` arrays (length == ``n_channels``) plus ``n_pixels`` (how many
    in-polygon pixels were averaged). A degenerate region (< 3 vertices) raises :class:`ValueError`
    (the API maps it to a graceful ``422``); a polygon that selects no pixel inside the map yields a
    spectrum of ``null`` with ``n_pixels == 0`` (honest emptiness, not a 500)."""
    if vertices is None or len(vertices) < 3:
        raise ValueError(
            f"a region needs at least 3 vertices to enclose an area, got "
            f"{0 if vertices is None else len(vertices)} (draw a polygon, not a point/line)"
        )
    nchan, ny, nx = cube.shape
    verts = np.asarray(vertices, dtype="float64")
    # Integer bounding box, clamped to the map — slice here so only these columns are pulled.
    x0 = max(0, int(math.floor(verts[:, 0].min())))
    y0 = max(0, int(math.floor(verts[:, 1].min())))
    x1 = min(nx, int(math.ceil(verts[:, 0].max())) + 1)
    y1 = min(ny, int(math.ceil(verts[:, 1].max())) + 1)
    velocity, vel_unit = _velocity_values(cube)
    if x1 <= x0 or y1 <= y0:
        # The polygon lies entirely outside the map: no pixels, an honest empty spectrum.
        return {
            "velocity": velocity,
            "velocity_unit": vel_unit,
            "value": [None] * nchan,
            "value_unit": str(cube.unit),
            "n_pixels": 0,
        }
    sub = cube[:, y0:y1, x0:x1]
    block = np.asarray(
        sub._data_quantity.value, dtype="float64"
    )  # (nchan, dy, dx), lazy slab
    # Mask in the bounding box's local frame (shift the polygon by the box origin).
    local = verts - np.array([x0, y0], dtype="float64")
    mask = _polygon_pixel_mask(local, y1 - y0, x1 - x0)
    n_pixels = int(mask.sum())
    if n_pixels == 0:
        return {
            "velocity": velocity,
            "velocity_unit": vel_unit,
            "value": [None] * nchan,
            "value_unit": str(cube.unit),
            "n_pixels": 0,
        }
    # Per-channel nanmean over the in-polygon spatial pixels; all-NaN channels → None (RuntimeWarning
    # suppressed because an all-blanked channel is a legitimate "unknown", not an error).
    region = block[:, mask]  # (nchan, n_pixels)
    with np.errstate(invalid="ignore"):
        all_nan = ~np.isfinite(region).any(axis=1)
        means = np.where(
            all_nan,
            np.nan,
            np.nanmean(np.where(np.isfinite(region), region, np.nan), axis=1),
        )
    return {
        "velocity": velocity,
        "velocity_unit": vel_unit,
        "value": [_finite_or_none(v) for v in means],
        "value_unit": str(cube.unit),
        "n_pixels": n_pixels,
    }


def _window_channel_range(velocity, vmin: float, vmax: float) -> tuple[int, int]:
    """The half-open channel index range ``[i0, i1)`` whose velocities fall within ``[vmin, vmax]``.

    *velocity* is the viewer's velocity-or-native axis (a list, possibly with ``None`` holes). The
    bounds are ordered defensively (a window dragged right-to-left is the same window). Raises
    :class:`ValueError` when no channel falls inside (an empty/degenerate window) so the API can
    return a graceful ``422`` rather than handing spectral-cube an empty slab."""
    lo, hi = (vmin, vmax) if vmin <= vmax else (vmax, vmin)
    selected = [i for i, v in enumerate(velocity) if v is not None and lo <= v <= hi]
    if not selected:
        raise ValueError(
            f"velocity window [{lo}, {hi}] selects no channel "
            f"(the axis spans the cube's own velocity range) — widen the window"
        )
    return selected[0], selected[-1] + 1


def windowed_moment0_map(cube, vmin: float, vmax: float) -> dict:
    """Moment-0 recomputed over **exactly** the channels in the velocity window ``[vmin, vmax]``.

    The velocity-window panel: select the channels whose (viewer) velocity is in ``[vmin, vmax]``,
    slice the cube to that contiguous channel range (``cube[i0:i1]`` — a public :class:`Cube` slice,
    so only the slab's channels are read and the integral stays lazy), and hand the slab to the same
    public :meth:`Cube.moment0` the full-map endpoint uses (dogfooding — the windowed map is the
    full map's integral restricted to the slab, not a reimplemented sum). Serialized exactly like
    :func:`moment0_map` (NaN→null, decimated if large, full-resolution sky extent) plus the
    ``vmin``/``vmax`` of the window and the channel span that realised it.

    An empty or invalid window (no channel inside) raises :class:`ValueError` (→ graceful ``422``)."""
    velocity, _ = _velocity_values(cube)
    i0, i1 = _window_channel_range(velocity, vmin, vmax)
    payload = moment0_map(cube[i0:i1])
    payload["vmin_window"] = _finite_or_none(min(vmin, vmax))
    payload["vmax_window"] = _finite_or_none(max(vmin, vmax))
    payload["channel_start"] = int(i0)
    payload["channel_stop"] = int(i1)  # exclusive
    return payload


__all__ = [
    "MAX_MAP_DIM",
    "cube_axes",
    "moment0_map",
    "channel_slice",
    "pixel_spectrum",
    "region_spectrum",
    "windowed_moment0_map",
]
