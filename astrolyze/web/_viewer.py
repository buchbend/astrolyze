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

The seams for #68 (region-averaged spectrum, velocity-window moment, linked-zoom) are deliberately
*not* built: the four basic panels need exactly these four reads. #68 will add endpoints alongside
these without reshaping them.
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
    JSON, and the same "blanked, not zero" meaning the validity machinery uses). ``object``-dtype
    via ``where`` keeps the nested-list shape while substituting ``None`` for the non-finite mask."""
    array = np.asarray(array, dtype="float64")
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


__all__ = [
    "MAX_MAP_DIM",
    "cube_axes",
    "moment0_map",
    "channel_slice",
    "pixel_spectrum",
]
