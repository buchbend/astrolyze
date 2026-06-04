# Coordinates and validity

*Embodies [ADR-0004](../adr/0004-object-model-cube-map-spectrum-wrappers.md) (thin wrappers
that surface, not re-derive) and
[ADR-0003](../adr/0003-unit-layer-pure-astropy-explicit-conventions.md) (no silent physics).*

## What it is

The core objects surface the WCS as first-class, unit-carrying value objects rather than bare
arrays. A `Cube` exposes two read-only properties:

- `Cube.coordinates` → an `AxisCoordinates` with the per-axis physical coordinate arrays, each
  a `Quantity` carrying its own unit:
  - `frequency` — the **authoritative** absolute frequency per channel;
  - `delta_v` — the per-line velocity offset of each channel from the rest frequency;
  - `longitude` / `latitude` — the sky coordinate maps;
  - `pixel_scale` — the celestial pixel scale.
- `Cube.validity` → a `Validity` descriptor: `data` (the values with blanked / edge /
  outside-coverage voxels exposed as `NaN`) plus `mask` (a boolean array, `True` exactly where
  the data is real).

`Map` exposes only the sky/pixel subset of `AxisCoordinates`; `Spectrum` only the spectral
subset. Absent fields are simply `None`.

## Why it works this way

- **Surfacing, not re-deriving (ADR-0004).** These properties read what the WCS / spectral
  axis astrolyze already parsed — no header reparse, no new geometry. The wrapper turns
  upstream-parsed axes into named, context-carrying values so callers never juggle a bare
  array or its unit.
- **The absolute frequency is authoritative; Δv is the matching inverse.** `frequency` is read
  from the spectral axis under the object's *stated* velocity convention and rest frequency;
  `delta_v` is derived from that same absolute frequency. If the spectral axis is a velocity
  and that context is absent, deriving an absolute frequency would require guessing the
  convention — so astrolyze **raises** `MissingContextError` rather than assuming one
  (ADR-0003). See [No silent physics](no-silent-physics.md).
- **The validity mask travels with the data.** Both `data` and `mask` are read from the same
  array, so slicing the wrapper slices both consistently: the mask of a subcube slice equals
  the slice of the mask. A consumer always knows where the data is real without re-deriving it.

## Usage

```python
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))

coords = cube.coordinates
coords.frequency      # absolute frequency per channel (a Quantity, e.g. in Hz)
coords.delta_v        # per-line velocity offset (a Quantity, e.g. in km/s)
coords.pixel_scale    # celestial pixel scale (angular Quantity, one per spatial axis)

valid = cube.validity
valid.mask            # boolean array, True where the data is finite
valid.data            # the values, with blanked voxels as NaN (carries the data's unit)
```

Because the mask is derived from the data, it slices for free:

```python
import numpy as np

subcube = cube[10:20, :, :]
np.array_equal(subcube.validity.mask, cube.validity.mask[10:20, :, :])  # True
```

If the cube's spectral axis is a velocity but the rest frequency / convention is absent,
asking for the absolute frequency raises rather than guessing:

```python
from astrolyze.units import MissingContextError

try:
    cube.coordinates.frequency
except MissingContextError as exc:
    ...  # the message names exactly what is missing
```

## See also

- [No silent physics](no-silent-physics.md) — the guard `coordinates` relies on, and how to
  probe before acting.
- [I/O and metadata](io-and-metadata.md) — where the velocity convention and rest frequency
  come from.
