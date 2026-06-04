# Beam and channel matching

*Embodies [ADR-0003](../adr/0003-unit-layer-pure-astropy-explicit-conventions.md) (no silent
physics) and [ADR-0004](../adr/0004-object-model-cube-map-spectrum-wrappers.md) (context-
carrying wrappers).*

## What it is

`Cube` offers four resolution-matching operations, each returning a new context-carrying
`Cube`:

- `convolve_to_beam(beam)` — smooth spatially to a **larger** beam.
- `spectral_bin(factor)` — bin the spectral axis by an integer `factor > 1` (fewer, wider
  channels).
- `spectral_smooth_to(width)` — smooth the spectral axis to a **broader** resolution `width`
  (a FWHM).
- `match_to(other, *, reproject=False)` — bring two cubes to a common beam (the smallest beam
  both can be smoothed to); returns `(matched_self, matched_other)`.

All of them only ever *lose* resolution. Asking for the gaining direction — a smaller beam,
finer channels — raises `LossyDirectionError`.

## Why it works this way

- **Never super-resolve silently (ADR-0003).** Deconvolving to a smaller beam or up-sampling to
  finer channels would invent structure the data never contained. astrolyze refuses it with
  `LossyDirectionError` (a `ValueError` — it reads as the programming error it is) rather than
  producing a plausible-looking but fabricated result. The lossy-direction check uses
  `radio_beam`'s deconvolve to decide whether a target beam is genuinely reachable by
  convolution.
- **astrolyze stays thin (ADR-0004).** `spectral-cube` does the convolution/binning maths and
  `radio_beam` decides which direction is resolution-losing. The value astrolyze adds is the
  guard plus the **context carry**: the new, larger beam is recorded on the returned cube's
  metadata, so downstream operations see the correct resolution.
- **Reprojection is explicit, never a side effect.** `match_to` brings both cubes to a common
  beam but only regrids onto a shared pixel grid when you pass `reproject=True`. Line-ratio
  work must opt in to regridding; it never happens implicitly.

## Usage

Spatial smoothing to a coarser beam:

```python
import astropy.units as u
import radio_beam
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))

target = radio_beam.Beam(major=15 * u.arcsec, minor=15 * u.arcsec, pa=0 * u.deg)
coarser = cube.convolve_to_beam(target)   # new Cube; carries `target` as its beam
```

Spectral coarsening — by integer factor, or to a broader FWHM:

```python
binned = cube.spectral_bin(2)                  # half as many channels, twice as wide
smoothed = cube.spectral_smooth_to(2 * u.km / u.s)   # smooth to a 2 km/s FWHM
```

The gaining direction is refused:

```python
from astrolyze.core import LossyDirectionError

try:
    cube.spectral_bin(1)        # a no-op / up-sample direction
except LossyDirectionError:
    ...   # astrolyze never invents a finer grid
```

Matching two cubes to a common beam (regridding only on opt-in):

```python
co21 = Cube.from_loaded(load("ngc0628_co21.fits"))
co10 = Cube.from_loaded(load("ngc0628_co10.fits"))

a, b = co21.match_to(co10)                    # common beam, spectral axes untouched
a, b = co21.match_to(co10, reproject=True)    # also share one pixel grid (for line ratios)
```

## See also

- [No silent physics](no-silent-physics.md) — the same refuse-rather-than-guess principle for
  unit conversions.
- [Coordinates and validity](coordinates-and-validity.md) — the channel grid and beam these
  operations preserve and carry.
