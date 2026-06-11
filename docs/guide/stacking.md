# Cutouts and stacking

## What it is

Two related capabilities for pulling a target out of a corpus: a single sky postage stamp, and a
container of stamps gathered from every cube that covers a target.

- `Cube.cutout(position, size)` — a spatial postage-stamp sub-`Cube` centred on a sky `SkyCoord`,
  keeping the **full** spectral axis (spectra are never cut). Lazy on a dask-backed cube.
- `Collection.stack(position_or_sources, size=…)` — gather a `Cube.cutout` from every cube that
  covers each position into a `Stack`: an aligned-cutout container that is **always safe to
  construct and browse**, however heterogeneous its members.

`Stack` is deliberately **two-stage**. Stage 1 is the browse-everything container above — no
homogeneity precondition, so a CO(2-1) single-dish stamp and an HI interferometer stamp sit side by
side and you can iterate, `filter`, `plot_grid`, and `map` over them freely. Stage 2 is the
**explicit, auditable** alignment-and-co-addition path: `to_common_beam()` → `to_velocity_grid()` →
`shift_to_rest(v_sys=…)`, then `coadd(weights="noise")`. Co-addition gates hard on physical
homogeneity and raises (`CoaddError`) rather than ever producing a meaningless average silently.

## Why it works this way

- **Cutting a stamp never cuts a spectrum.** `cutout` slices only the spatial axes; every channel
  survives, the sub-image WCS is correct, and the cube's context (beam / rest frequency / velocity
  convention) travels onto the result. `size` must be an **angular** `Quantity` (a scalar square or
  a `(height, width)` pair) — a bare number or a non-angular unit raises rather than guess
  pixels-vs-angle. An off-footprint centre, or a window that only partially overlaps under the
  default `partial="raise"`, raises rather than return a clipped or garbage stamp; opt into the
  trimmed stamp with `partial="trim"` (it is still a pure slice, never padded out to `size`). See
  [No silent physics](no-silent-physics.md).
- **Stage 1 is always safe; the physics is stage 2.** A `Stack` is a plain container — constructing,
  iterating, filtering, mapping, and plotting it have **no** homogeneity precondition, so you can
  look at everything a corpus holds for a target across surveys and lines before deciding what is
  comparable. The alignment steps are *separate, visible calls* precisely so every resampling is an
  auditable step in the analysis, not a silent transform buried inside `coadd`.
- **Each alignment step inherits a no-super-resolution guard.** `to_common_beam()` reuses
  `Cube.convolve_to_beam`, so it can only ever *lose* resolution; asking for a finer beam raises
  `LossyDirectionError`. With no target it convolves to the members' **common beam** — the smallest
  beam all of them can be smoothed to — so no member is super-resolved. `to_velocity_grid()` resamples
  onto the union range at the **coarsest** member channel width (the spectral analogue), filling
  unobserved target channels with `NaN`. See [Beam and channel matching](beam-and-channel-matching.md).
- **`v_sys` is required and resolved per member.** `shift_to_rest()` brings each member's line to
  rest (0 km/s) so the same channel samples the same rest velocity across sources. Each member's
  systemic velocity is resolved by a fixed precedence — the catalog's curated `v_sys` (the
  per-source authority), then the cube's own `Metadata.systemic_velocity`, then the value you pass.
  If a member has none at any tier it **raises**; a systemic velocity is never guessed. (A single
  scalar you pass only fills members the catalog/store did not curate — it cannot silently override
  a curated per-source value, which a scalar could never be correct for across a multi-source sample.)
- **`coadd` gates on homogeneity.** It raises `CoaddError` unless the members share one species, one
  transition, one brightness unit, **and** a compatible spatial + spectral grid; the message names
  the offending axis and the alignment step that fixes it. The only path to a valid coadd is the
  explicit one. `weights="noise"` (the default) inverse-variance weights from each member's noise
  companion; `weights="uniform"` weights equally.

## Usage

A single postage stamp, full spectral axis:

```python
import astropy.units as u
from astropy.coordinates import SkyCoord
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))
centre = SkyCoord("01h36m41.8s", "+15d47m00s")

stamp = cube.cutout(centre, 30 * u.arcsec)              # square stamp, every channel kept
stamp = cube.cutout(centre, (40, 60) * u.arcsec)        # (height, width)
stamp = cube.cutout(centre, 30 * u.arcsec, partial="trim")   # opt into a trimmed edge stamp
```

Stage 1 — gather everything that covers a target and browse it (heterogeneous, always safe):

```python
from astrolyze.collection import Collection

coll = Collection.open("/data/ism_corpus")
stack = coll.stack(centre, size=30 * u.arcsec)          # one cutout per covering cube

len(stack)                                              # member count
for member in stack:                                    # browse: identity + cube
    print(member.object, member.survey, member.species, member.transition)

stack.homogeneity_report().conflicts                    # e.g. ("species", "transition", "bunit")
fig, axes = stack.plot_grid()                           # one house-style panel per member
```

`stack()` also takes a **list** of catalog object names and/or `SkyCoord`s, building a multi-target
sample (names are resolved against *this collection's* catalog — there is no SIMBAD/Sesame in the
library; an unknown name raises `KeyError`):

```python
sample = coll.stack(["NGC3521", "NGC0628", centre], size=30 * u.arcsec)
```

Stage 2 — narrow to one kind of data, align explicitly, then co-add:

```python
co21 = stack.filter(species="CO", transition="2-1")     # the homogenising step
aligned = (
    co21
    .to_common_beam()          # smooth all to the common (largest) beam
    .to_velocity_grid()        # resample all onto one shared velocity grid
    .shift_to_rest(v_sys=850 * u.km / u.s)   # fills only members lacking a curated v_sys
)
coadded = aligned.coadd(weights="noise")     # inverse-variance from each member's noise companion
```

A heterogeneous (or unaligned) stack refuses to co-add, naming the fix:

```python
from astrolyze.collection import CoaddError

try:
    stack.coadd()              # mixed species / units, not on one grid
except CoaddError as exc:
    ...   # the message names the conflicting axis and the alignment step that fixes it
```

## A note on coadd noise propagation

`coadd(weights="noise")` uses each member's **as-published** noise companion σ. This is the exact
per-voxel σ when the members were not spatially smoothed differently on the way in, and the right
*relative* weighting (the inverse-variance mean depends on the σ ratios) otherwise. The known
limitation is that the companion is not yet re-propagated through a preceding `to_common_beam()`
convolution — so after heterogeneous smoothing the *absolute* combined σ is conservative. The cube
convolution machinery already supports analytic noise propagation; threading it through the stack is
tracked in [issue #82](https://github.com/buchbend/astrolyze/issues/82).

## See also

- [No silent physics](no-silent-physics.md) — the refuse-rather-than-guess principle behind the
  angular-`size` guard, the `partial="raise"` default, the required `v_sys`, and the `coadd` gate.
- [Beam and channel matching](beam-and-channel-matching.md) — the convolution / regrid machinery the
  alignment steps reuse, and the no-super-resolution guard they inherit.
- [Browsing a published corpus](collections.md) — `Collection.stack` / `covering`, and origin
  provenance on each member.
