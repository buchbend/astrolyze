# NEXT SESSION — start here

Single entry point for resuming the astrolyze build in a fresh session.

## Where you are

- **This repo (`~/git/astrolyze`, GitHub `buchbend/astrolyze`, public, BSD-3)** = the Core
  Toolkit. Scaffold + units + io + core are done; `pip install -e ".[dev]"` works, `pytest` is
  65/65 green (25 units, 25 io, 3 smoke, 12 core). Venv at `.venv/`.
- **Design lives in the PRIVATE repo `~/git/science/ppv-foundation`** (two-repo split,
  ADR-0007). Read these there before building:
  - `CONTEXT.md` — glossary / shared language.
  - `docs/adr/0001`–`0013` — the decisions that govern everything.
  - `docs/prd/0001-astrolyze-tracer-bullet.md` — the full PRD.
  - `docs/astrolyze-build-plan.md` — module layout + build order.
- The 2012-2016 original is preserved at **buchbend/astrolyze-legacy** (heritage; don't touch).

## What to build next — issue #6: `viz` (plotting engine + house style)

Tracer-bullet, **units-first**. GitHub issues (label `ready-for-agent`):

```
#1 PRD (tracking)
#2 scaffold ........ DONE (closed)
#3 units ........... DONE (merged, #3; 25 tests)                    ADR-0003
#4 io + schema ..... DONE (merged via PR; 25 tests)                 ADR-0006
#5 core ............ DONE (merged to main; 12 tests)                ADR-0004
#6 viz ............. NEXT  (needs #5)                                ADR-0005
#7 tracer + CLI .... needs #3-#6 (PRD acceptance)                   ADR-0011/0012
```

**#3 units — DONE** (`astrolyze/units/`, `tests/test_units.py`): pure-astropy, no I/O, no
import side effects. Aliases (`Tmb`, `Ta`, `Jy_beam`, `Jy_sr`, `MJy_sr`, `K_kms`);
equivalency builders (`brightness_temperature` with **both** Rayleigh-Jeans *and* an exact
**Planck** equivalency — astropy's built-in is RJ-only; `beam_angular_area`, `doppler`,
`spectral`); and `convert(quantity, target, *, rest_frequency, convention, beam,
temperature_scale)`. **Design notes for whoever picks this up:**
- Added `temperature_scale` (rayleigh_jeans|planck) to the mandated `convert` signature —
  RJ-vs-Planck is the #1 silent-error trap (ADR-0003), so it is explicit/mandatory like the
  velocity convention, never defaulted. Strings or enums accepted.
- Velocity-integrated intensity (K km/s ↔ Jy/beam km/s) is handled manually (astropy won't
  compose the equivalency across the km/s factor) and is **RJ-only by construction** —
  asking for Planck on an integrated quantity raises (B_ν is nonlinear; ∫B_ν(T)dv ≠ B_ν(∫T)).
- astropy does **not chain** two equivalencies in one `.to()`, so each conversion is built
  as exactly one registered pair (e.g. Jy/beam↔K folds the beam in, rather than
  beam_angular_area + surface-brightness).
- `#3` is **merged to `main`** (commit `Implement units module … (#3)`).

**#4 io + schema — DONE** (`astrolyze/io/`, `tests/test_io.py`; ADR-0006). The FITS header is
authoritative; the filename is a derived projection; loading is lazy. **Notes for #5:**
- `load(path) -> LoadedData(data, header, wcs, metadata, path)` — the seam #5 consumes to build
  `Cube`/`Map`/`Spectrum`. Byte I/O delegated to astropy; io adds the schema, not a reader.
- `Metadata` (`io/schema.py`) carries object, telescope, species, rest_frequency (Hz),
  velocity_convention, beam (`radio_beam.Beam`), bunit (`u.Unit`), distance (with unit),
  calibration_error, name_tag, + a `provenance=None` seam. `.is_complete` / `.missing` /
  `.ensure_complete()` — the latter raises the **shared** `MissingContextError` from
  `astrolyze.units` (reused, not re-defined). **#5 should call `metadata.ensure_complete()`
  at the door of any unit/velocity op.**
- Keyword schema: standard keys where they exist (`OBJECT`,`TELESCOP`,`BUNIT`,`RESTFRQ`,
  `BMAJ/BMIN/BPA`); everything astrolyze-specific under `HIERARCH ASTROLYZE …` (`VCONV`,
  `SPECIES`,`DISTANCE`+`DISTUNIT`,`CALERR`,`NAMETAG`,`SCHEMA`). **Incompleteness triggers =
  rest_frequency + velocity_convention** only (ADR-0006 ii); beam-ops are gated by units.
- Velocity convention is *read, not guessed*: `VCONV` keyword first, else WCS `CTYPE3`
  (`VRAD`→radio / `VOPT`→optical / `VELO`→relativistic).
- Filename projection (`io/naming.py`, legacy v1 grammar):
  `source_telescope_species_fluxunit_resolution.fits`; resolution from beam in arcsec
  (`12.00` / `12.00x10.00` / `…a30.0`); missing fields → `unknown`. `save(...)` writes a
  **derived file only** — raw inputs are never renamed.
- I added `coerce_velocity_convention` to the `astrolyze.units` public exports (io needs it).

**#5 core — DONE** (`astrolyze/core/`, `tests/test_core.py`; ADR-0004), merged to `main`
(commit `Implement core module … (#5)`; 12 tests). Thin `Cube`/`Map`/`Spectrum` that
**compose** (private handles: `Cube._sc` = SpectralCube, `Spectrum._spec` = specutils.Spectrum,
`Map._data` + `Map._wcs`) and carry the io `Metadata` as the single context source. **Notes for
#6 viz:**
- A shared `ContextCarrier` mixin (`core/_base.py`) gives all three `.to()`, `.plot()`, and the
  context shortcuts `.beam` / `.rest_frequency` / `.velocity_convention` / `.is_complete` /
  `.missing` (+ `require_complete()`, the door for velocity/spectral-axis ops).
- **`.plot()` is already wired** and is the seam #6 fills: it does `from astrolyze import viz`
  then `getattr(viz, self._viz_function)(self, **kwargs)`. The names it expects are **`plot_cube`**
  (Cube), **`plot_map`** (Map), **`plot_spectrum`** (Spectrum) — define those free functions in
  `astrolyze.viz` and `.plot()` lights up automatically; until they exist it raises a clear
  `NotImplementedError`.
- Free functions take the wrapper object (`plot_map(map, ax=...)`). Available for drawing:
  `map.wcs` (2D WCS), `map.data` (Quantity → colorbar unit), `map.beam` (`radio_beam.Beam` for
  the ellipse), `cube.spectral_axis`, `spectrum.spectral_axis` / `spectrum.flux`.
- Apply the house style **locally** (context manager); **never** mutate global rcParams on
  import (ADR-0005). Importing core already pulls matplotlib via spectral-cube but does NOT
  touch rcParams — keep it that way.
- `.to()` supplies beam/rest-freq/convention from the object but **never** defaults the
  RJ-vs-Planck `temperature_scale` (ADR-0003): pass it explicitly for brightness-temperature
  conversions. Type transitions: `Cube.moment0()` / `moment(order,axis)` → `Map`;
  `cube[:, y, x]` → `Spectrum`; `cube[k]` → `Map`; sub-cube slice → `Cube`. Build a cube with
  `Cube.from_loaded(io.load(path))`.

## Non-negotiable house rules (from the ADRs)

- Stay thin — delegate to spectral-cube/astropy; don't reimplement moments/convolution.
- No silent physics; no `SystemExit`/`print` in library code; never mutate global rcParams.
- Tests are the correctness obligation. **Verify claims against real artifacts** (incl. docs).
- Bash on this host needs `dangerouslyDisableSandbox: true` for real work (seccomp).

## Parked (not blocking)

- Stray duplicate `0005-ppv-velocity-2d-resolution-ladder.md` ADR in ppv-foundation (decide
  keep/renumber/remove).
- Peer bibliography in `~/private/w2/papers/` — themes F (surveys) & G (below-beam stats)
  still empty (arXiv rate-limited; needs a slow retry).
