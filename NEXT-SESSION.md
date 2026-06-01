# NEXT SESSION — start here

Single entry point for resuming the astrolyze build in a fresh session.

## Where you are

- **This repo (`~/git/astrolyze`, GitHub `buchbend/astrolyze`, public, BSD-3)** = the Core
  Toolkit. Scaffold is done and pushed; `pip install -e ".[dev]"` works, `pytest` is 3/3
  green (incl. the no-global-rcParams-mutation guard). Venv at `.venv/`.
- **Design lives in the PRIVATE repo `~/git/science/ppv-foundation`** (two-repo split,
  ADR-0007). Read these there before building:
  - `CONTEXT.md` — glossary / shared language.
  - `docs/adr/0001`–`0013` — the decisions that govern everything.
  - `docs/prd/0001-astrolyze-tracer-bullet.md` — the full PRD.
  - `docs/astrolyze-build-plan.md` — module layout + build order.
- The 2012-2016 original is preserved at **buchbend/astrolyze-legacy** (heritage; don't touch).

## What to build next — issue #4: `io + schema`

Tracer-bullet, **units-first**. GitHub issues (label `ready-for-agent`):

```
#1 PRD (tracking)
#2 scaffold ........ DONE (closed)
#3 units ........... DONE (28 tests green; not yet pushed/closed)   ADR-0003
#4 io + schema ..... NEXT  (was parallel with #3)                   ADR-0006
#5 core ............ needs #3 + #4                                  ADR-0004
#6 viz ............. needs #5                                       ADR-0005
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
- `#3` is committed locally on a branch; **not yet pushed / issue not yet closed** (awaiting
  go-ahead).

**#4 io + schema** (ADR-0006): header metadata schema is authoritative; lazy `load()` (loads
even if rest-freq/convention missing, but flags incomplete — contrast the strict Ingest
gate); filename projection synced from header. No reimplementing astropy/spectral-cube I/O.

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
