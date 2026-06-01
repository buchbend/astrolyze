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

## What to build next — issue #3: `units`

Tracer-bullet, **units-first**. GitHub issues (label `ready-for-agent`):

```
#1 PRD (tracking)
#2 scaffold ........ DONE (closed)
#3 units ........... NEXT  (deep module; tests first)   ADR-0003
#4 io + schema ..... parallel with #3                    ADR-0006
#5 core ............ needs #3 + #4                        ADR-0004
#6 viz ............. needs #5                             ADR-0005
#7 tracer + CLI .... needs #3-#6 (PRD acceptance)         ADR-0011/0012
```

**#3 units** (ADR-0003): pure-astropy, no I/O. Named radio aliases (Tmb, Jy/beam, K km/s,
MJy/sr); equivalency bundles (brightness temperature with correct **Planck** treatment, beam
angular area, spectral/Doppler); a converter. **Velocity convention (radio|optical|
relativistic) + rest frequency are explicit and mandatory** where needed — missing → raise,
never default. **Write tests first:** unit-zoo round-trips, one hand-checked Planck (non-RJ)
value, missing-context-raises.

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
