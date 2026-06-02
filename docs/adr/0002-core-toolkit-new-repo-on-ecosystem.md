# Core Toolkit: a new public package built on the modern ecosystem, with astrolyze2 as reference spec

The toolkit idea has been built three times (astrolyze v1, astrolyze3, astrolyze2)
without one iteration carrying across the finish line. v1/v3 hold the deep domain code
(CLASS/GILDAS I/O, LTE, the curated `cfg/` line/molecule/galaxy/calibration tables) on a
dated architecture; **astrolyze2** (2023–24) is the modern, clean-slate rewrite that
expresses Christof's current taste — unit-aware `FitsMap` carrying data+WCS+beam+distance,
new-object-returning map algebra, a hub-and-spoke unit converter (everything ↔ Jy/sr via
astropy equivalencies), centralized WCS+beam plotting, a `Stack` workflow, dynaconf config
— but it stalled at "works in my notebook": 2D-only (no Cube/PPV class), no test suite,
and it hand-rolls convolution/reproject/beam math the ecosystem now provides.

**Decision:** build the **Core Toolkit as a NEW, public open-source package from day one**,
re-based on the modern Astropy-affiliated ecosystem (`spectral-cube`, `astropy`,
`radio_beam`, `reproject`, `regions`). **astrolyze2 is the reference specification** (lift
its unit-hub, plotting, Stack, provenance and SED design deliberately — do not extend its
repo in place). **astrolyze v1/v3 are the domain-knowledge harvest source** (port the
curated tables + LTE physics + CLASS/GILDAS glue, re-derived against astropy units, with
tests). The from-scratch `Map`/convolution code is NOT carried over — depend on the
ecosystem instead.

## Considered options

- **(a) Continue astrolyze2 in place** — fastest, but inherits the WIP debt and the
  "extend the notebook" gravity, and keeps the from-scratch-Map path.
- **(c) Fork astrolyze2 → rename → refactor onto ecosystem** — keeps git history; rejected
  in favour of a clean slate that lets the ecosystem (esp. spectral-cube's Cube class) and
  the test-first discipline shape the architecture from the start.
- **(b) New repo, astrolyze2 as reference spec — CHOSEN.**

## Consequences

- Costs deliberately re-typing astrolyze2's good parts and porting v1/v3 tables; buys a
  clean architecture on spectral-cube, tests from line one, and no inherited footguns
  (the in-place `__getattr__` numpy mutation, `raise SystemExit`/`print` in library code,
  the latent unit bugs).
- Public from day one (BSD-class license, PyPI, citable) — supports astrolyze's open
  data-infrastructure-for-the-field positioning. The cost is day-one polish pressure on
  naming/licence/hygiene.
- Subject to the hard wall + promotion discipline that separates the core toolkit from
  research scaffolding (decided in the private research repo; see note in the ADR index).

## Name and license (resolved)

- **Name: `astrolyze`.** After exploring short fresh names (incl. sky/heaven-god themes),
  the decision is to reclaim the existing name. This new repo is the **canonical
  astrolyze** — the public, tested, ecosystem-based package that supersedes v1/astrolyze3
  (Harvest Source) and astrolyze2 (Reference Spec). Christof already owns the PyPI name and
  a decade of identity rides on it; the "fresh start" is in the *codebase*, not the brand.
- **License: BSD-3-Clause** — astrolyze has always been BSD, and it is the astropy/numpy
  ecosystem norm (permissive, max adoption, frictionless to depend on).
