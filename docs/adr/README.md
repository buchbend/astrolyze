# Architecture Decision Records

The design decisions that shape astrolyze, kept with the code so they travel with the
repository. Each ADR states a decision, the options considered, and the consequences.

These are **mirrored from the private research repository**, where the full design history
(PRD, build plan, and the strategy/topology ADRs) is canonical. The subset here is the set
that governs the **public toolkit's own architecture**; it is what the code in this repo
refers to (e.g. docstrings citing `ADR-0003`, `ADR-0006 ii`).

## The records

| ADR | Decision | In code today |
| --- | --- | --- |
| [0002](0002-core-toolkit-new-repo-on-ecosystem.md) | New public package on the modern Astropy ecosystem; astrolyze2 as reference spec | ✅ the repo itself |
| [0003](0003-unit-layer-pure-astropy-explicit-conventions.md) | Unit layer: pure-astropy substrate, object-carried context, **explicit** convention & rest frequency | ✅ `astrolyze/units/` |
| [0004](0004-object-model-cube-map-spectrum-wrappers.md) | Object model: thin `Cube`/`Map`/`Spectrum` that **wrap** (not subclass) the ecosystem | ✅ `astrolyze/core/` |
| [0005](0005-display-layer-engine-plus-sugar.md) | Display: free-function engine + thin method sugar; style applied **locally**, never global | ✅ `astrolyze/viz/` |
| [0006](0006-io-provenance-header-authoritative-filename-synced.md) | I/O: FITS header authoritative, filename a synced browsable projection, lazy load | ✅ `astrolyze/io/` (provenance seam reserved) |
| [0009](0009-data-layout-and-merciless-ingest.md) | Data layout: fixed experiment folder skeleton + merciless ingest gate + DB-backed manifest | ⬜ specced, not built |
| [0010](0010-experiment-log-always-runlog-optional-narrative.md) | Experiment log: always auto-capture *what ran*; narrative offered, never enforced | ⬜ specced, not built |
| [0011](0011-everything-through-astrolyze-formalize-before-use.md) | All data work goes through astrolyze; formalize-before-use; astrolyze is AI-independent | ◐ CLI path exists; discipline |
| [0012](0012-agent-integration-guide-cli-no-mcp.md) | Agent integration: usage guide + CLI/API, **no MCP** | ✅ `astrolyze/cli.py` |
| [0013](0013-traceable-correctable-results.md) | Traceable, correctable results — legibility over enforcement | ◐ partial; run-log pending (0010) |

## Why the numbering has gaps

ADRs **0001** (core/scaffolding hard wall), **0007** (two-repo topology), **0008** (agent-native
layer rationale), and a second **0005** (velocity-axis resolution ladder) are deliberately
**not** mirrored here — they describe the private research project and its strategy rather
than the public toolkit. The gaps are expected, and a few records below cite `ADR-0001` /
`ADR-0008`: those references are **external by design** and resolve in the private repo.

## Relationship to the originals

These copies are kept verbatim except for three small public-safe edits in 0002, 0006, and
0009 that removed a private document reference and internal strategy framing while preserving
all technical content. The private repo remains the canonical source; when an ADR changes
there, re-mirror it here.
