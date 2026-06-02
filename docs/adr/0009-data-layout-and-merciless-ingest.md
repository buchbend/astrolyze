# Data layout convention + merciless ingest gate (DB-backed manifest)

The AI Layer (ADR-0008) requires a cross-project data layout that is easy to find and
understand, and a guarantee that the crucial metadata is *always present* — the original
reason astrolyze had a filename convention. Real-world cubes are missing required header
entries far more often than one expects, so the guarantee must be enforced, not hoped for.

## Data layout — fixed skeleton + manifest (Q12 = c)

Every experiment uses the same skeleton:

```
<experiment>/
  data/
    raw/         # immutable inputs — SACRED, never written to or renamed
    interim/     # derived intermediates (harmonised cubes, zarr stores)
    processed/   # final analysis-ready products
  outputs/
    figures/     # house-style plots
    tables/      # ECSV/CSV results
  logs/          # the experiment log(s)
  config.toml    # dynaconf
```

- **`raw/` is sacred:** read-only, never modified in place, upstream filenames preserved
  (renaming breaks provenance back to the source). All derivation flows raw → interim →
  processed.
- **Filename projection (ADR-0006) applies to derived files only** (`interim/`,
  `processed/`) so astrolyze-written products are self-describing and browsable; raw files
  keep their upstream names, mapped via the manifest.

## Merciless ingest gate

A dedicated **ingest** step scans `raw/`, validates each file's header against the required
Metadata Schema (rest frequency, velocity convention, beam, distance, species, calibration
error, object — ADR-0003/0006), and is **merciless**: missing required entries are flagged
and the user is *required* to supply them before the dataset is accepted/registered. This
is the modern realization of the old filename convention's goal — crucial info guaranteed
present — enforced at the right moment.

**Two strictness levels, reconciled with ADR-0006:**
- **Ingest = merciless** (the gate; nothing registers incomplete).
- **Load = lazy/(ii)** (mid-analysis reads open incomplete, raise on ops needing missing
  context). Ingest makes lazy-load safe in practice: registered data is already complete.

## Manifest = DB-backed registry

The per-experiment manifest is a **database-backed registry** of datasets: identity,
location, full provenance (beam/channel/rest-freq — the dataset heterogeneity that makes
the registry scientifically valuable), and source DOI for reproducibility/citation.
Generated, not hand-maintained.

**Seam (not built now):** a DB-backed manifest lets astrolyze later serve as the **backbone
of a graphical UI frontend** — an explicit design goal, deliberately unbuilt (YAGNI). Shape
the manifest/DB so a UI *could* read it; build no UI yet.

## Consequences

- One predictable layout agents/skills can assume; raw provenance never lost.
- Ingest is where data-quality is enforced — expect to spend real effort filling missing
  headers (by design).
- DB choice (sqlite to start, per astrolyze v1 precedent) and schema are build-time details;
  keep it swappable and UI-readable.
