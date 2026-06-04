# Working with astro data in astrolyze (guide for humans and agents)

This is the house guide for using astrolyze. A coding agent uses astrolyze the **same way a
human does** — by writing astrolyze Python/CLI — producing ordinary, reviewable, reproducible
scripts. There is no AI-only interface (no MCP); the CLI and Python API *are* the interface.

## Principles

1. **Everything goes through astrolyze.** Don't hand-roll matplotlib, unit math, or FITS
   parsing in an analysis. If astrolyze can't do something yet, **extend astrolyze first**
   (with tests, by its conventions), then use it.
2. **Be explicit about physics.** Always supply the velocity convention and rest frequency
   where needed; astrolyze will raise if they're missing — that's intentional, not an error to
   work around.
3. **Use the house display.** Plot via `.plot()` / `astrolyze.viz`; don't restyle globally.
4. **Results are traceable.** Surface how a result was produced (operations + reasoning) so a
   reader can follow and correct it; cite the real artifact behind any claim.
5. **Stay thin.** Prefer delegating to `spectral-cube`/`astropy` over reimplementing.

## The core objects

`Cube` (PPV), `Map` (2D / moment), and `Spectrum` (1D) are thin wrappers that **compose**
spectral-cube / astropy / specutils and carry the physical context (beam + rest frequency +
velocity convention) parsed from the FITS header. Build one from a file through the `io` seam:

```python
from astrolyze.io import load
from astrolyze.core import Cube

cube = Cube.from_loaded(load("ngc0628_co21.fits"))   # <Cube NGC0628 (98, 1600, 1600) [K]>
```

Type transitions carry the context for you:

| operation                | result     |
| ------------------------ | ---------- |
| `cube.moment0()`         | `Map`      |
| `cube.moment(order, axis)` | `Map`    |
| `cube[k]` (channel)      | `Map`      |
| `cube[:, y, x]`          | `Spectrum` |
| `cube[:, y0:y1, x0:x1]`  | `Cube`     |

## The tracer recipe (load → moment0 → convert → plot)

```python
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))
integrated = cube.moment0().to("K km/s")   # beam / rest freq / convention come from the object
fig, ax = integrated.plot()                # cividis + WCS axes + beam ellipse + unit colorbar
fig.savefig("ngc0628_co21_moment0.png")
```

`.to(unit)` supplies the object's beam, rest frequency, and velocity convention to the unit
layer so you don't repeat them. It **never** defaults the genuinely-ambiguous
Rayleigh-Jeans-vs-Planck scale — a brightness-temperature conversion needs it stated:

```python
channel = cube[48]                                   # a Jy/beam or K channel Map
channel.to("K", temperature_scale="rayleigh_jeans")  # or "planck"; omitting it raises
```

If the header is missing a rest frequency or velocity convention, the file still **loads**
(so you can inspect archival data), but any operation that needs the missing context raises
`astrolyze.units.MissingContextError` rather than guessing.

## The CLI (same path from the shell)

The CLI exposes the identical spine; an agent reaches for it exactly as a human would.

```bash
astrolyze init  study/                            # scaffold an experiment (skeleton + config)
astrolyze ingest study/                           # validate data/raw/ + register what is complete
astrolyze manifest list study/ --object NGC0628   # query the registry (no grepping filenames)
astrolyze info  ngc0628_co21.fits                 # metadata schema + completeness (read-only)
astrolyze moment0 ngc0628_co21.fits -u "K km/s" -o ngc0628_mom0.png
astrolyze moment0 ngc0628_co21.fits --temperature-scale planck   # for a brightness conversion
astrolyze --help
```

`info` prints a table of the parsed schema (object, telescope, species, rest frequency,
velocity convention, beam, unit, distance, calibration error) and whether the file is
*complete*. `moment0` runs load → moment0 → `.to(unit)` → plot and writes a house-style PNG
(default `<input>_moment0.png`). A conversion that needs absent context exits non-zero with a
clear message — the CLI surfaces the library's refusal, it doesn't paper over it. `ingest` and
`manifest list` drive the merciless gate and registry (below).

## Organising an analysis: the experiment skeleton

An *experiment* is the fixed home of one analysis. `astrolyze init` scaffolds it so every
study has the same predictable, portable shape — you never hand-build folders, and an agent
assumes the identical layout (ADR-0009):

```bash
astrolyze init study/        # create the skeleton + a default config.toml (idempotent)
```

```
study/
  data/
    raw/         # immutable inputs — SACRED: astrolyze never writes to or renames these
    interim/     # derived intermediates (harmonised cubes, …)
    processed/   # analysis-ready products
  outputs/
    figures/     # house-style plots
    tables/      # ECSV/CSV results
  logs/          # experiment run log(s)
  config.toml    # dynaconf settings
```

- **`raw/` is sacred.** Derivation flows raw → interim → processed; raw files keep their
  upstream names so provenance back to the source is never broken. The header-derived filename
  projection (ADR-0006) applies to derived files only.
- **`init` is idempotent.** Re-running it adds nothing and never overwrites a `config.toml`
  you have edited.
- **`config.toml`** holds the configurable knobs, read through dynaconf — currently the
  manifest database URL and an optional override of the mandatory-context fields.

The same skeleton is a value object from Python:

```python
from astrolyze.experiment import Experiment, role_of

exp = Experiment.init("study")            # or Experiment("study") to resolve paths only
exp.raw, exp.interim, exp.figures         # resolved skeleton paths
role_of(exp, exp.raw / "archival.fits")   # -> Role.RAW  (interim / processed / output / None)
```

The **dataset manifest**, the merciless **`ingest`** gate, and the always-on **run log** (all
below) build on this skeleton. Saving derived products with header-derived names follows in a
later slice (ADR-0006).

### The dataset manifest

The manifest is the **generated registry of an experiment's datasets** — one row per dataset,
backed by a database (sqlite by default, one file inside the experiment). Each row holds:

- **identity** — the source path (relative to the experiment) and an optional source **DOI**;
- **full provenance**, copied from the FITS header via the `io` `Metadata` schema — object,
  telescope, species, rest frequency, velocity convention, beam, `bunit`, distance,
  calibration error;
- the **schema version** and an **ingested-at** timestamp.

**Why it exists.** `raw/` is sacred, so archival files keep their upstream names — the manifest
is what maps those names to identity and physical provenance, so you can *see what is in an
experiment at a glance* and query it (by object, species, …) instead of grepping filenames. A
DOI per dataset keeps inputs citable and reproducible. It is **generated and kept in sync by
ingest, never hand-edited**, so it never drifts from what is actually on disk. The backend is
deliberately **swappable** (sqlite to start) and the registry is shaped to be **UI-readable** —
a future graphical frontend can render it without rework (a seam, not built; ADR-0009).

**Opening it.** Resolve it from the experiment config (the config `db_url` is relative; this
anchors it to the experiment root, never the current working directory):

```python
from astrolyze.experiment import Experiment, Manifest

exp = Experiment.init("study")
manifest = Manifest.for_experiment(exp)        # uses config.toml [manifest] db_url
# …or point it at any backend URL directly (backend-swappable):
# manifest = Manifest("sqlite:////abs/path/registry.db")
```

**Registering** is normally driven by `ingest` (one call per accepted `raw/` file); done by
hand it looks like this. It is **idempotent on the source path** — re-registering the same path
updates the row in place (the id stays stable), never duplicating, so re-ingesting after fixing
a header is a normal step:

```python
from astrolyze.io import load

loaded = load(exp.raw / "ngc0628_co21.fits")
record = manifest.register(loaded.metadata, "data/raw/ngc0628_co21.fits", doi="10.3847/…")
```

**Reading it back** — `get` / `query` / `all` return `DatasetRecord` value objects, each
carrying a *reconstructed* `Metadata` (provenance is stored as primitives — Hz, degrees, unit
strings — and the astropy/`radio_beam` objects are rebuilt on read, so you round-trip the
value):

```python
manifest.get(record.id)                  # one dataset, or None
manifest.query(object="NGC0628")         # filter by object / species / telescope / …
manifest.query(species="CO21")
manifest.all()                           # every registered dataset

rec = manifest.query(object="NGC0628")[0]
rec.source_path, rec.doi, rec.ingested_at
rec.metadata.beam, rec.metadata.rest_frequency, rec.metadata.velocity_convention
```

The manifest itself is *generic storage*: it records whatever provenance is present (optional
fields as nulls). Guaranteeing the mandatory context is filled in is `ingest`'s job — the
merciless gate — not the manifest's. Don't hand-edit the database; let `ingest` generate it.

### Merciless ingest: the gate into the manifest

`ingest` is the **strict counterpart to the lazy `load`**. It walks `data/raw/`, parses each
file's header through the *same* `io` `Metadata` schema, and partitions every dataset:

- **accepted** — carries the mandatory physical context → **registered** in the manifest;
- **rejected** — missing at least one mandatory field → **not** registered, and the report
  **names the exact missing fields** so you know what to fix.

That is what *merciless* means: nothing incomplete is ever registered, so a dataset in the
manifest is always safe to compute on. `raw/` stays sacred — ingest never renames or writes
anything under it (and reads **only headers**, never the multi-GB cube body). Non-FITS files
(READMEs, sidecars) are ignored; an unreadable/corrupt file is reported, not allowed to abort
the pass.

```bash
astrolyze ingest study/        # rich report: accepted (+object/species), rejected (+reasons)
```

```
accepted — registered in the manifest      rejected — not registered
  data/raw/ngc0628_co21.fits  NGC0628 CO21    data/raw/archival_map.fits
                                                missing mandatory context:
                                                rest_frequency, velocity_convention
ingested study/: 1 accepted, 1 rejected
```

Fixing a header and re-running is the normal, iterative path — re-ingest **updates the same
manifest row** (stable id), it never duplicates, and the now-complete file is accepted.

From Python:

```python
from astrolyze.experiment import Experiment, ingest

report = ingest(Experiment("study"))     # -> IngestReport
[a.source_path for a in report.accepted]  # registered datasets (each .record is a DatasetRecord)
[(r.source_path, r.missing) for r in report.rejected]  # refused, with the named gaps
```

The mandatory-context set defaults to the schema baseline (`rest_frequency`,
`velocity_convention`) and is **overridable** — explicitly with `ingest(exp, required=(…))`, or
per-experiment via `[ingest] required_context` in `config.toml` (e.g. add `"beam"` for a study
that requires a resolved beam). An override naming a non-schema field is a clear error.

### The always-on run log

"What ran" is captured **automatically** — you never write a log by hand, and an agent never
has to *remember* to (ADR-0010). Open a run, and from then on every astrolyze operation appends
a record:

```python
from astrolyze.experiment import Experiment, RunLog
from astrolyze.core import Cube
from astrolyze.io import load

exp = Experiment.init("study")
with RunLog.open(exp) as run:                       # per-run append-only JSONL in logs/
    cube = Cube.from_loaded(load(exp.raw / "ngc0628_co21.fits"))
    mom0 = cube.moment0().to("K km/s")
    fig, ax = mom0.plot()
    print(run.path)                                 # study/logs/run-<id>.jsonl
```

Each `load` / `save` / `moment*` / `.to` / `.plot` writes one JSON line carrying the operation,
its params, inputs/outputs (paths/ids), the software versions (astrolyze + key deps), the data
version (the `io` schema version + any manifest ids), and a UTC timestamp:

```json
{"schema": 1, "run_id": "20260604T093737-dc781abb", "timestamp": "2026-06-04T09:37:37.057+00:00",
 "op": "moment", "params": {"order": 0, "axis": 0}, "inputs": [], "outputs": [],
 "software": {"astrolyze": "0.1.0.dev0", "astropy": "7.2.0", …}, "data": {"schema_version": 1, "manifest_ids": []}}
```

- **Per-run file/id.** Every `RunLog.open` gets its own id and `run-<id>.jsonl`, so separate
  sessions never tangle. The file is **append-only** (records are never overwritten) — read it
  back with `run.entries()` (a list of parsed records) or just parse the JSONL.
- **A no-op outside an experiment.** With no run open, the emit seam does nothing — the library
  (and the tracer spine) runs unchanged for quick one-off use. There is no "log" command and no
  flag to remember: opening a run is the only switch.
- **Narrative is offered, never enforced** (ADR-0010): the run log records *what* ran; the
  human "why / what it means" account is voluntary, and any narrative you write should point at
  real run-log/figure artifacts (the natural place to compose it is the future UI — not built).

```python
from astrolyze.experiment import RunLog, emit

with RunLog.open(exp) as run:
    ...                                  # operations emit automatically
    emit("note", params={"why": "first look at the inner ring"})   # add your own record too
records = run.entries()                  # [{op: "load", …}, {op: "moment", …}, …]
```

## Try it on real data

A 128×128×50 cutout of the PHANGS-ALMA NGC 628 CO(2-1) cube ships in
`tests/data/ngc0628_co21_cutout.fits.gz` (see `tests/data/PROVENANCE.md`):

```bash
astrolyze info tests/data/ngc0628_co21_cutout.fits.gz
python examples/tracer_ngc628.py tests/data/ngc0628_co21_cutout.fits.gz /tmp/mom0.png
```

Point `$ASTROLYZE_TRACER_CUBE` at a full cube to run the example on real survey data without
an argument. For a guided, interactive walk-through — slicing out spectra and channel maps,
converting units, composing multi-panel figures — open
[`examples/tutorial.ipynb`](examples/tutorial.ipynb) (`pip install -e ".[notebook]"`).
