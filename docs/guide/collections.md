# Browsing a published corpus

## What it is

A *collection* is a published corpus opened **read only** for analysis: a directory of
astrolyze-flavoured Zarr stores indexed by one `catalog.parquet` at its root. `Collection` is the
consumer-facing handle on it — open it, list what is in it, filter it, ask which cubes cover a sky
position, and open any record into a lazy `Cube`. It never writes to the corpus; lineage flows
*out* through each opened cube's origin `Metadata`, never back in.

- `Collection.open(path)` — open a corpus from a local directory, a `file://` URL, or an `s3://`
  URL with **no change in call shape** (path handling routes through fsspec from day one).
- `list()` — the **object-first** overview: one `ObjectSummary` per source (its surveys, species,
  store count, beam range).
- `describe(object, deep=…)` — drill one source down into per-store physical detail
  (`StoreDetail`), **catalog-row-first** (no store read unless you ask for `deep`).
- `query(**filters)` — slice the corpus along any catalog field; returns a **composable**
  sub-`Collection`.
- `covering(SkyCoord)` — which cubes cover a sky position, across surveys; returns a composable
  sub-`Collection`.
- `records` — the flat per-store view (`Record` each); `Record.open()` opens the store into a
  lazy, dask-backed `Cube` stamped with its corpus origin.

The catalog format is a versioned data-format contract specified in a separate producer repository;
astrolyze only ever *reads* it. The two codebases never import each other — the versioned spec is
the whole interface.

## Why it works this way

- **The catalog is an index; the WCS is the authority.** `covering()` is deliberately two-layered.
  The catalog's footprint columns — the always-present center+radius circle (`ra_deg` / `dec_deg` /
  `radius_deg`), and an exact-coverage MOC when one is present — only *prefilter* candidates
  cheaply, without opening a store. The final containment decision for every surviving candidate is
  made from that store's **own** celestial WCS at open time (the same half-open pixel-cell rule
  `Cube.cutout` uses). A stale or coarse index can never yield a wrong containment answer, and an
  unverifiable store is *not* asserted to cover the position rather than trusting the prefilter —
  refuse to guess. See [No silent physics](no-silent-physics.md).
- **`describe` is cheap by default.** Everything `describe(object)` returns is read straight off the
  typed catalog row, so describing a source on a remote corpus costs no store read. The three
  velocity-axis fields the catalog does not carry (native channel width, velocity coverage) are
  filled **only on request** (`deep=True`, which opens each store's 1-D spectral axis — still lazy,
  not the cube body); on a shallow describe they stay `None` — unknown, never guessed.
- **An unknown filter key raises.** `query(surveys=…)` (a typo for `survey=`) does not silently
  return an empty result; it raises `ValueError` naming the offending key and the valid axes. A
  silent no-op would hide the typo as a wrongly-empty result.
- **A catalog-less directory just works.** When the root carries no `catalog.parquet`, `open()`
  transparently falls back to the **scan-builder**: every astrolyze Zarr store is self-describing
  (it carries the `Metadata` schema and a verbatim WCS in its group attrs), so a catalog is
  *reconstructed* in memory by reading those attrs and computing each store's footprint from its own
  WCS. A bare directory of stores is browsable with zero extra tooling. The scan is read-only — it
  builds an in-memory catalog and never writes one into the corpus (persist explicitly with
  `astrolyze.collection.scan.write_catalog`). A store missing the celestial WCS a footprint needs is
  skipped with a `ScanWarning`, never catalogued with invented geometry.
- **Origin provenance bridges to a writable experiment.** Opening a record stamps the cube's
  `Metadata` with its origin — the source store URI and the catalog version — so a product later
  saved into an experiment traces back to the exact corpus snapshot it came from. This is the only
  thing the collection adds over a bare `Cube.from_zarr`. See [I/O and metadata](io-and-metadata.md).

## Usage

Open a corpus and see what is in it, object-first (the same call shape on disk or on S3):

```python
from astrolyze.collection import Collection

coll = Collection.open("/data/ism_corpus")            # local directory
coll = Collection.open("s3://my-bucket/ism_corpus")   # object store — same call

for summary in coll.list():
    print(summary.object, summary.surveys, summary.species, summary.n_stores)
    print("  beam range (arcsec):", summary.beam_range_arcsec)   # (coarsest, finest) major axis
```

Drill one source down to per-store detail — shallow (catalog only) or `deep` (opens each store's
spectral axis):

```python
for detail in coll.describe("NGC3521"):              # cheap: no store opened
    print(detail.survey, detail.species, detail.transition, detail.bunit)
    print("  beam:", detail.beam_major_arcsec, detail.beam_minor_arcsec, detail.beam_pa_deg)

for detail in coll.describe("NGC3521", deep=True):   # opens each store's 1-D axis (still lazy)
    print(detail.channel_width_kms, detail.velocity_min_kms, detail.velocity_max_kms)
```

Filter the corpus along its own fields; the result is a composable sub-`Collection`, so every facade
method chains on it:

```python
co = coll.query(species="CO", transition="2-1")      # conjunctive: must match all
co.list()                                            # works — query() returns a Collection
co.query(survey="HERACLES").describe("NGC3521")      # chain further

coll.query(surveys="HERACLES")                       # raises ValueError — unknown key (it's `survey`)
```

Ask which cubes cover a sky position, across surveys, then open and analyse them:

```python
import astropy.units as u
from astropy.coordinates import SkyCoord

pos = SkyCoord("11h05m48.6s", "-00d02m09s")
hits = coll.covering(pos)                            # a sub-Collection of covering cubes
for record in hits.records:
    cube = record.open()                             # lazy, dask-backed; carries origin provenance
    print(record.object, record.survey, cube.metadata.origin_store_uri)

coll.covering(pos).query(species="CO")               # composes like any sub-Collection
```

A position covered by nothing returns an **empty** collection — an honest "no coverage", never an
error.

## The CLI

The `astrolyze collection …` group mirrors these facade methods from the shell (read-only, like
`astrolyze manifest`). A local directory and an `s3://` URL share each command — fsspec resolves
both:

```bash
astrolyze collection list  /data/ism_corpus                       # object-first overview
astrolyze collection describe NGC3521 /data/ism_corpus            # per-store detail (catalog only)
astrolyze collection describe NGC3521 /data/ism_corpus --deep     # + channel width + velocity range
astrolyze collection query /data/ism_corpus --species CO --transition 2-1
astrolyze collection covering /data/ism_corpus --ra 166.45 --dec -0.036
astrolyze collection covering /data/ism_corpus --coord "11h05m48.6s -00d02m09s"
```

A catalog-less directory is scanned transparently, so a bare directory of Zarr stores lists too. An
unknown `catalog_schema_version` exits non-zero with a clear message — the CLI surfaces the
library's refusal rather than papering over it.

## Optional extras

A bare install browses local (and `file://` / `memory://`) corpora end to end. Two opt-in extras
light up extra paths, and their absence never breaks the bare path:

- `pip install "astrolyze[s3]"` — the fsspec S3 driver, so `Collection.open("s3://…")` works.
  Credentials and storage options are fsspec's, untouched (see [I/O and
  metadata](io-and-metadata.md)).
- `pip install "astrolyze[coverage]"` — `mocpy`, so `covering()` can prefilter candidates by the
  exact-footprint MOC instead of only the coarse center+radius circle. Without it, the center+radius
  prefilter plus the store-WCS authority decide containment regardless — the bare install runs
  `covering()` end to end.

## See also

- [No silent physics](no-silent-physics.md) — the refuse-rather-than-guess principle behind the
  catalog-as-index / WCS-as-authority split and the unknown-filter-key raise.
- [I/O and metadata](io-and-metadata.md) — origin provenance on `Metadata`, and remote stores via
  fsspec URLs.
- [Cutouts and stacking](stacking.md) — gathering a corpus into a `Stack` of postage stamps.
