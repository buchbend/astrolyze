# I/O and metadata

*Embodies [ADR-0006](../adr/0006-io-provenance-header-authoritative-filename-synced.md)
(header authoritative, lazy load) and
[ADR-0003](../adr/0003-unit-layer-pure-astropy-explicit-conventions.md) (no silent physics).*

## What it is

astrolyze treats a dataset's physical context — rest frequency, velocity convention, beam,
brightness unit, species, distance, calibration error — as a typed schema, `Metadata`. The
FITS header is the authoritative store of that context; `Metadata` is the parsed, in-memory
view of it.

There are now **two projections** of that one dataclass:

- `Metadata.to_header()` / `Metadata.from_header()` — the FITS header projection.
- `Metadata.to_attrs()` / `Metadata.from_attrs()` — a JSON-able dict projection, for
  backends that have no FITS header (e.g. a Zarr group's `.attrs`).

The dataclass stays the in-memory truth; the header and the `attrs` dict are two views of it.

`load()` returns a `LoadedData`: the array, the WCS, the parsed `Metadata`, and — so a
non-FITS backend can reconstruct the exact WCS — the **verbatim FITS-WCS header string**
(`header_string`). A live `astropy` `Header` is *optional* on `LoadedData`; only the FITS
loader fills it.

## Why it works this way

- **The header is authoritative, loading is lazy (ADR-0006).** A real archival header that is
  missing mandatory context (rest frequency or velocity convention) still opens. It is merely
  flagged incomplete; nothing is guessed and no file is refused. The gap only bites when an
  operation actually needs the absent context.
- **Two projections, one truth.** Keeping `Metadata` as the single source and treating the
  header and `attrs` as projections means a future non-FITS backend carries the *same*
  physical context verbatim without a second schema to drift out of sync. The `attrs`
  projection is JSON-serializable (`Quantity` → `{"value", "unit"}`, beam → `bmaj/bmin/bpa`,
  enums and units → strings), and absent fields are simply omitted — a partial schema
  round-trips while staying flagged incomplete (no invented context, ADR-0003).
- **`header_string` is the WCS round-trip vehicle.** A non-FITS backend has no live header,
  so the exact WCS travels as the verbatim header string; any backend rebuilds it with
  `WCS(fits.Header.fromstring(header_string))`. astropy still owns the reconstruction —
  astrolyze adds the schema, not a FITS reader.

## Usage

Load a cube and inspect completeness (nothing raises on an incomplete header):

```python
from astrolyze.io import load

loaded = load("ngc0628_co21.fits")
loaded.metadata.is_complete        # True if rest frequency + convention are present
loaded.metadata.missing            # e.g. [] or ["rest_frequency"]
loaded.header_string is not None   # the verbatim FITS-WCS header string is always present
```

Round-trip the metadata through the JSON-able `attrs` projection (#22):

```python
import json
from astrolyze.io import Metadata

attrs = loaded.metadata.to_attrs()   # a dict of JSON primitives
json.dumps(attrs)                    # succeeds — no Quantity/Beam objects left
Metadata.from_attrs(attrs) == loaded.metadata   # the inverse of to_attrs
```

Reconstruct the exact WCS from the verbatim header string (the path a non-FITS backend uses):

```python
from astropy.io import fits
from astropy.wcs import WCS

wcs = WCS(fits.Header.fromstring(loaded.header_string))
```

Writing only ever creates a derived file, named from the header, and never touches the raw
input (ADR-0006):

```python
from astrolyze.io import save

out_path = save(loaded.data, loaded.metadata, "derived/")
```

## Remote stores: fsspec URLs and S3 (#63)

The Zarr backend takes an **fsspec URL** anywhere it takes a path, so the same call shape serves
a local directory, an explicit `file://` URL, an in-memory `memory://` store, or an object store
on `s3://` — with no code change. A remote store opens through an fsspec mapper and reads stay
**lazy** (dask pulls only the chunks you touch, never the whole cube):

```python
from astrolyze.core import Cube
from astrolyze.collection import Collection

# Identical call shape — on disk or on S3:
cube = Cube.from_zarr("/data/ism_corpus/l1/heracles.zarr")     # local path
cube = Cube.from_zarr("file:///data/ism_corpus/l1/heracles.zarr")  # explicit file:// URL
cube = Cube.from_zarr("s3://my-bucket/ism_corpus/l1/heracles.zarr")  # object store

coll = Collection.open("s3://my-bucket/ism_corpus")            # a whole corpus on S3
```

`s3://` access needs the fsspec S3 driver, which is an **opt-in extra** (a bare install stays
fully functional on local / `file://` / `memory://` corpora):

```
pip install "astrolyze[s3]"
```

**Credentials and storage options are fsspec's, untouched.** A public/anonymous bucket and a
credentialed one differ only in how you configure fsspec, e.g. via standard AWS environment
variables (`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_PROFILE`) or, for a public bucket,
fsspec's anonymous mode (`{"anon": True}` in your fsspec/`s3fs` configuration). astrolyze never
inspects these — it hands them straight to fsspec.

**Caching is fsspec passthrough — there is no astrolyze cache layer.** Opt in to read caching with
fsspec's chained URLs, e.g. wrap a remote root in `simplecache::`:

```python
# Cache chunks locally on first read; subsequent reads hit the local cache:
cube = Cube.from_zarr("simplecache::s3://my-bucket/ism_corpus/l1/heracles.zarr")
```

This works because fsspec resolves the chained URL; astrolyze adds nothing. A missing catalog at a
URL raises a clear `FileNotFoundError` (the corpus root has no `catalog.parquet`), not an obscure
deferred read error.

## See also

- [No silent physics](no-silent-physics.md) — what happens when you act on incomplete metadata.
- [Coordinates and validity](coordinates-and-validity.md) — the WCS surfaced as physical arrays.
