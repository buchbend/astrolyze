# The web corpus explorer

## What it is

`astrolyze explore <corpus>` serves a published corpus as a browsable single-page web app — the GUI
sibling of the `astrolyze collection …` CLI group. It is a local FastAPI server whose endpoints are
**thin wrappers over the same read-only `Collection`/`Cube` public API**, with a Vue frontend on
top. The browser shows three views:

- the object-first **list** view (one row per source — surveys, species, store count, beam range;
  mirrors `Collection.list()`);
- the per-store **detail** view (one source expanded into its stores; mirrors
  `Collection.describe(object)`, shallow — no store opened);
- a **cube viewer** for one store, served lazily from the public `Cube` API.

It works on a local corpus and an `s3://` URL alike — it rides `Collection.open`'s fsspec layer, so
local and remote serve identically.

## The cube viewer

The viewer is modelled on the GO VIEW idea, served lazily from the public `Cube` API (each request
opens the store into a dask-backed cube and serves a **bounded JSON slice** — no full-cube load per
request, no server-side matplotlib). Its panels:

- an **integrated (moment-0) map** (`Cube.moment0`);
- a **channel map** with a velocity **slider** and **keyboard stepping** (click the viewer to focus
  it first; Arrow keys = ∓1 channel, Shift+Arrow or PageUp/PageDown = ∓5, Home/End jump to the
  first/last channel), the current velocity always shown;
- a **pixel spectrum** that updates when you click a position on either map;
- a **polygon region-averaged spectrum** — draw a region on the moment-0 map and its mean spectrum
  overlays the pixel spectrum;
- a **velocity-window moment** — brush a `[vmin, vmax]` window on the spectrum and the integrated map
  recomputes over *exactly* the channels in that window (replacing the moment-0 panel until cleared);
- **linked panels** — the clicked position is shared across panels (a crosshair on both maps, the
  spectrum's anchor), with a link/unlink control to share or split the two maps' pan/zoom.

## Why it works this way

- **The backend is a thin wrapper over the public API.** Every endpoint *dogfoods* the same public
  surface a Python user or the CLI uses — it never reaches into facade internals or the parquet
  reader. The list/detail endpoints serialize `Collection.list()` / `describe()` field-for-field;
  the viewer endpoints open a store with the public `Record.open()` and slice it with public `Cube`
  operations (`moment0`, `cube[index]`, `cube[:, y, x]`, spatial slicing + `nanmean` for a region,
  a `cube[i0:i1]` slab → `moment0` for a velocity window). The wire JSON mirrors the Python shape.
- **It stays lazy and bounded.** Opening a store reads its attrs/WCS, not its data plane, so the
  per-request open is cheap; each slice is bounded (a 2-D map, one channel, one spectrum, or one
  reduction), so a huge corpus cube never loads whole. The list/detail slice opens **no** store at
  all — `describe` is shallow, so the velocity-axis fields the catalog does not carry surface as
  `null`, the same honesty the CLI's shallow describe shows.
- **The bare install errors helpfully.** The web stack (fastapi + uvicorn) is the opt-in
  `astrolyze[web]` extra, not a hard dependency: a bare install never imports fastapi/uvicorn, the
  imports stay lazy, and `astrolyze explore` without the extra exits non-zero with a single
  actionable line (`pip install 'astrolyze[web]'`) rather than an opaque traceback — the same
  opt-in-extra contract as `astrolyze[s3]` and `astrolyze[coverage]`.

## Usage

```bash
pip install "astrolyze[web]"

astrolyze explore /data/ism_corpus                 # local corpus
astrolyze explore s3://my-bucket/ism_corpus        # object store — same call
astrolyze explore /data/ism_corpus --port 8001     # default 127.0.0.1:8000
```

It binds loopback by default (a local analysis tool); pass `--host 0.0.0.0` to expose it on the
network. Open the printed URL in a browser. On a bare install (no `[web]` extra) the command exits
with the install line above and nothing else.

### From a source checkout, build the frontend first

The published wheel ships the built Vue bundle, but a fresh `git clone` does not — there is no
`astrolyze/web/static/` until you build it. Without the build, `astrolyze explore` still serves the
JSON API, but the browser page is empty (the app degrades to the API alone, by design). Build the
bundle once:

```bash
cd astrolyze/web/frontend
npm install
npm run build      # emits the bundle into ../static, where the app mounts it at "/"
```

For frontend work, run the hot-reload dev server instead of rebuilding: keep `astrolyze explore
<corpus>` in one terminal and `npm run dev` in another — Vite proxies `/api` to the backend on
`127.0.0.1:8000`, so the SPA reloads on change against a real corpus (`npm run lint` runs eslint).
In production one process serves both the API and the built bundle; the proxy is dev-only.

## The JSON API

The endpoints that back the SPA are usable directly — handy for a quick `curl`, a notebook, or your
own client. All are read-only `GET`s under `/api` (plus one `POST` for a region polygon). A store is
addressed by an opaque `store_id` (the catalog `store_path`, base64-url-safe encoded so a slashed
path stays one URL segment); read it off the list/detail responses rather than building it.

| Endpoint | Returns | Public call |
| --- | --- | --- |
| `GET /api/collection` | corpus root + catalog version + object list | `Collection.list()` |
| `GET /api/objects/{object}` | per-store details for one object | `Collection.describe(object)` |
| `GET /api/stores/{store_id}` | one store's catalog row + resolved URI | a single `Record` |
| `GET /api/stores/{store_id}/cube` | axis metadata (shape, velocity axis, sky extent, units) | reads headers only |
| `GET /api/stores/{store_id}/moment0` *(opt. `?vmin=&vmax=`)* | moment-0 map, optionally over a velocity window | `Cube.moment0()` |
| `GET /api/stores/{store_id}/channel/{index}` | one channel slice + its velocity | `cube[index]` |
| `GET /api/stores/{store_id}/spectrum?x=&y=` | spectrum at a pixel | `cube[:, y, x]` |
| `POST /api/stores/{store_id}/region-spectrum` | region-averaged spectrum for a `{vertices}` polygon | spatial slice + `nanmean` |

An unknown object is `404`; an out-of-range channel/pixel or an empty velocity window is `422`,
naming the valid range. FastAPI also serves interactive API docs at `/docs`.

```bash
curl -s http://127.0.0.1:8000/api/collection | python -m json.tool
```

## See also

- [Browsing a published corpus](collections.md) — the `Collection` API the explorer wraps, and the
  `astrolyze collection …` CLI it is the GUI sibling of.
- [I/O and metadata](io-and-metadata.md) — remote corpora via fsspec URLs and the `astrolyze[s3]`
  extra the explorer reuses.
