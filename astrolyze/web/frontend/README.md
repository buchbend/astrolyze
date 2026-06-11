# astrolyze corpus explorer — frontend

A fresh Vue 3 + Vuex + Vue Router + D3 single-page app for the web corpus explorer (issue #66).
GO VIEW-inspired *conceptually only* — written from scratch, no yoda code or coupling.

Two views in this slice (phase-1):

- **List** (`/`) — object-first overview, mirroring `Collection.list()`: one row per source, its
  surveys / species / store count / beam range (with a small D3 span bar on a shared scale).
- **Detail** (`/objects/:object`) — per-store parameters, mirroring `Collection.describe(object)`:
  one card per store with beam, rest frequency, bunit, transition, provenance.

The four-panel **cube viewer** (integrated map / channel maps / spectra) is deferred to #67/#68.
A clean seam is left: a `store-viewer` route and a labelled placeholder card; no viewer logic.

## Develop

```bash
cd astrolyze/web/frontend
npm install
# In another terminal, serve a corpus so /api is live:
#   astrolyze explore /path/to/corpus
npm run dev        # hot-reload dev server, proxies /api -> http://127.0.0.1:8000
npm run lint       # eslint (vue3-recommended)
npm run build      # emits the production bundle into ../static (the wheel ships it)
```

## How it fits the backend

The store is a thin mirror of the public API the FastAPI backend serves:

| view   | store action      | endpoint                  | library method            |
| ------ | ----------------- | ------------------------- | ------------------------- |
| list   | `loadCollection`  | `GET /api/collection`     | `Collection.list()`       |
| detail | `loadObject`      | `GET /api/objects/{obj}`  | `Collection.describe(obj)`|

No business logic lives in the frontend — grouping and aggregation are the library's job; the app
fetches and renders what the API serialized. Production serves one process (FastAPI serves both the
API and this built bundle); `npm run dev` proxies to a running backend for hot-reload.
