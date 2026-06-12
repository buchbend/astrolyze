// The Vuex store: the single source of truth, mirroring the public Collection API endpoints.
//
// House convention is vuex (global CLAUDE.md). The store holds exactly what the backend endpoints
// return — the corpus overview + object-first list (GET /api/collection), the per-object store
// details (GET /api/objects/:object), and the cube viewer's slices (GET .../cube|moment0|channel|
// spectrum, issue #67) — with loading/error flags so the views stay declarative. No business logic
// lives here: grouping/aggregation are the library's job (Collection.list/.describe) and the slices
// are the library's lazy Cube reads (moment0/channel/spectrum); this just fetches and caches what
// the API serialized.
//
// The `viewer` slice is the cube viewer's shared state — the {store, channel, position} the three
// linked panels read so a click on either map updates the spectrum and both crosshairs. The #68
// interactions fill the once-reserved slots with logic: `region` is the drawn polygon + its
// region-averaged spectrum (GET region-spectrum), `velocityWindow` is the [v0, v1] brushed on the
// spectrum that recomputes the integrated map (GET moment0?vmin&vmax), and `linked` toggles whether
// the two map panels share one pan/zoom view (the link/unlink control).

import { createStore } from "vuex";

// All API calls share one tiny fetch helper so error handling is uniform (a non-2xx becomes a
// thrown Error carrying the backend's JSON `detail`, which the views surface verbatim).
async function getJson(url) {
  const response = await fetch(url);
  if (!response.ok) {
    let detail = `${response.status} ${response.statusText}`;
    try {
      const body = await response.json();
      if (body && body.detail) detail = body.detail;
    } catch {
      // Non-JSON error body — keep the status line.
    }
    throw new Error(detail);
  }
  return response.json();
}

export default createStore({
  state() {
    return {
      // GET /api/collection
      rootUri: null,
      catalogVersion: null,
      objects: [], // [{ object, surveys, species, n_stores, beam_min_arcsec, beam_max_arcsec }]
      // GET /api/objects/:object — keyed by object name so navigating back is instant (cached).
      detailsByObject: {},
      loading: false,
      error: null,
      // -- cube viewer (#67) --------------------------------------------------------------
      // The state the three linked panels share. `store` is the opened store_id; `axes` is the
      // GET .../cube response (shape, velocity axis, sky extent, units); `moment0` / `channel` are
      // the two heatmaps; `spectrum` is the pixel line. `channelIndex` drives the slider/keyboard;
      // `position` is the clicked {x, y} pixel shared as the crosshair on both maps and the spectrum
      // anchor. Loading/error are viewer-local so a slow slice does not blank the masthead.
      viewer: {
        store: null,
        axes: null,
        moment0: null,
        channel: null,
        channelIndex: 0,
        // Fetched channel slices keyed by channel index — stepping back to a seen channel is a
        // cache hit (no refetch). Reset per store so a switch never shows a stale cube's slices.
        channels: {},
        position: null, // { x, y } image pixel, or null until the first click
        spectrum: null,
        loading: false,
        error: null,
        // -- #68 interactions ----------------------------------------------------------
        // region: the drawn polygon (image-pixel [x, y] vertices) and its averaged spectrum.
        //   { vertices: [[x, y], …], spectrum: { velocity, value, n_pixels, … } | null }
        // velocityWindow: the [v0, v1] brushed on the spectrum (cube velocity unit), or null.
        // windowedMoment0: the moment-0 recomputed over that window (replaces the map panel), or null.
        // linked: whether the two map panels share one pan/zoom view (the link/unlink control).
        // mapView: the shared { k, x, y } zoom transform the linked maps render at.
        region: null,
        velocityWindow: null,
        windowedMoment0: null,
        linked: true,
        mapView: null,
      },
    };
  },
  getters: {
    // The list view sorts client-side too so a future column-sort is a one-line change; the API
    // already returns objects sorted by name, so this is a stable no-op until then.
    sortedObjects(state) {
      return [...state.objects].sort((a, b) =>
        String(a.object).localeCompare(String(b.object)),
      );
    },
    detailFor: (state) => (object) => state.detailsByObject[object] || null,
    // The current channel's velocity (km/s) for the slider/keyboard readout — read off the axes so
    // no extra request is needed. Falls back to the channel slice's own velocity, else null.
    currentVelocity(state) {
      const v = state.viewer;
      if (v.axes && Array.isArray(v.axes.velocity)) {
        return v.axes.velocity[v.channelIndex] ?? null;
      }
      return v.channel ? v.channel.velocity : null;
    },
    velocityUnit(state) {
      return state.viewer.axes ? state.viewer.axes.velocity_unit : null;
    },
  },
  mutations: {
    setLoading(state, value) {
      state.loading = value;
    },
    setError(state, message) {
      state.error = message;
    },
    setCollection(state, payload) {
      state.rootUri = payload.root_uri;
      state.catalogVersion = payload.catalog_version;
      state.objects = payload.objects;
    },
    setDetail(state, { object, stores }) {
      // Reassign the whole map so Vue reactivity tracks the new key.
      state.detailsByObject = { ...state.detailsByObject, [object]: stores };
    },
    // -- cube viewer (#67) ----------------------------------------------------------------
    setViewerLoading(state, value) {
      state.viewer.loading = value;
    },
    setViewerError(state, message) {
      state.viewer.error = message;
    },
    // Open a (new) store in the viewer: reset the per-store slices + selection so a switch never
    // shows the previous cube's panels. The #68 slots are reset to their reserved defaults too.
    resetViewer(state, storeId) {
      state.viewer.store = storeId;
      state.viewer.axes = null;
      state.viewer.moment0 = null;
      state.viewer.channel = null;
      state.viewer.channelIndex = 0;
      state.viewer.channels = {};
      state.viewer.position = null;
      state.viewer.spectrum = null;
      state.viewer.error = null;
      state.viewer.region = null;
      state.viewer.velocityWindow = null;
      state.viewer.windowedMoment0 = null;
      state.viewer.linked = true;
      state.viewer.mapView = null;
    },
    setViewerAxes(state, axes) {
      state.viewer.axes = axes;
    },
    setViewerMoment0(state, moment0) {
      state.viewer.moment0 = moment0;
    },
    setViewerChannel(state, channel) {
      state.viewer.channel = channel;
      if (channel && Number.isInteger(channel.index)) {
        state.viewer.channelIndex = channel.index;
        // Cache the slice so stepping back is a hit. Reassign the map so Vue tracks the new key
        // (same convention as detailsByObject).
        state.viewer.channels = {
          ...state.viewer.channels,
          [channel.index]: channel,
        };
      }
    },
    // The slider/keyboard set the index optimistically (snappy UI); the channel slice fetch follows.
    setChannelIndex(state, index) {
      state.viewer.channelIndex = index;
    },
    setViewerPosition(state, position) {
      state.viewer.position = position;
    },
    setViewerSpectrum(state, spectrum) {
      state.viewer.spectrum = spectrum;
    },
    // -- cube viewer interactions (#68) ---------------------------------------------------
    // The drawn polygon: keep the vertices immediately (so the overlay renders while the
    // region-averaged spectrum is in flight) and attach the spectrum when it lands.
    setViewerRegion(state, region) {
      state.viewer.region = region;
    },
    setViewerRegionSpectrum(state, spectrum) {
      if (state.viewer.region) {
        state.viewer.region = { ...state.viewer.region, spectrum };
      }
    },
    clearViewerRegion(state) {
      state.viewer.region = null;
    },
    // The velocity window brushed on the spectrum (or null to clear it and the windowed map).
    setViewerVelocityWindow(state, window) {
      state.viewer.velocityWindow = window;
      if (window == null) state.viewer.windowedMoment0 = null;
    },
    setViewerWindowedMoment0(state, moment0) {
      state.viewer.windowedMoment0 = moment0;
    },
    // Link/unlink the two map panels' pan/zoom; the shared view they render at.
    setViewerLinked(state, linked) {
      state.viewer.linked = linked;
    },
    setViewerMapView(state, view) {
      state.viewer.mapView = view;
    },
  },
  actions: {
    // The list view's loader: corpus overview + object-first list.
    async loadCollection({ commit }) {
      commit("setLoading", true);
      commit("setError", null);
      try {
        const data = await getJson("/api/collection");
        commit("setCollection", data);
      } catch (err) {
        commit("setError", err.message);
      } finally {
        commit("setLoading", false);
      }
    },
    // The detail view's loader: per-store details for one object (cached after first fetch).
    async loadObject({ commit, state }, object) {
      if (state.detailsByObject[object]) return; // already cached
      commit("setLoading", true);
      commit("setError", null);
      try {
        const data = await getJson(`/api/objects/${encodeURIComponent(object)}`);
        commit("setDetail", { object, stores: data.stores });
      } catch (err) {
        commit("setError", err.message);
      } finally {
        commit("setLoading", false);
      }
    },

    // -- cube viewer (#67) --------------------------------------------------------------
    // Open a store in the viewer: reset state, then fetch axes + moment-0 + the first channel in
    // parallel. The axes drive the slider/extent; moment-0 + channel-0 give the two heatmaps. The
    // pixel spectrum waits for the first click (no default pixel is invented).
    async openViewer({ commit, dispatch }, storeId) {
      commit("resetViewer", storeId);
      commit("setViewerLoading", true);
      commit("setViewerError", null);
      try {
        const base = `/api/stores/${storeId}`;
        const [axes, moment0] = await Promise.all([
          getJson(`${base}/cube`),
          getJson(`${base}/moment0`),
        ]);
        commit("setViewerAxes", axes);
        commit("setViewerMoment0", moment0);
        // Start the channel panel at the middle channel — the line is usually near the band centre.
        const mid = axes.n_channels ? Math.floor(axes.n_channels / 2) : 0;
        await dispatch("loadChannel", mid);
      } catch (err) {
        commit("setViewerError", err.message);
      } finally {
        commit("setViewerLoading", false);
      }
    },
    // Fetch one channel's 2-D slice (the velocity slider / keyboard stepping calls this). The index
    // is set optimistically first so the slider thumb tracks the input even before the slice lands.
    // A cached slice is committed immediately (no fetch); only a miss hits the server.
    async loadChannel({ commit, state }, index) {
      const storeId = state.viewer.store;
      if (storeId == null) return;
      commit("setChannelIndex", index);
      const cached = state.viewer.channels[index];
      if (cached) {
        commit("setViewerChannel", cached);
        return;
      }
      try {
        const channel = await getJson(
          `/api/stores/${storeId}/channel/${index}`,
        );
        // Ignore a stale response: the user may have stepped on while this was in flight.
        if (state.viewer.channelIndex === index) {
          commit("setViewerChannel", channel);
        }
      } catch (err) {
        commit("setViewerError", err.message);
      }
    },
    // Select a pixel (the click handler on either map): set the shared crosshair position, then
    // fetch that pixel's spectrum (a single lazy column server-side).
    async selectPixel({ commit, state }, { x, y }) {
      const storeId = state.viewer.store;
      if (storeId == null) return;
      commit("setViewerPosition", { x, y });
      try {
        const spectrum = await getJson(
          `/api/stores/${storeId}/spectrum?x=${x}&y=${y}`,
        );
        // Only apply if the position is still the one requested (clicks can race).
        const p = state.viewer.position;
        if (p && p.x === x && p.y === y) {
          commit("setViewerSpectrum", spectrum);
        }
      } catch (err) {
        commit("setViewerError", err.message);
      }
    },

    // -- cube viewer interactions (#68) -------------------------------------------------
    // Select a polygon region (the map's region-draw finishes): keep the vertices, then POST them
    // to region-spectrum and attach the region-averaged spectrum. A degenerate region (< 3) is
    // never sent (the map only emits a closed ≥3-vertex polygon), but the backend also guards it.
    async selectRegion({ commit, state }, vertices) {
      const storeId = state.viewer.store;
      if (storeId == null) return;
      if (!Array.isArray(vertices) || vertices.length < 3) return;
      commit("setViewerRegion", { vertices, spectrum: null });
      try {
        const response = await fetch(
          `/api/stores/${storeId}/region-spectrum`,
          {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ vertices }),
          },
        );
        if (!response.ok) {
          let detail = `${response.status} ${response.statusText}`;
          try {
            const body = await response.json();
            if (body && body.detail) detail = body.detail;
          } catch {
            // keep the status line
          }
          throw new Error(detail);
        }
        const spectrum = await response.json();
        // Ignore a stale response (the user may have redrawn while this was in flight).
        if (
          state.viewer.region &&
          state.viewer.region.vertices === vertices
        ) {
          commit("setViewerRegionSpectrum", spectrum);
        }
      } catch (err) {
        commit("setViewerError", err.message);
      }
    },
    clearRegion({ commit }) {
      commit("clearViewerRegion");
    },

    // Select a velocity window on the spectrum (the brush): set it, then recompute the integrated
    // map over EXACTLY that window (GET moment0?vmin&vmax) and show it in the map panel. A null
    // window clears both the window and the windowed map (back to the full-band moment-0).
    async selectVelocityWindow({ commit, state }, window) {
      const storeId = state.viewer.store;
      if (storeId == null) return;
      if (window == null) {
        commit("setViewerVelocityWindow", null);
        return;
      }
      const [v0, v1] = window;
      commit("setViewerVelocityWindow", [v0, v1]);
      try {
        const moment0 = await getJson(
          `/api/stores/${storeId}/moment0?vmin=${v0}&vmax=${v1}`,
        );
        // Apply only if this is still the active window (brush drags can race).
        const w = state.viewer.velocityWindow;
        if (w && w[0] === v0 && w[1] === v1) {
          commit("setViewerWindowedMoment0", moment0);
        }
      } catch (err) {
        commit("setViewerError", err.message);
      }
    },

    // Link/unlink the two map panels' pan/zoom (the control). Unlinking leaves each map at its
    // current view; re-linking adopts the shared view so they snap back together.
    setLinked({ commit }, linked) {
      commit("setViewerLinked", linked);
    },
    // The shared pan/zoom view the linked maps render at (one map's zoom drives both when linked).
    setMapView({ commit, state }, view) {
      if (!state.viewer.linked) return;
      commit("setViewerMapView", view);
    },
  },
});
