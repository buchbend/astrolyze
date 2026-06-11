// The Vuex store: the single source of truth, mirroring the public Collection API endpoints.
//
// House convention is vuex (global CLAUDE.md). The store holds exactly what the two backend
// endpoints return — the corpus overview + object-first list (GET /api/collection) and the
// per-object store details (GET /api/objects/:object) — with loading/error flags so the views
// stay declarative. No business logic lives here: grouping and aggregation are the library's job
// (Collection.list / .describe); this just fetches and caches what the API serialized.

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
  },
});
