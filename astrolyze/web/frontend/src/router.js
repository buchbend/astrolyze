// Client-side routes for the corpus explorer (issue #66).
//
// Two real routes — the object-first LIST view ("/") and the per-store-grouped DETAIL view
// ("/objects/:object") — plus ONE reserved seam for the cube viewer (#67/#68): a named
// placeholder route that today renders the detail view's "viewer coming in #67/#68" note rather
// than any viewer logic. Leaving the route name (`store-viewer`) in place means the follow-up
// slices fill a known hook instead of inventing navigation.
//
// History mode (createWebHistory) matches the backend SPA fallback: any non-/api path the browser
// opens directly is served index.html by FastAPI, and this router resolves it client-side. base
// "./" in vite.config keeps asset URLs relative so the app works served from "/".

import { createRouter, createWebHistory } from "vue-router";
import ListView from "./views/ListView.vue";
import DetailView from "./views/DetailView.vue";

const routes = [
  {
    path: "/",
    name: "list",
    component: ListView,
  },
  {
    path: "/objects/:object",
    name: "detail",
    component: DetailView,
    props: true,
  },
  {
    // Reserved cube-viewer seam (#67/#68). It resolves to the detail view today, which shows the
    // "viewer arrives in #67/#68" placeholder for the addressed store — no viewer logic ships in
    // this slice. The route NAME is the stable hook the viewer slices will repoint.
    path: "/objects/:object/stores/:storeId/viewer",
    name: "store-viewer",
    component: DetailView,
    props: true,
  },
];

const router = createRouter({
  history: createWebHistory(),
  routes,
});

export default router;
