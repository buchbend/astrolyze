// Client-side routes for the corpus explorer (issues #66 + #67).
//
// Three real routes — the object-first LIST view ("/"), the per-store-grouped DETAIL view
// ("/objects/:object"), and the cube VIEWER (#67) reached from a store card's "open viewer"
// affordance. The viewer route keeps the stable `store-viewer` name #66 reserved; #67 repoints it
// from the detail placeholder to the real ViewerView (the three linked panels). The #68 panels
// extend ViewerView in place — no new route.
//
// History mode (createWebHistory) matches the backend SPA fallback: any non-/api path the browser
// opens directly is served index.html by FastAPI, and this router resolves it client-side. base
// "./" in vite.config keeps asset URLs relative so the app works served from "/".

import { createRouter, createWebHistory } from "vue-router";
import ListView from "./views/ListView.vue";
import DetailView from "./views/DetailView.vue";
import ViewerView from "./views/ViewerView.vue";

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
    // The cube viewer (#67): the three linked panels for one store of an object. The route NAME
    // (`store-viewer`) is the stable hook #66 reserved; the component is now the real viewer.
    path: "/objects/:object/stores/:storeId/viewer",
    name: "store-viewer",
    component: ViewerView,
    props: true,
  },
];

const router = createRouter({
  history: createWebHistory(),
  routes,
});

export default router;
