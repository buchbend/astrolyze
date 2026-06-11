// App bootstrap for the corpus explorer (issue #66): wire Vue + Vuex + Vue Router and mount.
//
// The app is deliberately small — two views (object-first list, per-store detail) over a Vuex
// store that mirrors the public Collection API endpoints one-to-one. The cube viewer (#67/#68) is
// not here; the router leaves a clean seam (a documented placeholder route) but no viewer code.

import { createApp } from "vue";
import App from "./App.vue";
import router from "./router";
import store from "./store";
import "./styles.css";

createApp(App).use(store).use(router).mount("#app");
