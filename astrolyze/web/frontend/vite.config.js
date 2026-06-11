// Vite build config for the astrolyze corpus explorer (issue #66).
//
// `npm run build` emits the bundle straight into ../static — the directory the Python wheel ships
// as package data and the FastAPI app mounts at "/". Keeping the output INSIDE the package (rather
// than a top-level dist/) is what makes "the built frontend is in the wheel" a one-line packaging
// rule: setuptools package-data globs web/static/**/*.
//
// `base: "./"` makes the built index.html reference its assets by RELATIVE path, so the SPA works
// whether served from "/" or a sub-path, and the StaticFiles mount needs no rewrite.
//
// The dev server proxies /api to the running FastAPI backend (uvicorn on :8000) so `npm run dev`
// gives hot-reload against a real corpus: run `astrolyze explore <corpus>` in one terminal and
// `npm run dev` in another. In production one process serves both (the API and the built SPA), so
// the proxy is dev-only.

import { fileURLToPath, URL } from "node:url";
import { defineConfig } from "vite";
import vue from "@vitejs/plugin-vue";

export default defineConfig({
  plugins: [vue()],
  base: "./",
  resolve: {
    alias: {
      "@": fileURLToPath(new URL("./src", import.meta.url)),
    },
  },
  build: {
    outDir: "../static",
    emptyOutDir: true,
  },
  server: {
    proxy: {
      "/api": {
        target: "http://127.0.0.1:8000",
        changeOrigin: true,
      },
    },
  },
});
