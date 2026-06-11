/* ESLint config for the corpus explorer frontend (issue #66).
 * Vue 3 recommended + the browser/node globals the app and Vite config use. `npm run lint` runs
 * this in CI alongside the build. Kept minimal — house style, not a style war. */
module.exports = {
  root: true,
  env: { browser: true, es2022: true, node: true },
  extends: ["eslint:recommended", "plugin:vue/vue3-recommended"],
  parserOptions: { ecmaVersion: "latest", sourceType: "module" },
  rules: {
    // Single-word view/component names are clear in this tiny app (ListView, DetailView).
    "vue/multi-word-component-names": "off",
  },
};
