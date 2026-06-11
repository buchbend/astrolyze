<script setup>
// The app shell: a fixed masthead (corpus identity + breadcrumb) over the routed view. The
// masthead reads the corpus root + catalog version from the store so it is present on every view;
// the breadcrumb is the whole navigation model in this two-view app (list <-> one object).
import { computed } from "vue";
import { useStore } from "vuex";
import { useRoute } from "vue-router";

const store = useStore();
const route = useRoute();

const rootUri = computed(() => store.state.rootUri);
const catalogVersion = computed(() => store.state.catalogVersion);
// The current object (detail view) drives the breadcrumb tail; null on the list view.
const currentObject = computed(() => route.params.object || null);
</script>

<template>
  <div class="app">
    <header class="masthead">
      <div class="brand">
        <span class="brand-mark">astrolyze</span>
        <span class="brand-sub">corpus explorer</span>
      </div>
      <div
        v-if="rootUri"
        class="corpus-id"
      >
        <span
          class="corpus-root"
          :title="rootUri"
        >{{ rootUri }}</span>
        <span
          v-if="catalogVersion"
          class="corpus-version"
        >catalog {{ catalogVersion }}</span>
      </div>
    </header>

    <nav
      class="breadcrumb"
      aria-label="Breadcrumb"
    >
      <router-link :to="{ name: 'list' }">
        corpus
      </router-link>
      <template v-if="currentObject">
        <span class="sep">/</span>
        <span class="crumb-current">{{ currentObject }}</span>
      </template>
    </nav>

    <main class="content">
      <router-view />
    </main>

    <footer class="footer">
      <span>read-only view over the public Collection + Cube API — list · detail · cube viewer
        (issues #66, #67)</span>
    </footer>
  </div>
</template>
