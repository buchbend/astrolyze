<script setup>
// The object-first LIST view — the explorer's landing page (issue #66).
//
// Mirrors Collection.list(): one row per source object, aggregating its surveys, species, store
// count, and beam range. The row is the primary navigation affordance — click it to drill into the
// per-store DETAIL view for that object. The beam range renders both as text and as a small D3
// span bar (BeamRangeBar) so a researcher sees coarse-to-fine resolution at a glance across the
// whole corpus, on a shared scale.
import { computed, onMounted } from "vue";
import { useStore } from "vuex";
import { useRouter } from "vue-router";
import BeamRangeBar from "../components/BeamRangeBar.vue";

const store = useStore();
const router = useRouter();

onMounted(() => store.dispatch("loadCollection"));

const objects = computed(() => store.getters.sortedObjects);
const loading = computed(() => store.state.loading);
const error = computed(() => store.state.error);

// The shared beam scale: the corpus-wide [min, max] major axis, so every row's span bar reads on
// the same axis (a 1.5" interferometer beam and a 13" single-dish beam are comparable by length).
const beamDomain = computed(() => {
  const vals = [];
  for (const o of objects.value) {
    if (o.beam_min_arcsec != null) vals.push(o.beam_min_arcsec);
    if (o.beam_max_arcsec != null) vals.push(o.beam_max_arcsec);
  }
  if (vals.length === 0) return [0, 1];
  return [Math.min(...vals), Math.max(...vals)];
});

function openObject(object) {
  router.push({ name: "detail", params: { object } });
}
</script>

<template>
  <section>
    <div class="page-head">
      <h1 class="page-title">
        Sources
      </h1>
      <span
        v-if="objects.length"
        class="page-sub"
      >{{ objects.length }} object{{ objects.length === 1 ? "" : "s" }} in this
        corpus</span>
    </div>

    <div
      v-if="loading"
      class="state"
    >
      loading corpus…
    </div>
    <div
      v-else-if="error"
      class="state error"
    >
      {{ error }}
    </div>
    <div
      v-else-if="objects.length === 0"
      class="state"
    >
      no datasets in this corpus
    </div>

    <div
      v-else
      class="card"
    >
      <table class="catalog">
        <thead>
          <tr>
            <th>Object</th>
            <th>Surveys</th>
            <th>Species</th>
            <th>Stores</th>
            <th>Beam range (major)</th>
          </tr>
        </thead>
        <tbody>
          <tr
            v-for="o in objects"
            :key="o.object ?? '∅'"
            class="clickable"
            @click="openObject(o.object)"
          >
            <td class="object-name">
              <router-link :to="{ name: 'detail', params: { object: o.object } }">
                {{
                  o.object ?? "—"
                }}
              </router-link>
            </td>
            <td>
              <div
                v-if="o.surveys.length"
                class="chips"
              >
                <span
                  v-for="s in o.surveys"
                  :key="s"
                  class="chip"
                >{{ s }}</span>
              </div>
              <span
                v-else
                class="muted"
              >—</span>
            </td>
            <td>
              <div
                v-if="o.species.length"
                class="chips"
              >
                <span
                  v-for="s in o.species"
                  :key="s"
                  class="chip species"
                >{{
                  s
                }}</span>
              </div>
              <span
                v-else
                class="muted"
              >—</span>
            </td>
            <td class="num">
              {{ o.n_stores }}
            </td>
            <td>
              <BeamRangeBar
                :min="o.beam_min_arcsec"
                :max="o.beam_max_arcsec"
                :domain="beamDomain"
              />
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </section>
</template>
