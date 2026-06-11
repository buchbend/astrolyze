<script setup>
// The per-store DETAIL view (issue #66).
//
// Mirrors Collection.describe(object): one card per store of the source, showing its physical
// parameters — survey / telescope / species / transition, beam (major × minor @ PA), rest
// frequency, bunit, and provenance (store path + checksum). The describe is SHALLOW (the backend
// does not open stores in this slice), so the native channel width / velocity coverage show as "—"
// (not yet read), exactly as the CLI's shallow describe renders them. Each store carries a clearly
// labelled placeholder for the cube viewer — the reserved #67/#68 seam — but no viewer logic.
import { computed, watch, onMounted } from "vue";
import { useStore } from "vuex";

const props = defineProps({
  object: { type: String, required: true },
  // storeId is present only on the reserved viewer route (#67/#68); unused in this slice.
  storeId: { type: String, default: null },
});

const store = useStore();

function load() {
  store.dispatch("loadObject", props.object);
  // The detail view may be deep-linked before the list ever loaded — make sure the masthead has
  // the corpus identity too.
  if (!store.state.rootUri) store.dispatch("loadCollection");
}
onMounted(load);
watch(() => props.object, load);

const stores = computed(() => store.getters.detailFor(props.object));
const loading = computed(() => store.state.loading);
const error = computed(() => store.state.error);

function beamText(s) {
  if (s.beam_major_arcsec == null || s.beam_minor_arcsec == null) return "—";
  const pa = s.beam_pa_deg == null ? "" : ` @ ${s.beam_pa_deg.toFixed(1)}°`;
  return `${s.beam_major_arcsec.toFixed(2)}″ × ${s.beam_minor_arcsec.toFixed(2)}″${pa}`;
}
function restFreq(s) {
  if (s.rest_frequency_hz == null) return "—";
  return `${(s.rest_frequency_hz / 1e9).toFixed(6)} GHz`;
}
function velocityRange(s) {
  if (s.velocity_min_kms == null || s.velocity_max_kms == null) return "—";
  return `${s.velocity_min_kms.toFixed(2)}–${s.velocity_max_kms.toFixed(2)} km/s`;
}
function channelWidth(s) {
  return s.channel_width_kms == null ? "—" : `${s.channel_width_kms.toFixed(2)} km/s`;
}
</script>

<template>
  <section>
    <div class="page-head">
      <h1 class="page-title">
        {{ object }}
      </h1>
      <span
        v-if="stores && stores.length"
        class="page-sub"
      >{{ stores.length }} store{{ stores.length === 1 ? "" : "s" }}</span>
    </div>

    <div
      v-if="loading && !stores"
      class="state"
    >
      loading {{ object }}…
    </div>
    <div
      v-else-if="error"
      class="state error"
    >
      {{ error }}
    </div>

    <template v-else-if="stores">
      <div class="store-grid">
        <article
          v-for="s in stores"
          :key="s.store_id"
          class="card store-card"
        >
          <header class="store-head">
            <span class="store-survey">{{ s.survey ?? "—" }}</span>
            <span class="store-telescope muted">{{ s.telescope ?? "—" }}</span>
          </header>
          <dl class="params">
            <div class="param">
              <dt>Species</dt>
              <dd>{{ s.species ?? "—" }}</dd>
            </div>
            <div class="param">
              <dt>Transition</dt>
              <dd>{{ s.transition ?? "—" }}</dd>
            </div>
            <div class="param">
              <dt>Beam</dt>
              <dd class="num">
                {{ beamText(s) }}
              </dd>
            </div>
            <div class="param">
              <dt>Rest frequency</dt>
              <dd class="num">
                {{ restFreq(s) }}
              </dd>
            </div>
            <div class="param">
              <dt>Bunit</dt>
              <dd class="num">
                {{ s.bunit ?? "—" }}
              </dd>
            </div>
            <div class="param">
              <dt>Channel width</dt>
              <dd class="num">
                {{ channelWidth(s) }}
              </dd>
            </div>
            <div class="param">
              <dt>Velocity coverage</dt>
              <dd class="num">
                {{ velocityRange(s) }}
              </dd>
            </div>
            <div class="param wide">
              <dt>Store</dt>
              <dd
                class="num store-path"
                :title="s.store_uri"
              >
                {{ s.store_path }}
              </dd>
            </div>
            <div
              v-if="s.content_checksum"
              class="param wide"
            >
              <dt>Checksum</dt>
              <dd class="num muted">
                {{ s.content_checksum }}
              </dd>
            </div>
          </dl>

          <!-- Reserved cube-viewer seam (#67/#68). Labelled placeholder, NO viewer logic. -->
          <div class="viewer-seam">
            <strong>Cube viewer</strong><br>
            integrated map · channel maps · spectra — arrives in #67/#68
          </div>
        </article>
      </div>
    </template>
  </section>
</template>

<style scoped>
.store-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
  gap: 1.1rem;
}
.store-card {
  padding: 1.1rem 1.2rem 1.2rem;
}
.store-head {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  margin-bottom: 0.8rem;
  padding-bottom: 0.6rem;
  border-bottom: 1px solid var(--line);
}
.store-survey {
  font-weight: 700;
  font-size: 1.05rem;
}
.store-telescope {
  font-size: 0.85rem;
}
.params {
  margin: 0;
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 0.55rem 1rem;
}
.param.wide {
  grid-column: 1 / -1;
}
.param dt {
  font-size: 0.7rem;
  letter-spacing: 0.03em;
  text-transform: uppercase;
  color: var(--ink-faint);
}
.param dd {
  margin: 0.1rem 0 0;
  font-size: 0.92rem;
}
.store-path {
  word-break: break-all;
  font-size: 0.82rem;
}
</style>
