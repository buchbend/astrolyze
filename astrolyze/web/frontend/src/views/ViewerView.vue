<script setup>
// The cube viewer — three linked panels for one store of a source (issue #67).
//
// Mirrors the GO VIEW idea, served lazily from the public Cube API: an integrated (moment-0) map, a
// channel map with a velocity SLIDER + KEYBOARD stepping, and a pixel SPECTRUM that updates when the
// user clicks a position on either map. The selected position is SHARED across the panels (a
// crosshair on both maps, the spectrum's anchor) via the Vuex `viewer` slice — this view holds no
// data of its own, it dispatches actions and renders store state (single source of truth).
//
// Keyboard stepping (when the viewer has focus): ArrowLeft/ArrowRight = ∓1 channel; PageUp/PageDown
// and Shift+Arrow = ∓5 channels (the larger-skip modifier). The current velocity is always shown.
//
// The #68 panels (region-averaged spectrum, velocity-window-on-spectrum, link/unlink) are NOT built
// here — their Vuex slots (region, velocityWindow, linked) are reserved but unused, and a labelled
// seam notes where they land.
import { computed, onMounted, onBeforeUnmount, watch } from "vue";
import { useStore } from "vuex";
import HeatMap from "../components/HeatMap.vue";
import SpectrumPlot from "../components/SpectrumPlot.vue";

const props = defineProps({
  object: { type: String, required: true },
  storeId: { type: String, required: true },
});

const store = useStore();

// The larger-skip stride for PageUp/PageDown and Shift+Arrow (the "+N channels" modifier).
const SKIP = 5;

const viewer = computed(() => store.state.viewer);
const axes = computed(() => viewer.value.axes);
const moment0 = computed(() => viewer.value.moment0);
const channel = computed(() => viewer.value.channel);
const spectrum = computed(() => viewer.value.spectrum);
const position = computed(() => viewer.value.position);
const channelIndex = computed(() => viewer.value.channelIndex);
const nChannels = computed(() => (axes.value ? axes.value.n_channels : 0));
const loading = computed(() => viewer.value.loading);
const error = computed(() => viewer.value.error);

const currentVelocity = computed(() => store.getters.currentVelocity);
const velocityUnit = computed(() => store.getters.velocityUnit);

const velocityText = computed(() => {
  if (currentVelocity.value == null) return "—";
  return `${currentVelocity.value.toFixed(2)} ${velocityUnit.value || ""}`.trim();
});

const channelTitle = computed(() => {
  if (channel.value == null) return "Channel map";
  return `Channel ${channelIndex.value} · v = ${velocityText.value}`;
});

function load() {
  store.dispatch("openViewer", props.storeId);
  // Deep-linked viewer: make sure the masthead carries the corpus identity too.
  if (!store.state.rootUri) store.dispatch("loadCollection");
}
onMounted(load);
watch(() => props.storeId, load);

// -- channel stepping (slider + keyboard) ---------------------------------------------
function goToChannel(index) {
  if (nChannels.value === 0) return;
  const clamped = Math.max(0, Math.min(nChannels.value - 1, index));
  if (clamped === channelIndex.value && channel.value) return;
  store.dispatch("loadChannel", clamped);
}

function onSlider(event) {
  goToChannel(Number(event.target.value));
}

function onKey(event) {
  if (nChannels.value === 0) return;
  const big = event.shiftKey ? SKIP : 1;
  let handled = true;
  switch (event.key) {
    case "ArrowRight":
    case "ArrowUp":
      goToChannel(channelIndex.value + big);
      break;
    case "ArrowLeft":
    case "ArrowDown":
      goToChannel(channelIndex.value - big);
      break;
    case "PageUp":
      goToChannel(channelIndex.value + SKIP);
      break;
    case "PageDown":
      goToChannel(channelIndex.value - SKIP);
      break;
    case "Home":
      goToChannel(0);
      break;
    case "End":
      goToChannel(nChannels.value - 1);
      break;
    default:
      handled = false;
  }
  if (handled) event.preventDefault();
}

onMounted(() => window.addEventListener("keydown", onKey));
onBeforeUnmount(() => window.removeEventListener("keydown", onKey));

// -- shared position (click on either map) --------------------------------------------
function onSelect({ x, y }) {
  store.dispatch("selectPixel", { x, y });
}
</script>

<template>
  <section
    class="viewer"
    tabindex="0"
  >
    <div class="page-head">
      <h1 class="page-title">
        {{ object }}
        <span class="muted viewer-store">cube viewer</span>
      </h1>
      <router-link
        :to="{ name: 'detail', params: { object } }"
        class="back-link"
      >
        ← back to {{ object }}
      </router-link>
    </div>

    <div
      v-if="loading && !axes"
      class="state"
    >
      opening cube…
    </div>
    <div
      v-else-if="error"
      class="state error"
    >
      {{ error }}
    </div>

    <template v-else-if="axes">
      <div class="panels">
        <!-- integrated (moment-0) map -->
        <div class="card panel">
          <HeatMap
            v-if="moment0"
            title="Integrated (moment 0)"
            :data="moment0.data"
            :vmin="moment0.vmin"
            :vmax="moment0.vmax"
            :unit="moment0.unit"
            :extent="moment0.extent"
            :downsample="moment0.downsample"
            :position="position"
            @select="onSelect"
          />
        </div>

        <!-- channel map + velocity slider + keyboard stepping -->
        <div class="card panel">
          <HeatMap
            v-if="channel"
            :title="channelTitle"
            :data="channel.data"
            :vmin="axes.value_min_hint"
            :vmax="axes.value_max_hint"
            :unit="channel.unit"
            :extent="channel.extent"
            :downsample="channel.downsample"
            :position="position"
            @select="onSelect"
          />
          <div class="channel-controls">
            <input
              type="range"
              class="vel-slider"
              :min="0"
              :max="Math.max(0, nChannels - 1)"
              :value="channelIndex"
              :aria-label="`channel ${channelIndex} of ${nChannels}`"
              @input="onSlider"
            >
            <div class="vel-readout">
              <span class="num">ch {{ channelIndex }} / {{ nChannels - 1 }}</span>
              <span class="num vel">v = {{ velocityText }}</span>
            </div>
            <p class="kbd-hint muted">
              ← → step ±1 channel · Shift / PageUp-PageDown step ±{{ SKIP }} · Home/End jump
            </p>
          </div>
        </div>

        <!-- pixel spectrum (full-width row) -->
        <div class="card panel panel-wide">
          <SpectrumPlot
            :velocity="spectrum ? spectrum.velocity : null"
            :value="spectrum ? spectrum.value : null"
            :velocity-unit="spectrum ? spectrum.velocity_unit : velocityUnit"
            :value-unit="spectrum ? spectrum.value_unit : axes.bunit"
            :channel-index="channelIndex"
            :velocity-window="viewer.velocityWindow"
          />
          <p
            v-if="position"
            class="pixel-readout num muted"
          >
            pixel (x={{ position.x }}, y={{ position.y }})
          </p>
        </div>
      </div>

      <!-- reserved #68 seam: region-averaged spectrum · velocity-window moment · link/unlink -->
      <div class="viewer-seam">
        <strong>Coming in #68</strong><br>
        region selection &amp; region-averaged spectrum · select a velocity window on the spectrum to
        recompute the moment · link/unlink the maps' zoom
      </div>
    </template>
  </section>
</template>

<style scoped>
.viewer {
  outline: none;
}
.viewer-store {
  font-size: 0.85rem;
  font-weight: 400;
}
.back-link {
  font-size: 0.85rem;
}
.panels {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(440px, 1fr));
  gap: 1.1rem;
}
.panel {
  padding: 1rem 1.1rem 1.1rem;
}
.panel-wide {
  grid-column: 1 / -1;
}
.channel-controls {
  margin-top: 0.7rem;
}
.vel-slider {
  width: 100%;
  accent-color: var(--accent);
}
.vel-readout {
  display: flex;
  justify-content: space-between;
  font-size: 0.82rem;
  margin-top: 0.3rem;
}
.vel-readout .vel {
  color: var(--accent-ink);
  font-weight: 600;
}
.kbd-hint {
  font-size: 0.72rem;
  margin: 0.45rem 0 0;
}
.pixel-readout {
  font-size: 0.78rem;
  margin: 0.4rem 0 0;
}
</style>
