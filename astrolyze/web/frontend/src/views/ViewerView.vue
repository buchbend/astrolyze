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
// The #68 interactions are wired here: a REGION-DRAW toggle puts the moment-0 map into polygon mode
// (a drawn polygon → its region-averaged spectrum, overlaid on the pixel spectrum); a velocity-window
// BRUSH on the spectrum recomputes the integrated map over exactly that window (replacing the
// moment-0 panel until cleared); and a LINK/UNLINK control shares (or splits) the two map panels'
// pan/zoom. The state lives in the Vuex `viewer` slice (region, velocityWindow, windowedMoment0,
// linked, mapView).
import { computed, onMounted, onBeforeUnmount, ref, watch } from "vue";
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

// -- #68 interactions -----------------------------------------------------------------
const region = computed(() => viewer.value.region);
const velocityWindow = computed(() => viewer.value.velocityWindow);
const windowedMoment0 = computed(() => viewer.value.windowedMoment0);
const linked = computed(() => viewer.value.linked);
const mapView = computed(() => viewer.value.mapView);

// The moment-0 panel shows the windowed map when a velocity window is active, else the full band.
const momentPanel = computed(() => windowedMoment0.value || moment0.value);
const momentTitle = computed(() =>
  windowedMoment0.value
    ? `Integrated · v ∈ [${velocityWindow.value[0].toFixed(1)}, ${velocityWindow.value[1].toFixed(1)}] ${velocityUnit.value || ""}`.trim()
    : "Integrated (moment 0)",
);

// Region-draw mode: a local toggle that puts the moment-0 map into polygon mode.
const regionMode = ref(false);
function toggleRegionMode() {
  regionMode.value = !regionMode.value;
}
function onRegion(vertices) {
  regionMode.value = false; // one polygon at a time; drawing finished
  store.dispatch("selectRegion", vertices);
}
function clearRegion() {
  store.dispatch("clearRegion");
}

// Velocity-window brush on the spectrum.
function onWindow(window) {
  store.dispatch("selectVelocityWindow", window);
}
function onClearWindow() {
  store.dispatch("selectVelocityWindow", null);
}

// Per-map local view when UNLINKED (so each map keeps its own pan/zoom); shared view when LINKED.
const momentLocalView = ref(null);
const channelLocalView = ref(null);
function onMomentView(view) {
  momentLocalView.value = view;
  if (linked.value) store.dispatch("setMapView", view);
}
function onChannelView(view) {
  channelLocalView.value = view;
  if (linked.value) store.dispatch("setMapView", view);
}
const momentView = computed(() =>
  linked.value ? mapView.value : momentLocalView.value,
);
const channelView = computed(() =>
  linked.value ? mapView.value : channelLocalView.value,
);
function toggleLinked() {
  store.dispatch("setLinked", !linked.value);
}

// The region-averaged spectrum's values, overlaid on the pixel spectrum (aligned to its velocity).
const regionSpectrumValue = computed(() =>
  region.value && region.value.spectrum ? region.value.spectrum.value : null,
);
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
      <!-- #68 viewer interaction controls: region draw + link/unlink -->
      <div class="viewer-toolbar">
        <button
          type="button"
          class="tool-btn"
          :class="{ active: regionMode }"
          @click="toggleRegionMode"
        >
          {{ regionMode ? "drawing region… (double-click to close)" : "✏ draw region" }}
        </button>
        <button
          v-if="region"
          type="button"
          class="tool-btn"
          @click="clearRegion"
        >
          clear region
        </button>
        <button
          v-if="windowedMoment0"
          type="button"
          class="tool-btn"
          @click="onClearWindow"
        >
          clear velocity window
        </button>
        <span class="tool-spacer" />
        <button
          type="button"
          class="tool-btn link-btn"
          :class="{ active: linked }"
          @click="toggleLinked"
        >
          {{ linked ? "🔗 maps linked" : "⛓ maps unlinked" }}
        </button>
      </div>

      <div class="panels">
        <!-- integrated (moment-0) map — region-draw target; shows the windowed map when set -->
        <div class="card panel">
          <HeatMap
            v-if="momentPanel"
            :title="momentTitle"
            :data="momentPanel.data"
            :vmin="momentPanel.vmin"
            :vmax="momentPanel.vmax"
            :unit="momentPanel.unit"
            :extent="momentPanel.extent"
            :downsample="momentPanel.downsample"
            :position="position"
            :mode="regionMode ? 'region' : 'select'"
            :region="region ? region.vertices : null"
            :view="momentView"
            @select="onSelect"
            @region="onRegion"
            @viewchange="onMomentView"
          />
          <p
            v-if="region"
            class="region-readout num muted"
          >
            region: {{ region.vertices.length }} vertices<span v-if="region.spectrum">
              · {{ region.spectrum.n_pixels }} px averaged</span>
          </p>
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
            :view="channelView"
            @select="onSelect"
            @viewchange="onChannelView"
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

        <!-- pixel spectrum (full-width row) — brush a velocity window; overlays the region spectrum -->
        <div class="card panel panel-wide">
          <SpectrumPlot
            :velocity="spectrum ? spectrum.velocity : null"
            :value="spectrum ? spectrum.value : null"
            :velocity-unit="spectrum ? spectrum.velocity_unit : velocityUnit"
            :value-unit="spectrum ? spectrum.value_unit : axes.bunit"
            :channel-index="channelIndex"
            :velocity-window="velocityWindow"
            :region-value="regionSpectrumValue"
            @window="onWindow"
            @clear-window="onClearWindow"
          />
          <div class="spectrum-legend muted num">
            <span
              v-if="position"
              class="pixel-readout"
            >pixel (x={{ position.x }}, y={{ position.y }})</span>
            <span
              v-if="regionSpectrumValue"
              class="legend-region"
            >— region average ({{ region.spectrum.n_pixels }} px)</span>
            <span class="legend-hint">drag the spectrum to integrate a velocity window</span>
          </div>
        </div>
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
}
.viewer-toolbar {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  margin-bottom: 0.9rem;
  flex-wrap: wrap;
}
.tool-spacer {
  flex: 1;
}
.tool-btn {
  font-size: 0.78rem;
  padding: 0.32rem 0.7rem;
  border: 1px solid var(--line-strong);
  border-radius: 6px;
  background: var(--surface);
  color: var(--ink-soft);
  cursor: pointer;
}
.tool-btn:hover {
  border-color: var(--accent);
  color: var(--accent-ink);
}
.tool-btn.active {
  background: var(--accent-soft);
  border-color: var(--accent);
  color: var(--accent-ink);
  font-weight: 600;
}
.link-btn.active {
  background: var(--accent-soft);
}
.region-readout {
  font-size: 0.74rem;
  margin: 0.4rem 0 0;
}
.spectrum-legend {
  display: flex;
  gap: 1rem;
  flex-wrap: wrap;
  font-size: 0.74rem;
  margin: 0.4rem 0 0;
}
.legend-region {
  color: #ff5d3b;
}
.legend-hint {
  font-style: italic;
}
</style>
