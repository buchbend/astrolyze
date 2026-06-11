<script setup>
// A D3 line plot of a pixel spectrum — velocity (x) vs value (y) (issue #67).
//
// Draws the (velocity, value) arrays the spectrum endpoint serves as a line, with axes labelled in
// the cube's velocity unit and the data's bunit. Blanked channels (null) break the line into
// segments rather than spanning a gap with a straight line (honest, not interpolated). A vertical
// marker shows the current channel (the one the channel-map slider points at) so the two panels are
// visually linked.
//
// D3 owns the MATH (the linear x/y scales), Vue owns the DOM (the path + ticks are bound), matching
// the BeamRangeBar / HeatMap idiom. The #68 velocity-window selection is now wired: dragging across
// the plot brushes a [v0, v1] band (the `velocityWindow`), emitted so the parent recomputes the
// integrated map over exactly that window; a second click outside clears it. A region-averaged
// spectrum (`regionValue`) is overlaid as a second line so the region and pixel spectra compare.
import { computed, ref } from "vue";
import { scaleLinear, line as d3line } from "d3";

const props = defineProps({
  // Parallel arrays from GET .../spectrum: velocity[] and value[] (value may contain nulls).
  velocity: { type: Array, default: null },
  value: { type: Array, default: null },
  velocityUnit: { type: String, default: "" },
  valueUnit: { type: String, default: "" },
  // The current channel index (the slider's), marked with a vertical guide so the panels link.
  channelIndex: { type: Number, default: null },
  // The [v0, v1] velocity window brushed here (shaded); drives the windowed-moment recompute.
  velocityWindow: { type: Array, default: null },
  // A region-averaged spectrum's values (same velocity axis), overlaid as a second line, or null.
  regionValue: { type: Array, default: null },
});

const emit = defineEmits(["window", "clear-window"]);

const W = 460;
const H = 240;
const PAD = { top: 12, right: 16, bottom: 34, left: 52 };

const hasData = computed(
  () =>
    Array.isArray(props.velocity) &&
    Array.isArray(props.value) &&
    props.velocity.length > 0,
);

const xScale = computed(() => {
  if (!hasData.value) return scaleLinear().domain([0, 1]).range([PAD.left, W - PAD.right]);
  const vs = props.velocity.filter((v) => v != null && Number.isFinite(v));
  const lo = Math.min(...vs);
  const hi = Math.max(...vs);
  return scaleLinear()
    .domain([lo, hi === lo ? lo + 1 : hi])
    .range([PAD.left, W - PAD.right]);
});

const yScale = computed(() => {
  if (!hasData.value) return scaleLinear().domain([0, 1]).range([H - PAD.bottom, PAD.top]);
  // The y range spans BOTH the pixel spectrum and the overlaid region spectrum so neither clips.
  const vals = [
    ...props.value,
    ...(Array.isArray(props.regionValue) ? props.regionValue : []),
  ].filter((v) => v != null && Number.isFinite(v));
  let lo = vals.length ? Math.min(...vals) : 0;
  let hi = vals.length ? Math.max(...vals) : 1;
  if (hi === lo) {
    lo -= 1;
    hi += 1;
  }
  const margin = (hi - lo) * 0.08;
  return scaleLinear()
    .domain([lo - margin, hi + margin])
    .range([H - PAD.bottom, PAD.top]);
});

// Build the line path with gaps at null channels (defined() breaks the line, no fake interpolation).
const path = computed(() => {
  if (!hasData.value) return "";
  const gen = d3line()
    .defined((d) => d.v != null && Number.isFinite(d.v) && d.x != null)
    .x((d) => xScale.value(d.x))
    .y((d) => yScale.value(d.v));
  const pts = props.velocity.map((x, i) => ({ x, v: props.value[i] }));
  return gen(pts) || "";
});

const xTicks = computed(() =>
  xScale.value.ticks(5).map((t) => ({ t, x: xScale.value(t) })),
);
const yTicks = computed(() =>
  yScale.value.ticks(5).map((t) => ({ t, y: yScale.value(t) })),
);

// The current-channel guide (vertical line at the slider's velocity).
const channelGuideX = computed(() => {
  if (
    props.channelIndex == null ||
    !hasData.value ||
    props.channelIndex < 0 ||
    props.channelIndex >= props.velocity.length
  )
    return null;
  const v = props.velocity[props.channelIndex];
  return v == null ? null : xScale.value(v);
});

// The region-averaged spectrum overlaid as a second line (same x = velocity axis), gaps at nulls.
const regionPath = computed(() => {
  if (!hasData.value || !Array.isArray(props.regionValue)) return "";
  const gen = d3line()
    .defined((d) => d.v != null && Number.isFinite(d.v) && d.x != null)
    .x((d) => xScale.value(d.x))
    .y((d) => yScale.value(d.v));
  const pts = props.velocity.map((x, i) => ({ x, v: props.regionValue[i] }));
  return gen(pts) || "";
});

// Shade the velocity window if one is set (the band the moment was integrated over).
const windowBand = computed(() => {
  if (!Array.isArray(props.velocityWindow) || props.velocityWindow.length !== 2)
    return null;
  const [a, b] = props.velocityWindow.map((v) => xScale.value(v));
  return { x: Math.min(a, b), w: Math.abs(b - a) };
});

// -- velocity-window brush ------------------------------------------------------------
// Drag horizontally across the plot to select a [v0, v1] band; on release the window is emitted (the
// parent recomputes the windowed moment). A click without a drag clears the window. The brush works
// in screen x, converted to velocities via the x scale's inverse; the band is clamped to the data
// range so a drag past the edge selects the edge channel.
const brushing = ref(null); // { x0 } screen x where the drag started, else null
const brushPreview = ref(null); // { x, w } screen-space preview rect during a drag

function clampX(sx) {
  const [lo, hi] = [PAD.left, W - PAD.right];
  return Math.max(lo, Math.min(hi, sx));
}

function eventX(event) {
  const svg = event.currentTarget.ownerSVGElement || event.currentTarget;
  const rect = svg.getBoundingClientRect();
  return clampX((event.clientX - rect.left) * (W / rect.width));
}

function onBrushDown(event) {
  if (!hasData.value) return;
  brushing.value = { x0: eventX(event) };
  brushPreview.value = null;
}
function onBrushMove(event) {
  if (!brushing.value) return;
  const x1 = eventX(event);
  const x0 = brushing.value.x0;
  brushPreview.value = { x: Math.min(x0, x1), w: Math.abs(x1 - x0) };
}
function onBrushUp(event) {
  if (!brushing.value) return;
  const x1 = eventX(event);
  const x0 = brushing.value.x0;
  brushing.value = null;
  brushPreview.value = null;
  if (Math.abs(x1 - x0) < 3) {
    // A click (no real drag): clear any existing window.
    if (props.velocityWindow) emit("clear-window");
    return;
  }
  const v0 = xScale.value.invert(Math.min(x0, x1));
  const v1 = xScale.value.invert(Math.max(x0, x1));
  emit("window", [v0, v1]);
}
</script>

<template>
  <figure class="spectrum">
    <figcaption class="sp-title">
      Pixel spectrum
    </figcaption>
    <svg
      :viewBox="`0 0 ${W} ${H}`"
      class="sp-svg"
      role="img"
      aria-label="pixel spectrum"
      @pointerdown="onBrushDown"
      @pointermove="onBrushMove"
      @pointerup="onBrushUp"
      @pointerleave="onBrushUp"
    >
      <!-- the committed velocity-window band (the moment was integrated over this) -->
      <rect
        v-if="windowBand"
        :x="windowBand.x"
        :y="PAD.top"
        :width="windowBand.w"
        :height="H - PAD.bottom - PAD.top"
        class="sp-window"
      />

      <!-- live brush preview while dragging a window -->
      <rect
        v-if="brushPreview"
        :x="brushPreview.x"
        :y="PAD.top"
        :width="brushPreview.w"
        :height="H - PAD.bottom - PAD.top"
        class="sp-brush"
      />

      <!-- y grid + ticks -->
      <g
        v-for="tick in yTicks"
        :key="`y${tick.t}`"
      >
        <line
          :x1="PAD.left"
          :x2="W - PAD.right"
          :y1="tick.y"
          :y2="tick.y"
          class="sp-grid"
        />
        <text
          :x="PAD.left - 6"
          :y="tick.y + 3"
          class="sp-tick end"
        >{{ tick.t }}</text>
      </g>

      <!-- x ticks -->
      <g
        v-for="tick in xTicks"
        :key="`x${tick.t}`"
      >
        <line
          :x1="tick.x"
          :x2="tick.x"
          :y1="H - PAD.bottom"
          :y2="H - PAD.bottom + 4"
          class="sp-grid"
        />
        <text
          :x="tick.x"
          :y="H - PAD.bottom + 16"
          class="sp-tick mid"
        >{{ tick.t }}</text>
      </g>

      <!-- zero line for reference if the y-range straddles it -->
      <line
        v-if="yScale.domain()[0] < 0 && yScale.domain()[1] > 0"
        :x1="PAD.left"
        :x2="W - PAD.right"
        :y1="yScale(0)"
        :y2="yScale(0)"
        class="sp-zero"
      />

      <!-- the region-averaged spectrum (overlaid behind the pixel line for comparison) -->
      <path
        v-if="hasData && regionPath"
        :d="regionPath"
        class="sp-region-line"
      />

      <!-- the spectrum line -->
      <path
        v-if="hasData"
        :d="path"
        class="sp-line"
      />
      <text
        v-else
        :x="W / 2"
        :y="H / 2"
        class="sp-empty"
      >click a position on a map to see its spectrum</text>

      <!-- current-channel guide (links to the channel-map slider) -->
      <line
        v-if="channelGuideX != null"
        :x1="channelGuideX"
        :x2="channelGuideX"
        :y1="PAD.top"
        :y2="H - PAD.bottom"
        class="sp-channel-guide"
      />

      <!-- axis labels -->
      <text
        :x="(PAD.left + W - PAD.right) / 2"
        :y="H - 4"
        class="sp-axis-title mid"
      >velocity [{{ velocityUnit || "—" }}]</text>
      <text
        :x="-(PAD.top + H - PAD.bottom) / 2"
        :y="13"
        class="sp-axis-title mid"
        transform="rotate(-90)"
      >value [{{ valueUnit || "—" }}]</text>
    </svg>
  </figure>
</template>

<style scoped>
.spectrum {
  margin: 0;
}
.sp-title {
  font-size: 0.85rem;
  font-weight: 600;
  color: var(--ink-soft);
  margin-bottom: 0.4rem;
}
.sp-svg {
  max-width: 100%;
  height: auto;
  background: var(--surface);
  cursor: col-resize;
  touch-action: none;
}
.sp-line {
  fill: none;
  stroke: var(--accent);
  stroke-width: 1.5;
}
.sp-grid {
  stroke: var(--line);
  stroke-width: 1;
}
.sp-zero {
  stroke: var(--line-strong);
  stroke-width: 1;
  stroke-dasharray: 3 3;
}
.sp-channel-guide {
  stroke: #ff5d3b;
  stroke-width: 1.5;
  stroke-dasharray: 4 2;
  opacity: 0.85;
}
.sp-window {
  fill: var(--accent-soft);
  opacity: 0.6;
}
.sp-brush {
  fill: var(--accent);
  opacity: 0.18;
}
.sp-region-line {
  fill: none;
  stroke: #ff5d3b;
  stroke-width: 1.5;
  stroke-dasharray: 5 2;
  opacity: 0.9;
}
.sp-tick {
  font-family: var(--mono);
  font-size: 0.62rem;
  fill: var(--ink-soft);
}
.sp-tick.end {
  text-anchor: end;
}
.sp-tick.mid {
  text-anchor: middle;
}
.sp-axis-title {
  font-size: 0.7rem;
  fill: var(--ink-soft);
}
.sp-axis-title.mid {
  text-anchor: middle;
}
.sp-empty {
  text-anchor: middle;
  font-size: 0.78rem;
  fill: var(--ink-faint);
}
</style>
