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
// the BeamRangeBar / HeatMap idiom. The #68 velocity-window selection (brush on this axis to
// recompute the moment) is a reserved seam: the `velocityWindow` prop is accepted and, if set, the
// band is shaded — but NO brush interaction is wired here (that is #68).
import { computed } from "vue";
import { scaleLinear, line as d3line } from "d3";

const props = defineProps({
  // Parallel arrays from GET .../spectrum: velocity[] and value[] (value may contain nulls).
  velocity: { type: Array, default: null },
  value: { type: Array, default: null },
  velocityUnit: { type: String, default: "" },
  valueUnit: { type: String, default: "" },
  // The current channel index (the slider's), marked with a vertical guide so the panels link.
  channelIndex: { type: Number, default: null },
  // Reserved #68 seam: a [v0, v1] velocity window to shade. NO brush logic in this slice.
  velocityWindow: { type: Array, default: null },
});

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
  const vals = props.value.filter((v) => v != null && Number.isFinite(v));
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

// Reserved #68 seam: shade the velocity window if one is set. No interaction here.
const windowBand = computed(() => {
  if (!Array.isArray(props.velocityWindow) || props.velocityWindow.length !== 2)
    return null;
  const [a, b] = props.velocityWindow.map((v) => xScale.value(v));
  return { x: Math.min(a, b), w: Math.abs(b - a) };
});
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
    >
      <!-- reserved #68 velocity-window band (shaded only if a window is set; no brush) -->
      <rect
        v-if="windowBand"
        :x="windowBand.x"
        :y="PAD.top"
        :width="windowBand.w"
        :height="H - PAD.bottom - PAD.top"
        class="sp-window"
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
