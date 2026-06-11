<script setup>
// A D3 heatmap for a 2-D corpus map — the moment-0 map and the channel map both use it (issue #67).
//
// Draws a value array (rows × cols, with nulls for blanked pixels) as a colour grid, with a colour
// scale legend and RA/Dec extent labels. A click reports the PIXEL position (image x = column, y =
// row, in FULL-resolution coordinates even when the server decimated the array) so the parent can
// fetch that pixel's spectrum and share the crosshair across panels. The selected position renders
// as a crosshair marker so both maps and the spectrum agree on "where".
//
// D3 is used the minimal, idiomatic way (matching BeamRangeBar): d3 owns the MATH (the sequential
// colour scale, the pixel→screen geometry) and Vue owns the DOM (the <rect> grid is bound
// reactively), so there is one source of truth and the SVG stays inside Vue's reactivity. The grid
// is drawn as one <rect> per cell — fine for the bounded map sizes the backend serves (it decimates
// anything over MAX_MAP_DIM), and it keeps hit-testing trivial (no canvas pixel readback).
import { computed } from "vue";
import { scaleSequential, interpolateViridis } from "d3";

const props = defineProps({
  // The 2-D value array (rows of columns); null entries are blanked (drawn transparent).
  data: { type: Array, default: null },
  // Colour-scale bounds; when null they are derived from the data's finite min/max.
  vmin: { type: Number, default: null },
  vmax: { type: Number, default: null },
  // The value unit, shown on the legend.
  unit: { type: String, default: "" },
  // Sky bounding box { ra_min_deg, ra_max_deg, dec_min_deg, dec_max_deg } or null.
  extent: { type: Object, default: null },
  // The server's decimation factor: a clicked cell maps back to factor*col / factor*row in the
  // FULL-resolution pixel grid the spectrum/channel endpoints expect.
  downsample: { type: Number, default: 1 },
  // The shared crosshair position { x, y } in FULL-resolution pixels, or null.
  position: { type: Object, default: null },
  // A short title above the map (e.g. "Integrated (moment 0)" or "Channel 12 · v = 5.0 km/s").
  title: { type: String, default: "" },
});

const emit = defineEmits(["select"]);

// Layout: a fixed drawing box with room for axis labels and a legend strip on the right.
const PLOT = 360; // square map area in px
const PAD = { top: 8, right: 64, bottom: 28, left: 44 };
const LEGEND_W = 14;

const rows = computed(() => (props.data ? props.data.length : 0));
const cols = computed(() => (props.data && props.data.length ? props.data[0].length : 0));

const cellW = computed(() => (cols.value ? PLOT / cols.value : 0));
const cellH = computed(() => (rows.value ? PLOT / rows.value : 0));

// The colour scale: viridis is the perceptually-uniform default for scalar maps. The domain is the
// caller's [vmin, vmax] when given, else the data's finite extent (a flat map gets a degenerate
// domain nudged so the scale stays valid).
const domain = computed(() => {
  let lo = props.vmin;
  let hi = props.vmax;
  if (lo == null || hi == null) {
    let dmin = Infinity;
    let dmax = -Infinity;
    for (const row of props.data || []) {
      for (const v of row) {
        if (v == null || !Number.isFinite(v)) continue;
        if (v < dmin) dmin = v;
        if (v > dmax) dmax = v;
      }
    }
    if (dmin === Infinity) return [0, 1];
    lo = lo == null ? dmin : lo;
    hi = hi == null ? dmax : hi;
  }
  if (hi <= lo) hi = lo + 1;
  return [lo, hi];
});

const color = computed(() =>
  scaleSequential(interpolateViridis).domain(domain.value),
);

// Flatten the grid into bound cells: image y is the ROW index, drawn top-to-bottom. (Astronomical
// "up" handling is the WCS/extent labels' job; the array is served in the natural numpy row order.)
const cells = computed(() => {
  const out = [];
  const w = cellW.value;
  const h = cellH.value;
  for (let r = 0; r < rows.value; r++) {
    for (let c = 0; c < cols.value; c++) {
      const v = props.data[r][c];
      out.push({
        key: `${r}-${c}`,
        x: PAD.left + c * w,
        y: PAD.top + r * h,
        w,
        h,
        fill: v == null || !Number.isFinite(v) ? "none" : color.value(v),
        blank: v == null || !Number.isFinite(v),
      });
    }
  }
  return out;
});

// The crosshair, in screen coords: map a FULL-resolution pixel back to the (possibly decimated)
// displayed cell centre. A position outside the displayed grid is hidden.
const crosshair = computed(() => {
  if (!props.position || !cols.value) return null;
  const c = props.position.x / props.downsample;
  const r = props.position.y / props.downsample;
  if (c < 0 || c >= cols.value || r < 0 || r >= rows.value) return null;
  return {
    cx: PAD.left + (c + 0.5) * cellW.value,
    cy: PAD.top + (r + 0.5) * cellH.value,
  };
});

// Legend: a vertical gradient of the colour scale with min/max labels.
const legendStops = computed(() => {
  const [lo, hi] = domain.value;
  const n = 16;
  const stops = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    stops.push({ offset: `${t * 100}%`, color: color.value(lo + t * (hi - lo)) });
  }
  return stops;
});

const legendId = `hm-grad-${Math.random().toString(36).slice(2)}`;

function fmtDeg(v) {
  return v == null ? "—" : `${v.toFixed(3)}°`;
}

const totalW = PAD.left + PLOT + PAD.right;
const totalH = PAD.top + PLOT + PAD.bottom;

// Translate a click on the map area into a FULL-resolution pixel and emit it. Uses the SVG-local
// offset (the rect's geometry is known), so no getBoundingClientRect round-trip is needed.
function onClick(event) {
  if (!cols.value) return;
  const svg = event.currentTarget;
  const rect = svg.getBoundingClientRect();
  // The svg is rendered at its viewBox size (width/height attrs == viewBox), so client→viewBox is a
  // simple scale by the rendered/intrinsic ratio.
  const sx = (event.clientX - rect.left) * (totalW / rect.width);
  const sy = (event.clientY - rect.top) * (totalH / rect.height);
  const c = Math.floor((sx - PAD.left) / cellW.value);
  const r = Math.floor((sy - PAD.top) / cellH.value);
  if (c < 0 || c >= cols.value || r < 0 || r >= rows.value) return;
  // Map the displayed cell back to a full-resolution pixel (cell centre of the decimated block).
  const x = Math.round((c + 0.5) * props.downsample - 0.5);
  const y = Math.round((r + 0.5) * props.downsample - 0.5);
  emit("select", { x, y });
}
</script>

<template>
  <figure class="heatmap">
    <figcaption
      v-if="title"
      class="hm-title"
    >
      {{ title }}
    </figcaption>
    <svg
      :viewBox="`0 0 ${totalW} ${totalH}`"
      :width="totalW"
      :height="totalH"
      class="hm-svg"
      role="img"
      :aria-label="title || 'map'"
      @click="onClick"
    >
      <defs>
        <linearGradient
          :id="legendId"
          x1="0"
          y1="1"
          x2="0"
          y2="0"
        >
          <stop
            v-for="(s, i) in legendStops"
            :key="i"
            :offset="s.offset"
            :stop-color="s.color"
          />
        </linearGradient>
      </defs>

      <!-- the cell grid (Vue-bound rects; D3 supplied the colour) -->
      <rect
        v-for="cell in cells"
        :key="cell.key"
        :x="cell.x"
        :y="cell.y"
        :width="cell.w + 0.5"
        :height="cell.h + 0.5"
        :fill="cell.fill"
        :class="{ blank: cell.blank }"
      />

      <!-- map frame -->
      <rect
        :x="PAD.left"
        :y="PAD.top"
        :width="PLOT"
        :height="PLOT"
        class="hm-frame"
      />

      <!-- crosshair at the shared selected position -->
      <g
        v-if="crosshair"
        class="hm-cross"
      >
        <line
          :x1="crosshair.cx - 8"
          :x2="crosshair.cx + 8"
          :y1="crosshair.cy"
          :y2="crosshair.cy"
        />
        <line
          :x1="crosshair.cx"
          :x2="crosshair.cx"
          :y1="crosshair.cy - 8"
          :y2="crosshair.cy + 8"
        />
        <circle
          :cx="crosshair.cx"
          :cy="crosshair.cy"
          r="5"
        />
      </g>

      <!-- extent labels (RA along the bottom, Dec along the left) -->
      <text
        v-if="extent"
        :x="PAD.left"
        :y="totalH - 8"
        class="hm-axis-label start"
      >RA {{ fmtDeg(extent.ra_max_deg) }}</text>
      <text
        v-if="extent"
        :x="PAD.left + PLOT"
        :y="totalH - 8"
        class="hm-axis-label end"
      >{{ fmtDeg(extent.ra_min_deg) }}</text>
      <text
        v-if="extent"
        :x="6"
        :y="PAD.top + 10"
        class="hm-axis-label"
      >Dec {{ fmtDeg(extent.dec_max_deg) }}</text>
      <text
        v-if="extent"
        :x="6"
        :y="PAD.top + PLOT"
        class="hm-axis-label"
      >{{ fmtDeg(extent.dec_min_deg) }}</text>

      <!-- colour legend -->
      <rect
        :x="PAD.left + PLOT + 14"
        :y="PAD.top"
        :width="LEGEND_W"
        :height="PLOT"
        :fill="`url(#${legendId})`"
        class="hm-frame"
      />
      <text
        :x="PAD.left + PLOT + 14 + LEGEND_W + 4"
        :y="PAD.top + 8"
        class="hm-legend-label"
      >{{ domain[1].toPrecision(3) }}</text>
      <text
        :x="PAD.left + PLOT + 14 + LEGEND_W + 4"
        :y="PAD.top + PLOT"
        class="hm-legend-label"
      >{{ domain[0].toPrecision(3) }}</text>
      <text
        v-if="unit"
        :x="PAD.left + PLOT + 14 + LEGEND_W / 2"
        :y="PAD.top + PLOT + 18"
        class="hm-legend-unit"
      >{{ unit }}</text>
    </svg>
  </figure>
</template>

<style scoped>
.heatmap {
  margin: 0;
}
.hm-title {
  font-size: 0.85rem;
  font-weight: 600;
  color: var(--ink-soft);
  margin-bottom: 0.4rem;
}
.hm-svg {
  max-width: 100%;
  height: auto;
  cursor: crosshair;
  background: var(--surface);
}
.hm-frame {
  fill: none;
  stroke: var(--line-strong);
  stroke-width: 1;
}
.blank {
  fill: var(--paper);
}
.hm-cross line {
  stroke: #ff5d3b;
  stroke-width: 1.5;
}
.hm-cross circle {
  fill: none;
  stroke: #ff5d3b;
  stroke-width: 1.5;
}
.hm-axis-label {
  font-family: var(--mono);
  font-size: 0.62rem;
  fill: var(--ink-soft);
}
.hm-axis-label.end {
  text-anchor: end;
}
.hm-legend-label {
  font-family: var(--mono);
  font-size: 0.62rem;
  fill: var(--ink-soft);
}
.hm-legend-unit {
  font-family: var(--mono);
  font-size: 0.6rem;
  fill: var(--ink-faint);
  text-anchor: middle;
}
</style>
