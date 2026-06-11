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
// The #68 interactions extend this map: a REGION-DRAW mode (mode="region") collects clicked
// vertices into a polygon and emits it on finish (the region-averaged spectrum), and a shared
// pan/zoom VIEW (the `view` prop + `viewchange` emit) keeps the two map panels LINKED so one map's
// zoom drives both. The polygon is reported in the same full-resolution pixel coordinates as a
// click, so the region endpoint speaks the viewer's one coordinate convention.
import { computed, ref } from "vue";
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
  // Interaction mode: "select" (click → pixel) or "region" (clicks build a polygon).
  mode: { type: String, default: "select" },
  // The committed region polygon to overlay: full-resolution [x, y] vertices, or null.
  region: { type: Array, default: null },
  // The shared pan/zoom view { k, x, y } the linked maps render at (k = scale, x/y = translate).
  view: { type: Object, default: null },
});

const emit = defineEmits(["select", "region", "viewchange"]);

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
const clipId = `hm-clip-${Math.random().toString(36).slice(2)}`;

function fmtDeg(v) {
  return v == null ? "—" : `${v.toFixed(3)}°`;
}

const totalW = PAD.left + PLOT + PAD.right;
const totalH = PAD.top + PLOT + PAD.bottom;

// -- shared pan/zoom view -------------------------------------------------------------
// The plot group is transformed by the shared { k, x, y } view (identity when null). Linking is the
// parent's job: it forwards one map's `viewchange` to the other's `view` prop, so both render the
// same transform. Applying it as one SVG transform keeps every overlay (cells, crosshair, polygon)
// in sync for free — they all live inside the transformed group.
const viewTransform = computed(() => {
  const v = props.view;
  if (!v) return "translate(0,0) scale(1)";
  return `translate(${v.x || 0},${v.y || 0}) scale(${v.k || 1})`;
});

// A wheel zoom about the cursor and a drag pan, emitted as a { k, x, y } view for the parent to
// share. Kept minimal (no d3-zoom behaviour binding) so the transform stays a plain reactive value
// Vue owns; the maths (zoom-about-point) is the standard transform composition.
const ZOOM_MIN = 1;
const ZOOM_MAX = 8;
const dragging = ref(null); // { startX, startY, baseX, baseY } during a pan, else null

function currentView() {
  return props.view || { k: 1, x: 0, y: 0 };
}

function svgPoint(event) {
  const svg = event.currentTarget.ownerSVGElement || event.currentTarget;
  const rect = svg.getBoundingClientRect();
  return {
    sx: (event.clientX - rect.left) * (totalW / rect.width),
    sy: (event.clientY - rect.top) * (totalH / rect.height),
  };
}

function onWheel(event) {
  event.preventDefault();
  const v = currentView();
  const { sx, sy } = svgPoint(event);
  const factor = event.deltaY < 0 ? 1.15 : 1 / 1.15;
  const k = Math.max(ZOOM_MIN, Math.min(ZOOM_MAX, v.k * factor));
  if (k === v.k) return;
  // Keep the point under the cursor fixed: x' = sx - (sx - x) * (k / v.k).
  const ratio = k / v.k;
  const next = {
    k,
    x: sx - (sx - v.x) * ratio,
    y: sy - (sy - v.y) * ratio,
  };
  emit("viewchange", k === ZOOM_MIN ? { k: 1, x: 0, y: 0 } : next);
}

function onPointerDown(event) {
  if (props.mode === "region") return; // region-draw owns clicks; pan only in select mode
  const v = currentView();
  const { sx, sy } = svgPoint(event);
  dragging.value = { startX: sx, startY: sy, baseX: v.x, baseY: v.y };
}
function onPointerMove(event) {
  if (!dragging.value) return;
  const { sx, sy } = svgPoint(event);
  const d = dragging.value;
  emit("viewchange", {
    k: currentView().k,
    x: d.baseX + (sx - d.startX),
    y: d.baseY + (sy - d.startY),
  });
}
function onPointerUp() {
  dragging.value = null;
}

// -- region draw ----------------------------------------------------------------------
// In region mode, each click appends a vertex (full-resolution pixel); a double-click (or clicking
// near the first vertex) closes the polygon and emits it. A polygon needs ≥3 vertices to enclose an
// area — fewer is dropped on close (the backend also guards this, but the UI never sends a degenerate
// region). The in-progress vertices render as a dashed open path with dots.
const drawing = ref([]); // [[x, y], …] full-resolution vertices being placed

// Convert a screen event to the FULL-resolution pixel under the (possibly zoomed) cursor.
function eventToPixel(event) {
  const { sx, sy } = svgPoint(event);
  const v = currentView();
  // Undo the view transform to get plot coords, then plot → cell → full-resolution pixel.
  const px = (sx - v.x) / v.k;
  const py = (sy - v.y) / v.k;
  const c = (px - PAD.left) / cellW.value;
  const r = (py - PAD.top) / cellH.value;
  if (c < 0 || c >= cols.value || r < 0 || r >= rows.value) return null;
  return {
    x: c * props.downsample,
    y: r * props.downsample,
    col: c,
    row: r,
  };
}

// Translate a click on the map area into a FULL-resolution pixel and emit it (select mode), or
// append a polygon vertex (region mode). The view transform is inverted so a click on a zoomed map
// still lands on the right pixel.
function onClick(event) {
  if (!cols.value) return;
  const hit = eventToPixel(event);
  if (!hit) return;
  if (props.mode === "region") {
    // Close the polygon if the click lands near the first vertex (and we have ≥3 total).
    if (drawing.value.length >= 2) {
      const [fx, fy] = drawing.value[0];
      const near =
        Math.abs(hit.x - fx) < props.downsample &&
        Math.abs(hit.y - fy) < props.downsample;
      if (near && drawing.value.length >= 3) {
        finishRegion();
        return;
      }
    }
    drawing.value = [...drawing.value, [hit.x, hit.y]];
    return;
  }
  // select mode: snap to the cell centre (the full-resolution pixel the spectrum endpoint expects).
  const x = Math.round((hit.col + 0.5) * props.downsample - 0.5);
  const y = Math.round((hit.row + 0.5) * props.downsample - 0.5);
  emit("select", { x, y });
}

function onDblClick(event) {
  if (props.mode !== "region") return;
  event.preventDefault();
  finishRegion();
}

function finishRegion() {
  if (drawing.value.length >= 3) {
    emit("region", drawing.value);
  }
  drawing.value = [];
}

// The in-progress polyline (screen coords, inside the transformed group so the view applies).
const drawPath = computed(() => {
  if (!drawing.value.length) return null;
  const pts = drawing.value.map(([x, y]) => {
    const c = x / props.downsample;
    const r = y / props.downsample;
    return [PAD.left + c * cellW.value, PAD.top + r * cellH.value];
  });
  return {
    points: pts,
    line: pts.map((p, i) => `${i ? "L" : "M"}${p[0]},${p[1]}`).join(" "),
  };
});

// The committed region polygon (closed), rendered as a filled outline.
const regionPath = computed(() => {
  if (!Array.isArray(props.region) || props.region.length < 3) return null;
  const pts = props.region.map(([x, y]) => {
    const c = x / props.downsample;
    const r = y / props.downsample;
    return `${PAD.left + c * cellW.value},${PAD.top + r * cellH.value}`;
  });
  return `M${pts.join(" L")} Z`;
});
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
      :class="{ 'mode-region': mode === 'region' }"
      role="img"
      :aria-label="title || 'map'"
      @click="onClick"
      @dblclick="onDblClick"
      @wheel="onWheel"
      @pointerdown="onPointerDown"
      @pointermove="onPointerMove"
      @pointerup="onPointerUp"
      @pointerleave="onPointerUp"
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
        <!-- clip the zoomable content to the plot box so a panned map never spills over the axes -->
        <clipPath :id="clipId">
          <rect
            :x="PAD.left"
            :y="PAD.top"
            :width="PLOT"
            :height="PLOT"
          />
        </clipPath>
      </defs>

      <!-- zoomable content: cells + crosshair + region overlays share one pan/zoom transform -->
      <g :clip-path="`url(#${clipId})`">
        <g :transform="viewTransform">
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

          <!-- committed region polygon (the averaged region) -->
          <path
            v-if="regionPath"
            :d="regionPath"
            class="hm-region"
          />

          <!-- in-progress region draw (open dashed path + vertex dots) -->
          <g
            v-if="drawPath"
            class="hm-region-draw"
          >
            <path
              :d="drawPath.line"
              class="hm-region-draw-line"
            />
            <circle
              v-for="(p, i) in drawPath.points"
              :key="`dv${i}`"
              :cx="p[0]"
              :cy="p[1]"
              r="3"
            />
          </g>

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
        </g>
      </g>

      <!-- map frame (outside the transform so it stays a fixed border) -->
      <rect
        :x="PAD.left"
        :y="PAD.top"
        :width="PLOT"
        :height="PLOT"
        class="hm-frame"
      />

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
.hm-svg.mode-region {
  cursor: crosshair;
}
.hm-region {
  fill: rgba(255, 93, 59, 0.14);
  stroke: #ff5d3b;
  stroke-width: 1.5;
  vector-effect: non-scaling-stroke;
}
.hm-region-draw-line {
  fill: none;
  stroke: #ff5d3b;
  stroke-width: 1.5;
  stroke-dasharray: 4 3;
  vector-effect: non-scaling-stroke;
}
.hm-region-draw circle {
  fill: #ff5d3b;
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
