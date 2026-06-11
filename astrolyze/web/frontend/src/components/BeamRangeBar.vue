<script setup>
// A small D3 span bar for a source's beam range (issue #66).
//
// Draws the [min, max] major-axis beam span on a SHARED corpus-wide scale (the `domain` prop), so
// a glance down the list column compares angular resolution across sources: a short bar near the
// left is a fine (interferometer) beam, a long bar reaching right is a coarse (single-dish) one. A
// single-valued range (min == max) renders a dot. The numeric text accompanies the bar so the plot
// is a quick read, not the only read (the table stays the authority — D3 here is sugar, not the
// data path).
//
// D3 is used the minimal, idiomatic way: a d3.scaleLinear maps arcsec -> pixels; the SVG geometry
// is bound reactively (Vue owns the DOM, D3 owns the math) rather than d3 mutating the DOM, which
// keeps it inside Vue's reactivity and avoids a second source of truth.
import { computed } from "vue";
import { scaleLinear } from "d3";

const props = defineProps({
  min: { type: Number, default: null },
  max: { type: Number, default: null },
  domain: { type: Array, required: true }, // [lo, hi] arcsec, corpus-wide
});

const WIDTH = 120;
const HEIGHT = 18;

const scale = computed(() => {
  const [lo, hi] = props.domain;
  // Guard a degenerate domain (single distinct beam in the whole corpus) so the scale is valid.
  const dlo = lo;
  const dhi = hi > lo ? hi : lo + 1;
  return scaleLinear().domain([dlo, dhi]).range([2, WIDTH - 2]);
});

const hasBeam = computed(() => props.min != null && props.max != null);

const x1 = computed(() => (hasBeam.value ? scale.value(props.min) : 0));
const x2 = computed(() => (hasBeam.value ? scale.value(props.max) : 0));
const isPoint = computed(() => hasBeam.value && props.min === props.max);

const label = computed(() => {
  if (!hasBeam.value) return "—";
  if (props.min === props.max) return `${props.min.toFixed(2)}″`;
  return `${props.min.toFixed(2)}–${props.max.toFixed(2)}″`;
});
</script>

<template>
  <div class="beam-range">
    <svg
      :width="WIDTH"
      :height="HEIGHT"
      class="beam-svg"
      role="img"
      :aria-label="`beam major axis range ${label}`"
    >
      <line
        :x1="2"
        :x2="WIDTH - 2"
        :y1="HEIGHT / 2"
        :y2="HEIGHT / 2"
        class="beam-axis"
      />
      <template v-if="hasBeam">
        <line
          v-if="!isPoint"
          :x1="x1"
          :x2="x2"
          :y1="HEIGHT / 2"
          :y2="HEIGHT / 2"
          class="beam-span"
        />
        <circle
          :cx="x1"
          :cy="HEIGHT / 2"
          r="3"
          class="beam-cap"
        />
        <circle
          :cx="x2"
          :cy="HEIGHT / 2"
          r="3"
          class="beam-cap"
        />
      </template>
    </svg>
    <span class="num beam-label">{{ label }}</span>
  </div>
</template>

<style scoped>
.beam-range {
  display: flex;
  align-items: center;
  gap: 0.6rem;
}
.beam-svg {
  flex: none;
}
.beam-axis {
  stroke: var(--line-strong);
  stroke-width: 1;
}
.beam-span {
  stroke: var(--accent);
  stroke-width: 3;
  stroke-linecap: round;
}
.beam-cap {
  fill: var(--accent);
}
.beam-label {
  font-size: 0.8rem;
  color: var(--ink-soft);
  white-space: nowrap;
}
</style>
