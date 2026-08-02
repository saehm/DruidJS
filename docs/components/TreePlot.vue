<script setup lang="ts">
import { extent, interpolateViridis, scaleLinear, schemeObservable10 } from "d3";
import { computed } from "vue";

interface Props {
  /** N ⨯ 2 node positions. */
  data: number[][];
  /** Tree edges as `[u, v, weight]`. */
  edges?: [number, number, number][];
  /** Per-node categorical label, used for colour when `values` is absent. */
  labels?: any[];
  /** Per-node continuous value (e.g. curvature); overrides `labels` for colour. */
  values?: number[];
  /** Explicit per-node colour, taking precedence over `labels`. Needed when there are more classes
   *  than the default 10-colour scheme can distinguish. */
  colors?: string[];
  /** Node indices forming a path to highlight. */
  path?: number[];
  /** Draw edges whose endpoints have different labels in a contrasting colour. */
  showBridges?: boolean;
  radius?: number;
  size?: number;
  selected?: number[];
}

const props = withDefaults(defineProps<Props>(), {
  edges: () => [],
  radius: 4,
  size: 500,
  showBridges: false,
  path: () => [],
  selected: () => [],
});

const emit = defineEmits(["node-click"]);

const margin = 16;

const scales = computed(() => {
  if (!props.data?.length) {
    return { x: scaleLinear(), y: scaleLinear() };
  }
  let [x_min, x_max] = extent(props.data, (d: number[]) => d[0]);
  let [y_min, y_max] = extent(props.data, (d: number[]) => d[1]);
  const x_span = (x_max ?? 1) - (x_min ?? 0);
  const y_span = (y_max ?? 1) - (y_min ?? 0);
  // Square the domain so the tree is not stretched along whichever axis happens to be wider.
  const offset = Math.abs(x_span - y_span) / 2;
  if (x_span > y_span) {
    y_min = (y_min ?? 0) - offset;
    y_max = (y_max ?? 1) + offset;
  } else {
    x_min = (x_min ?? 0) - offset;
    x_max = (x_max ?? 1) + offset;
  }
  return {
    x: scaleLinear()
      .domain([x_min ?? 0, x_max ?? 1])
      .range([margin, props.size - margin]),
    y: scaleLinear()
      .domain([y_min ?? 0, y_max ?? 1])
      .range([props.size - margin, margin]),
  };
});

const valueScale = computed(() => {
  if (!props.values?.length) return null;
  const [lo, hi] = extent(props.values) as [number, number];
  return scaleLinear().domain([lo, hi]).range([0, 1]);
});

const labelOrder = computed(() =>
  props.labels ? Array.from(new Set(props.labels)).sort() : [],
);

function nodeColor(i: number) {
  const scale = valueScale.value;
  if (scale && props.values) return interpolateViridis(scale(props.values[i]));
  if (props.colors?.[i]) return props.colors[i];
  if (props.labels && props.labels[i] !== undefined) {
    const index = labelOrder.value.indexOf(props.labels[i]);
    return schemeObservable10[index % schemeObservable10.length];
  }
  return "var(--vp-c-text-1)";
}

/** Consecutive pairs of the highlighted path, for fast edge lookup. */
const pathEdges = computed(() => {
  const set = new Set<string>();
  for (let i = 0; i + 1 < props.path.length; ++i) {
    const [a, b] = [props.path[i], props.path[i + 1]];
    set.add(a < b ? `${a},${b}` : `${b},${a}`);
  }
  return set;
});

const renderedEdges = computed(() => {
  const { x, y } = scales.value;
  if (!props.data?.length) return [];
  return props.edges.map(([u, v]) => {
    const key = u < v ? `${u},${v}` : `${v},${u}`;
    const onPath = pathEdges.value.has(key);
    const bridge =
      props.showBridges &&
      props.labels &&
      props.labels[u] !== undefined &&
      props.labels[u] !== props.labels[v];
    return {
      key,
      x1: x(props.data[u][0]),
      y1: y(props.data[u][1]),
      x2: x(props.data[v][0]),
      y2: y(props.data[v][1]),
      onPath,
      bridge,
    };
  });
});

const pathSet = computed(() => new Set(props.path));
const selectedSet = computed(() => new Set(props.selected));
</script>

<template>
  <div class="tree-plot" :style="{ maxWidth: props.size + 'px' }">
    <svg :viewBox="`0 0 ${props.size} ${props.size}`" width="100%" height="100%">
      <g class="edges">
        <line
          v-for="edge in renderedEdges"
          :key="edge.key"
          :x1="edge.x1"
          :y1="edge.y1"
          :x2="edge.x2"
          :y2="edge.y2"
          :class="{ bridge: edge.bridge, 'on-path': edge.onPath }"
        />
      </g>
      <g class="nodes">
        <circle
          v-for="(point, i) in data"
          :key="`n-${i}`"
          :cx="scales.x(point[0])"
          :cy="scales.y(point[1])"
          :r="selectedSet.has(i) ? props.radius * 2 : pathSet.has(i) ? props.radius * 1.4 : props.radius"
          :fill="nodeColor(i)"
          :class="{ selected: selectedSet.has(i) }"
          @click="emit('node-click', i)"
        />
      </g>
    </svg>
  </div>
</template>

<style scoped>
.tree-plot {
  width: 100%;
  aspect-ratio: 1 / 1;
}

svg {
  display: block;
  max-width: 100%;
  height: auto;
}

.edges line {
  stroke: var(--vp-c-text-3);
  stroke-opacity: 0.45;
  stroke-width: 1;
}

.edges line.bridge {
  stroke: #e4572e;
  stroke-opacity: 0.95;
  stroke-width: 2.2;
}

.edges line.on-path {
  stroke: var(--vp-c-brand-1);
  stroke-opacity: 1;
  stroke-width: 3;
}

.nodes circle {
  cursor: pointer;
  stroke: var(--vp-c-bg);
  stroke-width: 0.75;
  transition: r 0.15s ease;
}

.nodes circle.selected {
  stroke: var(--vp-c-text-1);
  stroke-width: 2;
}
</style>
