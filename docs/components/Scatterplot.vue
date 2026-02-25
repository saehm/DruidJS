<script setup lang="ts">
import { computed } from "vue";

interface Props {
  /** 2D array of points: [[x1, y1], [x2, y2], ...] */
  data: number[][];
  /** Width of the SVG in pixels */
  width?: number;
  /** Height of the SVG in pixels */
  height?: number;
  /** Margin around the plot area */
  margin?: { top: number; right: number; bottom: number; left: number };
  /** Color of the points (can be a single color or array of colors) */
  color?: string | string[];
  /** Radius of each point in pixels */
  pointRadius?: number;
  /** Whether to show axes */
  showAxes?: boolean;
  /** Whether to show grid lines */
  showGrid?: boolean;
  /** X-axis label */
  xLabel?: string;
  /** Y-axis label */
  yLabel?: string;
  /** Title of the plot */
  title?: string;
}

const props = withDefaults(defineProps<Props>(), {
  width: 600,
  height: 400,
  margin: () => ({ top: 40, right: 40, bottom: 60, left: 60 }),
  color: "#3b82f6",
  pointRadius: 4,
  showAxes: true,
  showGrid: true,
  xLabel: "X",
  yLabel: "Y",
  title: "",
});

// Computed dimensions for the plot area (inside margins)
const innerWidth = computed(() => props.width - props.margin.left - props.margin.right);
const innerHeight = computed(() => props.height - props.margin.top - props.margin.bottom);

// Compute min and max values for scaling
const xValues = computed(() => props.data.map((d: number) => d[0]));
const yValues = computed(() => props.data.map((d: number) => d[1]));

const xMin = computed(() => Math.min(...xValues.value));
const xMax = computed(() => Math.max(...xValues.value));
const yMin = computed(() => Math.min(...yValues.value));
const yMax = computed(() => Math.max(...yValues.value));

// Add some padding to the domain (5%)
const xPadding = computed(() => (xMax.value - xMin.value) * 0.05 || 1);
const yPadding = computed(() => (yMax.value - yMin.value) * 0.05 || 1);

// Scale functions to map data values to SVG coordinates
const xScale = computed(() => {
  const domainMin = xMin.value - xPadding.value;
  const domainMax = xMax.value + xPadding.value;
  return (value: number) => {
    const t = (value - domainMin) / (domainMax - domainMin);
    return t * innerWidth.value;
  };
});

const yScale = computed(() => {
  const domainMin = yMin.value - yPadding.value;
  const domainMax = yMax.value + yPadding.value;
  return (value: number) => {
    const t = (value - domainMin) / (domainMax - domainMin);
    return innerHeight.value - t * innerHeight.value; // Invert Y axis
  };
});

// Generate tick values for axes
const xTicks = computed(() => {
  const count = 5;
  const step = (xMax.value - xMin.value + 2 * xPadding.value) / (count - 1);
  return Array.from({ length: count }, (_, i) => xMin.value - xPadding.value + i * step);
});

const yTicks = computed(() => {
  const count = 5;
  const step = (yMax.value - yMin.value + 2 * yPadding.value) / (count - 1);
  return Array.from({ length: count }, (_, i) => yMin.value - yPadding.value + i * step);
});

// Format number for display
const formatNumber = (n: number) => {
  if (Math.abs(n) >= 1000) return n.toExponential(1);
  if (Math.abs(n) < 0.01) return n.toExponential(1);
  return n.toFixed(2);
};

// Get color for a specific point
const getPointColor = (index: number) => {
  if (Array.isArray(props.color)) {
    return props.color[index % props.color.length];
  }
  return props.color;
};
</script>

<template>
  <div class="scatterplot-container">
    <svg :width="width" :height="height" class="scatterplot">
      <!-- Background -->
      <rect :width="width" :height="height" fill="#ffffff" />

      <!-- Plot area group with margin transform -->
      <g :transform="`translate(${margin.left}, ${margin.top})`">
        <!-- Grid lines (X) -->
        <g v-if="showGrid" class="grid x-grid">
          <line
            v-for="(tick, i) in xTicks"
            :key="`x-grid-${i}`"
            :x1="xScale(tick)"
            :y1="0"
            :x2="xScale(tick)"
            :y2="innerHeight"
            stroke="#e5e7eb"
            stroke-width="1"
          />
        </g>

        <!-- Grid lines (Y) -->
        <g v-if="showGrid" class="grid y-grid">
          <line
            v-for="(tick, i) in yTicks"
            :key="`y-grid-${i}`"
            :x1="0"
            :y1="yScale(tick)"
            :x2="innerWidth"
            :y2="yScale(tick)"
            stroke="#e5e7eb"
            stroke-width="1"
          />
        </g>

        <!-- X Axis -->
        <g v-if="showAxes" class="axis x-axis">
          <line
            :x1="0"
            :y1="innerHeight"
            :x2="innerWidth"
            :y2="innerHeight"
            stroke="#374151"
            stroke-width="2"
          />
          <g
            v-for="(tick, i) in xTicks"
            :key="`x-tick-${i}`"
            :transform="`translate(${xScale(tick)}, ${innerHeight})`"
          >
            <line y2="6" stroke="#374151" stroke-width="2" />
            <text y="20" text-anchor="middle" font-size="12" fill="#374151">
              {{ formatNumber(tick) }}
            </text>
          </g>
          <!-- X Axis Label -->
          <text
            :x="innerWidth / 2"
            :y="innerHeight + 45"
            text-anchor="middle"
            font-size="14"
            font-weight="bold"
            fill="#111827"
          >
            {{ xLabel }}
          </text>
        </g>

        <!-- Y Axis -->
        <g v-if="showAxes" class="axis y-axis">
          <line :x1="0" :y1="0" :x2="0" :y2="innerHeight" stroke="#374151" stroke-width="2" />
          <g
            v-for="(tick, i) in yTicks"
            :key="`y-tick-${i}`"
            :transform="`translate(0, ${yScale(tick)})`"
          >
            <line x2="-6" stroke="#374151" stroke-width="2" />
            <text x="-10" dy="0.32em" text-anchor="end" font-size="12" fill="#374151">
              {{ formatNumber(tick) }}
            </text>
          </g>
          <!-- Y Axis Label -->
          <text
            :x="-innerHeight / 2"
            :y="-45"
            transform="rotate(-90)"
            text-anchor="middle"
            font-size="14"
            font-weight="bold"
            fill="#111827"
          >
            {{ yLabel }}
          </text>
        </g>

        <!-- Data points -->
        <g class="points">
          <circle
            v-for="(point, i) in data"
            :key="`point-${i}`"
            :cx="xScale(point[0])"
            :cy="yScale(point[1])"
            :r="pointRadius"
            :fill="getPointColor(i)"
            stroke="#1f2937"
            stroke-width="1"
            class="data-point"
          >
            <title>
              Point {{ i }}: ({{ formatNumber(point[0]) }}, {{ formatNumber(point[1]) }})
            </title>
          </circle>
        </g>
      </g>

      <!-- Title -->
      <text
        v-if="title"
        :x="width / 2"
        :y="25"
        text-anchor="middle"
        font-size="16"
        font-weight="bold"
        fill="#111827"
      >
        {{ title }}
      </text>
    </svg>
  </div>
</template>

<style scoped>
.scatterplot-container {
  display: inline-block;
  border-radius: 8px;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

.scatterplot {
  display: block;
}

.data-point {
  transition: r 0.2s ease;
  cursor: pointer;
}

.data-point:hover {
  stroke-width: 2;
}
</style>
