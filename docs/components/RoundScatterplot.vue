<script setup lang="ts">
import { extent, scaleLinear, schemeObservable10 } from "d3";
import { computed, ref } from "vue";

interface Props {
  data: number[][];
  labels?: any[];
  colors?: string[];
  radius?: number;
  margin?: { top: number; right: number; bottom: number; left: number };
  size?: number;
}

const props = withDefaults(defineProps<Props>(), {
  radius: 5,
  margin: () => ({ top: 10, right: 10, bottom: 10, left: 10 }),
  size: 500,
});

const mid = computed(() => props.size / 2);

const dataExtent = computed(() => {
  if (!props.data || props.data.length === 0)
    return [
      [0, 1],
      [0, 1],
    ];
  let [x_min, x_max] = extent(props.data, (d: number[]) => d[0]);
  let [y_min, y_max] = extent(props.data, (d: number[]) => d[1]);

  const x_span = (x_max ?? 1) - (x_min ?? 0);
  const y_span = (y_max ?? 1) - (y_min ?? 0);
  const offset = Math.abs(x_span - y_span) / 2;
  if (x_span > y_span) {
    y_min = (y_min ?? 0) - offset;
    y_max = (y_max ?? 1) + offset;
  } else {
    x_min = (x_min ?? 0) - offset;
    x_max = (x_max ?? 1) + offset;
  }
  return [
    [x_min, x_max],
    [y_min, y_max],
  ];
});

const xScale = computed(() =>
  scaleLinear()
    .domain(dataExtent.value[0])
    .range([props.margin.left, props.size - props.margin.right]),
);
const yScale = computed(() =>
  scaleLinear()
    .domain(dataExtent.value[1])
    .range([props.size - props.margin.bottom, props.margin.top]),
);

const getFill = (i: number) => {
  if (props.colors && props.colors[i]) return props.colors[i];
  if (props.labels && props.labels[i] !== undefined) {
    // Find distinct labels to map them to colors
    const labelsSet = Array.from(new Set(props.labels)).sort();
    const colorIndex = labelsSet.indexOf(props.labels[i]);
    return schemeObservable10[colorIndex % schemeObservable10.length];
  }
  return "var(--vp-c-text-1)";
};

const emit = defineEmits(["pointer-move"]);

// Tooltip State
const hoveredPoint = ref<number | null>(null);
const tooltipStyle = ref<Record<string, any>>({ top: "0px", left: "0px", opacity: 0 });
const svgRef = ref<SVGSVGElement | null>(null);

const onMouseOver = (event: MouseEvent, index: number) => {
  hoveredPoint.value = index;
  updateTooltipPosition(event);
};

const onMouseMove = (event: MouseEvent) => {
  if (hoveredPoint.value !== null) {
    updateTooltipPosition(event);
  }

  if (svgRef.value) {
    const pt = svgRef.value.createSVGPoint();
    pt.x = event.clientX;
    pt.y = event.clientY;
    const ctm = svgRef.value.getScreenCTM();
    if (ctm) {
      const svgP = pt.matrixTransform(ctm.inverse());
      const xData = xScale.value.invert(svgP.x);
      const yData = yScale.value.invert(svgP.y);
      emit("pointer-move", [xData, yData]);
    }
  }
};

const onMouseLeave = () => {
  hoveredPoint.value = null;
  tooltipStyle.value.opacity = 0;
};

const updateTooltipPosition = (event: MouseEvent) => {
  // Add small offset to keep it comfortably next to the mouse
  tooltipStyle.value = {
    top: `${event.offsetY - 40}px`,
    left: `${event.offsetX + 15}px`,
    opacity: 1,
  };
};
</script>

<template>
  <div
    class="scatterplot-wrapper"
    :style="{ maxWidth: props.size + 'px', width: '100%', aspectRatio: '1 / 1' }"
  >
    <svg
      ref="svgRef"
      :viewBox="`0 0 ${props.size} ${props.size}`"
      width="100%"
      height="100%"
      @mousemove="onMouseMove"
    >
      <g class="points">
        <circle
          v-for="(point, i) in data"
          :key="`point-${i}`"
          :cx="xScale(point[0])"
          :cy="yScale(point[1])"
          :r="hoveredPoint === i ? props.radius * 1.5 : props.radius"
          :fill="getFill(i)"
          :stroke="hoveredPoint === i ? 'var(--vp-c-text-1)' : 'var(--vp-c-bg)'"
          :stroke-width="hoveredPoint === i ? 1.5 : 0.5"
          @mouseover="onMouseOver($event, i)"
          @mouseleave="onMouseLeave"
          style="cursor: crosshair"
        ></circle>
      </g>
    </svg>

    <div class="custom-tooltip" :style="tooltipStyle" v-show="hoveredPoint !== null">
      <template v-if="hoveredPoint !== null">
        <template v-if="labels && labels[hoveredPoint] !== undefined">
          <div class="tooltip-label">
            <span
              class="color-indicator"
              :style="{ backgroundColor: getFill(hoveredPoint) }"
            ></span>
            {{ labels[hoveredPoint] }}
          </div>
        </template>
        <template v-else>
          <div class="tooltip-label">Point {{ hoveredPoint }}</div>
          <div class="tooltip-coords">
            [{{ data[hoveredPoint][0].toFixed(2) }}, {{ data[hoveredPoint][1].toFixed(2) }}]
          </div>
        </template>
      </template>
    </div>
  </div>
</template>

<style scoped>
.scatterplot-wrapper {
  position: relative;
  display: flex;
  justify-content: center;
  align-items: center;
}

svg {
  display: block;
  max-width: 100%;
  height: auto;
}

circle {
  transition:
    r 0.2s cubic-bezier(0.175, 0.885, 0.32, 1.275),
    stroke-width 0.2s ease;
}

.custom-tooltip {
  position: absolute;
  background: var(--vp-c-bg-elv);
  border: 1px solid var(--vp-c-divider);
  padding: 8px 12px;
  border-radius: 6px;
  font-size: 0.85rem;
  pointer-events: none;
  box-shadow: 0 4px 16px rgba(0, 0, 0, 0.15);
  color: var(--vp-c-text-1);
  z-index: 100;
  transition: opacity 0.15s ease-out;
  white-space: nowrap;
}

.tooltip-label {
  font-weight: 600;
  display: flex;
  align-items: center;
  gap: 8px;
}

.color-indicator {
  display: inline-block;
  width: 10px;
  height: 10px;
  border-radius: 50%;
  border: 1px solid rgba(0, 0, 0, 0.1);
}

.tooltip-coords {
  font-family: var(--vp-font-family-mono);
  font-size: 0.9em;
  color: var(--vp-c-text-2);
  margin-top: 4px;
}
</style>
