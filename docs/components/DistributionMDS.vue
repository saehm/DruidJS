<script setup>
import { schemeObservable10 } from "d3";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { runInWorker } from "./workerPool.js";

const distributions = [
  { name: "Gaussian 1", mean: 0.3, std: 0.1 },
  { name: "Gaussian 2", mean: 0.4, std: 0.1 },
  { name: "Gaussian 3", mean: 0.7, std: 0.1 },
  { name: "Uniform", type: "uniform" },
  { name: "Bimodal", type: "bimodal", peaks: [0.2, 0.8] },
  { name: "Exponential", type: "exponential" },
  { name: "Wide Gaussian", mean: 0.5, std: 0.3 },
];

const bins = 50;
const results = ref(null);
const loading = ref(true);
const histograms = ref([]);

function generateHistogram(dist) {
  const data = new Float64Array(bins).fill(0);
  const samples = 1000;
  for (let i = 0; i < samples; i++) {
    let val;
    if (dist.type === "uniform") {
      val = Math.random();
    } else if (dist.type === "bimodal") {
      val =
        Math.random() > 0.5
          ? Math.random() * 0.2 + (dist.peaks[0] - 0.1)
          : Math.random() * 0.2 + (dist.peaks[1] - 0.1);
    } else if (dist.type === "exponential") {
      val = -Math.log(1 - Math.random()) * 0.2;
    } else {
      // Box-Muller
      const u1 = Math.random();
      const u2 = Math.random();
      const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      val = z0 * (dist.std || 0.1) + dist.mean;
    }
    const bin = Math.floor(Math.max(0, Math.min(0.99, val)) * bins);
    data[bin]++;
  }
  // Normalize
  for (let i = 0; i < bins; i++) data[i] /= samples;
  return data;
}

function getColor(index) {
  return schemeObservable10[index % schemeObservable10.length];
}

function getHistogramPolygon(hist) {
  if (!hist || hist.length === 0) return "";
  let max = Math.max(...hist);
  if (max === 0) max = 1; // avoid div by 0
  const width = 120;
  const height = 40;
  let points = `0,${height} `;
  for (let i = 0; i < hist.length; i++) {
    const x = (i / hist.length) * width;
    const y = height - (hist[i] / max) * height;
    points += `${x},${y} `;
  }
  points += `${width},${height}`;
  return points;
}

async function run() {
  loading.value = true;
  const rawHistograms = distributions.map((d) => generateHistogram(d));
  histograms.value = rawHistograms;

  try {
    const projection = await runInWorker("DR", "MDS", rawHistograms, {
      metric: "wasserstein",
      d: 2,
    });
    results.value = projection;
  } catch (e) {
    console.error("Error in Distribution MDS:", e);
  } finally {
    loading.value = false;
  }
}

onMounted(() => {
  run();
});
</script>

<template>
  <div class="distribution-mds">
    <div class="header">
      <p class="description">
        Compare stochastic distributions using the Wasserstein metric and MDS.
      </p>
      <button @click="run" :disabled="loading" class="btn">Regenerate Data</button>
    </div>

    <div class="content">
      <div class="plot-box">
        <div v-if="loading" class="loading-overlay">
          <div class="spinner"></div>
          Projected with Wasserstein...
        </div>
        <RoundScatterplot
          v-if="results"
          :data="results"
          :labels="distributions.map((d) => d.name)"
          :colors="distributions.map((d, i) => getColor(i))"
          :size="400"
          :radius="8"
          :style="{ opacity: loading ? 0.3 : 1, transition: 'opacity 0.3s' }"
        />
      </div>

      <div class="legend">
        <h3>Sampled Distributions</h3>
        <div class="legend-items">
          <div v-for="(dist, i) in distributions" :key="i" class="legend-item" :title="dist.name">
            <div class="color-swatch" :style="{ backgroundColor: getColor(i) }"></div>
            <div class="dist-info">
              <span class="dist-name">{{ dist.name }}</span>
              <svg width="120" height="40" class="mini-hist">
                <polygon
                  :fill="getColor(i)"
                  opacity="0.6"
                  :points="getHistogramPolygon(histograms[i])"
                />
                <polyline
                  fill="none"
                  :stroke="getColor(i)"
                  stroke-width="1.5"
                  :points="getHistogramPolygon(histograms[i])"
                />
              </svg>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.distribution-mds {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 24px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
  flex-wrap: wrap;
  gap: 10px;
}

.description {
  margin: 0;
  color: var(--vp-c-text-2);
  font-size: 0.95em;
}

.btn {
  padding: 6px 16px;
  border-radius: 6px;
  border: 1px solid var(--vp-c-divider);
  background: var(--vp-c-bg);
  cursor: pointer;
  font-weight: 500;
  transition: border-color 0.2s;
}

.btn:hover:not(:disabled) {
  border-color: var(--vp-c-brand);
}

.btn:disabled {
  opacity: 0.6;
  cursor: not-allowed;
}

.content {
  display: flex;
  flex-wrap: wrap;
  gap: 40px;
  align-items: flex-start;
}

.plot-box {
  position: relative;
  min-height: 400px;
  min-width: 400px;
  display: flex;
  justify-content: center;
  align-items: center;
  border-radius: 8px;
  background: var(--vp-c-bg);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
}

.loading-overlay {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 12px;
  color: var(--vp-c-text-2);
  z-index: 10;
}

.spinner {
  width: 30px;
  height: 30px;
  border: 3px solid var(--vp-c-divider);
  border-top-color: var(--vp-c-brand);
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

.legend {
  flex: 1;
  min-width: 150px;
  background: var(--vp-c-bg);
  padding: 16px;
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
}

.legend h3 {
  margin: 0 0 16px 0;
  font-size: 1.1em;
  font-weight: 600;
}

.legend-items {
  display: flex;
  flex-direction: column;
  gap: 12px;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 12px;
}

.color-swatch {
  width: 14px;
  height: 14px;
  border-radius: 50%;
  flex-shrink: 0;
}

.dist-info {
  display: flex;
  flex-direction: column;
  gap: 4px;
}

.dist-name {
  font-size: 0.85em;
  font-weight: 500;
  color: var(--vp-c-text-1);
}

.mini-hist {
  border-bottom: 1px solid var(--vp-c-divider);
  background: var(--vp-c-bg-soft);
  border-radius: 2px;
}

@keyframes spin {
  to {
    transform: rotate(360deg);
  }
}
</style>
