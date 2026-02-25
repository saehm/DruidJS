<script setup>
import * as mistle from "@saehrimnir/mistle";
import { computed, onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const datasetOptions = [
  { name: "Iris", key: "IRIS" },
  { name: "Wine", key: "WINE" },
  { name: "Penguins", key: "PENGUINS" },
];

const selectedDataset = ref("IRIS");
const labels = computed(() => mistle[selectedDataset.value].labels);

const metrics = [
  { name: "Euclidean", key: "euclidean" },
  { name: "Cosine", key: "cosine" },
  { name: "Manhattan", key: "manhattan" },
  { name: "Chebyshev", key: "chebyshev" },
];

const results = ref({});
const isLoading = ref(false);

const computeProjections = async () => {
  results.value = {};
  isLoading.value = true;
  const rawData = applyZScore(mistle[selectedDataset.value].values);
  const data = JSON.parse(JSON.stringify(rawData));

  metrics.forEach(async (metric) => {
    try {
      const result = await runInWorker("DR", "SAMMON", data, { metric: metric.key }, 100);
      results.value[metric.key] = result;
    } catch (e) {
      console.error(`Error with metric ${metric.name} in worker:`, e);
    }
  });
};

onMounted(() => {
  computeProjections();
});
</script>

<template>
  <div class="metric-comparison">
    <div class="controls">
      <label>Dataset: </label>
      <select v-model="selectedDataset" @change="computeProjections">
        <option v-for="opt in datasetOptions" :key="opt.key" :value="opt.key">
          {{ opt.name }}
        </option>
      </select>
    </div>

    <div class="description">
      <p>
        Compare how different distance measures (metrics) affect a
        <strong>Sammon Mapping</strong> projection. Each metric defines "distance" differently,
        leading to varied cluster shapes.
      </p>
    </div>

    <div class="grid">
      <div v-for="metric in metrics" :key="metric.key" class="metric-item">
        <h4>{{ metric.name }}</h4>
        <div class="plot-container">
          <RoundScatterplot
            v-if="results[metric.key]"
            :data="results[metric.key]"
            :labels="labels"
            :size="180"
            :radius="selectedDataset === 'IRIS' ? 3 : 2"
          />
          <div v-else class="loading">
            <div class="spinner"></div>
            Calculating...
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.metric-comparison {
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  padding: 20px;
  background-color: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.controls {
  margin-bottom: 15px;
}

select {
  padding: 6px 12px;
  border-radius: 4px;
  border: 1px solid var(--vp-c-divider);
  background: var(--vp-c-bg);
}

.description {
  font-size: 0.9em;
  color: var(--vp-c-text-2);
  margin-bottom: 20px;
}

.grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
  gap: 15px;
}

.metric-item {
  text-align: center;
}

.metric-item h4 {
  margin: 0 0 10px 0;
  font-size: 1em;
}

.plot-container {
  background: var(--vp-c-bg);
  border-radius: 4px;
  padding: 5px;
  display: flex;
  justify-content: center;
  align-items: center;
  height: 190px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
}

.loading {
  font-size: 0.8em;
  font-style: italic;
  color: var(--vp-c-text-2);
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 5px;
}

.spinner {
  width: 16px;
  height: 16px;
  border: 2px solid var(--vp-c-divider);
  border-top-color: var(--vp-c-brand);
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  to {
    transform: rotate(360deg);
  }
}
</style>
