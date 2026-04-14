<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const methods = [
  { name: "PCA", class: "PCA" },
  { name: "MDS", class: "MDS" },
  { name: "t-SNE", class: "TSNE", iterations: 1000 },
  { name: "UMAP", class: "UMAP", iterations: 350 },
  { name: "PaCMAP", class: "PaCMAP" },
  { name: "LocalMAP", class: "LocalMAP" },
  { name: "TriMap", class: "TriMap", iterations: 1000 },
  { name: "ISOMAP", class: "ISOMAP" },
  { name: "TopoMap", class: "TopoMap" },
  { name: "Sammon", class: "SAMMON", iterations: 1000 },
  { name: "LLE", class: "LLE" },
  { name: "FastMap", class: "FASTMAP" },
  { name: "SMACOF", class: "SMACOF", iterations: 300 },
  { name: "LDA", class: "LDA", needsLabels: true },
  { name: "LSP", class: "LSP" },
  { name: "LTSA", class: "LTSA" },
  { name: "SQDMDS", class: "SQDMDS", iterations: 1000 },
];

const results = ref({});
const labels = ref([]);

onMounted(async () => {
  const data = applyZScore(mistle.IRIS.values);
  labels.value = mistle.IRIS.labels;

  methods.forEach(async (method) => {
    try {
      const params = {};
      if (method.needsLabels) {
        params.labels = [...labels.value]; // Unwrap proxy to ensure serializability
      }
      const result = await runInWorker("DR", method.class, data, params, method.iterations);
      results.value[method.name] = result;
    } catch (e) {
      console.error(`Error calculating ${method.name} in worker:`, e);
    }
  });
});
</script>

<template>
  <div class="showcase-grid">
    <div v-for="method in methods" :key="method.name" class="showcase-item">
      <h3>{{ method.name }}</h3>
      <div class="plot-container">
        <RoundScatterplot
          v-if="results[method.name]"
          :data="results[method.name]"
          :labels="labels"
          :size="200"
          :radius="3"
        />
        <div v-else class="loading">
          <div class="spinner"></div>
          Calculating...
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.showcase-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
  gap: 20px;
  margin-top: 20px;
}

.showcase-item {
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  padding: 15px;
  text-align: center;
  background-color: var(--vp-c-bg-soft);
}

.showcase-item h3 {
  margin-top: 0;
  margin-bottom: 10px;
}

.plot-container {
  display: flex;
  justify-content: center;
  align-items: center;
  height: 200px;
  background: var(--vp-c-bg);
  border-radius: 8px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
}

.loading {
  color: var(--vp-c-text-2);
  font-style: italic;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 10px;
}

.spinner {
  width: 20px;
  height: 20px;
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
