<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { runInWorker } from "./workerPool.js";

const loading = ref(true);
const results = ref(null);
const projection = ref(null);
const kFound = ref(0);

onMounted(async () => {
  try {
    // Generate a dataset with 4 or 5 blobs
    const dataset = await mistle.blobs({ N: 600, clusters: 4, seed: 42 });
    const data = dataset.values; // applyZScore(dataset.values);

    // Run XMeans to find the clusters
    // And UMAP to project for visualization
    const [clusterList, projectResult] = await Promise.all([
      runInWorker("Clustering", "XMeans", data, { K_min: 2, K_max: 10 }),
      runInWorker("DR", "UMAP", data, { n_neighbors: 15 }, 350),
    ]);

    results.value = clusterList;
    projection.value = projectResult;
    kFound.value = new Set(clusterList).size;
  } catch (e) {
    console.error("XMeans showcase error:", e);
  } finally {
    loading.value = false;
  }
});
</script>

<template>
  <div class="xmeans-showcase">
    <div v-if="loading" class="loading">
      <div class="spinner"></div>
      Discovering clusters automatically...
    </div>
    <div v-else class="content">
      <div class="plot-box">
        <RoundScatterplot :data="projection" :labels="results" :size="400" :radius="4" />
      </div>
      <div class="info">
        <h3>Found {{ kFound }} Clusters</h3>
        <p>
          X-Means analyzed the dataset and determined that <b>{{ kFound }}</b> is the optimal number
          of clusters based on the Bayesian Information Criterion (BIC).
        </p>
        <p class="hint">
          No manual "K" input was required. The algorithm recursively split clusters until the data
          fit was optimal.
        </p>
        <div class="stat-grid">
          <div class="stat-item">
            <span class="val">2</span>
            <span class="lab">Min K</span>
          </div>
          <div class="stat-item">
            <span class="val">10</span>
            <span class="lab">Max K</span>
          </div>
          <div class="stat-item highlight">
            <span class="val">{{ kFound }}</span>
            <span class="lab">Optimal K</span>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.xmeans-showcase {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 24px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.loading {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 15px;
  height: 300px;
  color: var(--vp-c-text-2);
}

.content {
  display: flex;
  flex-wrap: wrap;
  gap: 40px;
  justify-content: center;
  align-items: center;
}

.info {
  flex: 1;
  min-width: 250px;
}

.hint {
  font-size: 0.9em;
  color: var(--vp-c-text-2);
  font-style: italic;
  margin: 15px 0;
}

.plot-box {
  background: var(--vp-c-bg);
  border-radius: 8px;
  padding: 10px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  display: inline-block;
}

.stat-grid {
  display: flex;
  gap: 20px;
  margin-top: 25px;
}

.stat-item {
  display: flex;
  flex-direction: column;
  align-items: center;
  background: var(--vp-c-bg);
  padding: 10px 15px;
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
  min-width: 70px;
}

.stat-item.highlight {
  border-color: var(--vp-c-brand);
  background: var(--vp-c-brand-soft);
}

.stat-item .val {
  font-size: 1.4em;
  font-weight: bold;
  color: var(--vp-c-brand);
}

.stat-item .lab {
  font-size: 0.7em;
  text-transform: uppercase;
  color: var(--vp-c-text-2);
}

.spinner {
  width: 30px;
  height: 30px;
  border: 3px solid var(--vp-c-divider);
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
