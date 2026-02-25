<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const loading = ref(true);
const results = ref({ tsne: null, topomap: null });
const labels = ref([]);

onMounted(async () => {
  try {
    // Use Moons dataset - it's disjoint and concave
    const dataset = await mistle.moons({ N: 600, open: 0 });
    labels.value = dataset.labels;
    const data = applyZScore(dataset.values);

    // Run both in parallel
    const [tsne, topomap] = await Promise.all([
      runInWorker("DR", "TSNE", data, { perplexity: 15 }, 1000),
      runInWorker("DR", "TopoMap", data, {}),
    ]);

    results.value.tsne = tsne;
    results.value.topomap = topomap;
  } catch (e) {
    console.error("Topology comparison error:", e);
  } finally {
    loading.value = false;
  }
});
</script>

<template>
  <div class="topology-comparison">
    <div v-if="loading" class="loading">
      <div class="spinner"></div>
      Calculating projections...
    </div>
    <div v-else class="content">
      <div class="side">
        <h3>t-SNE</h3>
        <p class="desc">Focuses on local similarities; often "shatters" complex manifolds.</p>
        <div class="plot-box">
          <RoundScatterplot :data="results.tsne" :labels="labels" :size="350" :radius="5" />
        </div>
      </div>
      <div class="side">
        <h3>TopoMap</h3>
        <p class="desc">Preserves 0-dimensional topology; maintains manifold continuity.</p>
        <div class="plot-box">
          <RoundScatterplot :data="results.topomap" :labels="labels" :size="350" :radius="5" />
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.topology-comparison {
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
  gap: 30px;
  justify-content: center;
}

.side {
  flex: 1;
  min-width: 250px;
  text-align: center;
}

.side h3 {
  margin-top: 0;
}

.desc {
  font-size: 0.85em;
  color: var(--vp-c-text-2);
  margin-bottom: 20px;
  min-height: 2.5em;
}

.plot-box {
  background: var(--vp-c-bg);
  border-radius: 8px;
  padding: 10px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  display: inline-block;
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
