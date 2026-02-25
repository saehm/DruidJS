<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const projection = ref(null);
const isLoading = ref(false);
const labels = mistle.IRIS.labels;

const runTopoMap = async () => {
  const data = applyZScore(mistle.IRIS.values);
  isLoading.value = true;
  try {
    const result = await runInWorker("DR", "TopoMap", data);
    projection.value = result;
  } catch (e) {
    console.error("Error in TopoMap worker:", e);
  } finally {
    isLoading.value = false;
  }
};

onMounted(() => {
  runTopoMap();
});
</script>

<template>
  <div class="topomap-showcase">
    <div class="plot-container">
      <RoundScatterplot
        v-if="projection"
        :data="projection"
        :labels="labels"
        :size="400"
        :radius="5"
      />
      <div v-if="isLoading" class="loading">
        <div class="spinner"></div>
        Calculating Topology...
      </div>
    </div>
    <div class="controls">
      <button @click="runTopoMap" class="btn-primary" :disabled="isLoading">Recalculate</button>
    </div>
  </div>
</template>

<style scoped>
.topomap-showcase {
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  padding: 20px;
  background-color: var(--vp-c-bg-soft);
  margin: 20px 0;
  display: flex;
  flex-direction: column;
  align-items: center;
}

.plot-container {
  display: flex;
  justify-content: center;
  align-items: center;
  background: var(--vp-c-bg);
  border-radius: 4px;
  padding: 10px;
  min-height: 420px;
  width: 100%;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
}

.loading {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 10px;
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

.controls {
  margin-top: 20px;
}

button {
  padding: 8px 16px;
  border-radius: 4px;
  background: var(--vp-c-brand);
  color: white;
  border: none;
  cursor: pointer;
}

button:disabled {
  opacity: 0.6;
}
</style>
