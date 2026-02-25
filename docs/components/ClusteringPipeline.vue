<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const datasetOptions = [
  { name: "Iris", key: "IRIS" },
  { name: "Wine", key: "WINE" },
  { name: "Penguins", key: "PENGUINS" },
  { name: "Blobs", key: "BLOBS" },
];

const selectedDataset = ref("IRIS");
const selectedMethod = ref("KMeans");
const methodParams = ref({
  KMeans: { K: 3 },
  KMedoids: { K: 3 },
  XMeans: { K_min: 2, K_max: 10, min_cluster_size: 25 },
  CURE: { K: 3, num_representatives: 5, shrink_factor: 0.5 },
  MeanShift: { bandwidth: 0.5 },
  OPTICS: { epsilon: 1.0, min_points: 4 },
  HierarchicalClustering: { linkage: "complete", cut_value: 5.0 },
});

const clusters = ref([]);
const projection = ref(null);
const isLoading = ref(false);

const runPipeline = async () => {
  const data = mistle[selectedDataset.value].values;
  isLoading.value = true;
  projection.value = null;
  clusters.value = [];

  try {
    // 1. Clustering in worker
    // Unwrap from Proxy to avoid DataCloneError
    const rawData = JSON.parse(JSON.stringify(applyZScore(data)));
    const rawParams = JSON.parse(JSON.stringify(methodParams.value[selectedMethod.value] || {}));

    // Some backend endpoints still map backwards-compatible K
    if (rawParams.K) rawParams.k = rawParams.K;

    const clusteringResult = await runInWorker(
      "Clustering",
      selectedMethod.value,
      rawData,
      rawParams,
    );
    clusters.value = clusteringResult;

    // 2. Projection in worker
    const projectionResult = await runInWorker("DR", "UMAP", rawData);
    projection.value = projectionResult;
  } catch (e) {
    console.error("Error in clustering pipeline worker:", e);
  } finally {
    isLoading.value = false;
  }
};

onMounted(() => {
  runPipeline();
});
</script>

<template>
  <div class="clustering-pipeline">
    <div class="controls-panel">
      <h3 class="panel-title">Pipeline Settings</h3>

      <div class="control-grid">
        <div class="control-item">
          <label>Dataset</label>
          <select
            v-model="selectedDataset"
            @change="runPipeline"
            :disabled="isLoading"
            class="modern-select"
          >
            <option v-for="opt in datasetOptions" :key="opt.key" :value="opt.key">
              {{ opt.name }}
            </option>
          </select>
        </div>
        <div class="control-item">
          <label>Algorithm</label>
          <select
            v-model="selectedMethod"
            @change="runPipeline"
            :disabled="isLoading"
            class="modern-select"
          >
            <option value="KMeans">K-Means</option>
            <option value="KMedoids">K-Medoids</option>
            <option value="XMeans">X-Means</option>
            <option value="CURE">CURE</option>
            <option value="MeanShift">Mean Shift</option>
            <option value="OPTICS">OPTICS</option>
            <option value="HierarchicalClustering">Hierarchical Clustering</option>
          </select>
        </div>
      </div>

      <div class="param-section">
        <div
          v-if="['KMeans', 'KMedoids'].includes(selectedMethod)"
          class="param-group animate-fade"
        >
          <div class="control-item slider-item">
            <label
              >Clusters (K): <span class="val">{{ methodParams[selectedMethod].K }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams[selectedMethod].K"
              min="2"
              max="15"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'XMeans'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Min Clusters: <span class="val">{{ methodParams.XMeans.K_min }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.XMeans.K_min"
              min="2"
              max="10"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Max Clusters: <span class="val">{{ methodParams.XMeans.K_max }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.XMeans.K_max"
              min="2"
              max="25"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Min Cluster Size:
              <span class="val">{{ methodParams.XMeans.min_cluster_size }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.XMeans.min_cluster_size"
              min="1"
              max="50"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'CURE'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Clusters (K): <span class="val">{{ methodParams.CURE.K }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.CURE.K"
              min="2"
              max="15"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Representatives:
              <span class="val">{{ methodParams.CURE.num_representatives }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.CURE.num_representatives"
              min="1"
              max="15"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Shrink Factor:
              <span class="val">{{ methodParams.CURE.shrink_factor.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.CURE.shrink_factor"
              min="0"
              max="1"
              step="0.05"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'MeanShift'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Bandwidth:
              <span class="val">{{ methodParams.MeanShift.bandwidth.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.MeanShift.bandwidth"
              min="0.1"
              max="5.0"
              step="0.1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'OPTICS'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Epsilon: <span class="val">{{ methodParams.OPTICS.epsilon.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.OPTICS.epsilon"
              min="0.1"
              max="10"
              step="0.1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Min Points: <span class="val">{{ methodParams.OPTICS.min_points }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.OPTICS.min_points"
              min="2"
              max="25"
              step="1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'HierarchicalClustering'" class="param-group animate-fade">
          <div class="control-item">
            <label>Linkage Strategy</label>
            <select
              v-model="methodParams.HierarchicalClustering.linkage"
              @change="runPipeline"
              :disabled="isLoading"
              class="modern-select"
              style="margin-top: 5px"
            >
              <option value="single">Single Link</option>
              <option value="average">Average Link</option>
              <option value="complete">Complete Link</option>
            </select>
          </div>
          <div class="control-item slider-item" style="margin-top: 10px">
            <label
              >Dendrogram Tree Cut (Distance):
              <span class="val">{{
                methodParams.HierarchicalClustering.cut_value.toFixed(2)
              }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.HierarchicalClustering.cut_value"
              min="0.1"
              max="30"
              step="0.1"
              @change="runPipeline"
              :disabled="isLoading"
            />
          </div>
        </div>

        <div class="description animate-fade">
          <p>
            This example runs <strong>{{ selectedMethod }}</strong> on the original high-dimensional
            data, assigning cluster labels locally, and then uses <strong>UMAP</strong> to visualize
            the manifold in 2D.
          </p>
        </div>
      </div>

      <div class="action-footer">
        <div class="action-buttons">
          <button @click="runPipeline" class="btn-action btn-start" :disabled="isLoading">
            <span class="btn-icon" v-if="!isLoading">⚡</span>
            <span class="btn-icon spinner-icon" v-else></span>
            {{ isLoading ? "Computing..." : "Re-run Pipeline" }}
          </button>
        </div>
      </div>
    </div>

    <div class="plot-container">
      <RoundScatterplot
        v-if="projection && clusters.length && !isLoading"
        :data="projection"
        :labels="clusters"
        :size="480"
        :radius="selectedDataset === 'IRIS' ? 5 : 3.5"
      />
      <div v-if="isLoading" class="loading-overlay">
        <div class="spinner-large"></div>
        <span>Running Pipeline...</span>
      </div>
    </div>
  </div>
</template>

<style scoped>
.clustering-pipeline {
  display: grid;
  grid-template-columns: 320px 1fr;
  gap: 20px;
  background-color: var(--vp-c-bg-soft);
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  overflow: hidden;
  margin: 24px 0;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.05);
  font-family: var(--vp-font-family-base);
}

@media (max-width: 900px) {
  .clustering-pipeline {
    grid-template-columns: 1fr;
    grid-template-rows: auto auto;
  }
}

.controls-panel {
  padding: 24px;
  background: var(--vp-c-bg);
  border-right: 1px solid var(--vp-c-divider);
  display: flex;
  flex-direction: column;
}

.panel-title {
  margin: 0 0 20px 0;
  font-size: 1.1em;
  font-weight: 600;
  color: var(--vp-c-text-1);
  border-bottom: 2px solid var(--vp-c-brand-soft);
  padding-bottom: 8px;
  display: inline-block;
}

.control-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: 12px;
  margin-bottom: 20px;
}

.control-item {
  display: flex;
  flex-direction: column;
  gap: 6px;
}

.control-item label {
  font-size: 0.85em;
  font-weight: 600;
  letter-spacing: 0.5px;
  color: var(--vp-c-text-2);
  text-transform: uppercase;
  display: flex;
  justify-content: space-between;
}

.control-item label .val {
  color: var(--vp-c-brand);
  font-family: var(--vp-font-family-mono);
  font-weight: bold;
}

.modern-select {
  padding: 8px 12px;
  border-radius: 6px;
  border: 1px solid var(--vp-c-divider);
  background-color: var(--vp-c-bg-soft);
  color: var(--vp-c-text-1);
  font-size: 0.9em;
  cursor: pointer;
  transition: all 0.2s ease;
  appearance: none;
  background-image: url("data:image/svg+xml;charset=UTF-8,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3e%3cpolyline points='6 9 12 15 18 9'%3e%3c/polyline%3e%3c/svg%3e");
  background-repeat: no-repeat;
  background-position: right 10px center;
  background-size: 14px;
}

.modern-select:hover:not(:disabled) {
  border-color: var(--vp-c-brand);
}

.param-section {
  display: flex;
  flex-direction: column;
  gap: 16px;
  margin-bottom: 24px;
  flex-grow: 1;
}

.param-group {
  display: flex;
  flex-direction: column;
  gap: 16px;
  padding: 12px;
  background: var(--vp-c-bg-soft);
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
}

.description {
  font-size: 0.9em;
  line-height: 1.5;
  color: var(--vp-c-text-2);
  padding: 0 4px;
}

.description strong {
  color: var(--vp-c-brand);
}

.animate-fade {
  animation: fadeIn 0.3s ease-out;
}

@keyframes fadeIn {
  from {
    opacity: 0;
    transform: translateY(-5px);
  }

  to {
    opacity: 1;
    transform: translateY(0);
  }
}

.slider-item input[type="range"] {
  width: 100%;
  cursor: pointer;
  accent-color: var(--vp-c-brand);
}

.action-footer {
  margin-top: auto;
  display: flex;
  flex-direction: column;
  gap: 16px;
}

.action-buttons {
  display: flex;
  gap: 10px;
}

button {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 8px;
  padding: 12px 16px;
  border-radius: 8px;
  font-weight: 600;
  font-size: 0.95em;
  cursor: pointer;
  transition: all 0.2s cubic-bezier(0.25, 0.8, 0.25, 1);
  border: 1px solid transparent;
}

button:disabled {
  opacity: 0.7;
  cursor: not-allowed;
  transform: none !important;
  box-shadow: none !important;
  filter: grayscale(0.5);
}

.btn-action {
  flex-grow: 1;
  color: white;
}

.btn-start {
  background: linear-gradient(
    135deg,
    var(--vp-c-brand-1, var(--vp-c-brand)) 0%,
    var(--vp-c-brand-3, var(--vp-c-brand)) 100%
  );
  color: #ffffff !important;
  box-shadow: 0 4px 12px rgba(var(--vp-c-brand-rgb, 0, 0, 0), 0.3);
}

.btn-start:hover:not(:disabled) {
  transform: translateY(-2px);
  box-shadow: 0 6px 16px rgba(var(--vp-c-brand-rgb, 0, 0, 0), 0.4);
}

.plot-container {
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 20px;
  min-height: 480px;
  position: relative;
  overflow: hidden;
  width: 100%;
}

.loading-overlay {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 16px;
  color: var(--vp-c-text-2);
  width: 100%;
  height: 100%;
  position: absolute;
  top: 0;
  left: 0;
  background: rgba(var(--vp-c-bg-rgb), 0.7);
  backdrop-filter: blur(4px);
  border-radius: 12px;
}

.spinner-large {
  width: 40px;
  height: 40px;
  border: 4px solid var(--vp-c-divider);
  border-top-color: var(--vp-c-brand);
  border-radius: 50%;
  animation: spin 1s cubic-bezier(0.25, 0.8, 0.25, 1) infinite;
}

.spinner-icon {
  display: inline-block;
  width: 14px;
  height: 14px;
  border: 2px solid rgba(255, 255, 255, 0.3);
  border-top-color: #fff;
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  to {
    transform: rotate(360deg);
  }
}
</style>
