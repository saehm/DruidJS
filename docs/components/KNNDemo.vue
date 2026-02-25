<script setup>
import { computed, onMounted, ref, watch } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { createPersistentInstance, usePersistentInstance } from "./workerPool.js";

const selectedMethod = ref("HNSW");
const numPoints = ref(500);
const queryK = ref(15);
const methodParams = ref({
  HNSW: { m: 16, ef_construction: 200 },
  NNDescent: { samples: 50, rho: 0.5, prng_seed: 1212 },
  KDTree: {},
  BallTree: {},
  Annoy: { numTrees: 10, maxPointsPerLeaf: 10 },
  LSH: { numHashFunctions: 5, numHashTables: 5 },
  NaiveKNN: {},
});

const labels = ref([]);
const projection = ref(null);
const isLoading = ref(false);
const buildTime = ref(0);
const searchTime = ref(0);

let currentInstance = null;
let currentRawData = null;
let isQuerying = false;
let lastPointerPos = null;

const pointColors = computed(() => {
  if (!labels.value.length) return [];
  return labels.value.map((l) => (l === "Neighbor" ? "#3b82f6" : "#cbd5e1"));
});

function generateUniformData(n) {
  const data = [];
  for (let i = 0; i < n; i++) {
    data.push([Math.random(), Math.random()]);
  }
  return data;
}

const runPipeline = async () => {
  isLoading.value = true;
  labels.value = [];

  try {
    currentRawData = generateUniformData(numPoints.value);

    // For uniform data, the projection is just the data itself!
    projection.value = currentRawData;
    labels.value = new Array(currentRawData.length).fill("Other");

    const rawParams = JSON.parse(JSON.stringify(methodParams.value[selectedMethod.value] || {}));

    // We don't terminate workers here because they are managed by a shared pool.
    // Terminating them would break subsequent calls to other components!
    currentInstance = null;

    const t0 = performance.now();
    currentInstance = await createPersistentInstance(
      selectedMethod.value,
      currentRawData,
      rawParams,
    );
    const t1 = performance.now();
    buildTime.value = (t1 - t0).toFixed(1);

    if (lastPointerPos) {
      await queryPoint(lastPointerPos);
    }
  } catch (e) {
    console.error("Error in KNN pipeline worker:", e);
  } finally {
    isLoading.value = false;
  }
};

const handlePointerMove = async (domainCoords) => {
  lastPointerPos = domainCoords;
  if (isQuerying || !currentInstance || !currentRawData) return;
  await queryPoint(domainCoords);
};

const queryPoint = async (query) => {
  if (!currentInstance) return;
  const instanceAtStart = currentInstance;
  isQuerying = true;
  try {
    const t0 = performance.now();
    const result = await usePersistentInstance(currentInstance, "Search", {
      query,
      k: queryK.value,
    });

    if (instanceAtStart !== currentInstance) return;

    const t1 = performance.now();
    searchTime.value = (t1 - t0).toFixed(2);

    const neighborIndices = result.map((n) =>
      n.index !== undefined ? n.index : n.i !== undefined ? n.i : n,
    );

    const newLabels = new Array(currentRawData.length).fill("Other");
    neighborIndices.forEach((idx) => {
      if (idx >= 0 && idx < newLabels.length) {
        newLabels[idx] = "Neighbor";
      }
    });

    labels.value = newLabels;
  } catch (e) {
    console.warn("KNN Query failed (likely instance changed):", e);
  } finally {
    isQuerying = false;
  }
};

// Debounce updates when slider or select changes
let pipelineTimeout = null;
watch(
  [numPoints, selectedMethod, methodParams],
  () => {
    if (pipelineTimeout) clearTimeout(pipelineTimeout);
    pipelineTimeout = setTimeout(runPipeline, 150);
  },
  { deep: true },
);

watch(queryK, () => {
  if (lastPointerPos) queryPoint(lastPointerPos);
});

onMounted(() => {
  runPipeline();
});
</script>

<template>
  <div class="knn-pipeline">
    <div class="controls-panel">
      <h3 class="panel-title">Pipeline Settings</h3>

      <div class="control-grid">
        <div class="control-item slider-item">
          <label
            >Uniform Points (N): <span class="val">{{ numPoints }}</span></label
          >
          <input type="range" v-model.number="numPoints" min="100" max="5000" step="100" />
        </div>
        <div class="control-item">
          <label>Algorithm</label>
          <select v-model="selectedMethod" class="modern-select">
            <option value="HNSW">HNSW</option>
            <option value="NNDescent">NN-Descent</option>
            <option value="Annoy">Annoy</option>
            <option value="KDTree">KD-Tree</option>
            <option value="BallTree">Ball Tree</option>
            <option value="LSH">Locality Sensitive Hashing (LSH)</option>
            <option value="NaiveKNN">Naive KNN (Exact)</option>
          </select>
        </div>
      </div>

      <div class="param-section">
        <div class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Neighbors (K): <span class="val">{{ queryK }}</span></label
            >
            <input type="range" v-model.number="queryK" min="1" max="100" step="1" />
          </div>
        </div>

        <!-- Algorithm Specific Params -->
        <div v-if="selectedMethod === 'HNSW'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >M (Links per node): <span class="val">{{ methodParams.HNSW.m }}</span></label
            >
            <input type="range" v-model.number="methodParams.HNSW.m" min="4" max="64" step="4" />
          </div>
          <div class="control-item slider-item">
            <label
              >EF Construction:
              <span class="val">{{ methodParams.HNSW.ef_construction }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.HNSW.ef_construction"
              min="10"
              max="500"
              step="10"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'NNDescent'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Samples: <span class="val">{{ methodParams.NNDescent.samples }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.NNDescent.samples"
              min="10"
              max="100"
              step="10"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Rho: <span class="val">{{ methodParams.NNDescent.rho.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.NNDescent.rho"
              min="0.1"
              max="1.0"
              step="0.1"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'Annoy'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Trees: <span class="val">{{ methodParams.Annoy.numTrees }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.Annoy.numTrees"
              min="1"
              max="50"
              step="1"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Max Pts / Leaf:
              <span class="val">{{ methodParams.Annoy.maxPointsPerLeaf }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.Annoy.maxPointsPerLeaf"
              min="2"
              max="50"
              step="1"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'LSH'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Hash Functions:
              <span class="val">{{ methodParams.LSH.numHashFunctions }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.LSH.numHashFunctions"
              min="2"
              max="16"
              step="1"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Hash Tables: <span class="val">{{ methodParams.LSH.numHashTables }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.LSH.numHashTables"
              min="1"
              max="20"
              step="1"
            />
          </div>
        </div>

        <div class="description animate-fade">
          <p>
            Index built in <strong>{{ buildTime }} ms</strong>.
          </p>
          <p>
            Hover query fetched in <strong>{{ searchTime }} ms</strong>.
          </p>
        </div>
      </div>
    </div>

    <div class="plot-container">
      <RoundScatterplot
        v-if="projection && labels.length"
        :data="projection"
        :labels="labels"
        :colors="pointColors"
        :size="480"
        :radius="5"
        @pointer-move="handlePointerMove"
      />
      <div v-if="isLoading" class="loading-overlay">
        <div class="spinner-large"></div>
        <span>Building Graph...</span>
      </div>
    </div>
  </div>
</template>

<style scoped>
.knn-pipeline {
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
  .knn-pipeline {
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
</style>
