<script setup>
import * as mistle from "@saehrimnir/mistle";
import { onMounted, ref } from "vue";
import { createPersistentInstance, usePersistentInstance } from "./workerPool.js";

const loading = ref(true);
const building = ref(false);
const dataset = ref(null);
const instance = ref(null);
const selectedIndex = ref(null);
const neighbors = ref([]);
const subsetIndices = ref([]);

const MNIST_SIZE = 3500;

onMounted(async () => {
  try {
    dataset.value = await mistle.fetch_mnist({ N: MNIST_SIZE });

    const counts = Array(10).fill(0);
    const indices = [];
    for (let i = 0; i < dataset.value.labels.length; i++) {
      const label = parseInt(dataset.value.labels[i]);
      if (!isNaN(label) && counts[label] < 15) {
        indices.push(i);
        counts[label]++;
      }
      if (indices.length === 150) break;
    }
    indices.sort((a, b) => parseInt(dataset.value.labels[a]) - parseInt(dataset.value.labels[b]));
    subsetIndices.value = indices;

    building.value = true;
    // Build HNSW index for the 784-dimensional vectors
    // Unwrap from Proxy to avoid DataCloneError
    const rawValues = JSON.parse(JSON.stringify(dataset.value.values));
    instance.value = await createPersistentInstance("HNSW", rawValues, {
      m: 16,
      ef_construction: 200,
    });
  } catch (e) {
    console.error("Error loading/building MNIST search:", e);
  } finally {
    building.value = false;
    loading.value = false;
  }
});

async function search(index) {
  if (!instance.value || building.value) return;
  selectedIndex.value = index;
  try {
    const query = JSON.parse(JSON.stringify(dataset.value.values[index]));
    const results = await usePersistentInstance(instance.value, "Search", {
      query,
      k: 14,
    });
    neighbors.value = results;
  } catch (e) {
    console.error("Search error:", e);
  }
}

// Helper to draw MNIST to canvas
function drawDigit(canvas, vector) {
  if (!canvas || !vector) return;
  const ctx = canvas.getContext("2d");
  const imageData = ctx.createImageData(28, 28);
  for (let i = 0; i < 784; i++) {
    const val = vector[i] * 255;
    const idx = i * 4;
    imageData.data[idx] = val;
    imageData.data[idx + 1] = val;
    imageData.data[idx + 2] = val;
    imageData.data[idx + 3] = 255;
  }
  ctx.putImageData(imageData, 0, 0);
}

// Custom directive for lazy drawing
const vDigit = {
  mounted(el, binding) {
    drawDigit(el, binding.value);
  },
  updated(el, binding) {
    drawDigit(el, binding.value);
  },
};
</script>

<template>
  <div class="image-search">
    <div v-if="loading" class="state">
      <div class="spinner"></div>
      Loading MNIST Dataset...
    </div>
    <div v-else-if="building" class="state">
      <div class="spinner"></div>
      Building HNSW Index (High-Dimensional Space)...
    </div>
    <div v-else class="search-layout">
      <div class="dataset-view">
        <h3>1. Select a Digit</h3>
        <p class="hint">Click an image to find similar ones in the 784-D space.</p>
        <div class="grid small-grid">
          <div
            v-for="i in subsetIndices"
            :key="i"
            class="digit-item"
            :class="{ active: selectedIndex === i }"
            @click="search(i)"
          >
            <canvas v-digit="dataset.values[i]" width="28" height="28"></canvas>
          </div>
        </div>
      </div>

      <div class="results-view">
        <div v-if="selectedIndex !== null">
          <h3>2. Nearest Neighbors (HNSW)</h3>
          <p class="hint">Top 14 matches found using approximate KNN.</p>
          <div class="grid large-grid">
            <div v-for="(neighbor, i) in neighbors" :key="i" class="neighbor-item">
              <canvas v-digit="neighbor.element" width="28" height="28"></canvas>
              <span class="dist">Dist: {{ neighbor.distance.toFixed(1) }}</span>
              <span class="label">True: {{ dataset.labels[neighbor.index] }}</span>
            </div>
          </div>
        </div>
        <div v-else class="empty-results">
          <p>Wait for results here...</p>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.image-search {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 24px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
  min-height: 600px;
}

.state {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 15px;
  height: 300px;
  color: var(--vp-c-text-2);
}

.search-layout {
  display: flex;
  flex-wrap: wrap;
  gap: 40px;
}

.dataset-view,
.results-view {
  flex: 1;
  min-width: 300px;
}

.grid {
  display: grid;
  gap: 8px;
}

.small-grid {
  grid-template-columns: repeat(auto-fill, minmax(32px, 1fr));
  max-height: 450px;
  overflow-y: auto;
  border: 1px solid var(--vp-c-divider);
  padding: 10px;
  border-radius: 8px;
  background: var(--vp-c-bg);
}

.large-grid {
  grid-template-columns: repeat(auto-fill, minmax(80px, 1fr));
}

.digit-item {
  cursor: pointer;
  padding: 2px;
  border-radius: 4px;
  transition: all 0.2s;
}

.digit-item:hover {
  background: var(--vp-c-brand);
}

.digit-item.active {
  background: var(--vp-c-brand);
  outline: 2px solid var(--vp-c-brand);
}

.neighbor-item {
  display: flex;
  flex-direction: column;
  align-items: center;
  background: var(--vp-c-bg);
  padding: 10px;
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
}

.neighbor-item canvas {
  width: 56px;
  height: 56px;
  image-rendering: pixelated;
}

.dist {
  font-size: 0.7em;
  color: var(--vp-c-text-2);
  margin-top: 5px;
}

.label {
  font-weight: bold;
  font-size: 0.8em;
}

.hint {
  font-size: 0.85em;
  color: var(--vp-c-text-2);
  margin-bottom: 15px;
}

.empty-results {
  height: 200px;
  display: flex;
  align-items: center;
  justify-content: center;
  border: 2px dashed var(--vp-c-divider);
  border-radius: 12px;
  color: var(--vp-c-text-3);
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
