<script setup>
import * as mistle from "@saehrimnir/mistle";
import { computed, onMounted, onUnmounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { applyZScore } from "./utils.js";
import { runStreamInWorker } from "./workerPool.js";

const datasetOptions = [
  { name: "Iris", key: "IRIS" },
  { name: "Swiss Roll", key: "SWISSROLL" },
  { name: "Moons", key: "MOONS" },
  { name: "Blobs", key: "BLOBS" },
  { name: "S-Shape", key: "SSHAPE" },
];

const selectedDataset = ref("IRIS");
const projection = ref(null);
const iteration = ref(0);
const isRunning = ref(false);
const selectedMethod = ref("TSNE");
const maxIterations = ref(200);

const methodParams = ref({
  TSNE: { perplexity: 30 },
  UMAP: { n_neighbors: 15, min_dist: 0.1 },
  PaCMAP: { n_neighbors: 10, MN_ratio: 0.5, FP_ratio: 2.0 },
  // low_dist_thres is an absolute embedding distance, and 10 -- the library default, tuned for
  // large datasets -- is too small for these, where it repels cluster-mates and breaks the
  // clusters apart. See docs/dimred/localmap.md, "Choosing low_dist_thres".
  LocalMAP: { n_neighbors: 10, MN_ratio: 0.5, FP_ratio: 2.0, low_dist_thres: 25 },
  TriMap: { n_inliers: 12, n_outliers: 4, n_random: 3 },
  SAMMON: { magic: 0.1 },
  SMACOF: {},
  StressMDS: { weights: -2 },
  SQDMDS: {},
});

/** Methods whose weight schedule is driven by `num_iters` rather than the iteration count. */
const phasedMethods = ["PaCMAP", "LocalMAP"];

/**
 * Methods whose `generator()` takes no argument and reads its budget from the `iterations`
 * parameter instead, so the slider has to be passed through rather than handed to the generator.
 */
const iterationParamMethods = ["SMACOF", "StressMDS"];

/**
 * PaCMAP splits its run into three phases and reads the boundaries from `num_iters`, not from how
 * many steps the caller asks for. Left at the default [100, 100, 250] the slider would only ever
 * animate the first two phases, so the phases are rescaled to it in the reference 100:100:250
 * ratio and the whole schedule stays visible at any length.
 */
const scaledPhases = (total) => {
  const p1 = Math.max(1, Math.round((total * 100) / 450));
  const p2 = Math.max(1, Math.round((total * 100) / 450));
  return [p1, p2, Math.max(1, total - p1 - p2)];
};

let currentStream = null;

const labels = computed(() => mistle[selectedDataset.value].labels);

/** Names the three exponents that coincide with a method the library already ships. */
const weightHint = computed(() => {
  const q = methodParams.value.StressMDS.weights;
  if (q === 0) return "raw stress, as SMACOF";
  if (q === -1) return "Sammon stress";
  if (q === -2) return "elastic / Kamada-Kawai";
  return "";
});

const startOptimization = async () => {
  if (isRunning.value) return;

  isRunning.value = true;
  iteration.value = 0;
  const data = mistle[selectedDataset.value].values;

  try {
    const rawData = JSON.parse(JSON.stringify(applyZScore(data)));
    const rawParams = JSON.parse(JSON.stringify(methodParams.value[selectedMethod.value] || {}));
    if (phasedMethods.includes(selectedMethod.value)) {
      rawParams.num_iters = scaledPhases(maxIterations.value);
    }
    if (iterationParamMethods.includes(selectedMethod.value)) {
      rawParams.iterations = maxIterations.value;
    }
    await runStreamInWorker(
      "AnimatedDR",
      selectedMethod.value,
      rawData,
      rawParams,
      maxIterations.value,
      (result, step) => {
        if (isRunning.value) {
          projection.value = result;
          iteration.value = step;
        }
      },
    );
  } catch (e) {
    console.error("Error in animated worker:", e);
  } finally {
    isRunning.value = false;
  }
};

const stopOptimization = () => {
  isRunning.value = false;
  // In a real implementation we'd probably want to tell the worker to stop
  // But for this showcase, just stopping the UI updates is okay for now,
  // or we could terminate the worker if we had a dedicated one.
};

const reset = () => {
  stopOptimization();
  projection.value = null;
  iteration.value = 0;
};

onMounted(() => {
  // We don't auto-start to avoid overwhelming on page load
});

onUnmounted(() => {
  stopOptimization();
});
</script>

<template>
  <div class="animated-projection">
    <div class="controls-panel">
      <h3 class="panel-title">Optimization Settings</h3>

      <div class="control-grid">
        <div class="control-item">
          <label>Dataset</label>
          <select
            v-model="selectedDataset"
            @change="reset"
            :disabled="isRunning"
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
            @change="reset"
            :disabled="isRunning"
            class="modern-select"
          >
            <option value="TSNE">t-SNE</option>
            <option value="UMAP">UMAP</option>
            <option value="PaCMAP">PaCMAP</option>
            <option value="LocalMAP">LocalMAP</option>
            <option value="TriMap">TriMap</option>
            <option value="SAMMON">Sammon</option>
            <option value="SMACOF">SMACOF</option>
            <option value="StressMDS">StressMDS</option>
            <option value="SQDMDS">SQDMDS</option>
          </select>
        </div>
      </div>

      <div class="param-section">
        <div class="control-item slider-item">
          <label
            >Iterations: <span class="val">{{ maxIterations }}</span></label
          >
          <input
            type="range"
            v-model.number="maxIterations"
            min="10"
            max="1000"
            step="10"
            :disabled="isRunning"
          />
        </div>

        <!-- Dependent parameters -->
        <div v-if="selectedMethod === 'TSNE'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Perplexity: <span class="val">{{ methodParams.TSNE.perplexity }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.TSNE.perplexity"
              min="5"
              max="100"
              step="1"
              :disabled="isRunning"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'UMAP'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Neighbors: <span class="val">{{ methodParams.UMAP.n_neighbors }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.UMAP.n_neighbors"
              min="2"
              max="100"
              step="1"
              :disabled="isRunning"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Min Dist: <span class="val">{{ methodParams.UMAP.min_dist.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.UMAP.min_dist"
              min="0"
              max="1"
              step="0.01"
              :disabled="isRunning"
            />
          </div>
        </div>

        <div
          v-if="selectedMethod === 'PaCMAP' || selectedMethod === 'LocalMAP'"
          class="param-group animate-fade"
        >
          <div class="control-item slider-item">
            <label
              >Neighbors:
              <span class="val">{{ methodParams[selectedMethod].n_neighbors }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams[selectedMethod].n_neighbors"
              min="2"
              max="50"
              step="1"
              :disabled="isRunning"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >MN Ratio:
              <span class="val">{{ methodParams[selectedMethod].MN_ratio.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams[selectedMethod].MN_ratio"
              min="0.1"
              max="2"
              step="0.1"
              :disabled="isRunning"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >FP Ratio:
              <span class="val">{{ methodParams[selectedMethod].FP_ratio.toFixed(1) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams[selectedMethod].FP_ratio"
              min="0.5"
              max="5"
              step="0.5"
              :disabled="isRunning"
            />
          </div>
          <div v-if="selectedMethod === 'LocalMAP'" class="control-item slider-item">
            <label
              >Low Dist Threshold:
              <span class="val">{{ methodParams.LocalMAP.low_dist_thres }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.LocalMAP.low_dist_thres"
              min="1"
              max="60"
              step="1"
              :disabled="isRunning"
            />
          </div>
          <p class="param-note">
            Phases scale to the iteration count: NN pairs attract, mid-near pairs pull the global
            layout together early, then fade out for local refinement.
          </p>
        </div>

        <div v-if="selectedMethod === 'TriMap'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Inliers: <span class="val">{{ methodParams.TriMap.n_inliers }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.TriMap.n_inliers"
              min="1"
              max="50"
              step="1"
              :disabled="isRunning"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Outliers: <span class="val">{{ methodParams.TriMap.n_outliers }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.TriMap.n_outliers"
              min="1"
              max="20"
              step="1"
              :disabled="isRunning"
            />
          </div>
          <div class="control-item slider-item">
            <label
              >Random: <span class="val">{{ methodParams.TriMap.n_random }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.TriMap.n_random"
              min="1"
              max="20"
              step="1"
              :disabled="isRunning"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'SAMMON'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Step Size (Magic):
              <span class="val">{{ methodParams.SAMMON.magic.toFixed(2) }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.SAMMON.magic"
              min="0.01"
              max="1.0"
              step="0.01"
              :disabled="isRunning"
            />
          </div>
        </div>

        <div v-if="selectedMethod === 'StressMDS'" class="param-group animate-fade">
          <div class="control-item slider-item">
            <label
              >Weight Exponent:
              <span class="val">{{ methodParams.StressMDS.weights.toFixed(1) }}</span>
              <span class="hint-inline">{{ weightHint }}</span></label
            >
            <input
              type="range"
              v-model.number="methodParams.StressMDS.weights"
              min="-3"
              max="0"
              step="0.5"
              :disabled="isRunning"
            />
          </div>
        </div>
      </div>

      <div class="action-footer">
        <div class="action-buttons">
          <button
            @click="isRunning ? stopOptimization() : startOptimization()"
            :class="['btn-action', isRunning ? 'btn-stop' : 'btn-start']"
          >
            <span class="btn-icon">{{ isRunning ? "⏹" : "▶" }}</span>
            {{ isRunning ? "Stop Animation" : "Start Animation" }}
          </button>
          <button @click="reset" class="btn-reset" :disabled="isRunning">Reset</button>
        </div>

        <div class="progress-container">
          <div class="progress-bar-bg">
            <div
              class="progress-bar-fill"
              :style="{ width: `${(iteration / maxIterations) * 100}%` }"
            ></div>
          </div>
          <span class="iteration-text">Iteration {{ iteration }} / {{ maxIterations }}</span>
        </div>
      </div>
    </div>

    <div class="plot-container">
      <RoundScatterplot
        v-if="projection"
        :data="projection"
        :labels="labels"
        :size="480"
        :radius="selectedDataset === 'IRIS' ? 5 : 3.5"
      />
      <div v-else class="placeholder">
        <div class="empty-state-icon">✨</div>
        <p>
          Select your parameters and click <strong>Start Animation</strong> to watch the structural
          manifold unfold dynamically.
        </p>
      </div>
    </div>
  </div>
</template>

<style scoped>
.animated-projection {
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
  .animated-projection {
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
  grid-template-columns: 1fr 1fr;
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

.control-item label .hint-inline {
  color: var(--vp-c-text-3);
  font-weight: 400;
  font-size: 0.9em;
  margin-left: 6px;
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

.animate-fade {
  animation: fadeIn 0.3s ease-out;
}

.param-note {
  margin: 0;
  font-size: 0.78em;
  line-height: 1.5;
  color: var(--vp-c-text-3);
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
  padding: 10px 16px;
  border-radius: 8px;
  font-weight: 600;
  font-size: 0.9em;
  cursor: pointer;
  transition: all 0.2s cubic-bezier(0.25, 0.8, 0.25, 1);
  border: 1px solid transparent;
}

button:disabled {
  opacity: 0.5;
  cursor: not-allowed;
  transform: none !important;
  box-shadow: none !important;
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
  box-shadow: 0 6px 16px rgba(var(--vp-c-brand-rgb), 0.4);
}

.btn-stop {
  background: linear-gradient(
    135deg,
    var(--vp-c-danger-1, #e74c3c) 0%,
    var(--vp-c-danger-2, #c0392b) 100%
  );
  color: #ffffff !important;
  box-shadow: 0 4px 12px rgba(231, 76, 60, 0.3);
}

.btn-stop:hover {
  transform: translateY(-2px);
  box-shadow: 0 6px 16px rgba(231, 76, 60, 0.4);
}

.btn-reset {
  background: var(--vp-c-bg-soft);
  color: var(--vp-c-text-1);
  border: 1px solid var(--vp-c-divider);
}

.btn-reset:hover:not(:disabled) {
  background: var(--vp-c-default-soft);
  border-color: var(--vp-c-text-3);
}

.progress-container {
  display: flex;
  flex-direction: column;
  gap: 6px;
}

.progress-bar-bg {
  height: 6px;
  background: var(--vp-c-divider);
  border-radius: 4px;
  overflow: hidden;
}

.progress-bar-fill {
  height: 100%;
  background: var(--vp-c-brand);
  border-radius: 4px;
  transition: width 0.1s linear;
}

.iteration-text {
  font-size: 0.8em;
  color: var(--vp-c-text-3);
  text-align: right;
  font-family: var(--vp-font-family-mono);
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

.placeholder {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  text-align: center;
  color: var(--vp-c-text-2);
  max-width: 300px;
  padding: 40px;
  border: 2px dashed var(--vp-c-divider);
  border-radius: 12px;
  background: rgba(var(--vp-c-bg-rgb), 0.5);
}

.empty-state-icon {
  font-size: 2.5em;
  margin-bottom: 12px;
  opacity: 0.8;
}

.placeholder strong {
  color: var(--vp-c-brand);
}
</style>
