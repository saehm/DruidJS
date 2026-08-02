<script setup>
import { interpolateViridis } from "d3";
import * as mistle from "@saehrimnir/mistle";
import { computed, onMounted, ref } from "vue";
import TreePlot from "./TreePlot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const loading = ref(true);
const error = ref(null);

const projection = ref([]);
const edges = ref([]);
const labels = ref([]);
const truth = ref([]);
const curvature = ref([]);
const beta = ref(0);

const colorBy = ref("cluster");
const showBridges = ref(true);
const selected = ref([]);

const DATASETS = {
  iris: {
    title: "Iris",
    clusters: 3,
    load: () => ({ values: mistle.IRIS.values, labels: mistle.IRIS.labels }),
  },
  blobs: {
    title: "Blobs",
    clusters: 5,
    load: async () => {
      // `centers`, not `clusters` — mistle.blobs ignores an unknown key and would quietly hand back
      // its default of 3 groups while MINFOTree was asked to find 5.
      const d = await mistle.blobs({ N: 400, centers: 5, seed: 42 });
      return { values: d.values, labels: d.labels };
    },
  },
  moons: {
    title: "Moons",
    clusters: 2,
    load: async () => {
      const d = await mistle.moons({ N: 400, open: 0 });
      return { values: d.values, labels: d.labels };
    },
  },
};
const dataset = ref("iris");

async function run() {
  loading.value = true;
  error.value = null;
  selected.value = [];
  try {
    const spec = DATASETS[dataset.value];
    const loaded = await spec.load();
    truth.value = loaded.labels;
    const data = applyZScore(loaded.values);

    const result = await runInWorker("Tree", "MINFOTree", data, { clusters: spec.clusters });
    projection.value = result.projection;
    edges.value = result.edges;
    labels.value = result.labels;
    curvature.value = result.curvature;
    beta.value = result.beta;
  } catch (e) {
    error.value = e.message;
  } finally {
    loading.value = false;
  }
}

onMounted(run);

/** Adjacency list over the tree, rebuilt whenever the tree changes. */
const adjacency = computed(() => {
  const adj = Array.from({ length: projection.value.length }, () => []);
  for (const [u, v] of edges.value) {
    adj[u].push(v);
    adj[v].push(u);
  }
  return adj;
});

/**
 * The unique path between the two selected nodes.
 *
 * This is the part a scatterplot cannot do. In a tree there is exactly one route between any two
 * points, so the path is not a choice among many — it is *the* sequence of intermediate observations
 * the data takes to get from one to the other.
 */
const path = computed(() => {
  const [a, b] = selected.value;
  if (a === undefined || b === undefined || a === b) return [];
  const previous = new Map([[a, -1]]);
  const queue = [a];
  while (queue.length) {
    const u = queue.shift();
    if (u === b) break;
    for (const v of adjacency.value[u]) {
      if (previous.has(v)) continue;
      previous.set(v, u);
      queue.push(v);
    }
  }
  if (!previous.has(b)) return [];
  const route = [];
  for (let at = b; at !== -1; at = previous.get(at)) route.push(at);
  return route.reverse();
});

function onNodeClick(i) {
  const current = selected.value;
  if (current.length >= 2) selected.value = [i];
  else if (current[0] === i) selected.value = [];
  else selected.value = [...current, i];
}

const bridgeCount = computed(
  () => edges.value.filter(([u, v]) => labels.value[u] !== labels.value[v]).length,
);

/** How many clusters the path crosses — its "transition count". */
const pathCrossings = computed(() => {
  let crossings = 0;
  for (let i = 0; i + 1 < path.value.length; ++i) {
    if (labels.value[path.value[i]] !== labels.value[path.value[i + 1]]) crossings++;
  }
  return crossings;
});

const values = computed(() => (colorBy.value === "curvature" ? curvature.value : undefined));

/** Viridis, sampled for the legend bar — the same scale {@link TreePlot} colours nodes with. */
const curvatureGradient = computed(() => {
  const stops = Array.from({ length: 9 }, (_, i) => interpolateViridis(i / 8));
  return `linear-gradient(to right, ${stops.join(", ")})`;
});

const curvatureRange = computed(() => {
  if (!curvature.value.length) return { lo: 0, hi: 1 };
  return { lo: Math.min(...curvature.value), hi: Math.max(...curvature.value) };
});

/**
 * Which end of the scale the interior points sit at, measured rather than assumed.
 *
 * The direction is not fixed: it depends on the estimated β, and hard-coding "high means interior"
 * would silently mislabel the legend on a dataset in the other regime. So the two groups are
 * compared on the tree that is actually on screen.
 */
const curvatureEnds = computed(() => {
  const c = curvature.value;
  const l = labels.value;
  if (!c.length || !edges.value.length) return null;

  /** @type {Set<number>[]} */
  const neighborhood = Array.from({ length: c.length }, () => new Set());
  for (const [u, v] of edges.value) {
    neighborhood[u].add(v);
    neighborhood[v].add(u);
  }
  let interior = 0;
  let interiorN = 0;
  let boundary = 0;
  let boundaryN = 0;
  for (let i = 0; i < c.length; ++i) {
    const mixes = [...neighborhood[i]].some((j) => l[j] !== l[i]);
    if (mixes) {
      boundary += c[i];
      boundaryN++;
    } else {
      interior += c[i];
      interiorN++;
    }
  }
  if (!interiorN || !boundaryN) return null;

  const { lo, hi } = curvatureRange.value;
  const span = hi - lo || 1;
  const at = (v) => `${Math.min(Math.max((v - lo) / span, 0), 1) * 100}%`;
  const interiorMean = interior / interiorN;
  const boundaryMean = boundary / boundaryN;

  return {
    interiorHigh: interiorMean > boundaryMean,
    // Marks on the bar rather than a min/max row: the curvature is normalised, so its endpoints read
    // 0.001 and 1.001 on every dataset and say nothing. Where the two populations *sit* is the part
    // that differs — far apart on well-separated blobs, nearly on top of each other on Iris.
    interior: { mean: interiorMean, left: at(interiorMean), n: interiorN },
    boundary: { mean: boundaryMean, left: at(boundaryMean), n: boundaryN },
  };
});
</script>

<template>
  <div class="minfo-showcase">
    <div class="controls">
      <label>
        Dataset
        <select v-model="dataset" @change="run">
          <option v-for="(spec, key) in DATASETS" :key="key" :value="key">{{ spec.title }}</option>
        </select>
      </label>
      <label>
        Colour by
        <select v-model="colorBy">
          <option value="cluster">Cluster label</option>
          <option value="curvature">Information curvature</option>
        </select>
      </label>
      <label class="check">
        <input type="checkbox" v-model="showBridges" />
        Highlight inter-cluster edges
      </label>
    </div>

    <div v-if="loading" class="state">
      <div class="spinner"></div>
      Building the information graph…
    </div>
    <div v-else-if="error" class="state error">{{ error }}</div>
    <div v-else class="content">
      <div class="plot">
        <TreePlot
          :data="projection"
          :edges="edges"
          :labels="labels"
          :values="values"
          :path="path"
          :selected="selected"
          :show-bridges="showBridges"
          :size="520"
          :radius="4"
          @node-click="onNodeClick"
        />

        <figcaption v-if="colorBy === 'curvature'" class="legend">
          <div class="legend-head">
            <span>Information curvature</span>
            <span class="legend-beta">β = {{ beta.toFixed(2) }}</span>
          </div>
          <div class="legend-bar" :style="{ background: curvatureGradient }">
            <template v-if="curvatureEnds">
              <span
                class="legend-mark boundary"
                :style="{ left: curvatureEnds.boundary.left }"
                :title="`mean curvature of boundary points: ${curvatureEnds.boundary.mean.toFixed(3)}`"
              ></span>
              <span
                class="legend-mark interior"
                :style="{ left: curvatureEnds.interior.left }"
                :title="`mean curvature of interior points: ${curvatureEnds.interior.mean.toFixed(3)}`"
              ></span>
            </template>
          </div>
          <div v-if="curvatureEnds" class="legend-ends">
            <span>{{ curvatureEnds.interiorHigh ? "boundary" : "interior" }}</span>
            <span>{{ curvatureEnds.interiorHigh ? "interior" : "boundary" }}</span>
          </div>
          <div v-if="curvatureEnds" class="legend-means">
            <span>
              <b class="swatch boundary"></b>
              boundary mean {{ curvatureEnds.boundary.mean.toFixed(2) }}
              <span class="muted">({{ curvatureEnds.boundary.n }} pts)</span>
            </span>
            <span>
              <b class="swatch interior"></b>
              interior mean {{ curvatureEnds.interior.mean.toFixed(2) }}
              <span class="muted">({{ curvatureEnds.interior.n }} pts)</span>
            </span>
          </div>
          <p class="legend-note">
            Derived from the Potts model: a point scores by how far its own label sits from what its
            neighbours predict, against how stable that neighbourhood is. Which end means
            <em>interior</em> depends on β and is measured here rather than assumed —
            <!-- Relative so it stays correct regardless of the site's `base`. -->
            <a href="../dimred/minfotree#a-note-on-the-curvature-ordering">why</a>.
          </p>
        </figcaption>
      </div>

      <aside class="panel">
        <h4>The tree</h4>
        <dl>
          <dt>Points</dt>
          <dd>{{ projection.length }}</dd>
          <dt>Tree edges</dt>
          <dd>{{ edges.length }}</dd>
          <dt>Inter-cluster edges</dt>
          <dd>{{ bridgeCount }} <span class="muted">(min {{ new Set(labels).size - 1 }})</span></dd>
          <dt>Potts β</dt>
          <dd>{{ beta.toFixed(3) }}</dd>
        </dl>

        <h4>Unique path</h4>
        <p v-if="selected.length < 2" class="hint">
          Click two points. A tree has exactly one route between them, so the path is the sequence of
          observations the data passes through — something a scatterplot cannot tell you.
        </p>
        <dl v-else>
          <dt>Length</dt>
          <dd>{{ path.length }} points, {{ Math.max(path.length - 1, 0) }} hops</dd>
          <dt>Cluster changes</dt>
          <dd>{{ pathCrossings }}</dd>
          <dt>Clusters visited</dt>
          <dd>{{ new Set(path.map((i) => labels[i])).size }}</dd>
        </dl>
        <button v-if="selected.length" class="clear" @click="selected = []">Clear selection</button>
      </aside>
    </div>
  </div>
</template>

<style scoped>
.minfo-showcase {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 20px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.controls {
  display: flex;
  flex-wrap: wrap;
  gap: 18px;
  align-items: center;
  margin-bottom: 18px;
  font-size: 0.85rem;
}

.controls label {
  display: flex;
  align-items: center;
  gap: 8px;
  color: var(--vp-c-text-2);
}

.controls select {
  background: var(--vp-c-bg);
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
  padding: 4px 8px;
  color: var(--vp-c-text-1);
}

.controls .check {
  gap: 6px;
}

.content {
  display: flex;
  flex-wrap: wrap;
  gap: 24px;
  align-items: flex-start;
}

.plot {
  flex: 2 1 340px;
  min-width: 300px;
  background: var(--vp-c-bg);
  border-radius: 8px;
  padding: 8px;
}

.legend {
  margin-top: 14px;
  padding: 0 4px 4px;
}

.legend-head {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  font-size: 0.78rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  color: var(--vp-c-text-2);
  margin-bottom: 6px;
}

.legend-beta {
  text-transform: none;
  letter-spacing: 0;
  font-family: var(--vp-font-family-mono);
  color: var(--vp-c-text-3);
}

.legend-bar {
  position: relative;
  height: 12px;
  border-radius: 3px;
  border: 1px solid var(--vp-c-divider);
}

/* Where each population sits on the scale. The gap between the two is the readable signal: wide
   means the clustering is cleanly separated, narrow means the groups blur into each other. */
.legend-mark {
  position: absolute;
  top: -3px;
  width: 2px;
  height: 18px;
  margin-left: -1px;
  border-radius: 1px;
  background: var(--vp-c-text-1);
  box-shadow: 0 0 0 1px var(--vp-c-bg);
}

.legend-mark.boundary {
  background: #e4572e;
}

.legend-ends,
.legend-means {
  display: flex;
  justify-content: space-between;
  gap: 12px;
}

.legend-ends {
  font-size: 0.72rem;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  color: var(--vp-c-text-3);
  margin-top: 7px;
}

.legend-means {
  font-size: 0.75rem;
  color: var(--vp-c-text-2);
  margin-top: 5px;
  flex-wrap: wrap;
}

.legend-ends span:last-child,
.legend-means span:last-child {
  text-align: right;
}

.swatch {
  display: inline-block;
  width: 8px;
  height: 8px;
  border-radius: 2px;
  margin-right: 3px;
  vertical-align: baseline;
}

.swatch.boundary {
  background: #e4572e;
}

.swatch.interior {
  background: var(--vp-c-text-1);
}

.legend-note {
  margin: 10px 0 0;
  font-size: 0.78rem;
  line-height: 1.5;
  color: var(--vp-c-text-2);
  text-align: left;
}

.panel {
  flex: 1 1 220px;
  min-width: 210px;
  font-size: 0.85rem;
}

.panel h4 {
  margin: 0 0 8px;
  font-size: 0.8rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  color: var(--vp-c-text-2);
}

.panel h4:not(:first-child) {
  margin-top: 22px;
}

.panel dl {
  display: grid;
  grid-template-columns: auto 1fr;
  gap: 4px 12px;
  margin: 0;
}

.panel dt {
  color: var(--vp-c-text-2);
}

.panel dd {
  margin: 0;
  font-variant-numeric: tabular-nums;
  font-weight: 600;
}

.muted {
  color: var(--vp-c-text-3);
  font-weight: 400;
}

.hint {
  margin: 0;
  color: var(--vp-c-text-2);
  line-height: 1.5;
}

.clear {
  margin-top: 14px;
  background: var(--vp-c-bg);
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
  padding: 5px 12px;
  cursor: pointer;
  color: var(--vp-c-text-1);
  font-size: 0.8rem;
}

.clear:hover {
  border-color: var(--vp-c-brand-1);
}

.state {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 14px;
  height: 320px;
  color: var(--vp-c-text-2);
}

.state.error {
  color: var(--vp-c-danger-1, #d33);
}

.spinner {
  width: 30px;
  height: 30px;
  border: 3px solid var(--vp-c-divider);
  border-top-color: var(--vp-c-brand-1);
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

@keyframes spin {
  to {
    transform: rotate(360deg);
  }
}
</style>
