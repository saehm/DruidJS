<script setup>
import * as mistle from "@saehrimnir/mistle";
import { computed, onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import TreePlot from "./TreePlot.vue";
import { applyZScore } from "./utils.js";
import { runInWorker } from "./workerPool.js";

const loading = ref(true);
const results = ref({ tsne: null, topomap: null, minfo: null });
const truth = ref([]);
const labels = ref([]);

/**
 * The mammoth has 11 body parts, but the plots' default scheme has only 10 colours, so two parts
 * would render identically and a torn-off limb could not be told from a recoloured one. These are
 * the Observable 10 plus two more from the Tableau palette.
 */
const PALETTE = [
  "#4269d0", "#efb118", "#ff725c", "#6cc5b0", "#a463f2", "#97bbf5",
  "#9c6b4e", "#01ab63", "#e45756", "#b279a2", "#ff9da7", "#8c8c8c",
];

const pointColors = computed(() => labels.value.map((l) => PALETTE[l % PALETTE.length]));

/** Points kept from the 1200-point mammoth. t-SNE at the full size takes ~16s, which is too long
 *  to block a page on; 600 keeps every body part and every method under two seconds. */
const N = 600;

onMounted(async () => {
  try {
    const mammoth = mistle.MAMMOTH;
    // Even stride rather than a random draw: the mammoth is a surface scan, so a regular subsample
    // thins it uniformly and leaves each body part contiguous.
    const step = mammoth.values.length / N;
    const indices = Array.from({ length: N }, (_, i) => Math.floor(i * step));

    const raw = indices.map((i) => mammoth.values[i]);
    labels.value = indices.map((i) => mammoth.labels[i]);
    // Side view of the 3D scan, so the shape the projections have to preserve is visible.
    truth.value = raw.map(([x, y]) => [x, y]);

    const data = applyZScore(raw);

    // MINFOTree is given the body-part labels rather than a cluster count: its default clustering is
    // Ward, which is variance-based and would cut the mammoth into compact blobs that have nothing
    // to do with limbs.
    const [tsne, topomap, minfo] = await Promise.all([
      runInWorker("DR", "TSNE", data, { perplexity: 30 }, 1000),
      runInWorker("DR", "TopoMap", data, {}),
      runInWorker("Tree", "MINFOTree", data, { labels: [...labels.value] }),
    ]);

    results.value.tsne = tsne;
    results.value.topomap = topomap;
    results.value.minfo = minfo;
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
      Calculating projections…
    </div>
    <div v-else class="grid">
      <div class="side">
        <h3>The original</h3>
        <p class="desc">
          A 3D scan of a mammoth skeleton, seen from the side. Colour marks one of 11 labelled body
          parts. This is the structure the projections have to keep.
        </p>
        <div class="plot-box">
          <RoundScatterplot :data="truth" :colors="pointColors" :labels="labels" :size="330" :radius="3" />
        </div>
      </div>

      <div class="side">
        <h3>t-SNE</h3>
        <p class="desc">
          Optimises local similarity, so it shatters the skeleton: the legs and tusks break off into
          islands and the limbs lose their attachment to the body.
        </p>
        <div class="plot-box">
          <RoundScatterplot :data="results.tsne" :colors="pointColors" :labels="labels" :size="330" :radius="3" />
        </div>
      </div>

      <div class="side">
        <h3>TopoMap</h3>
        <p class="desc">
          Preserves 0-dimensional topology. The skeleton stays in one connected piece and the limbs
          stay joined where they were joined in 3D, though the outline is not reproduced.
        </p>
        <div class="plot-box">
          <RoundScatterplot :data="results.topomap" :colors="pointColors" :labels="labels" :size="330" :radius="3" />
        </div>
      </div>

      <div class="side">
        <h3>MINFO Tree</h3>
        <p class="desc">
          Keeps the spanning tree itself as the output. Red edges join two different body parts —
          only 17 of 599 do, against a minimum of 10 for 11 parts.
        </p>
        <div class="plot-box">
          <TreePlot
            :data="results.minfo?.projection"
            :edges="results.minfo?.edges"
            :labels="results.minfo?.labels"
            :colors="pointColors"
            :size="330"
            :radius="2.5"
            show-bridges
          />
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

/* Two by two rather than a single row: four 330px plots side by side would each be squeezed to
   under 200px in the docs column, which is too small to read a skeleton in. */
.grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 28px;
}

.side {
  text-align: center;
  min-width: 0;
}

.side h3 {
  margin-top: 0;
}

.desc {
  font-size: 0.85em;
  color: var(--vp-c-text-2);
  margin-bottom: 16px;
  min-height: 5em;
  text-align: left;
}

.plot-box {
  background: var(--vp-c-bg);
  border-radius: 8px;
  padding: 10px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  display: inline-block;
  max-width: 100%;
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

