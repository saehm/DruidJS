<script setup>
import * as druid from "@saehrimnir/druidjs";
import * as mistle from "@saehrimnir/mistle";
import * as d3 from "d3";
import { onMounted, ref, watch } from "vue";
import { applyZScore } from "./utils.js";

const svgRef = ref(null);
const containerRef = ref(null);

const datasets = [
  { name: "Iris", key: "IRIS" },
  { name: "Wine", key: "WINE" },
  { name: "Penguins", key: "PENGUINS" },
];

const selectedDataset = ref("IRIS");
const selectedLinkage = ref("complete");
const svgWidth = ref(800);
const svgHeight = ref(500);

const colorScale = d3.scaleOrdinal(d3.schemeObservable10);

const renderDendrogram = () => {
  if (!svgRef.value) return;

  const ds = datasets.find((d) => d.key === selectedDataset.value);
  const rawData = mistle[selectedDataset.value].values;
  const labels = mistle[selectedDataset.value].labels;

  const data = applyZScore(rawData);
  svgHeight.value = Math.max(500, data.length * 15);

  // Run hierarchical clustering on the main thread
  const hc = new druid.HierarchicalClustering(data, { linkage: selectedLinkage.value });
  const rootNode = hc.root;

  // Convert DruidJS tree structure to D3 hierarchy format
  const toD3Hierarchy = (node) => {
    if (!node) return null;
    if (node.isLeaf) {
      // index could be an array or a single number depending on how the clustering handles it,
      // but the leaf nodes have a single index representing the row.
      const idx = Array.isArray(node.index) ? node.index[0] : node.index;
      return {
        name: labels[idx],
        originalIndex: idx,
        dist: 0,
      };
    }
    return {
      name: "",
      dist: node.dist,
      children: [toD3Hierarchy(node.left), toD3Hierarchy(node.right)].filter(Boolean),
    };
  };

  const d3Data = toD3Hierarchy(rootNode);

  // Clear previous SVG contents
  d3.select(svgRef.value).selectAll("*").remove();

  const width = svgWidth.value;
  const height = svgHeight.value;
  const margin = { top: 20, right: 120, bottom: 20, left: 40 };
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  const svg = d3
    .select(svgRef.value)
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Create cluster layout
  const cluster = d3.cluster().size([innerHeight, innerWidth]);
  const root = d3.hierarchy(d3Data);
  cluster(root);

  // If we want the y-position to map to distance rather than uniform depth, we can override root positions:
  // But for a simple dendrogram visual, d3.cluster() depth is fine. By default it spaces uniformly.
  // To use branch lengths proportional to `dist`, we can do:
  // Find max dist to create a scale.
  const maxDist = rootNode.dist;
  const xScale = d3.scaleLinear().domain([0, maxDist]).range([innerWidth, 0]);

  root.each((d) => {
    if (d.data.dist !== undefined && d.children) {
      d.y = xScale(d.data.dist);
    } else if (!d.children) {
      d.y = innerWidth; // leaves at the very right
    }
  });

  // Link generator with straight/orthogonal lines typical for dendrograms
  const linkGenerator = d3
    .linkHorizontal()
    .x((d) => d.y)
    .y((d) => d.x);

  // Draw links with right-angle elbows
  svg
    .selectAll("path.link")
    .data(root.links())
    .enter()
    .append("path")
    .attr("class", "link")
    .attr("d", (d) => {
      return `M${d.source.y},${d.source.x} 
                V${d.target.x} 
                H${d.target.y}`;
    })
    .attr("fill", "none")
    .attr("stroke", "var(--vp-c-divider)")
    .attr("stroke-width", 1.5)
    .attr("stroke-opacity", 0.8);

  // Draw nodes
  const node = svg
    .selectAll("g.node")
    .data(root.descendants())
    .enter()
    .append("g")
    .attr("class", "node")
    .attr("transform", (d) => `translate(${d.y},${d.x})`);

  node
    .append("circle")
    .attr("r", (d) => (d.children ? 3 : 5))
    .attr("fill", (d) => (d.children ? "var(--vp-c-text-3)" : colorScale(d.data.name)));

  // Text labels for leaves
  node
    .filter((d) => !d.children)
    .append("text")
    .attr("dx", 10)
    .attr("dy", 4)
    .text((d) => d.data.name)
    .attr("font-size", "12px")
    .attr("fill", "var(--vp-c-text-1)")
    .attr("font-family", "var(--vp-font-family-base)");
};

watch([selectedDataset, selectedLinkage], () => {
  renderDendrogram();
});

onMounted(() => {
  // Simple responsive handling
  if (containerRef.value) {
    svgWidth.value = containerRef.value.clientWidth;
  }
  renderDendrogram();

  window.addEventListener("resize", () => {
    if (containerRef.value) {
      svgWidth.value = containerRef.value.clientWidth;
      renderDendrogram();
    }
  });
});
</script>

<template>
  <div class="dendrogram-showcase">
    <div class="controls">
      <div class="control-item">
        <label>Dataset</label>
        <select v-model="selectedDataset" class="modern-select">
          <option v-for="opt in datasets" :key="opt.key" :value="opt.key">
            {{ opt.name }}
          </option>
        </select>
      </div>
      <div class="control-item">
        <label>Linkage Strategy</label>
        <select v-model="selectedLinkage" class="modern-select">
          <option value="single">Single Link</option>
          <option value="average">Average Link</option>
          <option value="complete">Complete Link</option>
        </select>
      </div>
    </div>

    <div class="description">
      <p>
        Hierarchical Clustering builds a tree representation (Dendrogram). The branch lengths
        represent the distance at which groups merge. Use the controls to see how different linkage
        criteria alter the clustering tree topology.
      </p>
    </div>

    <div class="plot-container" ref="containerRef">
      <svg ref="svgRef"></svg>
    </div>
  </div>
</template>

<style scoped>
.dendrogram-showcase {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 24px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.controls {
  display: flex;
  gap: 20px;
  margin-bottom: 20px;
  flex-wrap: wrap;
}

.control-item {
  display: flex;
  flex-direction: column;
  gap: 6px;
}

.control-item label {
  font-size: 0.85em;
  font-weight: 600;
  text-transform: uppercase;
  color: var(--vp-c-text-2);
}

.modern-select {
  padding: 8px 12px;
  border-radius: 6px;
  border: 1px solid var(--vp-c-divider);
  background: var(--vp-c-bg);
  min-width: 180px;
}

.description {
  margin-bottom: 20px;
  font-size: 0.9em;
  color: var(--vp-c-text-2);
}

.plot-container {
  background: var(--vp-c-bg);
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  overflow: auto;
  max-height: 600px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
  display: flex;
  justify-content: flex-start;
  align-items: flex-start;
}
</style>
