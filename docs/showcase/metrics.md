# Metric Sensitivity

The choice of distance metric (how "distance" is calculated) can drastically change the resulting projection. Each metric defines "distance" differently, leading to varied cluster shapes and separations.

<script setup>
import MetricComparison from '../components/MetricComparison.vue'
</script>

<MetricComparison />

## How-to (Code)

DruidJS provides several built-in metrics (Euclidean, Cosine, Manhattan, etc.). You can pass them as the `metric` parameter to any DR algorithm like [`UMAP`](/api/classes/UMAP) or [`TSNE`](/api/classes/TSNE).

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... data ... */
];

// Example: Using Cosine similarity with UMAP
const umap = new druid.UMAP(data, {
  metric: druid.cosine,
});
const projection = umap.transform();

// Example: Using Manhattan distance with TSNE
const tsne = new druid.TSNE(data, {
  metric: druid.manhattan,
});
const projection2 = tsne.transform();
```
