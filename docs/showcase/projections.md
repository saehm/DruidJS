# Standard Projections

The following grid shows several dimensionality reduction techniques applied to the classic Iris dataset. Each method attempts to project the 4D data into 2D while preserving different aspects of the original structure.

<script setup>
import ShowcaseGrid from '../components/ShowcaseGrid.vue'
</script>

<ShowcaseGrid />

## How-to (Code)

Most dimensionality reduction algorithms in DruidJS follow a standard interface. You instantiate the algorithm with your data and calling `transform()` to get the projection. You can find out more in the [API Reference](/api/index).

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
    /* ... multi-dimensional data ... */
];

// Example: PCA
const pca = new druid.PCA(data);
const projectionPCA = pca.transform();

// Example: t-SNE (Iterative)
const tsne = new druid.TSNE(data, { perplexity: 30 });
const projectionTSNE = tsne.transform(1000); // 1000 iterations
```
