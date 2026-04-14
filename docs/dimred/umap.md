<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.UMAP(mistle.IRIS.values);
    const umap = dr.transform()
</script>

# UMAP

Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualization similarly to t-SNE, but also for general non-linear dimension reduction.

## How It Works

Uniform Manifold Approximation and Projection (UMAP) constructs a high-dimensional graph representation of the data and optimizes a low-dimensional graph to be as structurally similar as possible, grounded in Riemannian geometry.

## Why or When to Use

An excellent, fast alternative to t-SNE that scales well to large datasets and tends to preserve both local and global data structures effectively. For a more recent alternative with explicit global structure control, see [PaCMAP](/dimred/pacmap) and [LocalMAP](/dimred/localmap).

## Example

<RoundScatterplot :data=umap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the iterative algorithm
const umap = new druid.UMAP(data, { n_neighbors: 15, min_dist: 0.1 });

// 2. Compute the projection (e.g. 500 iterations)
const projection = umap.transform(500);

// Alternatively, use a generator for animation:
// for (const proj of umap.generator(500)) { ... }
```
