<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.PaCMAP(mistle.IRIS.values);
    const pacmap = dr.transform()
</script>

# PaCMAP

Pairwise Controlled Manifold Approximation Projection (PaCMAP) is a dimensionality reduction technique that uses three explicit types of point pairs — nearest neighbor (NN), mid-near (MN), and further (FP) pairs — with a dynamic three-phase weight schedule to preserve both local and global structure.

## How It Works

PaCMAP constructs three categories of point pairs in the high-dimensional space: nearest neighbor pairs that capture local structure, mid-near pairs (the second-closest among random candidates) that capture intermediate-range relationships, and further pairs that act as repulsive anchors. A three-phase optimization schedule progressively shifts focus from global structure (via high MN weights) to local refinement (MN disabled), using Adam optimization throughout.

## Why or When to Use

A strong alternative to UMAP and t-SNE that explicitly controls global structure preservation through its MN pair mechanism. Often produces cleaner cluster separation and is less sensitive to hyperparameter choices than UMAP.

## Example

<RoundScatterplot :data=pacmap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const pacmap = new druid.PaCMAP(data, { n_neighbors: 10 });

// 2. Compute the projection (450 iterations across 3 phases by default)
const projection = pacmap.transform();

// Alternatively, use a generator for animation:
// for (const proj of pacmap.generator()) { ... }
```
