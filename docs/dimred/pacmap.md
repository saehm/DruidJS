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

PaCMAP constructs three categories of point pairs in the high-dimensional space: nearest neighbor pairs that capture local structure, mid-near pairs (the second-closest among six random candidates drawn from the whole dataset) that capture intermediate-range relationships, and further pairs that act as repulsive anchors. A three-phase optimization schedule progressively shifts focus from global structure (via high MN weights) to local refinement (MN disabled), using Adam optimization throughout.

Neighbors are not taken straight from the nearest-neighbor graph. Each candidate distance is first divided by a local scale — the mean distance to the 4th through 6th neighbor of both endpoints — and the pairs are chosen by that rescaled distance. This is what keeps dense regions from dominating the attractive term, and it is why the search asks for 50 more candidates than it finally keeps.

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

## Larger Datasets

By default PaCMAP runs an exact neighbor search whose selection happens inside a WASM SIMD kernel.
That is faster than a tree at the `n_neighbors + 50` candidates the rescaling needs, and it matches
the exact search the reference implementation performs — but it is still O(N²). Pass any KNN index
as `knn` to swap it for an approximate one:

```javascript
const knn = new druid.HNSW(data, { metric: druid.euclidean, ef: 100 });
const pacmap = new druid.PaCMAP(data, { n_neighbors: 10, knn });
```

[HNSW](/api/classes/HNSW), [Annoy](/api/classes/Annoy) and [NNDescent](/api/classes/NNDescent) all
work here. Note that the index must recall `n_neighbors + 50` candidates per point, so give it
enough search width — an index tuned for 10 neighbors will not do.

Inputs wider than 100 dimensions are reduced to 100 by PCA before the search and the
initialization; set `apply_pca: false` to keep the full width.
