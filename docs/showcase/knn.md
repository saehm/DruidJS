# K-Nearest Neighbors

K-Nearest Neighbors (KNN) algorithms are a fundamental building block for many complex Dimensionality Reduction and Clustering algorithms (such as UMAP, t-SNE, and OPTICS). DruidJS provides extremely fast web-native implementations of multiple exact and approximate KNN algorithms to find local neighborhoods efficiently.

This showcase allows you to interactively explore various KNN algorithms on synthetic datasets:

- **Uniform Data Generation**: Adjust the number of points (up to 5,000) to generate a 2D uniform distribution.
- **Real-time Pointer Queries**: Move your mouse over the scatterplot to perform instantaneous KNN searches at your cursor's position.
- **Algorithm Comparison**: Switch between different exact and approximate methods to observe their build and search times.

<script setup>
import KNNDemo from '../components/KNNDemo.vue'
</script>

<KNNDemo />

## How-to (Code)

DruidJS supports algorithms optimized for high dimensionality like [**HNSW**](/api/classes/HNSW), [**NNDescent**](/api/classes/NNDescent), [**Annoy**](/api/classes/Annoy), [**LSH**](/api/classes/LSH), and spatial partition trees like [**KD-Tree**](/api/classes/KDTree) and [**BallTree**](/api/classes/BallTree).

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional floats ... */
];

// 1. Build the index (e.g. NN-Descent)
const nndescent = new druid.NNDescent(data, {
  metric: druid.euclidean,
  samples: 50,
  rho: 0.8,
});

// 2. Perform a fast query (e.g. at a specific coordinate)
const queryVector = [0.5, 0.5];
const nearestNeighbors = nndescent.search(queryVector, 15);

// The list returned contains elements, their original index, and the computed distance
for (const match of nearestNeighbors) {
  console.log(`Neighbor ID: ${match.index}, Distance: ${match.distance}`);
}
```

## Querying by Index

Every index also answers queries by the position of an already-indexed element. `search_by_index`
returns the `k` nearest **other** elements — the queried point is never its own neighbor, which is
what building a neighborhood graph almost always wants:

```javascript
// The 5 points nearest to element 42, excluding element 42 itself
const neighbors = nndescent.search_by_index(42, 5);
```

Note the asymmetry with `search`: an arbitrary query point has no "self" to leave out, so there `k`
means "k results", while `search_by_index` means "k neighbors". An out-of-range index returns an
empty list.
