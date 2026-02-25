# KNN Image Search

This example demonstrates how DruidJS can be used for fast similarity search in high-dimensional datasets.

Instead of dimensionality reduction, we use a **Hierarchical Navigable Small World (HNSW)** index to find similar images in a dataset of handwritten digits (MNIST).

<script setup>
import ImageSearch from '../components/ImageSearch.vue'
</script>

<ImageSearch />

## Why use HNSW for Search?

- **High Dimensionality**: Each image is a $28 \times 28$ grid, resulting in a 784-dimensional vector. Traditional search is slow in such spaces (the "Curse of Dimensionality").
- **Efficient Retrieval**: HNSW is an approximate nearest neighbor algorithm that builds a multi-layer graph, allowing for sub-millisecond search times even in large datasets.
- **Off-Thread Performance**: The index is built and searched inside a **Web Worker**, ensuring the main UI thread never hitches while processing the 784-D comparisons.

## How to use

1. **Wait** for the MNIST dataset to load and the HNSW index to be built.
2. **Click** any digit in the left panel.
3. **View** the 12 most similar digits found by the algorithm in the right panel.

## How-to (Code)

To implement fast similarity search, use the [`HNSW`](/api/classes/HNSW) class. It builds an index once and allows for multiple high-speed queries.

```javascript
import * as druid from "@saehrimnir/druidjs";

// Your high-dimensional data (e.g. 784-D MNIST vectors)
const points = [
  /* ... thousands of points ... */
];

// 1. Build the HNSW index
const hnsw = new druid.HNSW(points, {
  m: 16, // Connections per element
  ef_construction: 200, // Construction quality
});

// 2. Define a query point
const query = points[0];

// 3. Search for the K nearest neighbors
const K = 12;
const neighbors = hnsw.search(query, K);

// Results are objects: { element, index, distance }
neighbors.forEach((n) => {
  console.log(`Found neighbor at index ${n.index} with distance ${n.distance}`);
});
```
