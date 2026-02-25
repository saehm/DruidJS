<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.TSNE(mistle.IRIS.values);
    const tsne = dr.transform()
</script>

# t-SNE

t-Distributed Stochastic Neighbor Embedding (t-SNE) is a statistical method for visualizing high-dimensional data by giving each datapoint a location in a two or three-dimensional map.

## How It Works

t-SNE converts pairwise Euclidean distances into conditional probabilities that represent similarities. It minimizes the Kullback-Leibler divergence between the probabilities of the low-dimensional embedding and the high-dimensional data, using a Student-t distribution for the lower dimension to alleviate the "crowding problem".

## Why or When to Use

Highly effective for exploring and visualizing high-dimensional data in 2D or 3D, especially for showing distinct clusters. It is mostly intended for visualization rather than downstream machine learning.

## Example

<RoundScatterplot :data=tsne :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the iterative algorithm
const tsne = new druid.TSNE(data, { perplexity: 30 });

// 2. Compute the projection (e.g. 1000 iterations)
const projection = tsne.transform(1000);

// Alternatively, use a generator for animation:
// for (const proj of tsne.generator(1000)) { ... }
```
