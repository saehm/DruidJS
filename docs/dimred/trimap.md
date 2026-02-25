<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.TriMap(mistle.IRIS.values);
    const trimap = dr.transform()
</script>

# TriMap

TriMap is a dimensionality reduction method that uses triplets of points to learn a low-dimensional embedding that preserves the global structure of the data.

## How It Works

TriMap samples triplets of points (an anchor, a positive point, and a negative point) to capture the relationships between them. It then optimizes a low-dimensional embedding to satisfy these triplet constraints.

## Why or When to Use

Offers a balance between preserving local and global structure. It often runs faster than t-SNE and preserves global structure better than UMAP or t-SNE.

## Example

<RoundScatterplot :data=trimap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const trimap = new druid.TriMap(data);

// 2. Compute the projection
const projection = trimap.transform();
```
