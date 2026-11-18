<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LSP(mistle.IRIS.values);
    const lsp = dr.transform()
</script>

# LSP

Least Square Projection (LSP) places a few representative points first and positions all other points relative to their neighbors, keeping local neighborhoods while staying fast on large datasets.

## How It Works

LSP picks control points (by default about √N of them, with k-medoids) and projects only those with classical MDS. Every other point should sit at the mean of its nearest neighbors, and every control point at its MDS position. Together these conditions form a sparse linear system, which LSP solves in the least-squares sense.

## Why or When to Use

Use LSP when the dataset is too large for a full distance-based method like MDS, but you still want a projection whose overall layout comes from distances and whose local neighborhoods are kept.

## Example

<RoundScatterplot :data=lsp :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const lsp = new druid.LSP(data);

// 2. Compute the projection
const projection = lsp.transform();
```
