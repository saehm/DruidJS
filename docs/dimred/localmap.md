<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LocalMAP(mistle.IRIS.values);
    const localmap = dr.transform()
</script>

# LocalMAP

LocalMAP is a variant of PaCMAP that improves local cluster separation by dynamically resampling further pairs (FP) in the third optimization phase using nearby points in the current low-dimensional embedding, rather than random non-neighbors.

## How It Works

LocalMAP runs identically to PaCMAP for the first two phases. In phase 3, instead of fixed random further pairs, it rebuilds the FP pairs from embedding-space neighbors and applies a distance-based weight scaling (`low_dist_thres / (2 × √d_ij)`) for pairs that are close in the current embedding. This pushes nearby clusters apart more aggressively, sharpening local boundaries.

## Why or When to Use

Use LocalMAP when PaCMAP produces clusters that are still somewhat merged or when fine-grained local separation is important. It adds negligible overhead over PaCMAP and generally produces crisper cluster boundaries.

## Example

<RoundScatterplot :data=localmap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const localmap = new druid.LocalMAP(data, {
    n_neighbors: 10,
    low_dist_thres: 10  // distance threshold for local FP resampling in phase 3
});

// 2. Compute the projection (450 iterations across 3 phases by default)
const projection = localmap.transform();

// Alternatively, use a generator for animation:
// for (const proj of localmap.generator()) { ... }
```
