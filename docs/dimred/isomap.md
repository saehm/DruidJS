<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.ISOMAP(mistle.IRIS.values, {neighbors: 58, project: "SMACOF"});
    const isomap = dr.transform()
</script>

# ISOMAP

Isomap is a non-linear dimensionality reduction method based on spectral theory which tries to preserve the geodesic distances in the lower dimension.

## How It Works

Isomap extends MDS by computing the shortest path (geodesic distance) between all pairs of points using a neighborhood graph (K-nearest neighbors), then applies MDS to these geodesic distances.

## Why or When to Use

Use Isomap when the data lies on a curved, non-linear manifold and preserving the global non-linear geometric structure is important.

## Example

<RoundScatterplot :data=isomap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const isomap = new druid.ISOMAP(data, { neighbors: 58, project: "SMACOF" });

// 2. Compute the projection
const projection = isomap.transform();
```
