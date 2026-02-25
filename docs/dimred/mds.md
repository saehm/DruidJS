<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.MDS(mistle.IRIS.values);
    const mds = dr.transform()
</script>

# MDS

Multidimensional Scaling (MDS) is a means of visualizing the level of similarity of individual cases of a dataset.

## How It Works

Multidimensional Scaling (MDS) takes a matrix of pairwise distances and finds a configuration of points in low-dimensional space such that the distances between the points are preserved as closely as possible.

## Why or When to Use

Use MDS when you only have distance/dissimilarity data between objects rather than feature vectors, and you want a simple global preservation of these spatial distances.

## Example

<RoundScatterplot :data=mds :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const mds = new druid.MDS(data);

// 2. Compute the projection
const projection = mds.transform();
```
