<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";

    const labels = mistle.IRIS.labels;

    const dr = new druid.FASTMAP(mistle.IRIS.values);
    const fastmap = dr.transform()
</script>

# FastMap

FastMap is an efficient algorithm that maps objects into a Euclidean space such that distances are preserved as much as possible, with a linear time complexity.

## How It Works

FastMap uses a heuristic approach to approximate MDS linearly. It iteratively projects data points onto an orthogonal line defined by two distant pivot points, extracting coordinates one dimension at a time.

## Why or When to Use

Use FastMap when you have a distance matrix and need a very fast (O(N)) linear dimensionality reduction, typically as a preprocessing step or for very large datasets where classical MDS is too slow.

## Example

<RoundScatterplot :data=fastmap 
        :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const fastmap = new druid.FASTMAP(data);

// 2. Compute the projection
const projection = fastmap.transform();
```
