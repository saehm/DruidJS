<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.TopoMap(mistle.IRIS.values);
    const topomap = dr.transform()
</script>

# TopoMap

TopoMap is a technique ensuring topological preservation while embedding data in lower dimensions.

## How It Works

TopoMap uses a minimum spanning tree to preserve topological features (like connectivity) of the data when projecting it into a lower dimension.

## Why or When to Use

Use TopoMap when it is crucial to guarantee that there are no topological intersections or when preserving the exact connectivity structure of the original data is paramount.

## Example

<RoundScatterplot :data=topomap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const topomap = new druid.TopoMap(data);

// 2. Compute the projection
const projection = topomap.transform();
```
