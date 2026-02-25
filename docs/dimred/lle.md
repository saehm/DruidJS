<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LLE(mistle.IRIS.values, {neighbors: 60});
    const lle = dr.transform()
</script>

# LLE

Locally Linear Embedding (LLE) seeks a lower-dimensional projection of the data which preserves distances within local neighborhoods.

## How It Works

Locally Linear Embedding (LLE) represents each data point as a linear combination of its nearest neighbors. It then finds a low-dimensional embedding that preserves these local linear relationships.

## Why or When to Use

Use LLE when the data lies on a non-linear manifold and you want to preserve local neighborhood distances (local geometry) while unfolding the manifold.

## Example

<RoundScatterplot :data=lle :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const lle = new druid.LLE(data);

// 2. Compute the projection
const projection = lle.transform();
```
