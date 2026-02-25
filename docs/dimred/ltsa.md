<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LTSA(mistle.IRIS.values, {neighbors: 51});
    const ltsa = dr.transform()
</script>

# LTSA

Local Tangent Space Alignment (LTSA) characterizes the local geometry of the data manifold using the tangent spaces at each data point, and then aligns these local tangent spaces to construct the global coordinate.

## How It Works

Local Tangent Space Alignment (LTSA) computes the tangent space at each data point using its nearest neighbors. It then aligns these local tangent spaces to construct a global low-dimensional coordinate system.

## Why or When to Use

Use LTSA for unfolding highly nonlinear and complex manifolds where the tangent space provides a good local approximation.

## Example

<RoundScatterplot :data=ltsa :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const ltsa = new druid.LTSA(data);

// 2. Compute the projection
const projection = ltsa.transform();
```
