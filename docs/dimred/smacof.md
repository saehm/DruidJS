<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.SMACOF(mistle.IRIS.values);
    const mds = dr.transform()
</script>

# SMACOF

SMACOF (Scaling by Majorizing a Complicated Function) is a method for metric multidimensional scaling.

## How It Works

SMACOF is an iterative majorization algorithm for solving metric multidimensional scaling problems. It aims to minimize the "stress" function, which measures the difference between the given true pairwise distances and the corresponding pairwise distances in the lower-dimensional embedded space. In each iterations, it bounds the complicated stress function with a simpler convex majorant and minimizes the surrogate, guaranteeing a monotonic decrease in stress until convergence.

## Why or When to Use

Use SMACOF as a more robust and flexible alternative to Classical MDS when you need accurate preservation of distances. It is advantageous because it can easily be modified to use weights (ignoring missing distances or discounting less reliable similarities) and generally yields a tighter fit to the original distances than classical methods.

## Example

<RoundScatterplot :data=mds :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const smacof = new druid.SMACOF(data);

// 2. Compute the projection
const projection = smacof.transform();
```
