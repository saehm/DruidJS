<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;  
    const dr = new druid.SQDMDS(mistle.IRIS.values);
    const sqdmds = dr.transform()
</script>

# SQDMDS

Sequential Quadratic Distance MDS (SQDMDS) is a method for multidimensional scaling that uses a sequential quadratic programming approach.

## How It Works

Sequential Quadratic Distance MDS (SQDMDS) formulates the MDS stress-minimization problem as a sequence of quadratic programming problems, offering efficient convergence.

## Why or When to Use

Use SQDMDS when classical MDS is too slow or you want an alternative mathematical optimization approach to minimize MDS stress more effectively.

## Example

<RoundScatterplot :data=sqdmds :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const sqdmds = new druid.SQDMDS(data);

// 2. Compute the projection
const projection = sqdmds.transform();
```
