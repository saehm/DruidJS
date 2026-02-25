<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.SAMMON(mistle.IRIS.values);
    const sammon = dr.transform();
</script>

# Sammon

Sammon Mapping is a non-linear approach for mapping higher-dimensional space to a space of lower dimensionality, which attempts to preserve structure.

## How It Works

Sammon Mapping is a variation of metric MDS that uses a specific cost function (Sammon's stress) which heavily penalizes errors in preserving smaller distances over larger ones.

## Why or When to Use

Use Sammon Mapping when preserving local distances (small distances between nearby points) is more important than global distances, often resulting in better cluster separation.

## Example

<RoundScatterplot :data=sammon :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const sammon = new druid.SAMMON(data);

// 2. Compute the projection
const projection = sammon.transform();
```
