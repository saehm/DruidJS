<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LSP(mistle.IRIS.values);
    const lsp = dr.transform()
</script>

# LSP

Local Space Preservation (LSP) aims to preserve the local geometry of the data points.

## How It Works

Local Space Preservation (LSP) focuses on preserving local geometric relationships by treating data points and their neighbors as a local patch, aiming to align these patches in the lower dimension.

## Why or When to Use

Use LSP when local structure is critical and computing a robust local patch alignment is desired over purely distance-based preserving methods.

## Example

<RoundScatterplot :data=lsp :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const lsp = new druid.LSP(data);

// 2. Compute the projection
const projection = lsp.transform();
```
