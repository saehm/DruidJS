<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";

    const labels = mistle.IRIS.labels;
    const kkmds = new druid.KKMDS(mistle.IRIS.values).transform();
</script>

# KKMDS

Kamada-Kawai multidimensional scaling: [StressMDS](./stressmds) fixed at `weights: -2`, so each pair
is weighted by the inverse square of its target distance.

$$\sigma(Y) = \sum_{i<j} \frac{\left( \lVert y_i - y_j \rVert - d_{ij} \right)^2}{d_{ij}^2}$$

## How It Works

That objective is the Kamada-Kawai energy from graph drawing, known in the MDS literature as *elastic
scaling*. Minimization is Jacobi-preconditioned gradient descent with a backtracking line search,
started from classical MDS on the same distances — the original paper's Newton-Raphson is replaced by
a warm start plus preconditioning, which avoids the O(N³) Hessian solves and most of the
initialization sensitivity Kamada-Kawai is known for.

## Why or When to Use

Its classic use is laying out a **graph** from its shortest-path distances: pass those as a
`"precomputed"` matrix. [MINFOTree](./minfotree) does exactly that with its spanning tree.

On ordinary feature data it behaves as the most locally-focused member of the stress family. Reach
for [StressMDS](./stressmds) directly if you want a different weighting; this class exists so the
Kamada-Kawai name resolves to the right objective without having to remember which exponent it is.

## Example

<RoundScatterplot :data="kkmds" :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

// On feature data
const projection = new druid.KKMDS(data).transform();

// On a graph, from its shortest-path distance matrix
const layout = new druid.KKMDS(graphDistances, { metric: "precomputed" }).transform();
```

## Parameters

Identical to [StressMDS](./stressmds#parameters), except that `weights` is fixed at `-2` and cannot
be overridden.

See also the [`KKMDS` API reference](/api/classes/KKMDS).
