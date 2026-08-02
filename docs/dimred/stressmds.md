<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";

    const labels = mistle.IRIS.labels;
    const raw = new druid.StressMDS(mistle.IRIS.values, { weights: 0 }).transform();
    const sammonish = new druid.StressMDS(mistle.IRIS.values, { weights: -1 }).transform();
    const elastic = new druid.StressMDS(mistle.IRIS.values, { weights: -2 }).transform();
</script>

# StressMDS

Weighted metric MDS. It minimizes

$$\sigma(Y) = \sum_{i<j} w_{ij} \left( \lVert y_i - y_j \rVert - d_{ij} \right)^2$$

for a weighting you choose, which turns a family of separately-named methods into one parameter.

## How It Works

Setting `weights` to an exponent $q$ gives $w_{ij} = d_{ij}^{\,q}$:

| `weights` | objective | also known as |
| --- | --- | --- |
| `0` | raw stress | the objective [SMACOF](./smacof) minimizes |
| `-1` | Sammon stress | the objective [Sammon](./sammon) minimizes |
| `-2` | elastic scaling | Kamada-Kawai energy, see [KKMDS](./kkmds) |

A more negative exponent concentrates the objective on short distances. Optimization is
Jacobi-preconditioned gradient descent with a backtracking line search, warm-started from classical
MDS on the same distances.

## Why or When to Use

StressMDS does **not** replace [SMACOF](./smacof) or [Sammon](./sammon). Those solve the same
objectives at `0` and `-1` with their own historical algorithms and reach different local minima, so
swapping them out would silently change existing layouts. Reach for StressMDS when you want
something they cannot express:

- **exponents between or beyond** the three named ones, as a continuous local↔global dial;
- **explicit weight matrices**, where a weight of zero drops the pair from the objective — the
  standard way to handle missing or untrusted dissimilarities, which nothing else here offers;
- **fewer iterations**: on a benchmark suite it reached lower Sammon stress in 52 iterations than
  Sammon's own optimizer reached in 200.

## Example

The same Iris data at three weightings. The difference is subtle on data this small and clean, which
is itself worth knowing.

<div style="display:flex;flex-wrap:wrap;gap:16px;justify-content:center">
  <div style="text-align:center">
    <strong style="font-size:0.85em">weights: 0</strong>
    <RoundScatterplot :data="raw" :labels="labels" :size="260" :radius="4" />
  </div>
  <div style="text-align:center">
    <strong style="font-size:0.85em">weights: -1</strong>
    <RoundScatterplot :data="sammonish" :labels="labels" :size="260" :radius="4" />
  </div>
  <div style="text-align:center">
    <strong style="font-size:0.85em">weights: -2</strong>
    <RoundScatterplot :data="elastic" :labels="labels" :size="260" :radius="4" />
  </div>
</div>

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// An exponent
const projection = new druid.StressMDS(data, { weights: -1 }).transform();

// Or the named constants
const elastic = new druid.StressMDS(data, { weights: druid.WEIGHTS_ELASTIC }).transform();

// Or a function of the target distance
const custom = new druid.StressMDS(data, { weights: (d) => 1 / (1 + d) }).transform();
```

Missing data — a zero weight removes the pair from the objective entirely:

```javascript
const W = new druid.Matrix(N, N, (i, j) => (wasObserved(i, j) ? 1 : 0));
const projection = new druid.StressMDS(D, { metric: "precomputed", weights: W }).transform();
```

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `weights` | `-2` | Exponent, matrix, or function. Zero weights drop the pair. |
| `init_DR` | `"MDS"` | Starting configuration. The objective is non-convex, and classical MDS on the same distances keeps the descent out of poor local minima. |
| `learning_rate` | `0.1` | Dimensionless — the gradient is preconditioned by weighted degree, so this needs no rescaling for the data or the weighting. |
| `iterations` | `300` | Maximum gradient steps. |
| `epsilon` | `1e-6` | Stop once the relative stress improvement falls below this. |

## What the exponent actually buys

Worth stating plainly, because the obvious guess is wrong. Holding the solver, initialization and
budget fixed and varying only the exponent, over 5 datasets × 5 seeds:

| `weights` | 10-NN preserved | distance correlation |
| --- | --- | --- |
| `0` | 56.3% | 0.8382 |
| `-1` | 57.9% | 0.8309 |
| `-2` | **59.1%** | 0.8144 |
| `-3` | 57.9% | 0.7962 |

Global fidelity falls **monotonically** as the exponent drops — that part is reliable. Local
neighborhood preservation, however, *peaks* near `-2` rather than improving without limit, and on
data with no manifold or cluster structure it moves the other way throughout (an isotropic
30-dimensional Gaussian ran 19.6% at `0` down to 13.1% at `-3`). So `-2` is a reasonable default
rather than a maximum, and on structureless data `0` is the better choice.

See also the [`StressMDS` API reference](/api/classes/StressMDS).
