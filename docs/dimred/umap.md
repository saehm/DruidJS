<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.UMAP(mistle.IRIS.values);
    const umap = dr.transform()
</script>

# UMAP

Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualization similarly to t-SNE, but also for general non-linear dimension reduction.

## How It Works

Uniform Manifold Approximation and Projection (UMAP) constructs a high-dimensional graph representation of the data and optimizes a low-dimensional graph to be as structurally similar as possible, grounded in Riemannian geometry.

## Why or When to Use

An excellent, fast alternative to t-SNE that scales well to large datasets and tends to preserve both local and global data structures effectively. For a more recent alternative with explicit global structure control, see [PaCMAP](/dimred/pacmap) and [LocalMAP](/dimred/localmap).

## Example

<RoundScatterplot :data=umap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the iterative algorithm
const umap = new druid.UMAP(data, { n_neighbors: 15, min_dist: 0.1 });

// 2. Compute the projection (e.g. 500 iterations)
const projection = umap.transform(500);

// Alternatively, use a generator for animation:
// for (const proj of umap.generator(500)) { ... }
```

## Reproducibility

Passing a `seed` pins the result **for a given engine and library build** — re-running the same code
in the same browser or Node version always produces the same embedding.

It does not pin it *across* environments. UMAP's stochastic gradient descent is chaotic: an epoch
chains thousands of edge updates, and the repulsive term's stiff `1 / (0.01 + d)` factor lets a
difference in the final bit grow into a visible one within two epochs. Such differences are
unavoidable, because ECMAScript specifies `Math.pow` and `Math.exp` as *implementation-approximated*
— two engines may legitimately disagree in the last bit, as may the WASM kernel and its JS fallback.

In practice this behaves like changing the seed, not like a loss of quality. On well-separated
clusters:

| comparison | shared 10-nearest-neighbours | cluster purity |
| :--- | ---: | ---: |
| WASM vs JS fallback, same seed | 65.8% | 100% / 100% |
| seed 1 vs seed 2, same path | 56.9% | 100% |

The two code paths differ *less* from each other than two seeds do, and the cluster structure is
preserved either way.

If you need a byte-identical picture, pin the engine and the library version, or store the resulting
coordinates rather than recomputing them. [t-SNE](./tsne) and [TriMap](./trimap) are not chaotic and
do reproduce bit-identically across paths; [SAMMON](./sammon) behaves like UMAP.
