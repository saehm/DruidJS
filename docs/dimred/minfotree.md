<script setup>
    import TreePlot from "../components/TreePlot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";

    const dr = new druid.MINFOTree(mistle.IRIS.values, { clusters: 3 });
    const projection = dr.transform();
    const edges = dr.edges;
    const labels = Array.from(dr.labels);
</script>

# MINFO Tree

Minimum Information Trees visualize the structure of **clustered** high-dimensional data by building
a tree over the points and laying that tree out, rather than fitting an embedding to the data
directly.

## How It Works

The points are read as one configuration of a q-state Potts model on a k-nearest-neighbor graph:
each point carries a cluster label, and neighbouring points prefer to agree. From that model the
method estimates an inverse temperature β by maximum pseudo-likelihood, then approximates a local
shape operator at every vertex from the first- and second-order Fisher information. Those curvatures
weight the graph's edges — shrunk by `alpha` when both ends share a label — and the minimum spanning
tree of the weighted graph is the Minimum Information Tree.

The layout is [KKMDS](./kkmds) over the tree's shortest-path distances, which is the Kamada-Kawai
algorithm the paper specifies.

## Why or When to Use

Use it when you already have a clustering and want to *inspect* it: which clusters touch, where the
boundaries are, and which points sit between groups. Because the output is a tree rather than a cloud
of coordinates, it answers questions a scatterplot cannot — most usefully, the unique path between
any two points. See [Cluster Diagnostics](/showcase/diagnostics) for that.

Do not reach for it as a general-purpose projection. The paper's own limitations section is blunt
about the dependency: the tree inherits whatever the clustering got wrong. That is also why
`clusters` (or `labels`) is required rather than defaulted.

## Example

Iris, clustered into three groups with Ward linkage. Points are coloured by the recovered cluster;
edges are the tree.

<TreePlot :data="projection" :edges="edges" :labels="labels" :size="500" :radius="4" show-bridges />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// `clusters` tells it how many groups to partition the data into for the labels field.
const minfo = new druid.MINFOTree(data, { clusters: 3 });

// The 2D layout, like any other DR method
const projection = minfo.transform();

// ...but the tree is the actual output
const edges = minfo.edges; // [[u, v, weight], ...], one per pair of adjacent points
const labels = minfo.labels; // the labels field, remapped to 0 … q-1
const curvature = minfo.curvature; // information curvature per point
const beta = minfo.beta; // the estimated Potts inverse temperature
```

Bring your own labels instead of clustering internally:

```javascript
const minfo = new druid.MINFOTree(data, { labels: myClassLabels });
```

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `clusters` | — | Number of clusters to partition the data into. Required unless `labels` is given. With the default hierarchical clustering the range is `2 … N-1`; `"kmeans"` goes finer. |
| `labels` | `null` | Precomputed labels, one per row. Bypasses the clustering step. |
| `clustering` | `"hierarchical"` | `"hierarchical"` uses Ward linkage, matching the paper. `"kmeans"` is also available. |
| `k` | `round(ln N)` | Neighbors in the k-NN graph. |
| `alpha` | `(√5−1)/2` | Shrinkage on intra-cluster edge weights. The paper picks the golden ratio conjugate for interpretability rather than by tuning. |
| `epsilon` | `1e-3` | Floor on the curvature denominator. See the note below. |
| `layout` | `"kamada_kawai"` | `"MDS"` stops at the classical-MDS warm start, which is much cheaper. |

## A note on the curvature ordering

The paper reads the curvature as "interior points low, boundary points high", but that only holds
while β is small. Both φ and ψ vanish for an interior point as β grows — a neighbourhood of a single
label makes the Gibbs weights degenerate — and by β ≈ 8 the second-order term has underflowed to
exactly zero, leaving the curvature pinned at the top of the range instead of the bottom. Measured on
clustered fixtures, the ordering flips somewhere around β ≈ 2, and well-separated clusters push β to
around 11.

The formulas are implemented as published in either regime, and the spanning tree stays faithful to
the clusters in both, because `alpha` and the sparsity of the k-NN graph govern which edges are even
available far more than the curvature ordering does. It is worth knowing before reading too much into
[`curvature`](/api/classes/MINFOTree) directly.

## Cost

All-pairs shortest paths over the tree needs O(N²) memory and the layout is O(N²) per iteration, so
this is meant for up to a few thousand points — the same ceiling as [ISOMAP](./isomap), and the one
the paper states for the original.

## Reference

Levada, A. L. M., *Minimum Information Trees for High Dimensional Data Visualization in Clustering*,
IEEE Access, 2025. [doi:10.1109/ACCESS.2025.3602730](https://doi.org/10.1109/ACCESS.2025.3602730)

See also the [`MINFOTree` API reference](/api/classes/MINFOTree).
