<script setup>
import MINFOTreeShowcase from '../components/MINFOTreeShowcase.vue'
</script>

# Cluster Diagnostics

You have run a clustering. Are the groups real, which of them actually touch, and what sits on the
boundary between them? A scatterplot of a projection is a poor witness here — it has already thrown
away structure to fit on the page, and clusters can look separated because the projection separated
them.

[MINFOTree](/dimred/minfotree) is built for this question. It returns a **spanning tree** over the
points, and the tree can be queried directly rather than eyeballed.

<MINFOTreeShowcase />

## Is the clustering clean?

A tree spanning *c* clusters needs at least *c − 1* edges to join them. Count how many it actually
uses: close to the minimum means the groups are genuinely separated in the data; many more means they
interleave. On Iris the tree uses exactly 2 for 3 clusters — as clean as it can be.

This is a verdict on the **clustering**, not on the picture, because it is computed on the tree
rather than on the 2D coordinates. Toggle **Highlight inter-cluster edges** to see where the groups
meet.

## What lies between two groups?

Click any two points. A tree has exactly one route between them, so the highlighted path is not a
choice among many — it is *the* sequence of observations the data passes through to get from one to
the other.

That answers questions a projection cannot. On Iris, a path from a *setosa* to a *virginica* is
forced through *versicolor*: the tree makes the ordering of the species explicit instead of leaving
you to infer it from a layout. Practical uses of the same idea:

- **Interpolation by example** — real data points between two states, not synthesized midpoints.
- **Interpreting a classifier** — a path crossing a boundary repeatedly marks genuinely interleaved
  labels, not merely adjacent ones on the page.
- **Finding transitional cases** — the points sitting on inter-cluster edges are the ambiguous ones.

## Which points are on the boundary?

Switch **Colour by** to *Information curvature*. Every point gets a score saying how comfortably it
sits inside its own cluster, and the legend under the plot gives the scale.

The score comes from reading the labels as a [Potts
model](https://en.wikipedia.org/wiki/Potts_model) — a physical model of neighbours that prefer to
agree, the same idea as magnetic spins aligning. Fitting it to your data yields one number, the
inverse temperature **β**, shown in the legend: high β means a point's label is almost completely
determined by its neighbours (clean, well-separated clusters), low β means the labels are only
loosely tied to the neighbourhood (interleaved or noisy ones). Each point is then scored on two
quantities:

- **How surprised the model is by this point** — the gap between the label it has and the label its
  neighbours predict.
- **How stable its neighbourhood is** — whether the surrounding labels are decisive or finely
  balanced.

The curvature is the ratio of the second to the first. Points deep inside a cluster and points on a
frontier land at opposite ends of it, which is exactly the split you want when hunting for
transitional cases.

::: warning Which end is which is not fixed
It depends on β, and it is **not** always "interior is low" — on all three datasets here, interior
points come out at the *top* of the scale. Rather than ask you to remember the rule, the legend
measures it: it compares the points whose tree neighbours all share their label against those whose
do not, and labels the two ends accordingly for the data on screen. The
[method page](/dimred/minfotree#a-note-on-the-curvature-ordering) explains why the direction moves,
and why that follows from the published formulas rather than contradicting them.
:::

The two marks on the legend bar are where each population actually sits: points whose tree
neighbours all share their label, and points whose do not. **The gap between them is the reading.**
On **Blobs**, where the groups are genuinely separate, the marks are far apart. On **Iris**, where
*versicolor* and *virginica* overlap, they nearly touch — the curvature is telling you the
clustering is less clear-cut, before you have looked at a single coordinate.

The bar's own endpoints are not worth reading: the curvature is rescaled to a fixed range on every
run, so the extremes always say the same thing. Where the populations fall between them does not.

## When to reach for it

Use MINFOTree when you already have a clustering and want to *inspect* it. Do not use it as a
general-purpose projection: the tree inherits whatever the clustering got wrong, which is why
`clusters` or `labels` is required rather than defaulted. For laying out data you have no labels for,
see [Standard Projections](/showcase/projections).

## Code

```javascript
import * as druid from "@saehrimnir/druidjs";

const minfo = new druid.MINFOTree(data, { clusters: 3 });
const projection = minfo.transform();
const edges = minfo.edges; // [[u, v, weight], ...]
const labels = minfo.labels;

// Build an adjacency list from the tree
const adjacency = Array.from({ length: data.length }, () => []);
for (const [u, v] of edges) {
  adjacency[u].push(v);
  adjacency[v].push(u);
}

// The unique path between two points is a plain breadth-first walk —
// no shortest-path weighting needed, because a tree has only one route.
function path(from, to) {
  const previous = new Map([[from, -1]]);
  const queue = [from];
  while (queue.length) {
    const u = queue.shift();
    if (u === to) break;
    for (const v of adjacency[u]) {
      if (previous.has(v)) continue;
      previous.set(v, u);
      queue.push(v);
    }
  }
  const route = [];
  for (let at = to; at !== -1; at = previous.get(at)) route.push(at);
  return route.reverse();
}

// How many cluster boundaries does the route cross?
const route = path(0, 100);
const crossings = route.filter((v, i) => i > 0 && labels[v] !== labels[route[i - 1]]).length;

// How separated is the clustering overall?
const bridges = edges.filter(([u, v]) => labels[u] !== labels[v]).length;
const minimum = new Set(labels).size - 1;
```

The spanning tree utility is exported on its own too, if you want to build one over your own weighted
graph:

```javascript
const tree = druid.minimum_spanning_tree(edges, N); // Kruskal, returns a forest if disconnected
```
