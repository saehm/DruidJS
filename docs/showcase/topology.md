<script setup>
import TopologyComparison from '../components/TopologyComparison.vue'
</script>

# Shape Preserver

This example illustrates the difference between algorithms that focus on **Local Probabilities** vs. **Topological Preservation**.

Most projections are hard to grade because you cannot see what the answer should be. This one you can: the data is a 3D scan of a **mammoth skeleton**, so the structure a projection is supposed to keep — four legs attached to a body, tusks attached to a head, a tail at the far end — is something you already recognise. Every point carries one of 11 body-part labels, which makes it obvious when a limb has been torn off.

<TopologyComparison />

## Comparison Details

- **t-SNE**: Excels at revealing fine-grained clusters, but it doesn't strictly preserve the global or topological shape. Here it shatters the skeleton — the legs float away as separate islands, and nothing tells you they were ever attached.
- **TopoMap**: A specialized algorithm that guarantees the preservation of **0-dimensional topology**. It uses Minimum Spanning Trees (MST) to ensure that the skeletal structure of the manifold stays intact, so the limbs stay joined to the body even though the silhouette is not reproduced.
- **[MINFO Tree](/dimred/minfotree)**: Also built on a spanning tree, but it *keeps* the tree as the output instead of using it only to place points. Its edges are weighted by an information-geometric curvature derived from a Potts model over the labels, so the tree runs along each body part and crosses between them as rarely as the anatomy allows.

Both TopoMap and MINFOTree use a minimum spanning tree, but for different ends. TopoMap uses it to guarantee no topological intersections while laying points out; MINFOTree treats the tree as the artefact you actually inspect. Choosing between them depends on whether you want a **projection** whose topology is trustworthy (TopoMap) or a **skeleton** you can query directly (MINFOTree) — see [Cluster Diagnostics](/showcase/diagnostics) for what the second buys you.

Choosing the right algorithm depends on whether you care more about **density/sub-clusters** (use t-SNE/UMAP) or the **underlying structure and continuity** of your data (use TopoMap or MINFOTree).

::: tip Reading the panels
The first panel is one axis-aligned view of the 3D scan, not a projection — it is there so the target is visible. A method is doing well when body parts stay contiguous and adjacent parts stay adjacent, **not** when the picture looks like a mammoth. None of these methods try to reproduce the outline.
:::

## How-to (Code)

`TopoMap` provides a simple interface consistent with other DR algorithms in DruidJS. It is particularly effective for high-dimensional data that forms continuous shapes. You can find out more in the [`TopoMap` API Reference](/api/classes/TopoMap).

```javascript
import * as druid from "@saehrimnir/druidjs";
import * as mistle from "@saehrimnir/mistle";

const { values, labels } = mistle.MAMMOTH; // 1200 points, 3D, 11 body parts

// Initialize TopoMap
const topomap = new druid.TopoMap(values);

// Project into 2D
// Topomap is non-iterative, so it returns the result directly
const projection = topomap.transform();

// MINFOTree, given the body parts as its labels field
const minfo = new druid.MINFOTree(values, { labels });
const layout = minfo.transform();
const edges = minfo.edges; // the spanning tree itself
```
