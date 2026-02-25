# Shape Preserver

This example illustrates the difference between algorithms that focus on **Local Probabilities** vs. **Topological Preservation**.

When visualizing complex, disjoint, or non-convex manifolds (like the "Two Moons" dataset), standard techniques like t-SNE can sometimes "shatter" the continuous structure as they try to flatten the points into 2D based on local similarity.

<script setup>
import TopologyComparison from '../components/TopologyComparison.vue'
</script>

<TopologyComparison />

## Comparison Details

- **t-SNE**: Excels at revealing fine-grained clusters, but it doesn't strictly preserve the global or topological shape. In the Moons dataset, it often separates the points into clusters that lose the continuous "arc" shape.
- **TopoMap**: A specialized algorithm that guarantees the preservation of **0-dimensional topology**. It uses Minimum Spanning Trees (MST) to ensure that the skeletal structure of the manifold stays intact, resulting in a cleaner preservation of the continuous arcs.

Choosing the right algorithm depends on whether you care more about **density/sub-clusters** (use t-SNE/UMAP) or the **underlying structure and continuity** of your data (use TopoMap).

## How-to (Code)

`TopoMap` provides a simple interface consistent with other DR algorithms in DruidJS. It is particularly effective for high-dimensional data that forms continuous shapes. You can find out more in the [`TopoMap` API Reference](/api/classes/TopoMap).

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// Initialize TopoMap
const topomap = new druid.TopoMap(data);

// Project into 2D
// Topomap is non-iterative, so it returns the result directly
const projection = topomap.transform();
```
