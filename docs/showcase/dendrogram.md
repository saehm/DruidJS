# Hierarchical Tree (Dendrogram)

Hierarchical clustering builds a recursive tree structure called a Dendrogram. The length of the branches represents the distance at which groups merge, and the leaves represent individual data points.

Unlike partitioning algorithms like K-Means, which just group points into $K$ buckets, hierarchical methods give you a full multiscale representation of cluster relationships.

<script setup>
import DendrogramShowcase from '../components/DendrogramShowcase.vue'
</script>

<DendrogramShowcase />

## How it works

1. **Initialization:** Every point starts as its own cluster (leaf node).
2. **Linkage Criterion:** The algorithm searches for the two closest clusters in the dataset. Distance between clusters can be measured using different **Linkage Strategies**:
   - **Single Link**: The distance between the two _closest_ elements of the clusters. Tends to "chain" points together.
   - **Complete Link**: The distance between the two _farthest_ elements of the clusters. Tends to form compact, spherical clusters.
   - **Average Link**: The average distance between all pairs of elements. A balanced approach.
3. **Merging:** The two closest clusters are merged into a new node. The distance where they merge is recorded as the branch length.
4. **Iterate:** Step 2 and 3 repeat until all points belong to a single root cluster.

## How-to (Code)

To get the full tree structure, access the `root` property of the `HierarchicalClustering` instance.

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize with your preferred linkage
const hc = new druid.HierarchicalClustering(data, {
  linkage: "complete",
});

// 2. Fetch the recursive tree structure
const rootNode = hc.root;

// Example parsing logic for tree traversal
function traverse(node) {
  if (node.isLeaf) {
    console.log(`Leaf point index: ${node.index}`);
  } else {
    console.log(`Merge distance: ${node.dist}`);
    traverse(node.left);
    traverse(node.right);
  }
}

traverse(rootNode);
```
