# Clustering Pipelines

A frequent task in data science is identifying clusters in high-dimensional space and then visualizing the learned structural manifolds. This showcase allows you to run multiple clustering algorithms—such as [**K-Means**](/api/classes/KMeans), [**OPTICS**](/api/classes/OPTICS), [**CURE**](/api/classes/CURE), and [**Hierarchical Clustering**](/api/classes/HierarchicalClustering)—on various datasets. It computes the assignments intrinsically and then uses [**UMAP**](/api/classes/UMAP), [**PaCMAP**](/api/classes/PaCMAP), or [**LocalMAP**](/api/classes/LocalMAP) for the 2D visual projection.

<script setup>
import ClusteringPipeline from '../components/ClusteringPipeline.vue'
</script>

# Clustering Pipelines

A frequent task in data science is identifying clusters in high-dimensional space and then visualizing the learned structural manifolds. This showcase allows you to run multiple clustering algorithms—such as [**K-Means**](/api/classes/KMeans), [**OPTICS**](/api/classes/OPTICS), [**CURE**](/api/classes/CURE), and [**Hierarchical Clustering**](/api/classes/HierarchicalClustering)—on various datasets. It computes the assignments intrinsically and then uses [**UMAP**](/api/classes/UMAP) for the 2D visual projection.

[**Hierarchical Clustering**](/api/classes/HierarchicalClustering) offers four linkage criteria. **Ward** minimizes within-cluster variance and tends to produce compact, similarly-sized groups, which suits the roughly spherical structure of Iris or Wine; **single** link chains along density and handles elongated shapes; **complete** and **average** sit between them. Ward's merge heights run on a larger scale than the others, so the tree cut needs a correspondingly larger value — around 10 for three clusters on Iris.

<ClusteringPipeline />

## How-to (Code)

Combining clustering and dimensionality reduction is a common pattern for exploratory data analysis. You can use any of DruidJS's clustering algorithms for the assignments and a technique like UMAP, PaCMAP, or LocalMAP for visualization.

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional float array ... */
];

// 1. Run a clustering algorithm (e.g., Hierarchical Clustering)
// linkage: "single" | "average" | "complete" | "ward"
const hc = new druid.HierarchicalClustering(data, { linkage: "ward" });
const clusters = hc.get_cluster_list(10.0); // Cut the dendrogram tree at distance 10.0

// Or use OPTICS
// const optics = new druid.OPTICS(data, { epsilon: 1.0, min_points: 4 });
// const clusters = optics.get_cluster_list();

// 2. Run UMAP for visualization
const umap = new druid.UMAP(data);
const projection = umap.transform();

// Combine `clusters` and `projection` for plotting...
```
