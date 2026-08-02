# DruidJS Showcases

Explore the capabilities of DruidJS through our interactive showcases. We've organized these examples into specialized categories to demonstrate different aspects of the library.

## Explore by Category

- [**Standard Projections**](/showcase/projections)  
  A gallery of classic dimensionality reduction algorithms on the Iris dataset. Include methods like [PCA](/dimred/pca), [t-SNE](/dimred/tsne), [UMAP](/dimred/umap), [PaCMAP](/dimred/pacmap), [LocalMAP](/dimred/localmap), [MDS](/dimred/mds), [TriMap](/dimred/trimap), [ISOMAP](/dimred/isomap), [TopoMap](/dimred/topomap), [Sammon](/dimred/sammon), [LLE](/dimred/lle), [FastMap](/dimred/fastmap), [SMACOF](/dimred/smacof), [StressMDS](/dimred/stressmds), [KKMDS](/dimred/kkmds), [LDA](/dimred/lda), [LSP](/dimred/lsp), [LTSA](/dimred/ltsa), [SQDMDS](/dimred/sqdmds).

- [**Clustering Pipelines**](/showcase/clustering)  
  Learn how to build pipelines using powerful clustering algorithms (including K-Means, OPTICS, CURE, and Hierarchical Clustering with single, average, complete or Ward linkage) combined with UMAP, PaCMAP, or LocalMAP to project high-dimensional structure.

- [**Hierarchical Dendrograms**](/showcase/dendrogram)  
  A visual interactive representation of the internal recursive nested-cluster logic powering the [**Hierarchical Clustering**](/api/classes/HierarchicalClustering) algorithm, rendered natively in D3.

- [**Cluster Diagnostics**](/showcase/diagnostics)  
  Interrogate a clustering rather than eyeball it. [MINFOTree](/dimred/minfotree) returns a spanning tree, so you can count how many edges cross between groups to judge whether they are really separated, and trace the *unique* path between any two points to see what lies on the boundary.

- [**Interactive Optimization**](/showcase/optimization)  
  Watch iterative algorithms (t-SNE, UMAP, PaCMAP, LocalMAP, TriMap, Sammon, SMACOF, StressMDS, SQDMDS) in action on synthetic and real datasets, and see how they optimize embeddings over time.

- [**Metric Sensitivity**](/showcase/metrics)  
  Discover how the choice of distance metrics (Euclidean, Cosine, Manhattan, Chebyshev) influences the final spatial layout in a Sammon projection.

- [**Topological Preservation**](/showcase/topology)  
  A 3D scan of a mammoth skeleton, projected to 2D. Because you can see what the answer should be, it is obvious when t-SNE tears the legs off — and equally obvious that [TopoMap](/dimred/topomap) and the [MINFO Tree](/dimred/minfotree) keep them attached.

- [**Earth Mover Analysis**](/showcase/distribution)  
  Compare non-negative distribution histograms (Gaussian, Uniform, Bimodal, Exponential) using the Wasserstein distance metric and MDS to project complex 1D distributions as 2D spatial coordinates.

- [**Global Projections**](/showcase/geospatial)  
  Geospatial analysis using the Haversine metric on world city coordinates, charting spherical data down to flat manifolds while retaining geographical distances.
- [**KNN Image Search**](/showcase/search)  
  Fast approximate similarity search (via an HNSW index graph) traversing a subset of the MNIST handwritten digits dataset in its native 784-dimensional space.

- [**Automatic Discovery**](/showcase/discovery)  
  Let X-Means automatically find the optimal number of clusters for a randomized blob dataset, dynamically tuning to the Bayesian Information Criterion without manual parameter nudging.

---

## About the Data

These showcases use various datasets provided by the [@saehrimnir/mistle](https://github.com/saehm/mistle) library, including:

- **Iris**: Classic botanical measurements of Iris flowers.
- **Wine**: Chemical analysis of wines from a specific region.
- **Penguins**: Body measurements of various penguin species.
- **MNIST**: High-dimensional representations of pixel-based handwritten digits.
- **Swiss Roll / Moons / Blobs**: Synthetic datasets designed to test specific nonlinear manifold structures and clustering characteristics.
- **Mammoth**: A 3D scan of a mammoth skeleton, labelled by body part — a rare case where the structure a projection ought to preserve is visible to the naked eye.
