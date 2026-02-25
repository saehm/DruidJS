# DruidJS Showcases

Explore the capabilities of DruidJS through our interactive showcases. We've organized these examples into specialized categories to demonstrate different aspects of the library.

## Explore by Category

- [**Standard Projections**](/showcase/projections)  
  A gallery of classic dimensionality reduction algorithms on the Iris dataset. Include methods like [PCA](/dimred/pca), [t-SNE](/dimred/tsne), [UMAP](/dimred/umap), [MDS](/dimred/mds), [TriMap](/dimred/trimap), [ISOMAP](/dimred/isomap), [TopoMap](/dimred/topomap), [Sammon](/dimred/sammon), [LLE](/dimred/lle), [FastMap](/dimred/fastmap), [SMACOF](/dimred/smacof), [LDA](/dimred/lda), [LSP](/dimred/lsp), [LTSA](/dimred/ltsa), [SQDMDS](/dimred/sqdmds).

- [**Clustering Pipelines**](/showcase/clustering)  
  Learn how to build pipelines using powerful clustering algorithms (including K-Means, OPTICS, CURE, and Hierarchical Clustering) combined with UMAP to project high-dimensional structure.
- [**Hierarchical Dendrograms**](/showcase/dendrogram)  
  A visual interactive representation of the internal recursive nested-cluster logic powering the [**Hierarchical Clustering**](/api/classes/HierarchicalClustering) algorithm, rendered natively in D3.

- [**Interactive Optimization**](/showcase/optimization)  
  Watch iterative algorithms (t-SNE, UMAP, TriMap, Sammon, SMACOF, SQDMDS) in action on synthetic and real datasets, and see how they optimize embeddings over time.

- [**Metric Sensitivity**](/showcase/metrics)  
  Discover how the choice of distance metrics (Euclidean, Cosine, Manhattan, Chebyshev) influences the final spatial layout in a Sammon projection.

- [**Topological Preservation**](/showcase/topology)  
  A deep dive into TopoMap and its unique approach to preserving data 0-dimensional topological connectivity on concave datasets like Moons.

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
