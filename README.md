## DruidJS — A JavaScript Library for Dimensionality Reduction.

<a href="#"><img src="https://raw.githubusercontent.com/saehm/DruidJS/refs/heads/master/icon.svg" width=80 align="left" hspace="10" vspace="6"></a>

Project high-dimensional data down to two or three dimensions while preserving the structure each method is designed to keep.
DruidJS gives you 20 dimensionality reduction methods, 7 nearest-neighbor indices and 7 clustering algorithms behind one API, with SIMD-accelerated WebAssembly kernels and a pure-JavaScript fallback everywhere.

<br>

![Codecov](https://img.shields.io/codecov/c/github/saehm/DruidJS)
![NPM Downloads](https://img.shields.io/npm/dw/%40saehrimnir%2Fdruidjs)
[![CI](https://img.shields.io/github/actions/workflow/status/saehm/DruidJS/ci.yml?branch=master&label=tests)](https://github.com/saehm/DruidJS/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/saehm/DruidJS)](https://raw.githubusercontent.com/saehm/DruidJS/refs/heads/master/LICENCE)
[![DOI:10.1109/VIS47514.2020.00029](https://img.shields.io/badge/DOI-10.1109%2FVIS47514.2020.00029-blue)](https://doi.org/10.1109/VIS47514.2020.00029)

### Installation

If you use npm, install with `npm install @saehrimnir/druidjs`, and use it with

```js
import * as druid from "@saehrimnir/druidjs";
```

Otherwise download the files [here](https://github.com/saehm/DruidJS/releases/latest), or use for instance [unpkg](https://unpkg.com/@saehrimnir/druidjs) this way:

```html
<script src="https://unpkg.com/@saehrimnir/druidjs"></script>
```

### Quick start

```js
import * as druid from "@saehrimnir/druidjs";

const X = druid.Matrix.from(data); // n ⨯ d
const Y = new druid.UMAP(X, { d: 2, seed: 1212 }).transform();

Y.to2dArray(); // [[x, y], ...] — ready for d3
```

Every method has the same shape: construct it with the data and a parameters object, then call `transform()`.
For a one-off projection there is a static shorthand:

```js
const Y = druid.PCA.transform(X, { d: 2 });
```

The iterative methods can be stepped instead of run to completion, so you can draw the embedding while it converges:

```js
for (const Y of new druid.TSNE(X, { d: 2, perplexity: 30 }).generator()) {
    draw(Y.to2dArray());
}
```

### Methods

**Dimensionality reduction** —
[PCA](https://saehm.github.io/DruidJS/dimred/pca) ·
[LDA](https://saehm.github.io/DruidJS/dimred/lda) ·
[MDS](https://saehm.github.io/DruidJS/dimred/mds) ·
[SMACOF](https://saehm.github.io/DruidJS/dimred/smacof) ·
[StressMDS](https://saehm.github.io/DruidJS/dimred/stressmds) ·
[KKMDS](https://saehm.github.io/DruidJS/dimred/kkmds) ·
[SQDMDS](https://saehm.github.io/DruidJS/dimred/sqdmds) ·
[t-SNE](https://saehm.github.io/DruidJS/dimred/tsne) ·
[UMAP](https://saehm.github.io/DruidJS/dimred/umap) ·
[TriMap](https://saehm.github.io/DruidJS/dimred/trimap) ·
[PaCMAP](https://saehm.github.io/DruidJS/dimred/pacmap) ·
[LocalMAP](https://saehm.github.io/DruidJS/dimred/localmap) ·
[Sammon](https://saehm.github.io/DruidJS/dimred/sammon) ·
[ISOMAP](https://saehm.github.io/DruidJS/dimred/isomap) ·
[LLE](https://saehm.github.io/DruidJS/dimred/lle) ·
[LTSA](https://saehm.github.io/DruidJS/dimred/ltsa) ·
[LSP](https://saehm.github.io/DruidJS/dimred/lsp) ·
[FastMap](https://saehm.github.io/DruidJS/dimred/fastmap) ·
[TopoMap](https://saehm.github.io/DruidJS/dimred/topomap) ·
[MINFOTree](https://saehm.github.io/DruidJS/dimred/minfotree)

**Nearest neighbors** — BallTree, KDTree, HNSW, Annoy, LSH, NNDescent, NaiveKNN

```js
const index = new druid.HNSW(points, { seed: 1212 });
index.search(points[0], 10); // the 10 nearest points
```

**Clustering** — KMeans, KMedoids, XMeans, OPTICS, CURE, MeanShift, HierarchicalClustering (single, complete, average and Ward linkage)

### Performance

The hot paths run through WebAssembly kernels once the input is large enough to pay for crossing the boundary, and fall back to JavaScript otherwise — automatically, with no difference in the API.
Measured between the published 0.8.0 and 0.9.0 packages: **2.35⨯ on dimensionality reduction, 2.24⨯ on nearest-neighbor search, 2.21⨯ on matrix operations.**

The kernels can be turned off, which is the easiest way to compare them against the fallback:

```js
druid.setWasmEnabled(false);
druid.isWasmAvailable();
```

### Matrix

DruidJS uses internally the Matrix class for storing data. You can use it by creating a `druid.Matrix` object for instance with the function `from`, in example:

```js
import * as druid from "@saehrimnir/druidjs";

let data = [[...], [...], ...];
let matrix = druid.Matrix.from(data);
```

You can create a `druid.Matrix` object programmatically by:

```js
let fn = (row, col) => (row == col ? 1 : 0);
let matrix = new druid.Matrix(rows, columns, fn);
```

If `rows == columns`, then `matrix` would be a identity matrix.
A shortcut for a identity matrix is:

```js
let matrix = new druid.Matrix(rows, columns, "I");
// or
let matrix = new druid.Matrix(rows, columns, "identity");
```

There are more shortcuts for creating matrices:

```js
let matrix = new druid.Matrix(3, 3, "zeros"); // matrix would be a 3x3 matrix with zeroes
let matrix = new druid.Matrix(3, 3, "center"); // matrix would be a 3x3 center matrix;
let number = 12;
let matrix = new druid.Matrix(3, 3, number); // matrix would be a 3x3 matrix filled with 'number'
```

If you want to use a `druid.Matrix` object, for instance, with [d3](https://d3js.org), you can use either the `to2dArray` method, the `iterate_rows` generator function, or just use the `druid.Matrix` object as an iterable (works with d3 since version 6).

```js
let data = await d3.csv("data.csv");
let matrix = druid.Matrix.from(data);
d3.selectAll("datapoints").data(matrix.to2dArray()); //...
d3.selectAll("datapoints").data(matrix.iterate_rows()); //...
d3.selectAll("datapoints").data(matrix); //...
```

### Showcases

Live, interactive examples in the documentation:

- [Projections](https://saehm.github.io/DruidJS/showcase/projections) — the DR methods side by side on the Iris dataset
- [Clustering](https://saehm.github.io/DruidJS/showcase/clustering) and [Cluster Diagnostics](https://saehm.github.io/DruidJS/showcase/diagnostics)
- [k-Nearest Neighbors](https://saehm.github.io/DruidJS/showcase/knn) — the indices compared
- [Topology](https://saehm.github.io/DruidJS/showcase/topology) · [Dendrograms](https://saehm.github.io/DruidJS/showcase/dendrogram) · [Geospatial](https://saehm.github.io/DruidJS/showcase/geospatial) · [Image Search](https://saehm.github.io/DruidJS/showcase/search)

[All showcases →](https://saehm.github.io/DruidJS/showcase) · [Every method, one page each →](https://saehm.github.io/DruidJS/dimred)

### Resources

```
@inproceedings{cutura2020druid,
  title={{DRUIDJS — A JavaScript Library for Dimensionality Reduction}},
  author={Cutura, Rene and Kralj, Christoph and Sedlmair, Michael},
  booktitle={2020 IEEE Visualization Conference (VIS)},
  pages={111--115},
  year={2020},
  organization={IEEE}
}
```

- [Documentation](https://saehm.github.io/DruidJS/)
- [API Reference](https://saehm.github.io/DruidJS/api/)
- [Demo](https://renecutura.eu/druid_demo)
- [Hello DruidJS](https://observablehq.com/@saehrimnir/hello-druidjs) on Observable
- [Conference Talk](https://youtu.be/bi6FfsWV_9k?t=2463) [IEEEVIS2020](http://ieeevis.org/year/2020/welcome)
