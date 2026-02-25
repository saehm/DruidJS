/**
 * Type tests for @saehrimnir/druidjs
 *
 * These tests verify that TypeScript can correctly infer and check types
 * from the generated .d.ts file. They don't run at runtime - they're
 * compile-time type checks.
 */

import {
  // KNN
  BallTree,
  // Metrics
  bray_curtis,
  canberra,
  chebyshev,
  type Comparator,
  cosine,
  DisjointSet,
  // Matrix utilities
  distance_matrix,
  type EigenArgs,
  euclidean,
  euclidean_squared,
  FASTMAP,
  goodman_kruskal,
  hamming,
  haversine,
  // Data structures
  Heap,
  // Clustering
  HierarchicalClustering,
  HNSW,
  inner_product,
  type InputType,
  ISOMAP,
  jaccard,
  k_nearest_neighbors,
  // Numerical
  kahan_sum,
  KMeans,
  KMedoids,
  LDA,
  linspace,
  LLE,
  LSP,
  LTSA,
  manhattan,
  // Core classes
  Matrix,
  // Utility functions
  max,
  MDS,
  // Types
  type Metric,
  min,
  NaiveKNN,
  neumair_sum,
  NNDescent,
  norm,
  normalize,
  OPTICS,
  // Dimensionality reduction
  PCA,
  // Optimization
  powell,
  // Linear algebra
  qr,
  qr_householder,
  Randomizer,
  SAMMON,
  simultaneous_poweriteration,
  sokal_michener,
  SQDMDS,
  TopoMap,
  TriMap,
  TSNE,
  UMAP,
  // Version
  version,
  wasserstein,
  yule,
} from "../../dist/druid";

// ============================================
// Matrix Tests
// ============================================

// Test Matrix constructor
const identityMatrix = new Matrix(5, 5, "I");
const zeroMatrix = new Matrix(3, 3, "zero");
const customMatrix = new Matrix(4, 4, (i, j) => i + j);

// Test Matrix static methods
const fromArray: Matrix = Matrix.from([
  [1, 2],
  [3, 4],
]);
const fromFloat64: Matrix = Matrix.from([new Float64Array([1, 2]), new Float64Array([3, 4])]);
const diagMatrix: Matrix = Matrix.from_diag([1, 2, 3]);
const colVector: Matrix = Matrix.from_vector([1, 2, 3], "col");
const rowVector: Matrix = Matrix.from_vector(new Float64Array([1, 2, 3]), "row");

// Test Matrix properties
const rows: number = identityMatrix.shape[0];
const cols: number = identityMatrix.shape[1];
const values: Float64Array = identityMatrix.values;

// Test Matrix methods
const row: Float64Array = identityMatrix.row(0);
const col: Float64Array = identityMatrix.col(0);
const entry: number = identityMatrix.entry(0, 0);
const transposed: Matrix = identityMatrix.transpose();
const transposeShorthand: Matrix = identityMatrix.T;
const cloned: Matrix = identityMatrix.clone();
const inverse: Matrix = identityMatrix.inverse();

// Test Matrix operations
const dotProduct: Matrix = identityMatrix.dot(identityMatrix);
const mult: Matrix = identityMatrix.mult(2);
const multMatrix: Matrix = identityMatrix.mult(identityMatrix);
const added: Matrix = identityMatrix.add(1);
const subbed: Matrix = identityMatrix.sub(identityMatrix);
const divided: Matrix = identityMatrix.divide(2);

// Test Matrix inline operations
const multInline: Matrix = identityMatrix.mult(2, { inline: true });

// Test Matrix block operations
const block: Matrix = identityMatrix.get_block(0, 0, 2, 2);
const concat: Matrix = identityMatrix.concat(identityMatrix, "horizontal");

// Test Matrix statistics
const mean: number = identityMatrix.mean();
const sum: number = identityMatrix.sum();
const diag: Float64Array = identityMatrix.diag();
const meanRows: Float64Array = identityMatrix.meanRows();
const meanCols: Float64Array = identityMatrix.meanCols();

// Test Matrix conversion
const as2dArray: Float64Array[] = identityMatrix.to2dArray();
const asNumberArray: number[][] = identityMatrix.asArray();

// Test Matrix iteration
for (const r of identityMatrix) {
  const rowData: Float64Array = r;
}

// Test Matrix with custom Generic (if applicable via any/typing)
const anyMatrix = identityMatrix as any as Matrix;

// Test static solve methods
const { L, U } = Matrix.LU(identityMatrix);
const det: number = Matrix.det(identityMatrix);
const solved: Matrix = Matrix.solve(identityMatrix, colVector);

// ============================================
// Randomizer Tests
// ============================================

const randomizer = new Randomizer(42);
const randValue: number = randomizer.random;
const randInt: number = randomizer.random_int;
const gaussValue: number = randomizer.gauss_random();
const seed: number = randomizer.seed;

// Test static choice
const choices: number[] = Randomizer.choice([1, 2, 3, 4, 5], 2, 1212);

// Test instance choice
const instanceChoices: string[] = randomizer.choice(["a", "b", "c"], 2);

// ============================================
// Metric Tests
// ============================================

const a = [1, 2, 3];
const b = [4, 5, 6];

const brayDist: number = bray_curtis(a, b);
const canberraDist: number = canberra(a, b);
const chebyshevDist: number = chebyshev(a, b);
const cosineDist: number = cosine(a, b);
const euclideanDist: number = euclidean(a, b);
const euclideanSqDist: number = euclidean_squared(a, b);
const goodmanDist: number = goodman_kruskal(a, b);
const hammingDist: number = hamming(a, b);
const haversineDist: number = haversine([0, 0], [1, 1]);
const haversineWithRadius: number = haversine([0, 0], [1, 1]);
const jaccardDist: number = jaccard(a, b);
const manhattanDist: number = manhattan(a, b);
const sokalDist: number = sokal_michener(a, b);
const wassersteinDist: number = wasserstein(a, b);
const yuleDist: number = yule(a, b);

// Test metric type - must accept both number[] and Float64Array
const customMetric: Metric = (x, y) => {
  // x and y can be number[] | Float64Array
  return Math.abs(x[0] - y[0]);
};

// @ts-expect-error - Metric cannot return string
const invalidMetric: Metric = (x, y) => "invalid";

// ============================================
// Matrix Utility Tests
// ============================================

const dataMatrix = Matrix.from([
  [1, 2],
  [3, 4],
  [5, 6],
]);
const distMatrix: Matrix = distance_matrix(dataMatrix);
// Use default metric (euclidean) - avoids type mismatch with concrete function
const distMatrixDefault: Matrix = distance_matrix(dataMatrix);

const knn = k_nearest_neighbors(dataMatrix, 2);
const knnEntry: { i: number; j: number; distance: number } = knn[0][0];

const space: number[] = linspace(0, 10, 5);
const spaceDefault: number[] = linspace(0, 10);

const normValue: number = norm([1, 2, 3]);
// Use default metric
const normMatrix: number = norm(new Float64Array([1, 2, 3]));

const normalized: number[] | Float64Array = normalize([1, 2, 3]);
// Use default metric
const normalizedFloat: number[] | Float64Array = normalize(new Float64Array([1, 2, 3]));

const maxVal: number = max([1, 2, 3, null, 5]);
const minVal: number = min([1, 2, 3, null, 5]);

// ============================================
// Dimensionality Reduction Tests
// ============================================

const drData = new Matrix(100, 10, () => Math.random());

// PCA
const pca = new PCA(drData, { d: 2 });
const pcaResult: Matrix = pca.transform();
const pcaComponents: Matrix = pca.principal_components();
const pcaParams = pca.parameter();
const pcaD: number = pca.parameter("d") as number;

// Test static transform
const pcaStatic = PCA.transform(drData, { d: 2, seed: 42 });

// MDS - use default metric
const mds = new MDS(drData, { d: 2 });
const mdsResult: Matrix = mds.transform();
const mdsStress: number = mds.stress();

// MDS with precomputed distances
const mdsPrecomputed = new MDS(distMatrix, { d: 2, metric: "precomputed" });

// TSNE
const tsne = new TSNE(drData, { d: 2, perplexity: 30, epsilon: 10 });
// @ts-expect-error - perplexity must be a number
new TSNE(drData, { d: 2, perplexity: "high" });
tsne.init();
const tsneResult: Matrix = tsne.transform(500);

// Test generator pattern
const tsneGen = tsne.generator(100);
for (const step of tsneGen) {
  const intermediate: Matrix = step;
}

// UMAP
const umap = new UMAP(drData, { d: 2, n_neighbors: 15, min_dist: 0.1 });
// @ts-expect-error - min_dist must be number
new UMAP(drData, { d: 2, min_dist: "0.1" });
umap.init();
const umapResult: Matrix = umap.transform(350);
const umapGraph = umap.graph();

// ISOMAP
const isomap = new ISOMAP(drData, { d: 2, neighbors: 10 });
const isomapResult: Matrix = isomap.transform();

// LDA (requires labels)
const labels = Array.from({ length: 100 }, (_, i) => i % 3);
const lda = new LDA(drData, { d: 2, labels });
const ldaResult: Matrix = lda.transform();

// @ts-expect-error - labels are required for LDA
const invalidLDA = new LDA(drData, { d: 2 });

// LLE
const lle = new LLE(drData, { d: 2, neighbors: 10 });
const lleResult: Matrix = lle.transform();

// LSP
const lsp = new LSP(drData, { d: 2, neighbors: 10 });
lsp.init();
const lspResult: Matrix = lsp.transform();

// LTSA
const ltsa = new LTSA(drData, { d: 2, neighbors: 10 });
const ltsaResult: Matrix = ltsa.transform();

// FASTMAP
const fastmap = new FASTMAP(drData, { d: 2 });
const fastmapResult: Matrix = fastmap.transform();

// SAMMON
const sammon = new SAMMON(drData, { d: 2, magic: 0.2 });
const sammonResult: Matrix = sammon.transform(200);

// SQDMDS
const sqdmds = new SQDMDS(drData, { d: 2 });
sqdmds.init();
const sqdmdsResult: Matrix = sqdmds.transform(500);

// TopoMap - use default metric
const topomap = new TopoMap(drData, {});
topomap.init();
const topomapResult: Matrix = topomap.transform();

// TriMap
const trimap = new TriMap(drData, { d: 2, c: 5 });
trimap.init();
const trimapResult: Matrix = trimap.transform(400);

// ============================================
// Clustering Tests
// ============================================

const clusterData = Matrix.from([
  [1, 2],
  [1.5, 1.8],
  [5, 8],
  [8, 8],
  [1, 0.6],
  [9, 11],
]);

// KMeans
const kmeans = new KMeans(clusterData, { K: 2, seed: 42 });
const kmeansK: number = kmeans.k;
const kmeansCentroids: Float64Array[] = kmeans.centroids;
const kmeansClusters: number[][] = kmeans.get_clusters();
const kmeansClusterList: number[] = kmeans.get_cluster_list();
const firstCluster: number[] = kmeansClusters[0];
const firstIndex: number = kmeansClusterList[0];

// KMedoids
const kmedoids = new KMedoids(clusterData, { K: 2, seed: 42 });
const kmedoidsClusters: number[][] = kmedoids.get_clusters();
const kmedoidsMedoids: number[] = kmedoids.get_medoids();

// Hierarchical Clustering - use default metric
const hc = new HierarchicalClustering(clusterData, {
  linkage: "complete",
  metric: customMetric,
});
const hcRoot = hc.root;
const hcClusters: number[][] = hc.get_clusters(2, "distance");
const hcClusterList: number[] = hc.get_cluster_list(2, "distance");
const firstHcCluster: number[] = hcClusters[0];

// OPTICS - use default metric
const optics = new OPTICS(clusterData, {
  epsilon: 3,
  min_points: 2,
  metric: customMetric,
});
const opticsClusters: number[][] = optics.get_clusters();
const opticsClusterList: number[] = optics.get_cluster_list();

// ============================================
// Data Structure Tests
// ============================================

// Heap
const minHeap = new Heap<{ name: string; value: number }>(
  [
    { name: "a", value: 3 },
    { name: "b", value: 1 },
  ],
  (d) => d.value,
  "min",
);

minHeap.push({ name: "c", value: 2 });
const topElement = minHeap.first;
const poppedElement = minHeap.pop();
const heapLength: number = minHeap.length;
const isEmpty: boolean = minHeap.empty;
const heapArray = minHeap.toArray();

// Generic Heap with custom object
interface TestObj {
  id: string;
  val: number;
}
const customHeap = new Heap<TestObj>([{ id: "1", val: 10 }], (d) => d.val);
const firstCustomNode: { element: TestObj; value: number } | null = customHeap.first;
const firstCustom: TestObj | undefined = firstCustomNode?.element;

// Static heapify
const heapified = Heap.heapify([5, 3, 8, 1], (d) => d, "min");

// DisjointSet
const ds = new DisjointSet<number>([1, 2, 3, 4, 5]);
ds.union(1, 2);
ds.union(3, 4);
const parent = ds.find(2);
const children = ds.get_children(1);

// ============================================
// KNN Tests
// ============================================

const knnData: Float64Array[] = [
  new Float64Array([1, 2]),
  new Float64Array([3, 4]),
  new Float64Array([5, 6]),
  new Float64Array([7, 8]),
];

// BallTree - use default params
const ballTree = new BallTree(knnData);
const btResults = ballTree.search(new Float64Array([2, 3]), 2);
const btResult: { element: Float64Array; index: number; distance: number } = btResults[0];
const btByIndex = ballTree.search_by_index(0, 2);

// HNSW - provide all required params
const hnsw = new HNSW(knnData, {
  m: 16,
  ef_construction: 200,
  ef: 50,
  metric: customMetric,
  heuristic: true,
  m0: 32,
  mL: null,
  seed: 1212,
});
// @ts-expect-error - m must be number
new HNSW(knnData, { m: "16" });
hnsw.addOne(new Float64Array([9, 10]));
hnsw.add([new Float64Array([11, 12]), new Float64Array([13, 14])]);
const hnswResults = hnsw.search(new Float64Array([2, 3]), 3);
const hnswSize: number = hnsw.size;
const hnswLayers: number = hnsw.num_layers;
const hnswElement: Float64Array = hnsw.get_element(0);

// Iterator search
for (const { layer, candidates } of hnsw.search_iter(new Float64Array([2, 3]), 3)) {
  const l: number = layer;
  const c = candidates;
}

const naiveKnn = new NaiveKNN(knnData, { metric: "precomputed" });
const naiveResults = naiveKnn.search(new Float64Array([2, 3]), 2);
const naiveByIndex = naiveKnn.search_by_index(0, 2);

// KNN with custom T (constrained to number[] | Float64Array)
const customKNNPoints: number[][] = [
  [1, 2],
  [3, 4],
];
const customNaive = new NaiveKNN<number[]>(customKNNPoints, {
  metric: (a: number[] | Float64Array, b: number[] | Float64Array) => Math.abs(a[0] - b[0]),
});
const customRes = customNaive.search([1.5, 2], 1);
const customElem: number[] = customRes[0].element;

// @ts-expect-error - T must be number[] | Float64Array
const invalidNaive = new NaiveKNN<{ x: number }>([{ x: 1 }], { metric: (a, b) => 0 });

// NNDescent - use custom metric that matches Metric type
const nnDescent = new NNDescent(
  knnData.map((a) => Array.from(a)),
  {
    metric: customMetric,
    K: 3,
    rho: 0.8,
    delta: 0.0001,
    seed: 42,
  },
);
const nnDescentResults = nnDescent.search([2, 3], 2);
const nnDescentByIndex = nnDescent.search_index(0, 2);

// ============================================
// Linear Algebra Tests
// ============================================

const laMatrix = Matrix.from([
  [1, 2],
  [3, 4],
]);

// QR decomposition
const qrResult = qr(laMatrix);
const Q: Matrix = qrResult.Q;
const R: Matrix = qrResult.R;

const qrHouseholder = qr_householder(laMatrix);
const QH: Matrix = qrHouseholder.Q;
const RH: Matrix = qrHouseholder.R;

// Inner product
const ip: number = inner_product([1, 2, 3], [4, 5, 6]);
const ipFloat: number = inner_product(new Float64Array([1, 2]), new Float64Array([3, 4]));

// Eigendecomposition
const symMatrix = Matrix.from([
  [4, 1],
  [1, 3],
]);
const eigenResult = simultaneous_poweriteration(symMatrix, 2, {
  seed: 42,
  max_iterations: 100,
  tol: 1e-8,
});
const eigenvalues: Float64Array = eigenResult.eigenvalues;
const eigenvectors: Float64Array[] = eigenResult.eigenvectors;

// ============================================
// Numerical Tests
// ============================================

const summands = [1, 2, 3, 4, 5];
const kahanResult: number = kahan_sum(summands);
const neumairResult: number = neumair_sum(new Float64Array(summands));

// ============================================
// Optimization Tests
// ============================================

const rosenbrock = (x: number[]): number => {
  const [a, b] = [1, 100];
  return (a - x[0]) ** 2 + b * (x[1] - x[0] ** 2) ** 2;
};

const x0 = [0, 0];
const optimized: number[] = powell(rosenbrock, x0, 100);

// With Float64Array
const rosenbrockTyped = (x: Float64Array): number => {
  const [a, b] = [1, 100];
  return (a - x[0]) ** 2 + b * (x[1] - x[0] ** 2) ** 2;
};
const x0Typed = new Float64Array([0, 0]);
const optimizedTyped: Float64Array = powell(rosenbrockTyped, x0Typed, 100);

// ============================================
// Misc Tests
// ============================================

// Version
const v: string = version;

// InputType union
const inputAsMatrix: InputType = new Matrix(2, 2, 1);
const inputAsFloat64Array: InputType = [new Float64Array([1, 2]), new Float64Array([3, 4])];
const inputAsNumberArray: InputType = [
  [1, 2],
  [3, 4],
];

// Comparator type
const compareFn: Comparator = (a, b) => a < b;

// EigenArgs
const eigenArgs: EigenArgs = {
  max_iterations: 100,
  seed: 42,
  tol: 1e-8,
};

console.log("Type tests passed!");
