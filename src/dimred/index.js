/**
 * @module dimred
 */
export { FASTMAP } from "./FASTMAP.js";
export { LocalMAP } from "./LocalMAP.js";
export { PaCMAP } from "./PaCMAP.js";
export { ISOMAP } from "./ISOMAP.js";
export { KKMDS } from "./KKMDS.js";
export { LDA } from "./LDA.js";
export { LLE } from "./LLE.js";
export { LocalMAP } from "./LocalMAP.js";
export { LSP } from "./LSP.js";
export { LTSA } from "./LTSA.js";
export { MDS } from "./MDS.js";
export { MINFOTree } from "./MINFOTree.js";
export { PaCMAP } from "./PaCMAP.js";
export { PCA } from "./PCA.js";
export { SAMMON } from "./SAMMON.js";
export { SMACOF } from "./SMACOF.js";
export { SQDMDS } from "./SQDMDS.js";
export { StressMDS, WEIGHTS_ELASTIC, WEIGHTS_SAMMON, WEIGHTS_UNIFORM } from "./StressMDS.js";
export { TopoMap } from "./TopoMap.js";
export { TriMap } from "./TriMap.js";
export { TSNE } from "./TSNE.js";
export { UMAP } from "./UMAP.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { Matrix } from "../matrix/index.js" */
/** @import { KNN } from "../knn/KNN.js" */
/** @import { EigenArgs } from "../linear_algebra/index.js" */
/** @import { ChooseDR } from "./SAMMON.js"; */

/**
 * @typedef {Object} ParametersLSP
 * @property {number} [neighbors] - number of neighbors to consider.
 * @property {number} [control_points] - number of controlpoints
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersFASTMAP
 * @property {number} [d=2] - The dimensionality of the projection
 * @property {Metric} [metric=euclidean] - The metric which defines the distance between two points.
 * @property {number} [seed=1212] - The seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersISOMAP
 * @property {number} [neighbors] - The number of neighbors ISOMAP should use to project the data.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {"MDS" | "SMACOF"} [project="MDS"] - Whether to use classical MDS or SMACOF for the final DR.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * How to weight each pair in the {@link StressMDS} objective.
 *
 * A number is an exponent `q`, giving `w_ij = d_ij^q` — `0` for raw stress, `-1` for Sammon stress,
 * `-2` for elastic scaling / Kamada-Kawai. A matrix supplies the weights directly. A function is
 * called per pair with the target distance and the two indices.
 *
 * A weight of zero, or any non-finite value, drops that pair from the objective — which is how
 * missing or untrusted dissimilarities are expressed.
 *
 * @typedef {number | Matrix | number[][] | ((d_ij: number, i: number, j: number) => number)} WeightSpec
 */

/**
 * @typedef {Object} ParametersStressMDS
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean] - the metric which defines the distance
 *   between two points. Pass graph shortest-path distances as `"precomputed"` for a graph layout.
 * @property {WeightSpec} [weights=-2] - Pair weighting, see {@link WeightSpec}.
 * @property {number} [iterations=300] - maximum number of gradient steps.
 * @property {number} [epsilon=1e-6] - stop once the relative stress improvement falls below this.
 * @property {number} [learning_rate=0.1] - initial step size. Dimensionless: the gradient is
 *   preconditioned by the weighted degree, so this needs no rescaling for the data or the weighting.
 *   Adapted by the line search, so it only sets where the search starts.
 * @property {"MDS" | "PCA" | "random"} [init_DR="MDS"] - starting configuration. `"MDS"` runs
 *   classical MDS on the same distances, which is what keeps the non-convex descent out of the poor
 *   local minima this objective is known for. `"PCA"` needs the original data, not a precomputed
 *   matrix.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * {@link ParametersStressMDS} without `weights`, which {@link KKMDS} fixes at `-2`.
 *
 * @typedef {Omit<ParametersStressMDS, "weights">} ParametersKKMDS
 */

/**
 * @typedef {Object} ParametersLDA
 * @property {any[] | Float64Array} labels - The labels / classes for each data point.
 * @property {number} [d=2] - The dimensionality of the projection.
 * @property {number} [seed=1212] - The seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersLLE
 * @property {number} [neighbors] - The number of neighbors for LLE.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {KNN<number[] | Float64Array, any> | null} [knn=null] - Index used to find the
 *   neighbors. If `null`, a KDTree or BallTree is built, depending on the metric. Pass an
 *   approximate index such as `HNSW`, `Annoy`, or `NNDescent` to avoid the exact O(N^2) search on
 *   large datasets. Default is `null`
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersLTSA
 * @property {number} [neighbors] - The number of neighbors for LTSA.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {KNN<number[] | Float64Array, any> | null} [knn=null] - Index used to find the
 *   neighbors. If `null`, a KDTree or BallTree is built, depending on the metric. Pass an
 *   approximate index such as `HNSW`, `Annoy`, or `NNDescent` to avoid the exact O(N^2) search on
 *   large datasets. Default is `null`
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersMDS
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersMINFOTree
 * @property {number} [k] - Neighbors in the k-NN graph. Defaults to `round(ln N)`, as in the paper.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [clusters] - Number of clusters to partition `X` into for the labels field.
 *   Required unless `labels` is given — the tree inherits whatever the clustering gets wrong, so
 *   there is no safe default. With the default hierarchical clustering the usable range is
 *   `2 … N-1`, the number of merges a dendrogram cut can separate; `"kmeans"` has no such limit.
 * @property {"hierarchical" | "kmeans"} [clustering="hierarchical"] - How to obtain the labels when
 *   `labels` is not given. `"hierarchical"` uses Ward linkage, matching the paper's experiments.
 * @property {any[] | Float64Array | Int32Array | null} [labels=null] - Precomputed labels, one per
 *   row of `X`. Bypasses the clustering step. Values may be of any type; they are remapped to
 *   `0 … q-1`.
 * @property {number} [alpha] - Shrinkage applied to intra-cluster edge weights. Defaults to the
 *   golden ratio conjugate `(√5-1)/2 ≈ 0.618`, which the paper picks for interpretability rather
 *   than by tuning.
 * @property {number} [epsilon=1e-3] - Floor on the curvature denominator, `S = -ψ / (φ + epsilon)`,
 *   needed because φ and ψ both vanish for interior points as β grows. Default matches the author's
 *   reference implementation. See `MINFOTree._information_curvature` on how the interior/boundary
 *   curvature ordering depends on β.
 * @property {"kamada_kawai" | "MDS"} [layout="kamada_kawai"] - How to lay the tree out. `"MDS"`
 *   stops after the classical-MDS warm start, which is far cheaper and often enough.
 * @property {number} [iterations=300] - Maximum Kamada-Kawai gradient steps.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersPCA
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {number} [seed=1212] - the seed for the random number generator.
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @template {keyof ChooseDR} K
 * @typedef {Object} ParametersSAMMON
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {K} [init_DR="random"] - Either "PCA" or "MDS", with which SAMMON initialiates the projection.
 * @property {ChooseDR[K]} [init_parameters] - Parameters for the "init"-DR method.
 * @property {number} [magic=0.1] - learning rate for gradient descent.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersSMACOF
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [iterations=300] - maximum number of iterations.
 * @property {number} [epsilon=1e-4] - tolerance for stress difference.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersSQDMDS
 * @property {number} [d=2]
 * @property {Metric | "precomputed"} [metric=euclidean]
 * @property {number} [decay_start=0.1] - Percentage of iterations using exaggeration phase.
 * @property {number} [decay_cte=0.34] - Controls the decay of the learning parameter.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef ParametersTopoMap
 * @property {Metric} metric = euclidean - The metric which defines the distance between
 *   two points.
 * @property {number} seed = 1212 - The seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersTriMap
 * @property {number} [weight_temp=0.5] - Temperature of the tempered log applied to the triplet
 *   weights. `1` recovers the ordinary logarithm; lower values compress large weights harder.
 * @property {number} [n_inliers=12] - number of inliers.
 * @property {number} [n_outliers=4] - number of outliers.
 * @property {number} [n_random=3] - number of random points.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {number} [lr=0.1] - learning rate of the delta-bar-delta optimizer.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersPaCMAP
 * @property {number} [n_neighbors=10] - Number of nearest neighbors forming the attractive pairs.
 * @property {number} [MN_ratio=0.5] - Mid-near pairs per point, as a fraction of `n_neighbors`.
 * @property {number} [FP_ratio=2.0] - Further pairs per point, as a multiple of `n_neighbors`.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [lr=1.0] - learning rate of the Adam optimizer.
 * @property {number[]} [num_iters=[100,100,250]] - Iterations in each of the three phases.
 * @property {KNN<number[] | Float64Array, any> | null} [knn=null] - Index used to find the
 *   neighbors. If `null`, an exact blocked search runs instead, which is faster than a tree at the
 *   `n_neighbors + 50` candidates the density rescaling needs. Pass an approximate index such as
 *   `HNSW`, `Annoy`, or `NNDescent` to avoid the O(N^2) search on large datasets. Default is `null`
 * @property {boolean} [apply_pca=true] - Reduce inputs wider than 100 dimensions to 100 via PCA
 *   before the search and the initialization.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersLocalMAP
 * @property {number} [n_neighbors=10] - Number of nearest neighbors forming the attractive pairs.
 * @property {number} [MN_ratio=0.5] - Mid-near pairs per point, as a fraction of `n_neighbors`.
 * @property {number} [FP_ratio=2.0] - Further pairs per point, as a multiple of `n_neighbors`.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [lr=1.0] - learning rate of the Adam optimizer.
 * @property {number[]} [num_iters=[100,100,250]] - Iterations in each of the three phases.
 * @property {number} [low_dist_thres=10] - Embedding distance below which a point may be redrawn as
 *   a further pair in phase 3. Also sets the local attraction scale, `low_dist_thres / 2`.
 * @property {KNN<number[] | Float64Array, any> | null} [knn=null] - Index used to find the
 *   neighbors. If `null`, an exact blocked search runs instead. Default is `null`
 * @property {boolean} [apply_pca=true] - Reduce inputs wider than 100 dimensions to 100 via PCA
 *   before the search and the initialization.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersTSNE
 * @property {number} [perplexity=50] - perplexity.
 * @property {number} [epsilon=10] - learning parameter.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean_squared] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersPaCMAP
 * @property {number} [n_neighbors=10] - Number of nearest neighbors for NN pairs.
 * @property {number} [MN_ratio=0.5] - Ratio of mid-near pairs to n_neighbors.
 * @property {number} [FP_ratio=2.0] - Ratio of further pairs to n_neighbors.
 * @property {number} [d=2] - The dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - The metric which defines the distance between two points.
 * @property {number} [lr=1.0] - Learning rate for the Adam optimizer.
 * @property {number[]} [num_iters=[100,100,250]] - Number of iterations for each of the three phases.
 * @property {"annoy" | "hnsw"} [knn_backend="annoy"] - KNN backend for nearest-neighbor search.
 * @property {Object} [knn_params={}] - Extra options forwarded to the KNN backend constructor (merged with defaults).
 * @property {boolean} [apply_pca=true] - If true and input has >100 dimensions, reduce to 100 dims via PCA before KNN and initialization.
 * @property {number} [seed=1212] - The seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersLocalMAP
 * @property {number} [n_neighbors=10] - Number of nearest neighbors for NN pairs.
 * @property {number} [MN_ratio=0.5] - Ratio of mid-near pairs to n_neighbors.
 * @property {number} [FP_ratio=2.0] - Ratio of further pairs to n_neighbors.
 * @property {number} [d=2] - The dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - The metric which defines the distance between two points.
 * @property {number} [lr=1.0] - Learning rate for the Adam optimizer.
 * @property {number[]} [num_iters=[100,100,250]] - Number of iterations for each of the three phases.
 * @property {number} [low_dist_thres=10] - Distance threshold for local FP pair resampling in phase 3.
 * @property {"annoy" | "hnsw"} [knn_backend="annoy"] - KNN backend for nearest-neighbor search.
 * @property {Object} [knn_params={}] - Extra options forwarded to the KNN backend constructor (merged with defaults).
 * @property {boolean} [apply_pca=true] - If true and input has >100 dimensions, reduce to 100 dims via PCA before KNN and initialization.
 * @property {number} [seed=1212] - The seed for the random number generator.
 */

/**
 * @typedef {Object} ParametersUMAP
 * @property {number} [n_neighbors=15] - size of the local neighborhood.
 * @property {number} [local_connectivity=1] - number of nearest neighbors connected in the local neighborhood.
 * @property {number} [min_dist=1] - controls how tightly points get packed together.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric | "precomputed"} [metric=euclidean] - the metric which defines the distance between two points in the high-dimensional space.
 * @property {number} [_spread=1] - The effective scale of embedded points.
 * @property {number} [_set_op_mix_ratio=1] - Interpolate between union and intersection.
 * @property {number} [_repulsion_strength=1] - Weighting applied to negative samples.
 * @property {number} [_negative_sample_rate=5] - The number of negative samples per positive sample.
 * @property {number} [_n_epochs=350] - The number of training epochs.
 * @property {number} [_initial_alpha=1] - The initial learning rate for the optimization.
 * @property {number} [seed=1212] - the seed for the random number generator.
 */
