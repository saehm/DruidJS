/**
 * @module dimred
 */
export { FASTMAP } from "./FASTMAP.js";
export { ISOMAP } from "./ISOMAP.js";
export { LDA } from "./LDA.js";
export { LLE } from "./LLE.js";
export { LSP } from "./LSP.js";
export { LTSA } from "./LTSA.js";
export { MDS } from "./MDS.js";
export { PCA } from "./PCA.js";
export { SAMMON } from "./SAMMON.js";
export { SMACOF } from "./SMACOF.js";
export { SQDMDS } from "./SQDMDS.js";
export { TopoMap } from "./TopoMap.js";
export { TriMap } from "./TriMap.js";
export { TSNE } from "./TSNE.js";
export { UMAP } from "./UMAP.js";

/** @import { Metric } from "../metrics/index.js" */
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
 * @property {Partial<EigenArgs>} [eig_args={}] - Parameters for the eigendecomposition algorithm.
 */

/**
 * @typedef {Object} ParametersLTSA
 * @property {number} [neighbors] - The number of neighbors for LTSA.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
 * @property {number} [seed=1212] - the seed for the random number generator.
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
 * @property {number} [weight_adj=500] - scaling factor.
 * @property {number} [n_inliers=5] - number of inliers.
 * @property {number} [n_outliers=5] - number of outliers.
 * @property {number} [n_random=5] - number of random points.
 * @property {number} [d=2] - the dimensionality of the projection.
 * @property {number} [tol=1e-8]
 * @property {Metric} [metric=euclidean] - the metric which defines the distance between two points.
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
