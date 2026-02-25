/**
 * @module clustering
 */
export { CURE } from "./CURE.js";
export { HierarchicalClustering } from "./Hierarchical_Clustering.js";
export { KMeans } from "./KMeans.js";
export { KMedoids } from "./KMedoids.js";
export { MeanShift } from "./MeanShift.js";
export { OPTICS } from "./OPTICS.js";
export { XMeans } from "./XMeans.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * @typedef ParametersHierarchicalClustering
 * @property {"single" | "complete" | "average"} linkage
 * @property {Metric | "precomputed"} metric
 */

/**
 * @typedef ParametersKMeans
 * @property {number} K
 * @property {Metric} metric Default is `euclidean`
 * @property {number} seed Default is `1212`
 * @property {number[][] | Float64Array[]} [initial_centroids] - Initial centroids. Default is `null`
 */

/** @typedef ParametersKMedoids
 * @property {number} K - Number of clusters
 * @property {number?} max_iter - Maximum number of iterations. Default is 10 * Math.log10(N). Default is `null`
 * @property {Metric} metric - Metric defining the dissimilarity. Default is `euclidean`
 * @property {number} seed - Seed value for random number generator. Default is `1212`
 */

/** @typedef ParametersOptics
 * @property {number} epsilon - The minimum distance which defines whether a point is a neighbor or not.
 * @property {number} min_points - The minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
 * @property {Metric} metric - The distance metric which defines the distance between two points of the points. Default is `euclidean`
 */

/**
 * @typedef ParametersXMeans
 * @property {number} K_min - Minimum number of clusters. Default is `2`
 * @property {number} K_max - Maximum number of clusters. Default is `10`
 * @property {Metric} metric - Distance metric function. Default is `euclidean`
 * @property {number} seed - Random seed. Default is `1212`
 * @property {number} min_cluster_size - Minimum points required to consider splitting a cluster. Default is `25`
 * @property {number} tolerance - Convergence tolerance for KMeans. Default is `0.001`
 */

/**
 * @typedef ParametersMeanShift
 * @property {number} bandwidth - bandwidth
 * @property {Metric} metric - Metric defining the dissimilarity. Default is `euclidean`
 * @property {number} seed - Seed value for random number generator. Default is `1212`
 * @property {"flat" | "gaussian" | ((dist: number) => number)} kernel - Kernel function. Default is `gaussian`
 * @property {number} [max_iter] - Maximum number of iterations. Default is `Math.max(10, Math.floor(10 * Math.log10(N)))`
 * @property {number} [tolerance] - Convergence tolerance. Default is `1e-3`
 */

/**
 * @typedef ParametersCURE
 * @property {number} K - Target number of clusters. Default is `2`
 * @property {number} num_representatives - Number of representative points per cluster. Default is `5`
 * @property {number} shrink_factor - Factor to shrink representatives toward centroid (0-1). Default is `0.5`
 * @property {Metric} metric - Distance metric function. Default is `euclidean`
 * @property {number} seed - Random seed. Default is `1212`
 */
