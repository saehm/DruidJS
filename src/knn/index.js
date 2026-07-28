/** @module knn */

/**
 * @module knn
 */
export { Annoy } from "./Annoy.js";
export { BallTree } from "./BallTree.js";
export { HNSW } from "./HNSW.js";
export { KDTree } from "./KDTree.js";
export { LSH } from "./LSH.js";
export { NaiveKNN } from "./NaiveKNN.js";
export { NNDescent } from "./NNDescent.js";
export { spatial_tree } from "./spatial_tree.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * @typedef {Object} ParametersAnnoy
 * @property {Metric} metric - Metric to use: (a, b) => distance. Default is `euclidean`
 * @property {number} numTrees - Number of random projection trees to build. Default is `10`
 * @property {number} maxPointsPerLeaf - Maximum points per leaf node. Default is `10`
 * @property {number} seed - Seed for random number generator. Default is `1212`
 */

/**
 * @typedef ParametersBallTree
 * @property {Metric} metric - Metric to use: (a, b) => distance. Default is `euclidean`
 * @property {number} leaf_size - Members below which a ball stops splitting and is scanned
 *   linearly. Default is `32`
 * @property {number} seed - Seed for random number generator. Default is `1212`
 */

/** @typedef {Object} ParametersHNSW
 * @property {Metric} metric - Metric to use: (a, b) => distance. Default is `euclidean`
 * @property {boolean} heuristic - Use heuristics or naive selection. Default is `true`
 * @property {number} m - Number of connections a newly inserted element creates on each layer it
 *   joins, layer 0 included. Default is `16`
 * @property {number} ef_construction - Size of candidate list during construction. Default is `200`
 * @property {number | null} m0 - Upper bound on connections at the ground layer (layer 0). Unlike
 *   `m` this is not a budget for the inserted element, it is the cap enforced once reverse
 *   connections from later insertions have accumulated. Default is `2 * m`
 * @property {number | null} mL - Normalization factor for level generation. Default is `1 / Math.log(m)`
 * @property {number} seed - Seed for random number generator. Default is `1212`
 * @property {number} ef - Size of candidate list during search. Default is `50`
 */

/**
 * @typedef {Object} ParametersKDTree
 * @property {Metric} metric - Metric to use: (a, b) => distance. Default is `euclidean`
 * @property {number} seed
 */

/**
 * @typedef {Object} ParametersLSH
 * @property {Metric} metric - Metric to use: (a, b) => distance. Default is `euclidean`
 * @property {number} numHashTables - Number of hash tables. Default is `10`
 * @property {number} numHashFunctions - Number of hash functions per table. Default is `10`
 * @property {number | null} bucketWidth - Quantization width w. Must match the scale of the data;
 *   estimated from the sampled distance distribution when `null`. Default is `null`
 * @property {number} seed - Seed for random number generator. Default is `1212`
 */

/**
 * @typedef {Object} ParametersNaiveKNN
 * @property {Metric | "precomputed"} [metric] Is either precomputed or a function to use: (a, b) => distance
 * @property {number} [seed]
 */

/** @typedef ParametersNNDescent
 * @property {Metric} metric - Called sigma in paper. Default is `euclidean`
 * @property {number} samples =10 - Number of samples. Default is `10`
 * @property {number} rho = 1 - Sample rate. Default is `1`
 * @property {number} delta = 0.001 - Precision parameter. Default is `0.001`
 * @property {number} seed = 1212 - Seed for the random number generator. Default is `1212`
 */
