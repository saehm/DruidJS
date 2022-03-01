/**
 * @class
 * @alias KMedoids
 */
export class KMedoids {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias KMedoids
     * @todo needs restructuring.
     * @param {Matrix} matrix - data matrix
     * @param {Number} K - number of clusters
     * @param {Number} [max_iter=null] - maximum number of iterations. Default is 10 * Math.log10(N)
     * @param {Function} [metric = euclidean] - metric defining the dissimilarity
     * @param {Number} [seed = 1212] - seed value for random number generator
     * @returns {KMedoids}
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(matrix: Matrix, K: number, max_iter?: number, metric?: Function, seed?: number);
    _metric: Function;
    _matrix: Matrix;
    _A: Float64Array[];
    _K: number;
    _N: any;
    _D: any;
    _max_iter: number;
    _distance_matrix: Matrix;
    _randomizer: Randomizer;
    _clusters: any[];
    _cluster_medoids: number[];
    _is_initialized: boolean;
    /**
     * @returns {Number[][]} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(): number[][];
    generator(): AsyncGenerator<number[][], void, unknown>;
    /**
     * Algorithm 1. FastPAM1: Improved SWAP algorithm
     */
    /** Algorithm 2. FastPAM2: SWAP with multiple candidates
     *
     */
    _iteration(): boolean;
    _get_distance(i: any, j: any, x_i?: any, x_j?: any): any;
    _nearest_medoid(x_j: any, j: any): {
        distance_nearest: any;
        index_nearest: any;
        distance_second: any;
        index_second: any;
    };
    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K: number, cluster_medoids: any): KMedoids;
    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @param {number} K - number of clusters
     *
     */
    _get_random_medoids(K: number): number[];
}
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
//# sourceMappingURL=KMedoids.d.ts.map