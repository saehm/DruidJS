/**
 * @class
 * @alias KMeans
 */
export class KMeans {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias KMeans
     * @todo needs restructuring.
     * @param {Matrix} matrix
     * @param {Number} K
     * @param {Function} [metric = euclidean]
     * @param {Number} [seed = 1212]
     * @param {Boolean} [init = true]
     * @returns {KMeans}
     */
    constructor(matrix: Matrix, K: number, metric?: Function, seed?: number, init?: boolean);
    _metric: Function;
    _matrix: Matrix;
    _K: number;
    _N: any;
    _D: any;
    _randomizer: Randomizer;
    _clusters: any[];
    _cluster_centroids: any[];
    /**
     * @returns {Number[][]} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(): number[][];
    /**
     * @private
     * @param {Array} points
     * @param {Array} candidates
     */
    private _furthest_point;
    _get_random_centroids(K: any): any[];
    _iteration(cluster_centroids: any): {
        clusters_changed: boolean;
        cluster_centroids: any;
    };
    _compute_centroid(cluster_centroids: any): void;
    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K: number, cluster_centroids: any): void;
}
import { Randomizer } from "../util/index.js";
//# sourceMappingURL=KMeans.d.ts.map