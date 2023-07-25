import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { Heap } from "../datastructure/index.js";
import { linspace } from "../matrix/index.js";

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
     * @param {Numbers} K 
     * @param {Function} [metric = euclidean] 
     * @param {Number} [seed = 1987]
     * @param {Boolean} [init = true]
     * @returns {KMeans}
     */
    constructor(matrix, K, metric = euclidean, seed=1987, init = true) {
        this._metric = metric;
        this._matrix = matrix;
        this._K = K;
        const [N, D] = matrix.shape;
        this._N = N;
        this._D = D;
        if (K > N) K = N;
        this._randomizer = new Randomizer(seed);
        this._clusters = new Array(N).fill(undefined);
        this._cluster_centroids = this._get_random_centroids(K);
        if (init) this.init(K, this._cluster_centroids);
        return this;
    }

    /**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */
    get_clusters() {
        const K = this._K;
        const clusters = this._clusters;
        const result = new Array(K).fill().map(() => new Array());
        clusters.forEach((c, i) => result[c].push(i));
        return result;
    }

    /**
     * @private
     * @param {Array} points 
     * @param {Array} candidates 
     */
    _furthest_point(points, candidates) {
        const A = this._matrix;
        const metric = this._metric;
        let i = points.length;
        let H = Heap.heapify(
            candidates, 
            (d) => {
                const Ad = A.row(d)
                let sum = 0;
                for (let j = 0; j < i; ++j) {
                    sum += metric(Ad, points[j])
                }
                return sum;
            }, 
            "max"
        )
        return H.pop().element;
    }

    _get_random_centroids(K) {
        const N = this._N;
        const randomizer = this._randomizer;
        const A = this._matrix;
        const cluster_centroids = new Array(K).fill()
        const indices = linspace(0, N - 1);
        const random_point = randomizer.random_int % (N - 1);
        cluster_centroids[0] = A.row(random_point);
        const init_points = [random_point];
        const sample_size = Math.floor((N - K) / K);// / K
        for (let i = 1; i < K; ++i) {
            // sampling + kmeans++ improvement?
            const sample = randomizer.choice(indices.filter(d => init_points.indexOf(d) == -1), sample_size);
            const furthest_point = this._furthest_point(cluster_centroids.slice(0, i), sample);
            init_points.push(furthest_point);
            cluster_centroids[i] = A.row(furthest_point);
        }
        return cluster_centroids;
    }

    _iteration(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const metric = this._metric;
        const clusters = this._clusters;
        let clusters_changed = false;
        // find nearest cluster centroid.
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i)
            let min_dist = Infinity;
            let min_cluster = null;
            for (let j = 0; j < K; ++j) {
                let d = metric(cluster_centroids[j], Ai);
                if (d < min_dist) {
                    min_dist = d;
                    min_cluster = j; 
                }
            }
            if (clusters[i] !== min_cluster) {
                clusters_changed = true;
            }
            clusters[i] = min_cluster;
        }
        // update cluster centroid
        // reset cluster centroids to 0
        for (let i = 0; i < K; ++i) {
            const centroid = cluster_centroids[i];
            for (let j = 0; j < D; ++j) {
                centroid[j] = 0;
            }
        }
        // compute centroid
        this._compute_centroid(cluster_centroids);

        return {   
            "clusters_changed": clusters_changed,
            "cluster_centroids": cluster_centroids
        };
    }

    _compute_centroid(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const clusters = this._clusters;
        const cluster_counter = new Array(K).fill(0);

        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            const ci = clusters[i];
            cluster_counter[ci]++;
            const centroid = cluster_centroids[ci];
            for (let j = 0; j < D; ++j) {
                centroid[j] += Ai[j];
            }
        }
        for (let i = 0; i < K; ++i) {
            const n = cluster_counter[i];
            cluster_centroids[i] = cluster_centroids[i].map(c => c / n);
        }
        
    }

    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K, cluster_centroids) {
        if (!K) K = this._K;
        if (!cluster_centroids) cluster_centroids = this._get_random_centroids(K);
        let clusters_changed = false;
        do {
            const iteration_result = this._iteration(cluster_centroids)
            cluster_centroids = iteration_result.cluster_centroids;
            clusters_changed = iteration_result.clusters_changed;
        } while (clusters_changed)
    }
    
}
