import { Heap } from "../datastructure/index.js";
import { linspace, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { Clustering } from "./Clustering.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersKMeans } from "./index.js" */
/**
 * K-Means Clustering
 *
 * A popular clustering algorithm that partitions data into K clusters where each point
 * belongs to the cluster with the nearest mean (centroid).
 *
 * @class
 * @extends Clustering<ParametersKMeans>
 * @category Clustering
 * @see {@link KMedoids} for a more robust alternative
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const points = [[1, 1], [1.5, 1.5], [5, 5], [5.5, 5.5]];
 * const kmeans = new druid.KMeans(points, { K: 2 });
 *
 * const clusters = kmeans.get_cluster_list(); // [0, 0, 1, 1]
 * const centroids = kmeans.centroids; // center points
 */
export class KMeans extends Clustering {
    /**
     * @param {InputType} points
     * @param {Partial<ParametersKMeans>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersKMeans} */ (Object.assign({ K: 4, metric: euclidean, seed: 1212 }, parameters)),
        );

        const K = this._parameters.K;
        const seed = parameters.seed;

        // Convert points to Matrix if needed
        if (points instanceof Matrix) {
            this._matrix = points;
        } else {
            this._matrix = Matrix.from(points);
        }

        const [N, D] = this._matrix.shape;
        this._N = N;
        this._D = D;

        this._K = K > N ? N : K;
        this._randomizer = new Randomizer(seed);

        /** @type {number[]} */
        this._clusters = new Array(N).fill(0);

        this._cluster_centroids = parameters.initial_centroids
            ? parameters.initial_centroids.map((c) => new Float64Array(c))
            : this._get_random_centroids(this._K);
        let cluster_centroids = this._cluster_centroids;
        let iterations = 0;
        const max_iterations = 300;
        let clusters_changed = true;

        while (clusters_changed && iterations < max_iterations) {
            const iteration_result = this._iteration(cluster_centroids);
            cluster_centroids = iteration_result.cluster_centroids;
            clusters_changed = iteration_result.clusters_changed;
            iterations++;
        }

        this._cluster_centroids = cluster_centroids;
    }

    /** @returns {number} The number of clusters */
    get k() {
        return this._K;
    }

    /** @returns {Float64Array[]} The cluster centroids */
    get centroids() {
        return this._cluster_centroids;
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        return this._clusters;
    }

    /** @returns {number[][]} An Array of clusters with the indices of the points. */
    get_clusters() {
        const K = this._K;
        const clusters = this._clusters;
        /** @type {number[][]} */
        const result = new Array(K).fill(0).map(() => []);
        clusters.forEach((c, i) => {
            if (c >= 0 && c < K) {
                result[c].push(i);
            }
        });
        return result;
    }

    /**
     * @private
     * @param {number[]} point_indices
     * @param {number[]} candidates
     * @returns {number}
     */
    _furthest_point(point_indices, candidates) {
        const A = this._matrix;
        const metric = this._parameters.metric;

        if (point_indices.length === 0 || candidates.length === 0) {
            return candidates[0] ?? 0;
        }

        const H = Heap.heapify(
            candidates,
            (d) => {
                const Ad = A.row(d);
                let sum = 0;
                for (let j = 0; j < point_indices.length; ++j) {
                    sum += metric(Ad, A.row(point_indices[j]));
                }
                return sum;
            },
            "max",
        );

        const furthest = H.pop();
        if (!furthest) throw new Error("Should not happen!");

        return furthest.element;
    }

    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    _get_random_centroids(K) {
        const N = this._N;
        const randomizer = this._randomizer;
        const A = this._matrix;
        /** @type {Float64Array[]} */
        const cluster_centroids = new Array(K);
        const indices = linspace(0, N - 1);

        // First centroid: random selection
        const random_point = randomizer.random_int % N;
        cluster_centroids[0] = A.row(random_point);
        const init_points = [random_point];

        const sample_size = Math.max(1, Math.floor((N - K) / K));

        for (let i = 1; i < K; ++i) {
            const remaining = indices.filter((d) => !init_points.includes(d));
            if (remaining.length === 0) break;

            const sample = randomizer.choice(remaining, Math.min(sample_size, remaining.length));
            const furthest_point = this._furthest_point(init_points, sample);

            init_points.push(furthest_point);
            cluster_centroids[i] = A.row(furthest_point);
        }

        return cluster_centroids;
    }

    /**
     * @private
     * @param {Float64Array[]} cluster_centroids
     * @returns {{ clusters_changed: boolean; cluster_centroids: Float64Array[] }}
     */
    _iteration(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const metric = this._parameters.metric;
        const A = this._matrix;
        const clusters = this._clusters;
        let clusters_changed = false;

        // Find nearest cluster centroid for each point
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            let min_dist = Infinity;
            let min_cluster = 0;

            for (let j = 0; j < K; ++j) {
                const d = metric(cluster_centroids[j], Ai);
                if (d < min_dist) {
                    min_dist = d;
                    min_cluster = j;
                }
            }

            if (clusters[i] !== min_cluster) {
                clusters_changed = true;
                clusters[i] = min_cluster;
            }
        }

        // Update cluster centroids
        const new_centroids = this._compute_centroid(K);

        return {
            clusters_changed: clusters_changed,
            cluster_centroids: new_centroids,
        };
    }

    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    _compute_centroid(K) {
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const clusters = this._clusters;

        // Initialize new centroids and counters
        /** @type {Float64Array[]} */
        const new_centroids = new Array(K);
        const cluster_counter = new Array(K).fill(0);

        for (let i = 0; i < K; ++i) {
            new_centroids[i] = new Float64Array(D);
        }

        // Sum up all points in each cluster
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            const ci = clusters[i];
            if (ci >= 0 && ci < K) {
                cluster_counter[ci]++;
                const centroid = new_centroids[ci];
                for (let j = 0; j < D; ++j) {
                    centroid[j] += Ai[j];
                }
            }
        }

        // Divide by count to get mean
        for (let i = 0; i < K; ++i) {
            const n = cluster_counter[i];
            if (n > 0) {
                const centroid = new_centroids[i];
                for (let j = 0; j < D; ++j) {
                    centroid[j] /= n;
                }
            }
        }

        return new_centroids;
    }
}
