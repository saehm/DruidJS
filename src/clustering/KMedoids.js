import { linspace, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { Clustering } from "./Clustering.js";

/** @import {InputType} from "../index.js" */
/** @import { ParametersKMedoids } from "./index.js" */

/**
 * K-Medoids (PAM - Partitioning Around Medoids)
 *
 * A robust clustering algorithm similar to K-Means, but uses actual data points (medoids)
 * as cluster centers and can work with any distance metric.
 *
 * @class
 * @extends Clustering<ParametersKMedoids>
 * @category Clustering
 * @see {@link KMeans} for a faster but less robust alternative
 */
export class KMedoids extends Clustering {
    /**
     * @param {InputType} points - Data matrix
     * @param {Partial<ParametersKMedoids>} parameters
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(points, parameters = {}) {
        super(points, Object.assign({ K: 4, max_iter: null, metric: euclidean, seed: 1212 }, parameters));
        this._A = this._matrix.to2dArray();
        let K = this._parameters.K;
        const N = this._N;
        this._max_iter = this._parameters.max_iter ?? 10 * Math.log10(N);
        this._distance_matrix = new Matrix(N, N, "zeros");

        if (K > N) {
            this._parameters.K = K = N;
        }
        this._randomizer = new Randomizer(this._parameters.seed);
        this._clusters = new Array(N).fill(-1);
        this._cluster_medoids = this._get_random_medoids(K);
        this._is_initialized = false;
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        if (!this._is_initialized) {
            this.get_clusters();
        }
        return this._clusters;
    }

    /** @returns {number[][]} - Array of clusters with the indices of the rows in given points. */
    get_clusters() {
        const K = this._parameters.K;
        const A = this._A;
        const N = this._N;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        /** @type {number[][]} */
        const result = new Array(K).fill(0).map(() => []);
        for (let j = 0; j < N; j++) {
            const nearest = this._nearest_medoid(A[j], j);
            const cluster_idx = nearest.index_nearest;
            result[cluster_idx].push(j);
            this._clusters[j] = cluster_idx;
        }
        return result;
    }

    /** @returns {number} */
    get k() {
        return this._parameters.K;
    }

    /** @returns {number[]} */
    get medoids() {
        return this.get_medoids();
    }

    /** @returns {number[]} */
    get_medoids() {
        const K = this._parameters.K;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        return this._cluster_medoids;
    }

    async *generator() {
        const max_iter = this._max_iter;
        if (!this._is_initialized) {
            this.get_clusters();
        }
        yield this.get_clusters();
        let i = 0;
        while (i < max_iter) {
            const finish = this._iteration();
            this._update_clusters();
            yield this.get_clusters();
            if (finish) break;
            i++;
        }
    }

    /** Algorithm 1. FastPAM1: Improved SWAP algorithm */
    /* _iteration_1() {
        const A = this._A;
        const N = this._N;
        const K = this._K;
        const medoids = this._cluster_medoids;
        let DeltaTD = 0;
        let m0 = null;
        let x0 = null;
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const nearest_medoid = this._nearest_medoid(x_j, j);
                const d_j = nearest_medoid.distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    // disance to new medoid
                    const d_oj = this._get_distance(o, j, x_o, x_j);
                    const {
                        "index_nearest": n,
                        "distance_nearest": d_n,
                        "distance_second": d_s,
                    } = this._nearest_medoid(x_o, o);
                    this._clusters[o] = n; // cached values
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change
                    if (d_oj < d_n) { // reassignment check
                        deltaTD.forEach((d_i, i) => {
                            if (n !== i) {
                                deltaTD[i] = d_i + d_oj - d_n; // update loss change
                            }
                        });
                    }
                });
                // choose best medoid i;
                const i = deltaTD
                    .map((d, i) => [d, i])
                    .sort((d1, d2) => d1[0] - d2[0])[0][1];
                const deltaTD_i = deltaTD[i];
                // store
                if (deltaTD_i < DeltaTD) {
                    DeltaTD = deltaTD_i;
                    m0 = i;
                    x0 = j;
                }
            }
        });

        if (DeltaTD >= 0) {
            return true // break loop if DeltaTD >= 0
        }
        // swap roles of medoid m and non-medoid x;
        medoids[m0] = x0;
        this._cluster_medoids = medoids;
        return false
    } */

    /**
     * FastPAM1: One best swap per iteration
     * @private
     * @returns {boolean}
     */
    _iteration() {
        const A = this._A;
        const K = this._parameters.K;
        const medoids = this._cluster_medoids;
        const N = this._N;

        // Precompute nearest and second nearest medoid for all points
        const cache = new Array(N);
        for (let i = 0; i < N; i++) {
            cache[i] = this._nearest_medoid(A[i], i);
        }

        let best_delta = 0;
        let best_swap = null; // { m_idx: index in medoids, x_idx: index in A }

        // For each non-medoid point j, evaluate swapping it with each medoid i
        const medoid_set = new Set(medoids);
        for (let j = 0; j < N; j++) {
            if (medoid_set.has(j)) continue;

            const x_j = A[j];
            const d_j = cache[j].distance_nearest;

            // deltaTD[i] will store the change in total distance if we swap medoid[i] with j
            const deltaTD = new Array(K).fill(-d_j);

            for (let o = 0; o < N; o++) {
                if (o === j) continue;
                const dist_o_j = this._get_distance(o, j, A[o], x_j);
                const { index_nearest: n, distance_nearest: d_n, distance_second: d_s } = cache[o];

                // If o is assigned to the current medoid being swapped out (n)
                deltaTD[n] += Math.min(dist_o_j, d_s) - d_n;

                // For all other medoids i != n, if j is closer to o than its current medoid
                if (dist_o_j < d_n) {
                    for (let i = 0; i < K; i++) {
                        if (i !== n) {
                            deltaTD[i] += dist_o_j - d_n;
                        }
                    }
                }
            }

            // Find best medoid to swap with j
            for (let i = 0; i < K; i++) {
                if (deltaTD[i] < best_delta) {
                    best_delta = deltaTD[i];
                    best_swap = { m_idx: i, x_idx: j };
                }
            }
        }

        if (best_swap && best_delta < 0) {
            medoids[best_swap.m_idx] = best_swap.x_idx;
            this._cluster_medoids = medoids;
            return false; // not finished
        }

        return true; // finished
    }

    /**
     * @private
     * Get distance between two points
     * @param {number} i
     * @param {number} j
     * @param {Float64Array?} x_i
     * @param {Float64Array?} x_j
     * @returns {number}
     */
    _get_distance(i, j, x_i = null, x_j = null) {
        if (i === j) return 0;
        const D = this._distance_matrix;
        const A = this._A;
        const metric = this._parameters.metric;
        let d_ij = D.entry(i, j);
        if (d_ij === 0) {
            d_ij = metric(x_i || A[i], x_j || A[j]);
            D.set_entry(i, j, d_ij);
            D.set_entry(j, i, d_ij);
        }
        return d_ij;
    }

    /**
     * @private
     * @param {Float64Array} x_j
     * @param {number} j
     * @returns
     */
    _nearest_medoid(x_j, j) {
        const medoids = this._cluster_medoids;
        const A = this._A;
        if (medoids.length === 0) {
            throw new Error("No medoids available. Initialization failed.");
        }

        let d_n = Infinity;
        let n = -1;
        let d_s = Infinity;
        let s = -1;

        for (let i = 0; i < medoids.length; i++) {
            const m = medoids[i];
            const d = this._get_distance(j, m, x_j, A[m]);
            if (d < d_n) {
                d_s = d_n;
                s = n;
                d_n = d;
                n = i;
            } else if (d < d_s) {
                d_s = d;
                s = i;
            }
        }

        if (s === -1) s = n;

        return {
            distance_nearest: d_n,
            index_nearest: n,
            distance_second: d_s,
            index_second: s,
        };
    }

    /**
     * @private
     */
    _update_clusters() {
        const N = this._N;
        const A = this._A;
        for (let j = 0; j < N; j++) {
            const nearest = this._nearest_medoid(A[j], j);
            this._clusters[j] = nearest.index_nearest;
        }
    }

    /**
     * Computes `K` clusters out of the `matrix`.
     * @param {number} K - Number of clusters.
     * @param {number[]} cluster_medoids
     */
    init(K, cluster_medoids) {
        if (!K) K = this._parameters.K;
        if (!cluster_medoids) cluster_medoids = this._get_random_medoids(K);
        this._cluster_medoids = cluster_medoids;
        const max_iter = this._max_iter;
        let finish = false;
        let i = 0;
        do {
            finish = this._iteration();
        } while (!finish && ++i < max_iter);
        this._update_clusters();
        this._is_initialized = true;
        return this;
    }

    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @private
     * @param {number} K - Number of clusters
     * @returns {number[]}
     */
    _get_random_medoids(K) {
        const N = this._N;
        const A = this._A;
        const indices = linspace(0, N - 1);
        const randomizer = this._randomizer;
        const n = Math.min(N, 10 + Math.ceil(Math.sqrt(N)));

        // Handle case where K >= N
        if (K >= N) {
            return indices.slice(0, N);
        }

        /** @type {number[]} */
        const medoids = [];

        // first medoid: select from a random sample of size n
        let best_j = -1;
        let min_td = Infinity;
        let S = randomizer.choice(indices, n);
        for (let j = 0; j < S.length; ++j) {
            let td = 0;
            const S_j = S[j];
            const x_j = A[S_j];
            for (let o = 0; o < S.length; ++o) {
                if (o === j) continue;
                td += this._get_distance(S_j, S[o], x_j, A[S[o]]);
            }
            if (td < min_td) {
                min_td = td;
                best_j = S_j;
            }
        }
        medoids.push(best_j);

        // other medoids: greedy additive selection (Algorithm LAB)
        for (let i = 1; i < K; ++i) {
            let best_idx = -1;
            let best_delta = Infinity;

            const remainingIndices = indices.filter((idx) => !medoids.includes(idx));
            if (remainingIndices.length === 0) break;

            S = randomizer.choice(remainingIndices, Math.min(n, remainingIndices.length));
            for (let j = 0; j < S.length; ++j) {
                let deltaTD = 0;
                const S_j = S[j];
                const x_j = A[S_j];

                // Estimate TD reduction on the sample S
                for (let o = 0; o < S.length; ++o) {
                    if (o === j) continue;
                    const S_o = S[o];
                    const x_o = A[S_o];

                    // Closest distance to current medoids
                    let min_d_existing = Infinity;
                    for (let m = 0; m < medoids.length; m++) {
                        const d = this._get_distance(S_o, medoids[m], x_o, A[medoids[m]]);
                        if (d < min_d_existing) min_d_existing = d;
                    }

                    const delta = this._get_distance(S_j, S_o, x_j, x_o) - min_d_existing;
                    if (delta < 0) {
                        deltaTD += delta;
                    }
                }

                if (deltaTD < best_delta) {
                    best_delta = deltaTD;
                    best_idx = S_j;
                }
            }
            if (best_idx !== -1) {
                medoids.push(best_idx);
            }
        }
        return medoids;
    }
}
