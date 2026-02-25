import { euclidean, euclidean_squared } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { Clustering } from "./Clustering.js";
import { KMeans } from "./KMeans.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersXMeans } from "./index.js" */

/**
 * @typedef SplitResult
 * @property {number} index - Index of the cluster being split
 * @property {number} bic_parent - BIC score of the parent cluster
 * @property {number} bic_children - BIC score of the split children
 * @property {number[][]} child_clusters - Clusters after splitting
 * @property {Float64Array[]} child_centroids - Centroids of child clusters
 */

/**
 * @typedef CandidateResult
 * @property {KMeans} kmeans - The KMeans instance for this K
 * @property {number} score - BIC score
 */

/**
 * X-Means Clustering
 *
 * An extension of K-Means that automatically determines the number of clusters (K)
 * using the Bayesian Information Criterion (BIC).
 *
 * @class
 * @extends Clustering<ParametersXMeans>
 * @category Clustering
 */
export class XMeans extends Clustering {
    /**
     * XMeans clustering algorithm that automatically determines the optimal number of clusters.
     *
     * X-Means extends K-Means by starting with a minimum number of clusters and iteratively
     * splitting clusters to improve the Bayesian Information Criterion (BIC).
     *
     * Algorithm:
     * 1. Start with K_min clusters using KMeans
     * 2. For each cluster, try splitting it into 2 sub-clusters
     * 3. If BIC improves after splitting, keep the split
     * 4. Run KMeans again with all (old + new) centroids
     * 5. Repeat until K_max is reached or no more improvements
     *
     * @param {InputType} points - The data points to cluster
     * @param {Partial<ParametersXMeans>} [parameters={}] - Configuration parameters
     * @see {@link https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf}
     * @see {@link https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py}
     * @see {@link https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java}
     */
    constructor(points, parameters = {}) {
        const defaults = {
            K_max: 10,
            K_min: 2,
            metric: euclidean,
            seed: 1212,
            min_cluster_size: 35,
            tolerance: 0.001,
        };
        super(points, /** @type {ParametersXMeans} */ (Object.assign(defaults, parameters)));
        this._randomizer = new Randomizer(this._parameters.seed);

        /** @type {KMeans | null} */
        this._best_kmeans = null;

        // Run XMeans algorithm
        this._run();
    }

    /**
     * Run the XMeans algorithm
     *
     * @private
     */
    _run() {
        /** @type {Map<number, CandidateResult>} */
        const candidates = new Map();
        const A = this._matrix;

        // Initialize with K_min clusters
        let current_kmeans = new KMeans(this._points, {
            K: this._parameters.K_min,
            metric: this._parameters.metric,
            seed: this._parameters.seed,
        });

        let K = this._parameters.K_min;

        candidates.set(K, {
            kmeans: current_kmeans,
            score: -Infinity,
        });

        // Iteratively improve clustering
        while (K < this._parameters.K_max) {
            const clusters = current_kmeans.get_clusters();
            const centroids = current_kmeans.centroids;

            // Try splitting each cluster
            /** @type {SplitResult[]} */
            const split_results = [];

            for (let j = 0; j < clusters.length; ++j) {
                const cluster = clusters[j];

                // Skip small clusters - need enough points for reliable BIC
                if (cluster.length < this._parameters.min_cluster_size) {
                    continue;
                }

                // Get subset data for this cluster
                /** @type {number[][]} */
                const subset_points = cluster.map((idx) => {
                    const row = A.row(idx);
                    return Array.from(row);
                });

                // Calculate BIC for parent (single cluster)
                const parent_bic = this._bic([cluster], [centroids[j]]);

                // Run KMeans with K=2 on subset
                const subset_kmeans = new KMeans(subset_points, {
                    K: 2,
                    metric: this._parameters.metric,
                    seed: this._randomizer.seed,
                });

                const child_clusters_local = subset_kmeans.get_clusters();
                const child_centroids = subset_kmeans.centroids;

                // Map local indices back to global indices
                /** @type {number[][]} */
                const child_clusters_global = child_clusters_local.map((local_cluster) =>
                    local_cluster.map((local_idx) => cluster[local_idx]),
                );

                // Calculate BIC for children (split into 2 clusters)
                const children_bic = this._bic(child_clusters_global, child_centroids);

                split_results.push({
                    index: j,
                    bic_parent: parent_bic,
                    bic_children: children_bic,
                    child_clusters: child_clusters_global,
                    child_centroids: child_centroids,
                });
            }

            // Keep all splits that improve BIC (BIC_children > BIC_parent)
            /** @type {SplitResult[]} */
            const accepted_splits = split_results.filter((result) => result.bic_children > result.bic_parent);

            // If no splits improve BIC, we're done
            if (accepted_splits.length === 0) {
                break;
            }

            // Build new centroids array: keep non-split centroids + add split centroids
            /** @type {Float64Array[]} */
            const new_centroids = [];
            const split_indices = new Set();

            // Sort accepted splits by improvement (descending)
            accepted_splits.sort((a, b) => b.bic_children - b.bic_parent - (a.bic_children - a.bic_parent));

            for (const split of accepted_splits) {
                if (centroids.length + split_indices.size + 1 <= this._parameters.K_max) {
                    split_indices.add(split.index);
                } else {
                    break;
                }
            }

            for (let i = 0; i < centroids.length; ++i) {
                if (split_indices.has(i)) {
                    // This cluster was split - add both child centroids
                    const split_result = accepted_splits.find((s) => s.index === i);
                    if (split_result) {
                        new_centroids.push(...split_result.child_centroids);
                    }
                } else {
                    // This cluster wasn't split - keep its centroid
                    new_centroids.push(centroids[i]);
                }
            }

            // Run KMeans on full dataset with new centroids as initialization
            // This is crucial - we need to reassign all points to all clusters
            const newK = new_centroids.length;

            // Create a new KMeans instance with K set to new number of clusters
            current_kmeans = new KMeans(this._matrix, {
                K: newK,
                metric: this._parameters.metric,
                seed: this._randomizer.seed,
                initial_centroids: new_centroids,
            });

            // Store the candidate with the BIC of the FULL dataset
            candidates.set(newK, {
                kmeans: current_kmeans,
                score: this._bic(current_kmeans.get_clusters(), current_kmeans.centroids),
            });

            K = newK;
        }

        // Select best candidate based on BIC score
        this._best_kmeans = this._select_best_candidate(candidates);
    }

    /**
     * Select the best candidate based on BIC score
     *
     * @private
     * @param {Map<number, CandidateResult>} candidates
     * @returns {KMeans}
     */
    _select_best_candidate(candidates) {
        if (candidates.size === 0) {
            throw new Error("No candidates found");
        }

        const first_candidate = candidates.get(this._parameters.K_min);
        if (!first_candidate) {
            throw new Error("Missing initial candidate");
        }

        let best_score = first_candidate.score;
        /** @type {KMeans} */
        let best_kmeans = first_candidate.kmeans;

        for (const candidate of candidates.values()) {
            if (candidate.score > best_score) {
                best_score = candidate.score;
                best_kmeans = candidate.kmeans;
            }
        }

        return best_kmeans;
    }

    /**
     * Calculate Bayesian Information Criterion for a set of clusters.
     *
     * Uses Kass's formula for BIC calculation:
     * BIC(θ) = L(D) - 0.5 * p * ln(N)
     *
     * Where:
     * - L(D) is the log-likelihood of the data
     * - p is the number of free parameters: (K-1) + D*K + 1
     * - N is the total number of points
     *
     * @private
     * @param {number[][]} clusters - Array of clusters with point indices
     * @param {Float64Array[]} centroids - Array of centroids
     * @returns {number} BIC score (higher is better)
     */
    _bic(clusters, centroids) {
        const A = this._matrix;
        const D = this._D;
        const K = centroids.length;

        let total_variance = 0;
        let N = 0;

        // Calculate total variance (sum of squared distances)
        for (let i = 0; i < K; ++i) {
            const cluster = clusters[i];
            const centroid = centroids[i];
            N += cluster.length;

            for (let j = 0; j < cluster.length; ++j) {
                const point_idx = cluster[j];
                const point = A.row(point_idx);
                // Sum of squared distances (variance term)
                total_variance += euclidean_squared(centroid, point);
            }
        }

        // Not enough points for meaningful BIC
        if (N <= K) {
            return -Infinity;
        }

        // Estimate variance (ML estimate)
        const variance = total_variance / (N - K);

        // Handle case of zero variance (all points identical)
        if (variance <= 0) {
            return -Infinity;
        }

        // Number of free parameters: (K-1) cluster weights + K*D centroid coordinates + 1 variance
        const p = K - 1 + D * K + 1;

        // Calculate log-likelihood
        let log_likelihood = 0;
        const log_2pi = Math.log(2 * Math.PI);

        for (let i = 0; i < K; ++i) {
            const n = clusters[i].length;
            if (n <= 1) continue;

            // Log-likelihood for cluster i
            const cluster_log_likelihood =
                n * Math.log(n / N) - 0.5 * n * log_2pi - 0.5 * n * D * Math.log(variance) - 0.5 * (n - 1);

            log_likelihood += cluster_log_likelihood;
        }

        // BIC = log_likelihood - 0.5 * p * ln(N)
        return log_likelihood - 0.5 * p * Math.log(N);
    }

    /**
     * Get the computed clusters
     *
     * @returns {number[][]} Array of clusters, each containing indices of points
     */
    get_clusters() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.get_clusters();
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.get_cluster_list();
    }

    /**
     * Get the final centroids
     *
     * @returns {Float64Array[]} Array of centroids
     */
    get centroids() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.centroids;
    }

    /**
     * Get the optimal number of clusters found
     *
     * @returns {number} The number of clusters
     */
    get k() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.k;
    }
}
