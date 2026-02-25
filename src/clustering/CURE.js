import { euclidean } from "../metrics/index.js";
import { Clustering } from "./Clustering.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersCURE } from "./index.js" */

/**
 * CURE (Clustering Using REpresentatives)
 *
 * An efficient clustering algorithm for large databases that is robust to outliers
 * and identifies clusters with non-spherical shapes and wide variances in size.
 *
 * @class
 * @extends Clustering<ParametersCURE>
 * @category Clustering
 */
export class CURE extends Clustering {
    /** @type {number} */
    _K;
    /** @type {number} */
    _num_representatives;
    /** @type {number} */
    _shrink_factor;
    /**
     * @private
     * @type {CURECluster[]}
     */
    _clusters = [];
    /** @type {number[]} */
    _cluster_ids = [];

    /**
     * @param {InputType} points
     * @param {Partial<ParametersCURE>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersCURE} */ (
                Object.assign(
                    { K: 2, num_representatives: 5, shrink_factor: 0.5, metric: euclidean, seed: 1212 },
                    parameters,
                )
            ),
        );

        this._K = this._parameters.K ?? 2;
        this._num_representatives = this._parameters.num_representatives ?? 5;
        this._shrink_factor = this._parameters.shrink_factor ?? 0.5;

        // Initialize clusters
        this._initialize_clusters();
        // Run CURE algorithm
        this._cure();
    }

    /**
     * Initialize each point as its own cluster
     * @private
     */
    _initialize_clusters() {
        const N = this._N;
        //const D = this._D;
        this._clusters = [];

        for (let i = 0; i < N; ++i) {
            const point = this._matrix.row(i);
            const centroid = new Float64Array(point);
            // For single point, representative is the point itself
            const representatives = [new Float64Array(point)];

            this._clusters.push(new CURECluster([i], centroid, representatives));
        }
    }

    /**
     * Compute distance between two clusters using representative points
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {number}
     */
    _cluster_distance(cluster1, cluster2) {
        const reps1 = cluster1.representatives;
        const reps2 = cluster2.representatives;
        const metric = this._parameters.metric;

        let min_dist = Infinity;
        for (const r1 of reps1) {
            for (const r2 of reps2) {
                const dist = metric(r1, r2);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
        }
        return min_dist;
    }

    /**
     * Find the closest pair of clusters
     * @private
     * @returns {[number, number, number]} [index1, index2, distance]
     */
    _find_closest_clusters() {
        let min_dist = Infinity;
        let min_i = 0;
        let min_j = 1;

        for (let i = 0; i < this._clusters.length; ++i) {
            for (let j = i + 1; j < this._clusters.length; ++j) {
                const dist = this._cluster_distance(this._clusters[i], this._clusters[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        return [min_i, min_j, min_dist];
    }

    /**
     * Merge two clusters
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {CURECluster}
     */
    _merge_clusters(cluster1, cluster2) {
        // Merge indices
        const merged_indices = [...cluster1.indices, ...cluster2.indices];

        // Calculate new centroid
        const size1 = cluster1.indices.length;
        const size2 = cluster2.indices.length;
        const total_size = size1 + size2;
        const D = this._D;
        const new_centroid = new Float64Array(D);

        for (let d = 0; d < D; ++d) {
            new_centroid[d] = (size1 * cluster1.centroid[d] + size2 * cluster2.centroid[d]) / total_size;
        }

        // Collect all points from both clusters
        /** @type {{index: number, point: Float64Array}[]} */
        const all_points = [];
        for (const idx of cluster1.indices) {
            all_points.push({ index: idx, point: this._matrix.row(idx) });
        }
        for (const idx of cluster2.indices) {
            all_points.push({ index: idx, point: this._matrix.row(idx) });
        }

        // Select representative points - pick points farthest from centroid
        const num_reps = Math.min(this._num_representatives, all_points.length);
        const metric = this._parameters.metric;

        // Calculate distances from centroid for all points
        const distances = all_points.map(({ point }) => metric(point, new_centroid));

        // Select num_reps points with maximum distance (farthest from centroid)
        const selected_indices = [];
        const used = new Set();

        for (let r = 0; r < num_reps; ++r) {
            let max_dist = -1;
            let max_idx = -1;

            for (let i = 0; i < distances.length; ++i) {
                if (!used.has(i) && distances[i] > max_dist) {
                    max_dist = distances[i];
                    max_idx = i;
                }
            }

            if (max_idx >= 0) {
                used.add(max_idx);
                selected_indices.push(max_idx);
            }
        }

        // Shrink representative points toward centroid
        const new_representatives = selected_indices.map((idx) => {
            const point = all_points[idx].point;
            const shrunk = new Float64Array(D);
            const alpha = this._shrink_factor;

            for (let d = 0; d < D; ++d) {
                shrunk[d] = point[d] + alpha * (new_centroid[d] - point[d]);
            }

            return shrunk;
        });

        return new CURECluster(merged_indices, new_centroid, new_representatives);
    }

    /**
     * Run CURE clustering algorithm
     * @private
     */
    _cure() {
        // Merge clusters until we have K clusters
        while (this._clusters.length > this._K) {
            const [i, j] = this._find_closest_clusters();

            // Merge clusters i and j
            const merged = this._merge_clusters(this._clusters[i], this._clusters[j]);

            // Remove the old clusters and add the merged one
            // Remove larger index first to maintain correct indices
            // min_i < min_j is always true from _find_closest_clusters
            this._clusters.splice(j, 1);
            this._clusters.splice(i, 1);

            this._clusters.push(merged);
        }

        // Build cluster list for get_cluster_list
        this._build_cluster_ids();
    }

    /**
     * Build the cluster list (point -> cluster assignment)
     * @private
     */
    _build_cluster_ids() {
        const N = this._N;
        this._cluster_ids = new Array(N).fill(-1);

        for (let c = 0; c < this._clusters.length; ++c) {
            for (const idx of this._clusters[c].indices) {
                this._cluster_ids[idx] = c;
            }
        }
    }

    /**
     * @returns {number[][]}
     */
    get_clusters() {
        return this._clusters.map((cluster) => cluster.indices);
    }

    /**
     * @returns {number[]}
     */
    get_cluster_list() {
        return this._cluster_ids;
    }
}

/**
 * @private
 * Represents a cluster in CURE algorithm
 */
class CURECluster {
    /**
     * @param {number[]} indices - Indices of points in the cluster
     * @param {Float64Array} centroid - Centroid of the cluster
     * @param {Float64Array[]} representatives - Representative points (shrunk toward centroid)
     */
    constructor(indices, centroid, representatives) {
        /** @type {number[]} */
        this.indices = indices;
        /** @type {Float64Array} */
        this.centroid = centroid;
        /** @type {Float64Array[]} */
        this.representatives = representatives;
    }
}
