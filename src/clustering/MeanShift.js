import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { Clustering } from "./Clustering.js";

/** @import { ParametersMeanShift } from "./index.js" */
/** @import { InputType } from "../index.js" */

/**
 * Mean Shift Clustering
 *
 * A non-parametric clustering technique that does not require prior knowledge of the
 * number of clusters. It identifies centers of density in the data.
 *
 * @class
 * @extends Clustering<ParametersMeanShift>
 * @category Clustering
 */
export class MeanShift extends Clustering {
    /** @type {number} */
    _bandwidth;
    /** @type {number} */
    _max_iter;
    /** @type {number} */
    _tolerance;
    /** @type {(dist: number) => number} */
    _kernel;
    /** @type {Matrix} */
    _points;
    /** @type {number[] | undefined} */
    _clusters;
    /** @type {number[][] | undefined} */
    _cluster_list;
    /**
     *
     * @param {InputType} points
     * @param {Partial<ParametersMeanShift>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersMeanShift} */ (
                Object.assign({ seed: 1212, metric: euclidean, bandwidth: 0, kernel: "gaussian" }, parameters)
            ),
        );

        // Ensure bandwidth is positive
        this._bandwidth = parameters.bandwidth ?? this._compute_bandwidth(this._matrix);
        this._max_iter = parameters.max_iter ?? Math.max(10, Math.floor(10 * Math.log10(this._N)));
        this._tolerance = parameters.tolerance ?? 1e-3;
        const kernel_param = parameters.kernel ?? "gaussian";
        // If kernel is a string, map to function
        if (typeof kernel_param === "string") {
            if (kernel_param === "flat") {
                this._kernel = (dist) => (dist <= this._bandwidth ? 1 : 0);
            } else {
                // gaussian (default)
                this._kernel = (dist) => Math.exp(-(dist * dist) / (2 * this._bandwidth * this._bandwidth));
            }
        } else {
            // custom function
            this._kernel = kernel_param;
        }

        // Copy points to a mutable matrix
        this._points = this._matrix.clone();

        this._mean_shift();
        this._assign_clusters();
    }

    // Helper to compute bandwidth if not provided
    /**
     * @param {Matrix} matrix
     * @returns {number}
     */
    _compute_bandwidth(matrix) {
        const N = matrix.shape[0];
        //const D = matrix.shape[1];
        // Compute average pairwise distance
        let totalDist = 0;
        for (let i = 0; i < N; ++i) {
            const row_i = matrix.row(i);
            for (let j = i + 1; j < N; ++j) {
                const row_j = matrix.row(j);
                const dist = this._parameters.metric(row_i, row_j);
                totalDist += dist;
            }
        }
        const avgDist = totalDist / ((N * (N - 1)) / 2);
        // Use a fraction of avgDist as bandwidth
        return avgDist / 2;
    }

    // Compute kernel weight
    /**
     * @param {number} dist
     * @returns {number}
     */
    _kernel_weight(dist) {
        return this._kernel(dist);
    }

    // Perform mean shift iterations
    _mean_shift() {
        const N = this._N;
        const D = this._D;
        const points = this._points;
        const metric = this._parameters.metric;
        //const bandwidth = this._bandwidth;
        const kernel = this._kernel_weight.bind(this);
        const tolerance = this._tolerance;

        for (let iter = 0; iter < this._max_iter; ++iter) {
            let max_shift = 0;
            // For each point compute shift
            for (let i = 0; i < N; ++i) {
                const row_i = points.row(i);
                let sum_weights = 0;
                const weighted_sum = new Float64Array(D);
                for (let j = 0; j < N; ++j) {
                    const row_j = points.row(j);
                    const dist = metric(row_i, row_j);
                    const weight = kernel(dist);
                    sum_weights += weight;
                    for (let d = 0; d < D; ++d) {
                        weighted_sum[d] += weight * row_j[d];
                    }
                }
                if (sum_weights === 0) {
                    // No neighbors within kernel, shift is zero
                    //const shift = new Float64Array(D);
                    // Compute shift magnitude
                    const shift_norm = Math.sqrt(weighted_sum.reduce((acc, v) => acc + v * v, 0));
                    max_shift = Math.max(max_shift, shift_norm);
                } else {
                    const shift = new Float64Array(D);
                    for (let d = 0; d < D; ++d) {
                        shift[d] = weighted_sum[d] / sum_weights - row_i[d];
                    }
                    const shift_norm = Math.sqrt(shift.reduce((acc, v) => acc + v * v, 0));
                    max_shift = Math.max(max_shift, shift_norm);
                    // Update point
                    for (let d = 0; d < D; ++d) {
                        row_i[d] += shift[d];
                    }
                }
            }
            if (max_shift < tolerance) {
                // Converged
                break;
            }
        }
    }

    // After convergence, assign clusters based on nearest mode
    _assign_clusters() {
        const N = this._N;
        const metric = this._parameters.metric;
        const bandwidth = this._bandwidth;

        // Group points that converged to the same mode
        // Two points are in the same mode if they're within bandwidth/2 of each other
        const mode_threshold = bandwidth * 0.5;
        /** @type {number[][]} */
        const modes = []; // Each mode contains indices of points in that mode
        const point_to_mode = new Array(N).fill(-1);

        for (let i = 0; i < N; ++i) {
            if (point_to_mode[i] !== -1) continue; // Already assigned to a mode

            const row_i = this._points.row(i);
            const mode = [i];
            point_to_mode[i] = modes.length;

            // Find all points close to this mode
            for (let j = i + 1; j < N; ++j) {
                if (point_to_mode[j] !== -1) continue;

                const row_j = this._points.row(j);
                const dist = metric(row_i, row_j);

                if (dist < mode_threshold) {
                    mode.push(j);
                    point_to_mode[j] = modes.length;
                }
            }

            modes.push(mode);
        }

        // Build final clusters - each mode becomes a cluster
        /** @type {number[][]} */
        const clusters = [];
        const cluster_ids = new Array(N).fill(-1);

        for (let mode_idx = 0; mode_idx < modes.length; ++mode_idx) {
            const mode = modes[mode_idx];
            clusters.push([...mode]);
            for (const point_idx of mode) {
                cluster_ids[point_idx] = mode_idx;
            }
        }

        this._clusters = cluster_ids;
        this._cluster_list = clusters;
    }

    /**
     * @returns {number[][]}
     */
    get_clusters() {
        // Ensure algorithm has been run
        if (!this._cluster_list) {
            this._mean_shift();
            this._assign_clusters();
        }
        return /** @type {number[][]} */ (this._cluster_list);
    }

    /**
     *
     * @returns {number[]}
     */
    get_cluster_list() {
        if (!this._clusters) {
            this._mean_shift();
            this._assign_clusters();
        }
        return /** @type {number[]} */ (this._clusters);
    }
}
