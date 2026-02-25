import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { Clustering } from "./Clustering.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersOptics } from "./index.js" */

/** @typedef {Object} DBEntry
 * @property {Float64Array} element
 * @property {number} index
 * @property {number} [reachability_distance]
 * @property {boolean} processed
 * @property {DBEntry[]} [neighbors]
 */

/**
 * OPTICS (Ordering Points To Identify the Clustering Structure)
 *
 * A density-based clustering algorithm that extends DBSCAN. It handles clusters of varying
 * densities and produces a reachability plot that can be used to extract clusters.
 *
 * @class
 * @extends Clustering<ParametersOptics>
 * @category Clustering
 */
export class OPTICS extends Clustering {
    /**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     *
     * @param {InputType} points - The data.
     * @param {Partial<ParametersOptics>} [parameters={}]
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersOptics} */ (
                Object.assign({ epsilon: 1, min_points: 4, metric: euclidean }, parameters)
            ),
        );
        const matrix = this._matrix;
        /**
         * @private
         * @type {DBEntry[]}
         */
        this._ordered_list = [];
        const ordered_list = this._ordered_list;
        /** @type {number[][]} */
        this._clusters = [];
        const clusters = this._clusters;

        const N = this._N;

        /**
         * @private
         * @type {DBEntry[]}
         */
        this._DB = new Array(N).fill(0).map((_, i) => {
            return {
                element: matrix.row(i),
                index: i,
                reachability_distance: undefined,
                processed: false,
            };
        });
        const DB = this._DB;

        this._cluster_index = 0;
        let cluster_index = this._cluster_index;

        for (const p of DB) {
            if (p.processed) continue;
            p.neighbors = this._get_neighbors(p);
            p.processed = true;
            clusters.push([p.index]);
            cluster_index = clusters.length - 1;
            ordered_list.push(p);
            if (this._core_distance(p) !== undefined) {
                const seeds = new Heap(null, (d) => d.reachability_distance, "min");
                this._update(p, seeds);
                this._expand_cluster(seeds, clusters[cluster_index]);
            }
        }
    }

    /**
     * @private
     * @param {DBEntry} p - A point of the data.
     * @returns {DBEntry[]} An array consisting of the `epsilon`-neighborhood of `p`.
     */
    _get_neighbors(p) {
        if (p?.neighbors) return p.neighbors;
        const DB = this._DB;
        const metric = this._parameters.metric;
        const epsilon = this._parameters.epsilon;
        const neighbors = [];
        for (const q of DB) {
            if (q.index === p.index) continue;
            if (metric(p.element, q.element) <= epsilon) {
                neighbors.push(q);
            }
        }
        return neighbors;
    }

    /**
     * @private
     * @param {DBEntry} p - A point of `matrix`.
     * @returns {number|undefined} The distance to the `min_points`-th nearest point of `p`, or undefined if the
     *   `epsilon`-neighborhood has fewer elements than `min_points`.
     */
    _core_distance(p) {
        const min_points = this._parameters.min_points;
        const metric = this._parameters.metric;
        // Need min_points - 1 other points plus the point itself
        if (!p.neighbors || p.neighbors.length < min_points - 1) {
            return undefined;
        }
        // Sort neighbors by distance to find the MinPts-th closest
        const sortedNeighbors = p.neighbors.toSorted(
            (a, b) => metric(p.element, a.element) - metric(p.element, b.element),
        );
        // MinPts-th closest is at index min_points - 2 (0-indexed, excluding p itself)
        return metric(p.element, sortedNeighbors[min_points - 2].element);
    }

    /**
     * Updates the reachability distance of the points.
     *
     * @private
     * @param {DBEntry} p
     * @param {Heap<DBEntry>} seeds
     */
    _update(p, seeds) {
        const metric = this._parameters.metric;
        const core_distance = this._core_distance(p);
        // If p is not a core point, don't update seeds
        if (core_distance === undefined) {
            return;
        }
        const neighbors = this._get_neighbors(p); //p.neighbors;
        for (const q of neighbors) {
            if (q.processed) continue;
            const new_reachability_distance = Math.max(core_distance, metric(p.element, q.element));
            //if (q.reachability_distance == undefined) { // q is not in seeds
            if (seeds.raw_data().findIndex((d) => d.element === q) < 0) {
                q.reachability_distance = new_reachability_distance;
                seeds.push(q);
            } else {
                // q is in seeds
                if (new_reachability_distance < (q.reachability_distance ?? Infinity)) {
                    q.reachability_distance = new_reachability_distance;
                    seeds = Heap.heapify(seeds.data(), (d) => d.reachability_distance ?? Infinity, "min"); // seeds change key =/
                }
            }
        }
    }

    /**
     * Expands the `cluster` with points in `seeds`.
     *
     * @private
     * @param {Heap<DBEntry>} seeds
     * @param {number[]} cluster
     */
    _expand_cluster(seeds, cluster) {
        const ordered_list = this._ordered_list;
        while (!seeds.empty) {
            const q = /** @type {{ element: DBEntry, value: number}} */ (seeds.pop()).element;
            q.neighbors = this._get_neighbors(q);
            q.processed = true;
            cluster.push(q.index);
            ordered_list.push(q);
            if (this._core_distance(q) !== undefined) {
                this._update(q, seeds);
                // Recursive call removed - while loop handles iteration correctly
            }
        }
    }

    /**
     * Returns an array of clusters.
     *
     * @returns {number[][]} Array of clusters with the indices of the rows in given `matrix`.
     */
    get_clusters() {
        const clusters = [];
        const outliers = [];
        const min_points = this._parameters.min_points;
        for (const cluster of this._clusters) {
            if (cluster.length < min_points) {
                outliers.push(...cluster);
            } else {
                clusters.push(cluster);
            }
        }
        clusters.push(outliers);
        return clusters;
    }

    /**
     * @returns {number[]} Returns an array, where the ith entry defines the cluster affirmation of the ith point of
     *   given data. (-1 stands for outlier)
     */
    get_cluster_list() {
        const N = this._matrix.shape[0];
        /** @type {number[]} */
        const result = new Array(N).fill(0);
        const clusters = this.get_clusters();
        for (let i = 0, n = clusters.length; i < n; ++i) {
            const cluster = clusters[i];
            for (const index of cluster) {
                result[index] = i < n - 1 ? i : -1;
            }
        }
        return result;
    }
}
