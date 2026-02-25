import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { Clustering } from "./Clustering.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersHierarchicalClustering } from "./index.js" */

/**
 * Hierarchical Clustering
 *
 * A bottom-up approach (agglomerative) to clustering that builds a tree of clusters (dendrogram).
 * Supports different linkage criteria: single, complete, and average.
 *
 * @class
 * @extends Clustering<ParametersHierarchicalClustering>
 * @category Clustering
 */
export class HierarchicalClustering extends Clustering {
    /** @type {Cluster | null} */
    root = null;

    /**
     * @param {InputType} points - Data or distance matrix if metric is 'precomputed'
     * @param {Partial<ParametersHierarchicalClustering>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersHierarchicalClustering} */ (
                Object.assign({ linkage: "complete", metric: euclidean }, parameters)
            ),
        );
        this._id = 0;
        if (this._parameters.metric === "precomputed" && this._matrix.shape[0] !== this._matrix.shape[1]) {
            throw new Error("If metric is 'precomputed', then matrix has to be square!");
        }

        const metric = this._parameters.metric;
        const A = this._matrix;
        const N = this._N;
        this._d_min = new Float64Array(N);
        const d_min = this._d_min;
        let distance_matrix;
        if (metric !== "precomputed") {
            distance_matrix = new Matrix(N, N, Infinity);
            for (let i = 0; i < N; ++i) {
                distance_matrix.set_entry(i, i, 0);
                d_min[i] = i; // temporary
                const Ai = A.row(i);
                for (let j = i + 1; j < N; ++j) {
                    const dist = metric(Ai, A.row(j));
                    distance_matrix.set_entry(i, j, dist);
                    distance_matrix.set_entry(j, i, dist);
                }
            }
            for (let i = 0; i < N; i++) {
                let min_j = 0;
                let min_d = Infinity;
                for (let j = 0; j < N; j++) {
                    if (i === j) continue;
                    const d = distance_matrix.entry(i, j);
                    if (d < min_d) {
                        min_d = d;
                        min_j = j;
                    }
                }
                d_min[i] = min_j;
            }
        } else {
            distance_matrix = this._matrix.clone();
            for (let i = 0; i < N; ++i) {
                distance_matrix.set_entry(i, i, 0);
                d_min[i] = i === 0 ? 1 : 0;
                for (let j = 0; j < N; ++j) {
                    if (i === j) continue;
                    if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        this._distance_matrix = distance_matrix;
        this._clusters = new Array(N);
        const clusters = this._clusters;
        this._c_size = new Uint16Array(N);
        const c_size = this._c_size;
        for (let i = 0; i < N; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        const D = this._distance_matrix;
        const linkage = this._parameters.linkage;
        const p_max = N - 1;
        for (let p = 0; p < p_max; ++p) {
            let c1 = -1;
            let min_dist = Infinity;
            for (let i = 0; i < N; ++i) {
                if (D.entry(i, i) === Infinity) continue;
                const dist = D.entry(i, d_min[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    c1 = i;
                }
            }
            if (c1 === -1) break;

            const c2 = d_min[c1];
            const c1_cluster = clusters[c1][0];
            const c2_cluster = clusters[c2][0];
            const c1_cluster_indices = c1_cluster.isLeaf ? [c1_cluster.index] : c1_cluster.index;
            const c2_cluster_indices = c2_cluster.isLeaf ? [c2_cluster.index] : c2_cluster.index;
            const indices = c1_cluster_indices.concat(c2_cluster_indices);
            const new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D.entry(c1, c2), null, indices);
            c1_cluster.parent = new_cluster;
            c2_cluster.parent = new_cluster;
            clusters[c1].unshift(new_cluster);

            const size1 = c_size[c1];
            const size2 = c_size[c2];
            c_size[c1] += size2;

            for (let j = 0; j < N; ++j) {
                if (j === c1 || j === c2 || D.entry(j, j) === Infinity) continue;
                const D_c1_j = D.entry(c1, j);
                const D_c2_j = D.entry(c2, j);
                let value;
                switch (linkage) {
                    case "single":
                        value = Math.min(D_c1_j, D_c2_j);
                        break;
                    case "complete":
                        value = Math.max(D_c1_j, D_c2_j);
                        break;
                    case "average":
                        value = (size1 * D_c1_j + size2 * D_c2_j) / (size1 + size2);
                        break;
                }
                D.set_entry(j, c1, value);
                D.set_entry(c1, j, value);
            }

            D.set_entry(c2, c2, Infinity);
            for (let i = 0; i < N; ++i) {
                D.set_entry(i, c2, Infinity);
                D.set_entry(c2, i, Infinity);
            }

            // Update d_min for all rows
            for (let i = 0; i < N; i++) {
                if (D.entry(i, i) === Infinity) continue;
                if (d_min[i] === c1 || d_min[i] === c2 || i === c1) {
                    let min_j = 0;
                    let min_d = Infinity;
                    for (let j = 0; j < N; j++) {
                        if (i === j || D.entry(j, j) === Infinity) continue;
                        const d = D.entry(i, j);
                        if (d < min_d) {
                            min_d = d;
                            min_j = j;
                        }
                    }
                    d_min[i] = min_j;
                } else {
                    if (D.entry(i, c1) < D.entry(i, d_min[i])) {
                        d_min[i] = c1;
                    }
                }
            }

            this.root = new_cluster;
        }
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {Cluster[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters_raw(value, type = "distance") {
        /** @type {Cluster[][]} */
        const clusters = [];
        /** @type {(d: {dist: number, depth: number}) => number} */
        let accessor;
        switch (type) {
            case "distance":
                accessor = (d) => d.dist;
                break;
            case "depth":
                accessor = (d) => d.depth;
                break;
            default:
                throw new Error("invalid type");
        }
        this._traverse(/** @type {Cluster} */ (this.root), accessor, value, clusters);
        return clusters;
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters(value, type = "distance") {
        /** @type {Cluster[][]} */
        const clusters = [];
        /** @type {(d: {dist: number, depth: number}) => number} */
        let accessor;
        switch (type) {
            case "distance":
                accessor = (d) => d.dist;
                break;
            case "depth":
                accessor = (d) => d.depth;
                break;
            default:
                throw new Error("invalid type");
        }
        if (this.root) this._traverse(this.root, accessor, value, clusters);
        return clusters.map((cluster) => cluster.map((d) => d.index));
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[]} - Array of clusters with the indices of the rows in given points.
     */
    get_cluster_list(value, type = "distance") {
        const clusters = this.get_clusters(value, type);
        /** @type {number[]} */
        const list = new Array(this._N).fill(0);
        for (let i = 0; i < clusters.length; ++i) {
            const cluster = clusters[i];
            for (let j = 0; j < cluster.length; ++j) {
                const index = cluster[j];
                list[index] = i;
            }
        }
        return list;
    }

    /**
     * @private
     * @param {Cluster} node
     * @param {(d: {dist: number, depth: number}) => number} f
     * @param {number} value
     * @param {Cluster[][]} result
     */
    _traverse(node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        } else {
            if (node.left) this._traverse(node.left, f, value, result);
            if (node.right) this._traverse(node.right, f, value, result);
        }
    }
}

/** @private */
class Cluster {
    /**@type {number} */
    size;
    /**@type {number} */
    depth;
    /**@type {Cluster | null} */
    parent;

    /**
     *
     * @param {number} id
     * @param {Cluster?} left
     * @param {Cluster?} right
     * @param {number} dist
     * @param {Float64Array?} centroid
     * @param {number} index
     * @param {number} [size]
     * @param {number} [depth]
     */
    constructor(id, left, right, dist, centroid, index, size, depth) {
        this.id = id;
        this.left = left;
        this.right = right;
        this.dist = dist;
        this.index = index;
        if (size) {
            this.size = size;
        } else {
            if (!left || !right) throw new Error("If size is not given, left & right cannot be null!");
            this.size = left.size + right.size;
        }

        if (depth !== undefined) {
            this.depth = depth;
        } else {
            if (!left || !right) throw new Error("If depth is not given, left & right cannot be null!");
            this.depth = Math.max(left.depth, right.depth) + 1;
        }

        if (centroid !== undefined && centroid !== null) {
            this.centroid = centroid;
        } else {
            if (!left || !right) throw new Error("If centroid is not given, left & right cannot be null!");

            this.centroid = this._calculate_centroid(left, right);
        }

        this.parent = null;
    }

    /**
     *
     * @param {Cluster} left
     * @param {Cluster} right
     * @returns {Float64Array}
     */
    _calculate_centroid(left, right) {
        const l_size = left.size;
        const r_size = right.size;
        const l_centroid = left.centroid;
        const r_centroid = right.centroid;
        const size = this.size;
        const n = left.centroid.length;
        const new_centroid = new Float64Array(n);
        for (let i = 0; i < n; ++i) {
            new_centroid[i] = (l_size * l_centroid[i] + r_size * r_centroid[i]) / size;
        }
        return new_centroid;
    }

    get isLeaf() {
        return this.depth === 0;
    }

    /**
     *
     * @returns {Cluster[]}
     */
    leaves() {
        if (this.isLeaf) return [this];
        const left = this.left;
        const right = this.right;
        return (left ? (left.isLeaf ? [left] : left.leaves()) : []).concat(
            right ? (right.isLeaf ? [right] : right.leaves()) : [],
        );
    }

    /**
     *
     * @returns {Cluster[]}
     */
    descendants() {
        if (this.isLeaf) return [this];
        const left_descendants = this.left ? this.left.descendants() : [];
        const right_descendants = this.right ? this.right.descendants() : [];
        return left_descendants.concat(right_descendants).concat([this]);
    }
}
