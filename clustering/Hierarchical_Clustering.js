import { euclidean } from "../metrics/index.js";
import { Matrix } from "../matrix/index.js";
/**
 * @class
 * @alias Hierarchical_Clustering
 */
export class Hierarchical_Clustering {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias Hierarchical_Clustering
     * @todo needs restructuring.
     * @param {Matrix} - Data or distance matrix if metric is 'precomputed'
     * @param {("single"|"complete"|"average")} [linkage = "complete"]
     * @param {Function|"precomputed"} [metric = euclidean]
     * @returns {Hierarchical_Clustering}
     */
    constructor(matrix, linkage = "complete", metric = euclidean) {
        this._id = 0;
        this._matrix = matrix instanceof Matrix ? matrix : Matrix.from(matrix);
        this._metric = metric;
        this._linkage = linkage;
        if (metric === "precomputed" && this._matrix.shape[0] !== this._matrix.shape[1]) {
            throw new Error("If metric is 'precomputed', then matrix has to be square!");
        }
        this.init();
        this.root = this.do();
        return this;
    }

    /**
     *
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(value, type = "distance") {
        let clusters = [];
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
        this._traverse(this.root, accessor, value, clusters);
        return clusters;
    }

    /**
     * @private
     * @param {} node
     * @param {*} f
     * @param {*} value
     * @param {*} result
     */
    _traverse(node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        } else {
            this._traverse(node.left, f, value, result);
            this._traverse(node.right, f, value, result);
        }
    }

    /**
     * computes the tree.
     */
    init() {
        const metric = this._metric;
        const A = this._matrix;
        const n = (this._n = A.shape[0]);
        const d_min = (this._d_min = new Float64Array(n));
        let distance_matrix;
        if (metric !== "precomputed") {
            distance_matrix = new Matrix(n, n, 0); //new Array(n);
            for (let i = 0; i < n; ++i) {
                d_min[i] = 0;
                //distance_matrix[i] = new Float64Array(n);
                for (let j = 0; j < n; ++j) {
                    distance_matrix.set_entry(i, j, i === j ? Infinity : metric(A.row(i), A.row(j)));
                    if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        } else {
            distance_matrix = this._matrix.clone();
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    if (i === j) {
                        distance_matrix.set_entry(i, j, Infinity);
                    } else if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        this._distance_matrix = distance_matrix;
        const clusters = (this._clusters = new Array(n));
        const c_size = (this._c_size = new Uint16Array(n));
        for (let i = 0; i < n; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        return this;
    }

    /**
     * computes the tree.
     */
    do() {
        const n = this._n;
        const d_min = this._d_min;
        const D = this._distance_matrix;
        const clusters = this._clusters;
        const c_size = this._c_size;
        const linkage = this._linkage;
        let root = null;
        for (let p = 0, p_max = n - 1; p < p_max; ++p) {
            let c1 = 0;
            for (let i = 0; i < n; ++i) {
                let D_i_min = D.entry(i, d_min[i]);
                for (let j = i + 1; j < n; ++j) {
                    if (D_i_min > D.entry(i, j)) {
                        d_min[i] = j;
                        D_i_min = D.entry(i, d_min[i]);
                    }
                }
            }
            for (let i = 0; i < n; ++i) {
                if (D.entry(i, d_min[i]) < D.entry(c1, d_min[c1])) {
                    c1 = i;
                }
            }
            let c2 = d_min[c1];
            let c1_cluster = clusters[c1][0];
            let c2_cluster = clusters[c2][0];
            let c1_cluster_indices = c1_cluster.isLeaf ? [c1_cluster.index] : c1_cluster.index;
            let c2_cluster_indices = c2_cluster.isLeaf ? [c2_cluster.index] : c2_cluster.index;
            let indices = c1_cluster_indices.concat(c2_cluster_indices);
            let new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D.entry(c1, c2), null, indices);
            c1_cluster.parent = new_cluster;
            c2_cluster.parent = new_cluster;
            clusters[c1].unshift(new_cluster);
            c_size[c1] += c_size[c2];
            for (let j = 0; j < n; ++j) {
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
                        value = (c_size[c1] * D_c1_j + c_size[c2] * D_c2_j) / (c_size[c1] + c_size[j]);
                        break;
                }
                D.set_entry(j, c1, value);
                D.set_entry(c1, j, value);
            }

            D.set_entry(c1, c1, Infinity);
            for (let i = 0; i < n; ++i) {
                D.set_entry(i, c2, Infinity);
                D.set_entry(c2, i, Infinity);
            }

            /* for (let j = 0; j < n; ++j) {
                if (d_min[j] === c2) {
                    d_min[j] = c1;
                }
                if (D.entry(c1, j) < D.entry(c1, d_min[c1])) {
                    d_min[c1] = j;
                }
            } */
            root = new_cluster;
        }
        return root;
    }
}

class Cluster {
    constructor(id, left, right, dist, centroid, index, size, depth) {
        this.id = id;
        this.left = left;
        this.right = right;
        this.dist = dist;
        this.index = index;
        this.size = size ?? left.size + right.size;
        this.depth = depth ?? 1 + Math.max(left.depth, right.depth);
        this.centroid = centroid ?? this._calculate_centroid(left, right);
        this.parent = null;
        return this;
    }

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

    leaves() {
        if (this.isLeaf) return [this];
        const left = this.left;
        const right = this.right;
        return (left.isLeaf ? [left] : left.leaves()).concat(right.isLeaf ? [right] : right.leaves());
    }

    descendants() {
        if (this.isLeaf) return [this];
        const left_descendants = this.left.descendants();
        const right_descendants = this.right.descendants();
        return left_descendants.concat(right_descendants).concat([this]);
    }
}
