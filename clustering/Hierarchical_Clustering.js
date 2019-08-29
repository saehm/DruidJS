import { euclidean } from "../metrics/index";

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
     * @param {Matrix} matrix 
     * @param {("single"|"complete"|"average")} [linkage = "single"] 
     * @param {Function} [metric = euclidean] 
     * @returns {Hierarchical_Clustering}
     */
    constructor(matrix, linkage="single", metric=euclidean) {
        this._id = 0;
        this._matrix = matrix;
        this._metric = metric;
        this._linkage = linkage;
        this.init();
        this.root = this.do();
        return this;
    }

    /**
     * 
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(value, type="distance") {
        let clusters = [];
        let accessor;
        switch (type) {
            case "distance":
                accessor = d => d.dist;
                break;
            case "depth":
                accessor = d => d.depth;
                break;
            default:
                throw "invalid type";
        }
        this._traverse(this.root, accessor, value, clusters)
        return clusters
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
            result.push(node.leaves())
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
        const n = this._n = A.shape[0];
        const d_min = this._d_min = new Float64Array(n);
        const distance_matrix = this._distance_matrix = new Array(n);
        for (let i = 0; i < n; ++i) {
            d_min[i] = 0;
            distance_matrix[i] = new Float64Array(n);
            for (let j = 0; j < n; ++j) {
                distance_matrix[i][j] = i === j ? Infinity : metric(A.row(i), A.row(j))
                if (distance_matrix[i][d_min[i]] > distance_matrix[i][j]) {
                    d_min[i] = j;
                }
            }
        }
        const clusters = this._clusters = new Array(n);
        const c_size = this._c_size = new Uint16Array(n);
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
                if (D[i][d_min[i]] < D[c1][d_min[c1]]) {
                    c1 = i;
                }
            }
            let c2 = d_min[c1];
            let c1_cluster = clusters[c1][0];
            let c2_cluster = clusters[c2][0];
            let new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D[c1][c2]);
            clusters[c1].unshift(new_cluster);
            c_size[c1] += c_size[c2];
            for (let j = 0; j < n; ++j) {
                switch(linkage) {
                    case "single":
                        if (D[c1][j] > D[c2][j]) {
                            D[j][c1] = D[c1][j] = D[c2][j];
                        }
                        break;
                    case "complete":
                        if (D[c1][j] < D[c2][j]) {
                            D[j][c1] = D[c1][j] = D[c2][j];
                        }
                        break;
                    case "average":
                        D[j][c1] = D[c1][j] = (c_size[c1] * D[c1][j] + c_size[c2] * D[c2][j]) / (c_size[c1] + c_size[j]);
                        break;
                }
            }
            D[c1][c1] = Infinity;
            for (let i = 0; i < n; ++i) {
                D[i][c2] = D[c2][i] = Infinity;
            }
            for (let j = 0; j < n; ++j) {
                if (d_min[j] === c2) {
                    d_min[j] = c1;
                }
                if (D[c1][j] < D[c1][d_min[c1]]) {
                    d_min[c1] = j;
                }
            }
            root = new_cluster
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
        this.size = size != null ? size : left.size + right.size;
        this.depth = depth != null ? depth : 1 + Math.max(left.depth, right.depth);
        this.centroid = centroid != null ? centroid : this._calculate_centroid(left, right);
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
        if (this.isLeaf) return [this.index];
        const left = this.left;
        const right = this.right;
        return (left.isLeaf ? [left.index] : left.leaves())
            .concat(right.isLeaf ? [right.index] : right.leaves())
    }
}