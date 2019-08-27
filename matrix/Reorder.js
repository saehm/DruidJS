import { euclidean } from "../metrics/index";
import { Hierarchical_Clustering } from "../clustering/index";
import { Matrix } from "./Matrix";
import { simultaneous_poweriteration } from "../linear_algebra/index";
import { linspace } from "../matrix/index";

// https://github.com/jdfekete/reorder.js/blob/master/reorder.v1.js
export class Reorder{
    constructor(A) {
        this._A = A;
        this._optimal_leaf_order = null;
    }

    get available() {
        return [
            "optimal_leaf_order",
            "spectral_order",
        ]
    }

    reorder(type="optimal_leaf_order", metric=euclidean) {
        let result = null;
        switch (type) {
            case "optimal_leaf_order":
                this._optimal_leaf_order = new Optimal_Leaf_Order(this._A, metric);
                result = this._optimal_leaf_order.ordering;
                break;

            case "spectral_order":
                this._spectral_order = new Spectral_Order(this._A, metric);
                result = this._spectral_order.ordering;
                break;
            case "barycenter_order":
                break;
        }
        return result;
    }
}

//class Barycenter_Order{    
//}

class Spectral_Order{
    constructor(A, metric=euclidean) {
        this._A = A;
        this._metric = metric;
        this._N = A.shape[0];
        const fiedler_vector = this._fiedler_vector(A);
        this._ordering = linspace(0, this._N - 1).sort((a, b) => fiedler_vector[a] - fiedler_vector[b]);
        return this;
    }

    get ordering() {
        return this._ordering;
    }

    _fiedler_vector(B) {
        const g = this._gershgorin_bound(B);
        const N = B.shape[0];
        const B_hat = new Matrix(N, N, (i, j) => i === j ? g - B.entry(i, j) : -B.entry(i, j))
        const eig = simultaneous_poweriteration(B_hat, 2);
        return eig.eigenvectors[1]
    }

    _gershgorin_bound(B) {
        let max = 0;
        let N = B.shape[0];
        for (let i = 0; i < N; ++i) {
            let t = B.entry(i, i);
            for (let j = 0; j < N; ++j) {
                if (i !== j) {
                    t += Math.abs(B.entry(i, j));
                }
            }
            max = max > t ? max : t;
        }
        return max;
    }
}

class Optimal_Leaf_Order{
    constructor(A, metric=euclidean) {
        this._A = A;
        const N = A.shape[0];
        const hclust = this._hclust = new Hierarchical_Clustering(A, "complete", metric)
        const distance_matrix = this._distance_matrix = new Array(N);
        for (let i = 0; i < N; ++i) {
            distance_matrix[i] = new Float64Array(N);
            for (let j = 0; j < N; ++j) {
                distance_matrix[i][j] = i === j ? Infinity : metric(A.row(i), A.row(j))
            }
        }
        this._order_map = new Map();
        let min = Infinity;
        this._optimal_order = null;
        let left = hclust.root.left.leaves();
        let right = hclust.root.right.leaves();

        for (let i = 0, n = left.length; i < n; ++i) {
            for (let j = 0, m = right.length; j < m; ++j) {
                let so = this.order(hclust.root, left[i], right[j]);
                if (so[0] < min) {
                    min = so[0];
                    this._optimal_order = so[1];
                }
            }
        }
        return this;
    }

    get ordering() {
        return this._optimal_order;
    }

    order(v, i, j) {
        const order_map = this._order_map;
        const key = `k${v.id}-${i}-${j}`; // ugly key
        /*if (key in order_map) 
            return order_map[key];
        return (order_map[key] = this._order(v, i, j))*/
        /*if (order_map.has(v)) {
            const v_map = order_map.get(v)
            if (v_map.has(`${i},${j}`)) {
                return v_map.get(`${i},${j}`)
            } else {
                let value = this._order(v, i, j);
                v_map.set(`${i},${j}`, value);
                return value;
            }
        } else {
            let value = this._order(v, i, j);
            const v_map = new Map();
            v_map.set(`${i},${j}`, value);
            order_map.set(v, v_map);
            return value;
        }*/
        if (order_map.has(key)) {
            return order_map.get(key);
        } else {
            let value = this._order(v, i, j);
            order_map.set(key, value);
            return value;
        }
    }

    _order(v, i, j) {
        if (v.isLeaf) {
            return [0, [v.index]];
        }
        const D = this._distance_matrix;
        let l = v.left;
        let r = v.right;
        let L = l ? l.leaves() : [];
        let R = r ? r.leaves() : [];
        let w;
        let x;
        if (L.indexOf(i) !== -1 && R.indexOf(j) !== -1) {
            w = l; 
            x = r;
        } else if (R.indexOf(i) !== -1 && L.indexOf(j) !== -1) {
            w = r;
            x = l;
        } else {
            throw "Node is not common ancestor of i and j";
        }

        let Wl = w.left ? w.left.leaves() : [];
        let Wr = w.right ? w.right.leaves() : [];
        let Ks = Wr.indexOf(i) != -1 ? Wl : Wr;
        if (Ks.length === 0) { 
            Ks = [i];
        }

        let Xl = x.left ? x.left.leaves() : [];
        let Xr = x.right ? x.right.leaves() : [];
        let Ls = Xr.indexOf(j) != -1 ? Xl : Xr;
        if (Ls.length === 0) {
            Ls = [j];
        }

        let min = Infinity;
        let optimal_order = [];
        for (let k = 0, Ks_length = Ks.length; k < Ks_length; ++k) {
            let w_min = this.order(w, i, Ks[k]);
            for (let m = 0, Ls_length = Ls.length; m < Ls_length; ++m) {
                let x_min = this.order(x, Ls[m], j);
                let dist = w_min[0] + D[Ks[k]][Ls[m]] + x_min[0];
                if (dist < min) {
                    min = dist;
                    optimal_order = w_min[1].concat(x_min[1]);
                }
            }
        }

        return [min, optimal_order];
    }
}