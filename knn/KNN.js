import { euclidean } from "../metrics/index.js";
import { Heap } from "../datastructure/index.js";
import { distance_matrix, Matrix } from "../matrix/index.js";

/**
 * @class
 * @alias KNN
 */
export class KNN {
    /**
     * Generates a KNN list with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias KNN
     * @param {Array=} elements - Elements which should be added to the KNN list
     * @param {Function|"precomputed"} [metric = euclidean] metric is either precomputed or a function to use: (a, b) => distance
     * @returns {KNN}
     */
    constructor(elements=null, metric=euclidean) {
        this._metric = metric;
        this._elements = elements instanceof Matrix ? elements : Matrix.from(elements);
        const N = this._elements.shape[0];
        if (metric === "precomputed") {
            this._D = this._elements.clone();
        } else {
            this._D = distance_matrix(this._elements, metric);
        }
        this.KNN = [];
        for (let row = 0; row < N; ++row) {
            const distances = this._D.row(row);
            const H = new Heap(null, d => d.value, "min");
            for (let j = 0; j < N; ++j) {
                H.push({
                    value: distances[j],
                    index: j,
                });
            }
            this.KNN.push(H);
        }
    }

    /**
     * 
     * @param {Array|Number} t - query element or index.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t, k = 5) {
        const metric = this._metric;
        const KNN = this.KNN;
        let H;
        if (Array.isArray(t)) {
            if (this._metric == "precomputed") {
                throw "Search by query element is only possible when not using a precomputed distance matrix!"
            } 
            const elements = this._elements;
            const N = KNN.length;
            let nearest_element_index = null;
            let nearest_dist = Infinity;
            for (let i = 0; i < N; ++i) {
                const element = elements.row(i);
                const dist = metric(t, element);
                if (dist < nearest_dist) {
                    nearest_element_index = i;
                    nearest_dist = dist;
                }
            }
            H = KNN[nearest_element_index];
        } else if (Number.isInteger(t)) {
            H = KNN[t]
        }

        let result = []
        for (let i = 0; i < k; ++i) {
            result.push(H.pop())
        }
        result.forEach(res => H.push(res.element))
        return result
    }    
}
