import { Heap } from "../datastructure/index.js";
import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { KNN } from "./KNN.js";

/** @import { ParametersNaiveKNN } from "./index.js" */

/**
 * Naive KNN implementation using a distance matrix.
 *
 * This implementation pre-computes the entire distance matrix and performs
 * an exhaustive search. Best suited for small datasets or when a distance
 * matrix is already available.
 *
 * @template {number[] | Float64Array} T
 * @category KNN
 * @class
 * @extends KNN<T, ParametersNaiveKNN>
 */
export class NaiveKNN extends KNN {
    /**
     * Generates a KNN list with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KNN list
     * @param {ParametersNaiveKNN} parameters
     */
    constructor(elements, parameters = {}) {
        const params = Object.assign({ metric: euclidean, seed: 1212 }, parameters);
        super(elements, params);
        const N =
            this._elements instanceof Matrix ? /** @type {any} */ (this._elements).shape[0] : this._elements.length;
        if (this._parameters.metric === "precomputed") {
            this._D = Matrix.from(/** @type {number[][] | Float64Array[]} */ (/** @type {any} */ (this._elements)));
        } else {
            this._D = distance_matrix(
                /** @type {number[][] | Float64Array[]} */ (this._elements),
                this._parameters.metric,
            );
        }

        /** @type {Heap<{ value: number; index: number }>[]} */
        this.KNN = [];
        for (let row = 0; row < N; ++row) {
            const distances = this._D.row(row);
            /** @type {Heap<{ value: number; index: number }>} */
            const H = new Heap(null, (d) => d.value, "min");
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
     * @param {number} i
     * @param {number} k
     */
    search_by_index(i, k = 5) {
        if (this._parameters.metric === "precomputed") {
            const H = this.KNN[i];
            /** @type {{ element: T; index: number; distance: number }[]} */
            const result = [];
            const data = H.toArray(); // Get array representation
            const temp_heap = new Heap(data, (d) => d.value, "min");
            const N =
                this._elements instanceof Matrix ? /** @type {any} */ (this._elements).shape[0] : this._elements.length;
            for (let j = 0; j < Math.min(k, N); ++j) {
                const node = temp_heap.pop();
                if (!node) break;
                result.push({
                    element: /** @type {T} */ (
                        this._elements instanceof Matrix
                            ? /** @type {any} */ (this._elements).row(node.element.index)
                            : this._elements[node.element.index]
                    ),
                    index: /** @type {number} */ (node.element.index),
                    distance: /** @type {number} */ (node.value),
                });
            }
            return result;
        }
        return this.search(
            /** @type {T} */ (
                this._elements instanceof Matrix ? /** @type {any} */ (this._elements).row(i) : this._elements[i]
            ),
            k,
        );
    }

    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        if (this._parameters.metric === "precomputed") {
            throw new Error("Search by query element is only possible when not using a precomputed distance matrix!");
        }
        /** @type {import("../metrics/index.js").Metric} */
        const metric = /** @type {any} */ (this._parameters.metric);

        const isMatrix = this._elements instanceof Matrix;
        const elementsAny = /** @type {any} */ (this._elements);
        const N = isMatrix ? elementsAny.shape[0] : this._elements.length;

        // Compute distances from query to ALL points
        const distances = [];
        for (let i = 0; i < N; i++) {
            const element = /** @type {T} */ (isMatrix ? elementsAny.row(i) : this._elements[i]);
            distances.push({
                element: element,
                index: i,
                distance: metric(t, element),
            });
        }

        // Sort by distance and return k nearest
        distances.sort((a, b) => a.distance - b.distance);
        return distances.slice(0, k);
    }
}
