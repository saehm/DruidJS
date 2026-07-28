import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { quickselect, Randomizer } from "../util/index.js";
import { KNN } from "./KNN.js";

/** @import { ParametersNaiveKNN } from "./index.js" */
/** @import { Metric } from "../metrics/index.js" */

/**
 * Naive KNN implementation performing an exhaustive scan.
 *
 * Every query measures the distance to all `N` elements and then selects the `k` smallest with
 * {@link quickselect}, which is O(N) on average — cheaper than sorting all `N` and far cheaper than
 * the N heaps of N entries this class used to build up front.
 *
 * The N x N distance matrix is only materialized when it actually pays off: a `"precomputed"` index
 * is handed one directly, and {@link NaiveKNN#search_by_index} builds one lazily on first use so
 * that building a whole kNN graph costs a single pass over the pairs instead of one per query.
 * Callers that only ever use {@link NaiveKNN#search} never allocate it.
 *
 * Best suited for small datasets, or when a distance matrix is already available. For larger `N`
 * prefer an approximate index such as {@link HNSW}, {@link Annoy}, or {@link NNDescent}.
 *
 * @template {number[] | Float64Array} T
 * @category KNN
 * @class
 * @extends KNN<T, ParametersNaiveKNN>
 */
export class NaiveKNN extends KNN {
    /**
     * Number of indexed elements.
     *
     * @type {number}
     */
    _N;
    /**
     * Pairwise distances. Supplied by the caller when `metric` is `"precomputed"`, built on demand
     * by {@link NaiveKNN#search_by_index} otherwise, and `null` until then.
     *
     * @type {Matrix | null}
     */
    _D;

    /**
     * Generates a KNN list with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KNN list
     * @param {ParametersNaiveKNN} parameters
     */
    constructor(elements, parameters = {}) {
        const params = Object.assign({ metric: euclidean, seed: 1212 }, parameters);
        super(elements, params);
        // Seeded pivots for `quickselect`, so a selection among tied distances is reproducible.
        /** @type {Randomizer} */
        this._randomizer = new Randomizer(params.seed);
        const elements_any = /** @type {any} */ (this._elements);
        this._N = elements_any instanceof Matrix ? elements_any.shape[0] : this._elements.length;
        this._D =
            this._parameters.metric === "precomputed"
                ? Matrix.from(/** @type {number[][] | Float64Array[]} */ (elements_any))
                : null;
    }

    /**
     * Reads the element stored at `i`, which may live in a `Matrix` or a plain array.
     *
     * @private
     * @param {number} i
     * @returns {T}
     */
    _element_at(i) {
        const elements = /** @type {any} */ (this._elements);
        return /** @type {T} */ (elements instanceof Matrix ? elements.row(i) : elements[i]);
    }

    /**
     * Selects the `k` elements closest to a query from its distances to every element.
     *
     * QuickSelect partitions in O(N) average time, after which only the `k` selected entries are
     * sorted. Ties break on the element index so the result never depends on the pivots drawn.
     *
     * @private
     * @param {ArrayLike<number>} distances - Distance from the query to each element, by index.
     * @param {number} k - Number of neighbors to return.
     * @returns {{ element: T; index: number; distance: number }[]} The `k` nearest, closest first.
     */
    _k_smallest(distances, k) {
        const N = this._N;
        const size = Math.min(Math.max(Math.floor(k), 0), N);
        if (size === 0) return [];

        /** @type {(a: number, b: number) => number} */
        const compare = (a, b) => distances[a] - distances[b] || a - b;
        /** @type {number[]} */
        const indices = new Array(N);
        for (let i = 0; i < N; ++i) indices[i] = i;

        if (size < N) {
            // Leaves the `size` smallest in indices[0..size-1], in no particular order.
            quickselect(indices, this._randomizer, size - 1, compare);
            indices.length = size;
        }
        indices.sort(compare);

        return indices.map((index) => ({
            element: this._element_at(index),
            index,
            distance: distances[index],
        }));
    }

    /**
     * Returns the pairwise distance matrix, computing it once on first use.
     *
     * @private
     * @returns {Matrix}
     */
    _distance_matrix() {
        if (this._D === null) {
            this._D = distance_matrix(
                /** @type {number[][] | Float64Array[]} */ (/** @type {any} */ (this._elements)),
                /** @type {Metric} */ (this._parameters.metric),
            );
        }
        return this._D;
    }

    /**
     * @param {number} i - Index of the query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest
     *   neighbors, closest first. `i` itself is included, at distance 0.
     */
    search_by_index(i, k = 5) {
        return this._k_smallest(this._distance_matrix().row(i), k);
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
        const metric = /** @type {Metric} */ (this._parameters.metric);
        const N = this._N;
        const distances = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            distances[i] = metric(t, this._element_at(i));
        }
        return this._k_smallest(distances, k);
    }
}
