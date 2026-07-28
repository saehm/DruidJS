import { Randomizer } from "../util/index.js";

/**
 * Base class for all K-Nearest Neighbors (KNN) search algorithms.
 *
 * Provides a common interface for elements management and search operations.
 *
 * @abstract
 * @category KNN
 * @template {number[] | Float64Array} T - Type of elements
 * @template {Object} Para - Type of parameters
 * @class
 */
export class KNN {
    /** @type {T[]} */
    _elements;
    /** @type {Para} */
    _parameters;
    /** @type {"typed" | "array"} */
    _type;
    /**
     * Seeded source of randomness shared by every index. Construction is randomized — the trees
     * pick quickselect pivots from it — so the `seed` parameter is what makes a built index, and
     * therefore its query results, reproducible.
     *
     * @type {Randomizer}
     */
    _randomizer;

    /**
     * @param {T[]} elements
     * @param {Para} parameters
     */
    constructor(elements, parameters) {
        if (elements.length === 0) throw new Error("Elements needs to contain at least one element!");
        if (elements[0] instanceof Float64Array) {
            this._type = "typed";
        } else {
            this._type = "array";
        }
        this._parameters = parameters;
        this._elements = elements;
        this._randomizer = new Randomizer(/** @type {{ seed?: number } | undefined} */ (parameters)?.seed ?? 1212);
    }

    /**
     * @abstract
     * @param {T} t
     * @param {number} k
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(t, k) {
        t;
        k;
        throw new Error("The function search must be implemented!");
    }

    /**
     * Searches the `k` nearest neighbors of the element stored at index `i`.
     *
     * The queried element is never part of the result. It is trivially its own closest neighbor at
     * distance 0, which is never what a caller asking "what is this point near?" wants, so every
     * caller used to strip it back out — each in its own, subtly different way. Note the asymmetry
     * with `search`: an arbitrary query point has no "self" to exclude, so there `k` means
     * "k results", while here it means "k neighbors".
     *
     * The self match is removed **by index** — not by position, and not by looking for a zero
     * distance. Position is wrong because an approximate index may order ties differently or miss
     * the element altogether (one extra candidate is requested to cover that), and a zero distance
     * is wrong because genuine duplicate points share it and must survive.
     *
     * @param {number} i - Index of the query element.
     * @param {number} [k=5] - Number of neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} The `k` nearest *other* elements,
     *   closest first. Empty when `i` is out of range.
     */
    search_by_index(i, k = 5) {
        const elements = this._elements;
        if (i < 0 || i >= elements.length) return [];
        const element = elements[i];
        if (!element) return [];
        return this.search(element, k + 1)
            .filter((neighbor) => neighbor.index !== i)
            .slice(0, k);
    }
}
