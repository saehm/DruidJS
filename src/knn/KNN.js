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
     * @abstract
     * @param {number} i
     * @param {number} k
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_by_index(i, k) {
        i;
        k;
        throw new Error("The function search_by_index must be implemented!");
    }
}
