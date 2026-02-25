/** @import { Comparator } from "./index.js" */

/**
 * @template T
 * @class
 * @category Data Structures
 */
export class Heap {
    /** @type {{ element: T; value: number }[]} */
    _container;

    /** @type {Comparator} */
    _comparator;

    /**
     * A heap is a datastructure holding its elements in a specific way, so that the top element would be the first
     * entry of an ordered list.
     *
     * @param {T[]?} elements - Contains the elements for the Heap. `elements` can be null.
     * @param {(d: T) => number} accessor - Function returns the value of the element.
     * @param {"min" | "max" | Comparator} [comparator="min"] - Function returning true or false
     *   defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
     *   Max_heap). Default is `"min"`
     * @see {@link https://en.wikipedia.org/wiki/Binary_heap}
     */
    constructor(elements = null, accessor, comparator = "min") {
        /** @type {(d: T) => number} */
        this._accessor = accessor;
        this._container = [];
        if (comparator === "min") {
            this._comparator = (a, b) => a < b;
        } else if (comparator === "max") {
            this._comparator = (a, b) => a > b;
        } else {
            this._comparator = comparator;
        }
        if (elements) {
            this._container = [];
            for (const e of elements) {
                this._container.push({
                    element: e,
                    value: accessor(e),
                });
            }
            for (let i = Math.floor(elements.length / 2 - 1); i >= 0; --i) {
                this._heapify_down(i);
            }
        }
    }

    /**
     * Creates a Heap from an Array
     *
     * @template T
     * @param {T[]} elements - Contains the elements for the Heap.
     * @param {(d: T) => number} accessor - Function returns the value of the element.
     * @param {"min" | "max" | Comparator} [comparator="min"] - Function returning true or false
     *   defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
     *   Max_heap). Default is `"min"`
     * @returns {Heap<T>}
     */
    static heapify(elements, accessor, comparator = "min") {
        const heap = new Heap(null, accessor, comparator);
        const container = heap._container;
        for (const e of elements) {
            container.push({
                element: e,
                value: accessor(e),
            });
        }
        for (let i = Math.floor(elements.length / 2 - 1); i >= 0; --i) {
            heap._heapify_down(i);
        }
        return heap;
    }

    /**
     * Swaps elements of container array.
     *
     * @private
     * @param {number} index_a
     * @param {number} index_b
     */
    _swap(index_a, index_b) {
        const container = this._container;
        [container[index_b], container[index_a]] = [container[index_a], container[index_b]];
        return;
    }

    /** @private */
    _heapify_up() {
        const container = this._container;
        let index = container.length - 1;
        while (index > 0) {
            const parentIndex = Math.floor((index - 1) / 2);
            if (!this._comparator(container[index].value, container[parentIndex].value)) {
                break;
            } else {
                this._swap(parentIndex, index);
                index = parentIndex;
            }
        }
    }

    /**
     * Pushes the element to the heap.
     *
     * @param {T} element
     * @returns {Heap<T>}
     */
    push(element) {
        const value = this._accessor(element);
        //const node = new Node(element, value);
        const node = { element: element, value: value };
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    /**
     * @private
     * @param {Number} [start_index=0] Default is `0`
     */
    _heapify_down(start_index = 0) {
        const container = this._container;
        const comparator = this._comparator;
        const length = container.length;
        const left = 2 * start_index + 1;
        const right = 2 * start_index + 2;
        let index = start_index;
        if (index >= length) throw "index higher than length";
        if (left < length && comparator(container[left].value, container[index].value)) {
            index = left;
        }
        if (right < length && comparator(container[right].value, container[index].value)) {
            index = right;
        }
        if (index !== start_index) {
            this._swap(start_index, index);
            this._heapify_down(index);
        }
    }

    /**
     * Removes and returns the top entry of the heap.
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`}).
     */
    pop() {
        const container = this._container;
        if (container.length === 0) {
            return null;
        } else if (container.length === 1) {
            const item = container.pop();
            if (!item) throw new Error("Cannot happen!");
            return item;
        }
        this._swap(0, container.length - 1);
        const item = container.pop();
        this._heapify_down();
        return item ?? null;
    }

    /**
     * Returns the top entry of the heap without removing it.
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`).
     */
    get first() {
        return this._container.length > 0 ? this._container[0] : null;
    }

    /**
     * Yields the raw data
     *
     * @yields {T} Object consists of the element and its value (computed by `accessor`}).
     */
    *iterate() {
        for (let i = 0, n = this._container.length; i < n; ++i) {
            yield this._container[i].element;
        }
    }

    /**
     * Returns the heap as ordered array.
     *
     * @returns {T[]} Array consisting the elements ordered by `comparator`.
     */
    toArray() {
        return this._container.sort((a, b) => (this._comparator(a.value, b.value) ? -1 : 1)).map((d) => d.element);
    }

    /**
     * Returns elements of container array.
     *
     * @returns {T[]} Array consisting the elements.
     */
    data() {
        return this._container.map((d) => d.element);
    }

    /**
     * Returns the container array.
     *
     * @returns {{ element: T; value: number }[]} The container array.
     */
    raw_data() {
        return this._container;
    }

    /**
     * The size of the heap.
     *
     * @returns {number}
     */
    get length() {
        return this._container.length;
    }

    /**
     * Returns false if the the heap has entries, true if the heap has no entries.
     *
     * @returns {boolean}
     */
    get empty() {
        return this.length === 0;
    }
}
