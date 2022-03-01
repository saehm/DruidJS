/**
 * @class
 * @alias Heap
 */
export class Heap {
    /**
     * Creates a Heap from an Array
     * @param {Array|Set} elements - Contains the elements for the Heap.
     * @param {Function=} [accessor = (d) => d] - Function returns the value of the element.
     * @param {(String=|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     */
    static heapify(elements: any[] | Set<any>, accessor?: Function | undefined, comparator?: (String?: number) => any): Heap;
    /**
     * A heap is a datastructure holding its elements in a specific way, so that the top element would be the first entry of an ordered list.
     * @constructor
     * @memberof module:datastructure
     * @alias Heap
     * @param {Array=} elements - Contains the elements for the Heap. {@link elements} can be null.
     * @param {Function} [accessor = (d) => d] - Function returns the value of the element.
     * @param {("min"|"max"|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     * @see {@link https://en.wikipedia.org/wiki/Binary_heap}
     */
    constructor(elements?: any[] | undefined, accessor?: Function, comparator?: ("min" | "max" | Function));
    _accessor: Function;
    _container: any[];
    _comparator: Function;
    /**
     * Swaps elements of container array.
     * @private
     * @param {Number} index_a
     * @param {Number} index_b
     */
    private _swap;
    /**
     * @private
     */
    private _heapify_up;
    /**
     * Pushes the element to the heap.
     * @param {} element
     * @returns {Heap}
     */
    push(element: any): Heap;
    /**
     * @private
     * @param {Number} [start_index = 0]
     */
    private _heapify_down;
    /**
     * Removes and returns the top entry of the heap.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    pop(): any;
    /**
     * Returns the top entry of the heap without removing it.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    get first(): any;
    /**
     * Yields the raw data
     * @yields {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    iterate(): Generator<any, void, unknown>;
    /**
     * Returns the heap as ordered array.
     * @returns {Array} Array consisting the elements ordered by {@link comparator}.
     */
    toArray(): any[];
    /**
     * Returns elements of container array.
     * @returns {Array} Array consisting the elements.
     */
    data(): any[];
    /**
     * Returns the container array.
     * @returns {Array} The container array.
     */
    raw_data(): any[];
    /**
     * The size of the heap.
     * @returns {Number}
     */
    get length(): number;
    /**
     * Returns false if the the heap has entries, true if the heap has no entries.
     * @returns {Boolean}
     */
    get empty(): boolean;
}
//# sourceMappingURL=Heap.d.ts.map