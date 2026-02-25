import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { KNN } from "./KNN.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersBallTree } from "./index.js" */

/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} ElementWithIndex
 * @property {number} index
 * @property {T} element
 */

/**
 * Ball Tree for efficient nearest neighbor search.
 *
 * A Ball Tree is a metric tree that partitions points into a nested set of
 * hyperspheres (balls). It is particularly effective for high-dimensional
 * data and supports any valid metric.
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersBallTree>
 */
export class BallTree extends KNN {
    /**
     * Generates a BallTree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the BallTree
     * @param {ParametersBallTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     */
    constructor(elements, parameters = { metric: euclidean, seed: 1212 }) {
        super(elements, Object.assign({ seed: 1212 }, parameters));
        /**
         * @private
         * @type {BallTreeNode<T> | BallTreeLeaf<T>}
         */
        this._root = this._construct(elements.map((element, index) => ({ index, element })));
    }

    /** @returns {Metric} */
    get _metric() {
        return this._parameters.metric;
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @returns {BallTreeNode<T> | BallTreeLeaf<T>} Root of balltree.
     */
    _construct(elements) {
        if (elements.length === 1) {
            return new BallTreeLeaf(elements);
        } else {
            const c = this._greatest_spread(elements);
            const sorted_elements = elements.sort((a, b) => a.element[c] - b.element[c]);
            const n = sorted_elements.length;
            const p_index = Math.floor(n / 2);
            const p = sorted_elements[p_index];
            const L = sorted_elements.slice(0, p_index);
            const R = sorted_elements.slice(p_index, n);
            const radius = Math.max(...elements.map((d) => this._metric(p.element, d.element)));
            let B;
            if (L.length > 0 && R.length > 0) {
                B = new BallTreeNode(p, this._construct(L), this._construct(R), radius);
            } else {
                B = new BallTreeLeaf(elements);
            }
            return B;
        }
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} B
     * @returns {number}
     */
    _greatest_spread(B) {
        const d = B[0].element.length;
        const start = new Array(d);

        for (let i = 0; i < d; ++i) {
            start[i] = [Infinity, -Infinity];
        }

        let spread = B.reduce((acc, current) => {
            for (let i = 0; i < d; ++i) {
                acc[i][0] = Math.min(acc[i][0], current.element[i]);
                acc[i][1] = Math.max(acc[i][1], current.element[i]);
            }
            return acc;
        }, start);
        spread = spread.map((d) => d[1] - d[0]);

        let c = 0;
        for (let i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c;
    }

    /**
     * @param {number} i
     * @param {number} k
     */
    search_by_index(i, k = 5) {
        return this.search(this._elements[i], k);
    }

    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        /** @type {Heap<ElementWithIndex<T>>} */
        const heap = new Heap(null, (d) => this._metric(d.element, t), "max");
        this._search(t, k, heap, this._root);

        // Convert heap to result array
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (heap.length > 0) {
            const item = /** @type {{ element: ElementWithIndex<T>; value: number }} */ (heap.pop());
            result.push({
                element: item.element.element,
                index: item.element.index,
                distance: item.value,
            });
        }
        return result.reverse(); // Reverse to get closest first
    }

    /**
     * @private
     * @param {T} t - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {Heap<ElementWithIndex<T>>} Q - Heap consists of the currently found `k` nearest neighbors.
     * @param {BallTreeNode<T> | BallTreeLeaf<T>} B
     */
    _search(t, k, Q, B) {
        if (!B) return;

        if (B instanceof BallTreeNode) {
            const dist_to_pivot = this._metric(t, B.pivot.element);
            if (Q.length >= k && dist_to_pivot - B.radius >= (Q.first?.value ?? -Infinity)) {
                return;
            }

            const c1 = B.child1;
            const c2 = B.child2;

            let d1 = Infinity;
            let d2 = Infinity;

            if (c1 instanceof BallTreeNode) d1 = this._metric(t, c1.pivot.element);
            else if (c1 instanceof BallTreeLeaf) d1 = this._metric(t, c1.points[0].element);

            if (c2 instanceof BallTreeNode) d2 = this._metric(t, c2.pivot.element);
            else if (c2 instanceof BallTreeLeaf) d2 = this._metric(t, c2.points[0].element);

            if (d1 < d2) {
                if (c1) this._search(t, k, Q, c1);
                if (c2) this._search(t, k, Q, c2);
            } else {
                if (c2) this._search(t, k, Q, c2);
                if (c1) this._search(t, k, Q, c1);
            }
        } else if (B instanceof BallTreeLeaf) {
            for (let i = 0, n = B.points.length; i < n; ++i) {
                const p = B.points[i];
                const dist = this._metric(p.element, t);
                if (Q.length < k) {
                    Q.push(p);
                } else if (dist < (Q.first?.value ?? Infinity)) {
                    Q.pop();
                    Q.push(p);
                }
            }
        }
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeNode {
    /**
     * @param {ElementWithIndex<T>} pivot
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child1
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child2
     * @param {number} radius
     */
    constructor(pivot, child1 = null, child2 = null, radius = 0) {
        this.pivot = pivot;
        this.child1 = child1;
        this.child2 = child2;
        this.radius = radius;
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeLeaf {
    /** @param {ElementWithIndex<T>[]} points */
    constructor(points) {
        this.points = points;
    }
}
