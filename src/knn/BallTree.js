import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { max, quickselectByAxis } from "../util/index.js";
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
 * Every ball stores a center and the distance from it to its furthest member, which bounds the
 * subtree from below by `d(t, center) - radius`. That bound holds for any metric by the triangle
 * inequality, so the center does not need to be one of the indexed points — it is the centroid,
 * which keeps radii tighter than an arbitrary member would.
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
     * @param {Partial<ParametersBallTree>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     */
    constructor(elements, parameters = {}) {
        super(elements, Object.assign({ metric: euclidean, leaf_size: 32, seed: 1212 }, parameters));
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
     * @returns {Float64Array} Componentwise mean of `elements`.
     */
    _centroid(elements) {
        const d = elements[0].element.length;
        const n = elements.length;
        const center = new Float64Array(d);
        for (let i = 0; i < n; ++i) {
            const element = elements[i].element;
            for (let j = 0; j < d; ++j) center[j] += element[j];
        }
        for (let j = 0; j < d; ++j) center[j] /= n;
        return center;
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @returns {BallTreeNode<T> | BallTreeLeaf<T>} Root of balltree.
     */
    _construct(elements) {
        const center = this._centroid(elements);
        // `max` from utils rather than `Math.max(...array)`: spreading the array blew the argument
        // limit and threw `Maximum call stack size exceeded` from roughly 125k elements upward.
        const radius = max(elements.map((d) => this._metric(center, d.element)));

        if (elements.length <= this._parameters.leaf_size) {
            return new BallTreeLeaf(center, radius, elements);
        }

        const c = this._greatest_spread(elements);
        const p_index = elements.length >> 1;
        quickselectByAxis(elements, this._randomizer, p_index, c);
        // Split at `p_index`, so each element lands on exactly one side. Keeping a separate pivot
        // element and leaving it in the right half stored every interior point twice, which made a
        // full traversal cost 2N distance evaluations.
        const L = elements.slice(0, p_index);
        const R = elements.slice(p_index);
        return new BallTreeNode(center, radius, this._construct(L), this._construct(R));
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
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        // The accessor reads the distance the traversal already computed. Recomputing the metric
        // here made every heap insertion cost a second distance evaluation of the same pair.
        /** @type {Heap<{ point: ElementWithIndex<T>; distance: number }>} */
        const heap = new Heap(null, (d) => d.distance, "max");
        this._search(t, k, heap, this._root);

        // Convert heap to result array
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (heap.length > 0) {
            const item = /** @type {{ element: { point: ElementWithIndex<T> }; value: number }} */ (heap.pop());
            result.push({
                element: item.element.point.element,
                index: item.element.point.index,
                distance: item.value,
            });
        }
        return result.reverse(); // Reverse to get closest first
    }

    /**
     * @private
     * @param {T} t - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {Heap<{ point: ElementWithIndex<T>; distance: number }>} Q - Heap consists of the currently found `k`
     *   nearest neighbors.
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} B
     * @param {number} [known_distance] - Distance from `t` to this ball's center, when the caller
     *   already measured it while ordering its children.
     */
    _search(t, k, Q, B, known_distance) {
        if (!B) return;

        // Applies to leaves too, which previously scanned every member unconditionally.
        const distance = known_distance ?? this._metric(t, B.center);
        if (Q.length >= k && distance - B.radius >= (Q.first?.value ?? -Infinity)) {
            return;
        }

        if (B instanceof BallTreeLeaf) {
            const points = B.points;
            for (let i = 0, n = points.length; i < n; ++i) {
                const p = points[i];
                const dist = this._metric(p.element, t);
                if (Q.length < k) {
                    Q.push({ point: p, distance: dist });
                } else if (dist < (Q.first?.value ?? Infinity)) {
                    Q.pop();
                    Q.push({ point: p, distance: dist });
                }
            }
            return;
        }

        const c1 = B.child1;
        const c2 = B.child2;

        // Measured once here and handed to the recursive call, which would otherwise repeat it.
        const d1 = c1 ? this._metric(t, c1.center) : Infinity;
        const d2 = c2 ? this._metric(t, c2.center) : Infinity;

        if (d1 < d2) {
            this._search(t, k, Q, c1, d1);
            this._search(t, k, Q, c2, d2);
        } else {
            this._search(t, k, Q, c2, d2);
            this._search(t, k, Q, c1, d1);
        }
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeNode {
    /**
     * @param {Float64Array} center
     * @param {number} radius
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child1
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child2
     */
    constructor(center, radius, child1 = null, child2 = null) {
        this.center = center;
        this.radius = radius;
        this.child1 = child1;
        this.child2 = child2;
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeLeaf {
    /**
     * @param {Float64Array} center
     * @param {number} radius
     * @param {ElementWithIndex<T>[]} points
     */
    constructor(center, radius, points) {
        this.center = center;
        this.radius = radius;
        this.points = points;
    }
}
