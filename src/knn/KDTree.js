import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { KNN } from "./KNN.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersKDTree } from "./index.js" */

/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} ElementWithIndex
 * @property {number} index
 * @property {T} element
 */

/**
 * KD-Tree (K-dimensional Tree) for efficient nearest neighbor search.
 *
 * KD-Trees partition k-dimensional space by recursively splitting along coordinate axes.
 * At each level, the tree splits points based on the median of the coordinate with the largest spread.
 * This creates a balanced binary tree structure that enables efficient O(log n) search on average.
 *
 * Best suited for:
 * - Low to moderate dimensional data (d < 20-30)
 * - When exact nearest neighbors are needed
 * - When dimensionality is not too high
 *
 * Performance degrades in high dimensions (curse of dimensionality) where approximate
 * methods like HNSW or LSH become more effective.
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersKDTree>
 * @see {@link https://en.wikipedia.org/wiki/K-d_tree}
 */
export class KDTree extends KNN {
    /**
     * Generates a KD-Tree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KD-Tree
     * @param {ParametersKDTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
     */
    constructor(elements, parameters = { metric: euclidean, seed: 1212 }) {
        super(elements, Object.assign({ seed: 1212 }, parameters));
        /**
         * @private
         * @type {KDTreeNode<T> | KDTreeLeaf<T> | null}
         */
        this._root = this._construct(
            elements.map((element, index) => ({ index, element })),
            0,
        );
    }

    /** @returns {Metric} */
    get _metric() {
        return this._parameters.metric;
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @param {number} depth - Current depth in the tree (determines splitting axis)
     * @returns {KDTreeNode<T> | KDTreeLeaf<T> | null} Root of KD-Tree.
     */
    _construct(elements, depth) {
        if (elements.length === 0) {
            return null;
        }

        if (elements.length === 1) {
            return new KDTreeLeaf(elements[0]);
        }

        const k = elements[0].element.length;
        const axis = depth % k;

        // Sort by the splitting axis and find median
        elements.sort((a, b) => a.element[axis] - b.element[axis]);
        const medianIndex = Math.floor(elements.length / 2);
        const medianPoint = elements[medianIndex];

        // Recursively build left and right subtrees
        const leftElements = elements.slice(0, medianIndex);
        const rightElements = elements.slice(medianIndex + 1);

        const left = this._construct(leftElements, depth + 1);
        const right = this._construct(rightElements, depth + 1);

        return new KDTreeNode(medianPoint, axis, left, right);
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
        /** @type {Heap<{ point: ElementWithIndex<T>; distance: number }>} */
        const best = new Heap(null, (d) => d.distance, "max");

        this._search_recursive(t, k, this._root, best);

        // Convert heap to result array (closest first)
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (best.length > 0) {
            const item = /** @type {{ element: { point: ElementWithIndex<T>; distance: number }; value: number }} */ (
                best.pop()
            );
            result.push({
                element: item.element.point.element,
                index: item.element.point.index,
                distance: item.value,
            });
        }
        return result.reverse();
    }

    /**
     * @private
     * @param {T} target - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} node - Current node.
     * @param {Heap<{ point: ElementWithIndex<T>; distance: number }>} best - Heap of k best found so far.
     */
    _search_recursive(target, k, node, best) {
        if (node === null) return;

        if (node instanceof KDTreeLeaf) {
            const dist = this._metric(target, node.point.element);
            if (best.length < k) {
                best.push({ point: node.point, distance: dist });
            } else if (dist < (best.first?.value ?? Infinity)) {
                best.pop();
                best.push({ point: node.point, distance: dist });
            }
            return;
        }

        // Node is an internal node
        const axis = node.axis;
        const point = node.point;
        const pointValue = point.element[axis];
        const targetValue = target[axis];

        // Determine which subtree to search first
        const firstSubtree = targetValue < pointValue ? node.left : node.right;
        const secondSubtree = targetValue < pointValue ? node.right : node.left;

        // Search the nearer subtree
        this._search_recursive(target, k, firstSubtree, best);

        // Check if we need to search the other subtree
        // The hyperplane could contain closer points
        const distToHyperplane = Math.abs(targetValue - pointValue);
        const currentMaxDist = best.first?.value ?? Infinity;

        // Calculate distance to current point
        const distToPoint = this._metric(target, point.element);
        if (best.length < k) {
            best.push({ point: point, distance: distToPoint });
        } else if (distToPoint < currentMaxDist) {
            best.pop();
            best.push({ point: point, distance: distToPoint });
        }

        // Check if we need to explore the other side of the hyperplane
        if (best.length < k || distToHyperplane < (best.first?.value ?? Infinity)) {
            this._search_recursive(target, k, secondSubtree, best);
        }
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class KDTreeNode {
    /**
     * @param {ElementWithIndex<T>} point
     * @param {number} axis - The splitting axis
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} left
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} right
     */
    constructor(point, axis, left = null, right = null) {
        this.point = point;
        this.axis = axis;
        this.left = left;
        this.right = right;
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class KDTreeLeaf {
    /**
     * @param {ElementWithIndex<T>} point
     */
    constructor(point) {
        this.point = point;
    }
}
