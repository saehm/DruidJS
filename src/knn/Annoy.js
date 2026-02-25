import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { KNN } from "./KNN.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersAnnoy } from "./index.js" */

/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} AnnoyNode
 * @property {boolean} isLeaf - Whether this is a leaf node
 * @property {number[]} indices - Indices of points in this node (leaf) or children (internal)
 * @property {number[]} normal - Hyperplane normal vector (internal nodes only)
 * @property {number} offset - Hyperplane offset (internal nodes only)
 * @property {AnnoyNode<T> | null} left - Left child (internal nodes only)
 * @property {AnnoyNode<T> | null} right - Right child (internal nodes only)
 */

/**
 * Annoy-style (Approximate Nearest Neighbors Oh Yeah) implementation using Random Projection Trees.
 *
 * This implementation builds multiple random projection trees where each tree randomly selects
 * two points and splits the space based on a hyperplane equidistant between them.
 *
 * Key features:
 * - Multiple random projection trees for better recall
 * - Each tree uses random hyperplanes for splitting
 * - Priority queue search for better recall
 * - Combines results from all trees
 *
 * Best suited for:
 * - High-dimensional data
 * - Approximate nearest neighbor search
 * - Large datasets
 * - When high recall is needed with approximate methods
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersAnnoy>
 * @see {@link https://github.com/spotify/annoy}
 * @see {@link https://erikbern.com/2015/09/24/nearest-neighbors-and-vector-models-epilogue-curse-of-dimensionality.html}
 */
export class Annoy extends KNN {
    /**
     * Creates a new Annoy-style index with random projection trees.
     *
     * @param {T[]} elements - Elements to index
     * @param {ParametersAnnoy} [parameters={}] - Configuration parameters
     */
    constructor(
        elements,
        parameters = {
            metric: euclidean,
            numTrees: 10,
            maxPointsPerLeaf: 10,
            seed: 1212,
        },
    ) {
        // Handle empty initialization - use dummy element
        const hasElements = elements && elements.length > 0;
        const firstElement = /** @type {T} */ (hasElements ? elements[0] : new Float64Array([0]));

        super([firstElement], parameters);

        this._metric = this._parameters.metric ?? euclidean;
        this._numTrees = this._parameters.numTrees ?? 10;
        this._maxPointsPerLeaf = this._parameters.maxPointsPerLeaf ?? 10;
        this._seed = this._parameters.seed ?? 1212;
        this._randomizer = new Randomizer(this._seed);

        /**
         * @private
         * @type {AnnoyNode<T>[]}
         */
        this._trees = [];

        // Build trees
        if (hasElements) {
            // Reset elements and rebuild properly
            /** @type {T[]} */
            this._elements = [];
            this._trees = [];
            this.add(elements);
        }
    }

    /**
     * Get the number of trees in the index.
     * @returns {number}
     */
    get num_trees() {
        return this._trees.length;
    }

    /**
     * Get the total number of nodes in all trees.
     * @returns {number}
     */
    get num_nodes() {
        let total = 0;
        for (const tree of this._trees) {
            total += this._countNodes(tree);
        }
        return total;
    }

    /**
     * @private
     * @param {any} node
     * @returns {number}
     */
    _countNodes(node) {
        if (!node) return 0;
        return 1 + this._countNodes(node.left) + this._countNodes(node.right);
    }

    /**
     * Add elements to the Annoy index.
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements) {
        // Extend elements array
        this._elements = this._elements.concat(elements);

        // Rebuild all trees with new elements
        this._trees = [];
        this._buildTrees();

        return this;
    }

    /**
     * Build all random projection trees.
     * @private
     */
    _buildTrees() {
        const elements = this._elements;
        const n = elements.length;

        for (let t = 0; t < this._numTrees; t++) {
            // Create index array for this tree
            const indices = Array.from({ length: n }, (_, i) => i);
            const tree = this._buildTreeRecursive(indices);
            this._trees.push(tree);
        }
    }

    /**
     * Recursively build a random projection tree.
     * @private
     * @param {number[]} indices - Indices of elements to include
     * @returns {AnnoyNode<T>}
     */
    _buildTreeRecursive(indices) {
        const elements = this._elements;

        // Base case: small enough to be a leaf
        if (indices.length <= this._maxPointsPerLeaf) {
            return {
                isLeaf: true,
                indices: indices,
                normal: [],
                offset: 0,
                left: null,
                right: null,
            };
        }

        // Select two random points to define the splitting hyperplane
        const idx1 = indices[Math.floor(this._randomizer.random * indices.length)];
        const idx2 = indices[Math.floor(this._randomizer.random * indices.length)];

        const point1 = elements[idx1];
        const point2 = elements[idx2];

        // Compute normal vector (point2 - point1)
        const dim = point1.length;
        /** @type {number[]} */
        const normal = new Array(dim);
        for (let i = 0; i < dim; i++) {
            normal[i] = point2[i] - point1[i];
        }

        // Normalize
        let norm = 0;
        for (let i = 0; i < dim; i++) {
            norm += normal[i] * normal[i];
        }
        norm = Math.sqrt(norm);

        if (norm > 1e-10) {
            for (let i = 0; i < dim; i++) {
                normal[i] /= norm;
            }
        }

        // Compute midpoint and offset
        /** @type {number[]} */
        const midpoint = new Array(dim);
        for (let i = 0; i < dim; i++) {
            midpoint[i] = (point1[i] + point2[i]) / 2;
        }

        // Compute offset: dot(normal, midpoint)
        let offset = 0;
        for (let i = 0; i < dim; i++) {
            offset += normal[i] * midpoint[i];
        }

        // Split points based on which side of hyperplane they fall
        const leftIndices = [];
        const rightIndices = [];

        for (const idx of indices) {
            const point = elements[idx];
            let dot = 0;
            for (let i = 0; i < dim; i++) {
                dot += normal[i] * point[i];
            }

            if (dot < offset) {
                leftIndices.push(idx);
            } else {
                rightIndices.push(idx);
            }
        }

        // Handle edge case where all points fall on one side
        if (leftIndices.length === 0 || rightIndices.length === 0) {
            return {
                isLeaf: true,
                indices: indices,
                normal: [],
                offset: 0,
                left: null,
                right: null,
            };
        }

        // Recursively build subtrees
        const left = this._buildTreeRecursive(leftIndices);
        const right = this._buildTreeRecursive(rightIndices);

        return {
            isLeaf: false,
            indices: [],
            normal: normal,
            offset: offset,
            left: left,
            right: right,
        };
    }

    /**
     * Compute distance from point to hyperplane.
     * @private
     * @param {T} point
     * @param {number[]} normal
     * @param {number} offset
     * @returns {number} Signed distance (positive = right side, negative = left side)
     */
    _distanceToHyperplane(point, normal, offset) {
        let dot = 0;
        for (let i = 0; i < point.length; i++) {
            dot += normal[i] * point[i];
        }
        return dot - offset;
    }

    /**
     * Search for k approximate nearest neighbors.
     * @param {T} query
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(query, k = 5) {
        const metric = this._metric;
        const elements = this._elements;

        if (elements.length === 0) return [];

        // Collect candidates from all trees using priority queue
        const candidates = new Set();

        // Collect more candidates for better recall
        // Search at least k * numTrees * 2 candidates
        const minCandidates = Math.min(k * this._numTrees * 3, elements.length);

        for (const tree of this._trees) {
            this._searchTreePriority(tree, query, candidates, minCandidates);
        }

        // Compute exact distances for all candidates
        /** @type {Heap<{ index: number; distance: number }>} */
        const best = new Heap(null, (d) => d.distance, "max");

        for (const idx of candidates) {
            const element = elements[idx];
            if (!element || element.length !== query.length) continue;

            const dist = metric(query, element);

            if (best.length < k) {
                best.push({ index: idx, distance: dist });
            } else if (dist < (best.first?.value ?? Infinity)) {
                best.pop();
                best.push({ index: idx, distance: dist });
            }
        }

        // If we still don't have enough candidates, do a linear scan fallback
        if (best.length < k) {
            for (let i = 0; i < elements.length && best.length < k; i++) {
                if (candidates.has(i)) continue;

                const element = elements[i];
                if (!element || element.length !== query.length) continue;

                const dist = metric(query, element);
                best.push({ index: i, distance: dist });
            }
        }

        // Convert to result format
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (best.length > 0) {
            const item = /** @type {{ element: { index: number; distance: number }; value: number }} */ (best.pop());
            result.push({
                element: elements[item.element.index],
                index: item.element.index,
                distance: item.value,
            });
        }

        return result.reverse();
    }

    /**
     * Search tree using priority queue for better recall.
     * Explores nodes in order of distance to hyperplane.
     * @private
     * @param {AnnoyNode<T>} node
     * @param {T} query
     * @param {Set<number>} candidates
     * @param {number} maxCandidates
     */
    _searchTreePriority(node, query, candidates, maxCandidates) {
        if (!node) return;

        // Priority queue entry: { node, distance }
        /** @type {Heap<{ node: AnnoyNode<T>; dist: number }>} */
        const pq = new Heap(null, (d) => d.dist, "min");
        pq.push({ node: node, dist: 0 });

        while (!pq.empty && candidates.size < maxCandidates) {
            const entry = pq.pop();
            if (!entry) continue;

            const currentNode = entry.element.node;

            // Leaf node: add all points
            if (currentNode.isLeaf) {
                for (const idx of currentNode.indices) {
                    candidates.add(idx);
                    if (candidates.size >= maxCandidates) return;
                }
                continue;
            }

            // Internal node: compute distance to hyperplane
            const dist = this._distanceToHyperplane(query, currentNode.normal, currentNode.offset);

            // Determine which side is closer
            const closerSide = dist < 0 ? currentNode.left : currentNode.right;
            const fartherSide = dist < 0 ? currentNode.right : currentNode.left;

            // Add closer side with priority 0 (explore first)
            if (closerSide) {
                pq.push({ node: closerSide, dist: 0 });
            }

            // Add farther side with priority = |dist| (explore later if needed)
            if (fartherSide && candidates.size < maxCandidates) {
                pq.push({ node: fartherSide, dist: Math.abs(dist) });
            }
        }
    }

    /**
     * @param {number} i
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_by_index(i, k = 5) {
        if (i < 0 || i >= this._elements.length) return [];
        return this.search(this._elements[i], k);
    }

    /**
     * Alias for search_by_index for backward compatibility.
     *
     * @param {number} i - Index of the query element
     * @param {number} [k=5] - Number of nearest neighbors to return
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_index(i, k = 5) {
        return this.search_by_index(i, k);
    }
}
