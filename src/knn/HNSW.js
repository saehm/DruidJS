import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { KNN } from "./KNN.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersHNSW } from "./index.js" */

/**
 * @typedef {Object} Layer
 * @property {number} l_c - Layer number
 * @property {number[]} point_indices - Global indices of points in this layer
 * @property {Map<number, number[]>} edges - Global index -> array of connected global indices
 */

/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} Candidate
 * @property {T} element - The actual data point
 * @property {number} index - Global index in the dataset
 * @property {number} distance - Distance from query
 */

/**
 * Hierarchical Navigable Small World (HNSW) graph for approximate nearest neighbor search.
 *
 * HNSW builds a multi-layer graph structure where each layer is a navigable small world graph.
 * The top layers serve as "highways" for fast traversal, while lower layers provide accuracy.
 * Each element is assigned to a random level, allowing logarithmic search complexity.
 *
 * Key parameters:
 * - `m`: Controls the number of connections per element (affects accuracy/memory)
 * - `ef_construction`: Controls the quality of the graph during construction (higher = better but slower)
 * - `ef`: Controls the quality of search (higher = better recall but slower)
 *
 * Based on:
 * - "Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs"
 *   by Malkov & Yashunin (2016)
 * - "Approximate Nearest Neighbor Search on High Dimensional Data"
 *   by Li et al. (2019)
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersHNSW>
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const points = [[1, 2], [3, 4], [5, 6], [7, 8]];
 * const hnsw = new druid.HNSW(points, {
 *     metric: druid.euclidean,
 *     m: 16,
 *     ef_construction: 200
 * });
 *
 * const query = [2, 3];
 * const neighbors = hnsw.search(query, 2);
 * // [{ element: [1, 2], index: 0, distance: 1.41 }, ...]
 */
export class HNSW extends KNN {
    /**
     * Creates a new HNSW index.
     *
     * @param {T[]} points - Initial points to add to the index
     * @param {ParametersHNSW} [parameters={}] - Configuration parameters
     */
    constructor(
        points,
        parameters = {
            metric: euclidean,
            heuristic: true,
            m: 16,
            ef_construction: 200,
            m0: null,
            mL: null,
            seed: 1212,
            ef: 50,
        },
    ) {
        // Handle empty initialization - use dummy element
        const hasElements = points && points.length > 0;
        let firstElement = /** @type {T} */ (hasElements ? points[0] : new Float64Array([0]));

        // Validate all points have consistent dimensions
        if (hasElements) {
            const expected_dim = firstElement.length;
            for (let i = 1; i < points.length; i++) {
                if (!points[i] || points[i].length !== expected_dim) {
                    console.warn(
                        `HNSW: Point ${i} has inconsistent dimensions (expected ${expected_dim}, got ${points[i]?.length})`,
                    );
                    // Remove invalid points
                    points = points.filter((_, idx) => idx === 0 || points[idx]?.length === expected_dim);
                    firstElement = points[0];
                }
            }
        }

        super([firstElement], parameters);

        // Store reference to elements before clearing
        const elementsToAdd = hasElements ? [...points] : [];
        /** @type {T[]} */
        this._elements = [];

        /** @type {Metric} */
        this._metric = this._parameters.metric || euclidean;

        /** @type {Function} */
        this._select = this._parameters.heuristic ? this._select_heuristic.bind(this) : this._select_simple.bind(this);

        /**
         * @private
         * @type {Map<number, Layer>}
         */
        this._graph = new Map();

        /** @type {number} */
        this._next_index = 0;

        // Validate and set parameters
        const m_param = this._parameters.m ?? 16;
        if (m_param <= 0 || !Number.isInteger(m_param)) {
            throw new Error("HNSW: parameter 'm' must be a positive integer");
        }
        /** @type {number} */
        this._m = Math.max(2, m_param);

        const ef_construction_param = this._parameters.ef_construction ?? 200;
        if (ef_construction_param <= 0 || !Number.isInteger(ef_construction_param)) {
            throw new Error("HNSW: parameter 'ef_construction' must be a positive integer");
        }
        /** @type {number} */
        this._ef_construction = ef_construction_param;

        const ef_param = this._parameters.ef ?? 50;
        if (ef_param <= 0 || !Number.isInteger(ef_param)) {
            throw new Error("HNSW: parameter 'ef' must be a positive integer");
        }
        /** @type {number} */
        this._ef = ef_param;

        const m0_param = this._parameters.m0 ?? 2 * this._m;
        if (m0_param <= 0 || !Number.isInteger(m0_param)) {
            throw new Error("HNSW: parameter 'm0' must be a positive integer");
        }
        /** @type {number} */
        this._m0 = m0_param;

        /** @type {number} */
        this._mL = this._parameters.mL ?? 1 / Math.log(this._m);

        /** @type {Randomizer} */
        this._randomizer = new Randomizer(this._parameters.seed);

        /** @type {number} - Current maximum layer in the graph */
        this._L = -1;

        /** @type {number[] | null} - Entry point indices for search */
        this._ep = null;

        // Add initial points
        if (elementsToAdd && elementsToAdd.length > 0) {
            this.add(elementsToAdd);
        }
    }

    /**
     * Add a single element to the index.
     *
     * @param {T} element - Element to add
     * @returns {HNSW<T>} This instance for chaining
     */
    addOne(element) {
        return this.add([element]);
    }

    /**
     * Add multiple elements to the index.
     *
     * @param {T[]} new_elements - Elements to add
     * @returns {HNSW<T>} This instance for chaining
     */
    add(new_elements) {
        // Handle empty array
        if (!new_elements || new_elements.length === 0) {
            return this;
        }

        const m = this._m;
        const ef_construction = this._ef_construction;
        const m0 = this._m0;
        const mL = this._mL;
        const randomizer = this._randomizer;
        const graph = this._graph;

        // Ensure _elements is a proper array that supports push
        if (!Array.isArray(this._elements)) {
            this._elements = Array.from(this._elements);
        }
        const elements = this._elements;

        // Get expected dimension from first existing element or first new element
        const expected_dim = elements.length > 0 ? elements[0].length : new_elements[0]?.length;

        for (const element of new_elements) {
            // Validate element
            if (!element || (!Array.isArray(element) && !(element instanceof Float64Array))) {
                console.warn("HNSW: Skipping invalid element (null, undefined, or not an array)");
                continue;
            }

            // Validate dimensions
            if (element.length !== expected_dim) {
                console.warn(
                    `HNSW: Skipping element with wrong dimensions (expected ${expected_dim}, got ${element.length})`,
                );
                continue;
            }

            elements.push(element);
            const global_index = elements.length - 1;

            // Assign random level to the element
            // Level is drawn from exponential distribution: l = floor(-ln(uniform(0,1)) * mL)
            const rand = Math.max(randomizer.random, 1e-10); // Avoid log(0)
            const l = Math.min(31, Math.floor(-Math.log(rand) * mL));

            let ep_indices = this._ep ? [...this._ep] : null;
            const L = this._L;

            if (L >= 0) {
                // Search from top layer down to min(L, l) + 1
                // These are the layers where element will NOT be inserted
                for (let l_c = L; l_c > l; --l_c) {
                    const search_result = this._search_layer(element, ep_indices, 1, l_c);
                    if (search_result.length > 0) {
                        ep_indices = [search_result[0].index];
                    }
                }

                // Insert element into layers l down to 0
                for (let l_c = Math.min(L, l); l_c >= 0; --l_c) {
                    const layer = graph.get(l_c);
                    if (!layer) continue;

                    layer.point_indices.push(global_index);

                    // Search for ef_construction nearest neighbors
                    let W = this._search_layer(element, ep_indices, ef_construction, l_c);

                    // If graph search returns no results (e.g., graph is empty or disconnected),
                    // fall back to linear search over all existing elements
                    if (W.length === 0 && elements.length > 1) {
                        const fallbackCandidates = [];
                        for (let i = 0; i < elements.length - 1; i++) {
                            const elem = elements[i];
                            if (elem && elem.length === element.length) {
                                fallbackCandidates.push({
                                    element: elem,
                                    index: i,
                                    distance: this._metric(element, elem),
                                });
                            }
                        }
                        fallbackCandidates.sort((a, b) => a.distance - b.distance);
                        W = fallbackCandidates.slice(0, ef_construction);
                        // Update ep_indices for next layer based on fallback results
                        if (l_c === Math.min(L, l)) {
                            ep_indices = W.map((c) => c.index);
                        }
                    }

                    // Select neighbors using heuristic or simple approach (respect heuristic setting on all layers)
                    const neighbor_indices = this._select(element, W, l_c === 0 ? m0 : m, l_c);

                    // Add bidirectional connections
                    for (const neighbor_idx of neighbor_indices) {
                        if (neighbor_idx === global_index) continue;

                        // Add connection from element to neighbor
                        if (!layer.edges.has(global_index)) {
                            layer.edges.set(global_index, []);
                        }
                        layer.edges.get(global_index)?.push(neighbor_idx);

                        // Add connection from neighbor to element
                        if (!layer.edges.has(neighbor_idx)) {
                            layer.edges.set(neighbor_idx, []);
                        }
                        const neighbor_edge_list = layer.edges.get(neighbor_idx);
                        if (neighbor_edge_list && !neighbor_edge_list.includes(global_index)) {
                            neighbor_edge_list.push(global_index);
                        }

                        // Prune connections if too many
                        const max_conn = l_c === 0 ? m0 : m;
                        const neighbor_edges = layer.edges.get(neighbor_idx);
                        if (neighbor_edges && neighbor_edges.length > max_conn) {
                            const neighbor_element = elements[neighbor_idx];
                            // Filter out self-connections before pruning
                            const valid_neighbor_edges = neighbor_edges.filter((idx) => idx !== neighbor_idx);
                            const neighbor_candidates = valid_neighbor_edges.map((idx) => ({
                                element: elements[idx],
                                index: idx,
                                distance: this._metric(neighbor_element, elements[idx]),
                            }));
                            const pruned =
                                l_c === 0
                                    ? this._select_simple(neighbor_element, neighbor_candidates, max_conn)
                                    : this._select(neighbor_element, neighbor_candidates, max_conn, l_c);
                            layer.edges.set(neighbor_idx, pruned);
                        }
                    }

                    // Use closest neighbor as entry point for next layer (following HNSW paper)
                    if (W.length > 0) {
                        ep_indices = [W[0].index];
                    }
                }
            }

            // If element's level is higher than current max, create new layers
            if (l > L) {
                for (let i = L + 1; i <= l; ++i) {
                    graph.set(i, {
                        l_c: i,
                        point_indices: [global_index],
                        edges: new Map(),
                    });
                }
                // Element becomes the new entry point
                this._ep = [global_index];
                this._L = l;
            }

            // Special case: if this is the first element (L was -1),
            // we need to ensure layer 0 has proper structure for future insertions
            if (L === -1) {
                if (!graph.has(0)) {
                    graph.set(0, {
                        l_c: 0,
                        point_indices: [global_index],
                        edges: new Map(),
                    });
                }
                const layer0 = graph.get(0);
                if (layer0 && !layer0.edges.has(global_index)) {
                    layer0.edges.set(global_index, []);
                }
            }
        }

        return this;
    }

    /**
     * Select neighbors using the heuristic approach.
     *
     * The heuristic extends candidates with their neighbors and selects
     * points that are closer to the query than to already selected points.
     * This maintains graph connectivity better than simple selection.
     *
     * @private
     * @param {T} q - Query element
     * @param {Candidate<T>[]} candidates - Candidate elements with distances
     * @param {number} M - Maximum number of neighbors to return
     * @param {number} l_c - Layer number
     * @param {boolean} [extend_candidates=true] - Whether to extend candidates with their neighbors
     * @param {boolean} [keep_pruned_connections=true] - Whether to add pruned connections back if needed
     * @returns {number[]} Selected neighbor indices
     */
    _select_heuristic(q, candidates, M, l_c, extend_candidates = true, keep_pruned_connections = true) {
        if (l_c > this._L) {
            return candidates.map((c) => c.index);
        }

        const metric = this._metric;
        const layer = this._graph.get(l_c);
        const elements = this._elements;

        // Extend candidate set with neighbors of candidates
        const W_set = new Set(candidates.map((c) => c.index));
        if (extend_candidates) {
            for (const c of candidates) {
                const edges = layer?.edges.get(c.index);
                if (edges) {
                    for (const neighbor_idx of edges) {
                        W_set.add(neighbor_idx);
                    }
                }
            }
        }

        // Create extended candidates with distances
        const W = [...W_set]
            .map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            }))
            .sort((a, b) => a.distance - b.distance);

        const R = [];
        const W_discarded = [];

        // Select neighbors: prefer points closer to query than to already selected points
        for (const e of W) {
            if (R.length >= M) break;

            let should_add = true;

            // Check if e is closer to query than to any already selected point
            for (const r of R) {
                const dist_er = metric(e.element, r.element);
                if (dist_er < e.distance) {
                    should_add = false;
                    break;
                }
            }

            if (should_add) {
                R.push(e);
            } else {
                W_discarded.push(e);
            }
        }

        // Add discarded connections if we need more
        if (keep_pruned_connections && R.length < M) {
            for (const e of W_discarded) {
                if (R.length >= M) break;
                R.push(e);
            }
        }

        return R.map((c) => c.index);
    }

    /**
     * Select neighbors using simple distance-based selection.
     *
     * Simply returns the M closest candidates to the query.
     *
     * @private
     * @param {T} q - Query element
     * @param {Candidate<T>[]} C - Candidate elements with distances
     * @param {number} M - Maximum number of neighbors to return
     * @returns {number[]} M nearest candidate indices
     */
    _select_simple(q, C, M) {
        if (C.length <= M) return C.map((c) => c.index);

        // Candidates already have distance computed, use it directly
        return C.slice()
            .sort((a, b) => a.distance - b.distance)
            .slice(0, M)
            .map((c) => c.index);
    }

    /**
     * Search a single layer for nearest neighbors.
     *
     * Implements the greedy search algorithm: start from entry points,
     * always expand the closest unvisited candidate, maintain a list
     * of the ef closest found neighbors.
     *
     * @private
     * @param {T} q - Query element
     * @param {number[] | null} ep_indices - Entry point indices
     * @param {number} ef - Number of nearest neighbors to find
     * @param {number} l_c - Layer number to search
     * @returns {Candidate<T>[]} ef nearest neighbors found with their distances
     */
    _search_layer(q, ep_indices, ef, l_c) {
        const metric = this._metric;
        const layer = this._graph.get(l_c);
        const elements = this._elements;

        if (!layer || layer.edges.size === 0 || !ep_indices || ep_indices.length === 0) {
            return [];
        }

        // Filter out invalid indices
        const valid_ep_indices = ep_indices.filter((idx) => elements[idx] !== undefined);
        if (valid_ep_indices.length === 0) {
            return [];
        }

        // Visited set to avoid cycles
        const visited = new Set(valid_ep_indices);

        // Candidate set (min-heap): closest unvisited candidates to expand
        const C = new Heap(
            valid_ep_indices.map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            })),
            (item) => item.distance,
            "min",
        );

        // Result set (max-heap): ef closest found neighbors
        const W = new Heap(
            valid_ep_indices.map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            })),
            (item) => item.distance,
            "max",
        );

        // Algorithm 2 stops when the distance from query to the next candidate is greater
        // than the distance to the furthest element in the result set W.
        while (!C.empty) {
            const c = C.pop();
            if (!c) break;
            const furthest_dist = W.first?.value ?? Infinity;

            // Stop if current candidate is farther than furthest result
            if (c.value > furthest_dist) {
                break;
            }

            const edges = layer.edges.get(c.element.index);
            if (!edges) continue;

            for (const neighbor_idx of edges) {
                if (!visited.has(neighbor_idx)) {
                    const neighbor_element = elements[neighbor_idx];
                    // Skip invalid elements or elements with different dimensions
                    if (!neighbor_element || neighbor_element.length !== q.length) continue;

                    // Skip self-connections
                    if (neighbor_idx === c.element.index) continue;

                    visited.add(neighbor_idx);
                    const dist_e = metric(neighbor_element, q);

                    const current_furthest = W.first?.value ?? Infinity;
                    if (dist_e < current_furthest || W.length < ef) {
                        C.push({
                            element: neighbor_element,
                            index: neighbor_idx,
                            distance: dist_e,
                        });
                        W.push({
                            element: neighbor_element,
                            index: neighbor_idx,
                            distance: dist_e,
                        });

                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }

        // Return sorted results for consistent entry point selection
        return W.data().sort((a, b) => a.distance - b.distance);
    }

    /**
     * Searches for the K nearest neighbors to a query element in the HNSW graph.
     *
     * Performs a multi-layer search starting from the entry point and traversing
     * each layer as entry points for the next.
     *
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @returns {Candidate<T>[]} K nearest neighbors with their distances
     */
    search(q, K) {
        // Validate K
        if (!Number.isInteger(K) || K <= 0) {
            throw new Error("HNSW: parameter 'K' must be a positive integer");
        }

        // Validate query dimensions
        if (!q || (!Array.isArray(q) && !(q instanceof Float64Array))) {
            throw new Error("HNSW: query must be an array");
        }

        const search_ef = this._ef;

        // Fallback to linear search if graph is not properly initialized
        if (this._L < 0 || !this._ep || this._elements.length === 0) {
            return this._linear_search(q, K);
        }

        let ep_indices = [...this._ep];

        // Search from top layer down to layer 1
        for (let l_c = this._L; l_c > 0; --l_c) {
            const result = this._search_layer(q, ep_indices, 1, l_c);
            if (result.length > 0) {
                ep_indices = [result[0].index];
            }
        }

        // Search layer 0 with ef candidates
        const result = this._search_layer(q, ep_indices, Math.max(search_ef, K), 0);

        // If graph search returns no results, fallback to linear search
        if (result.length === 0) {
            return this._linear_search(q, K);
        }

        // Return K closest
        return result.slice(0, K);
    }

    /**
     * Fallback linear search when graph search fails
     * @private
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @returns {Candidate<T>[]}
     */
    _linear_search(q, K) {
        const metric = this._metric;
        const elements = this._elements;
        const N = elements.length;

        if (N === 0) return [];

        /** @type {Candidate<T>[]} */
        const candidates = [];
        for (let i = 0; i < N; i++) {
            const element = elements[i];
            // Skip elements with different dimensions (can happen with inconsistent data)
            if (!element || element.length !== q.length) continue;

            candidates.push({
                element: element,
                index: i,
                distance: metric(q, element),
            });
        }

        candidates.sort((a, b) => a.distance - b.distance);
        return candidates.slice(0, K);
    }

    /**
     * Iterator for searching the HNSW graph layer by layer.
     *
     * Yields intermediate results at each layer for debugging or visualization.
     *
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @param {number?} [ef] - Size of dynamic candidate list
     * @yields {{layer: number, candidates: Candidate[]}}
     */
    *search_iter(q, K, ef = null) {
        const search_ef = ef ?? this._ef;

        if (this._L < 0 || !this._ep) {
            return;
        }

        let ep_indices = [...this._ep];

        // Yield entry points at top layer instead of query itself
        const top_layer = this._graph.get(this._L);
        if (top_layer && this._ep && this._ep.length > 0) {
            const entry_candidates = this._ep
                .filter((idx) => this._elements[idx] !== undefined)
                .map((idx) => ({
                    element: this._elements[idx],
                    index: idx,
                    distance: this._metric(this._elements[idx], q),
                }));
            yield {
                layer: this._L,
                candidates: entry_candidates,
            };
        }

        for (let l_c = this._L; l_c > 0; --l_c) {
            const result = this._search_layer(q, ep_indices, 1, l_c);
            yield { layer: l_c, candidates: result };
            // Use closest candidate as entry point for next layer (following HNSW paper)
            ep_indices = result.length > 0 ? [result[0].index] : ep_indices;
        }

        const result = this._search_layer(q, ep_indices, Math.max(search_ef, K), 0);
        yield { layer: 0, candidates: result };
    }

    /**
     * Get the number of elements in the index.
     *
     * @returns {number} Number of elements
     */
    get size() {
        return this._elements?.length ?? 0;
    }

    /**
     * Get the number of layers in the graph.
     *
     * @returns {number} Number of layers
     */
    get num_layers() {
        return this._L + 1;
    }

    /**
     * Get an element by its index.
     *
     * @param {number} index - Element index
     * @returns {T} The element at the given index
     */
    get_element(index) {
        return this._elements[index];
    }

    /**
     * Search for nearest neighbors using an element index as the query.
     *
     * @param {number} i - Index of the query element
     * @param {number} [K=5] - Number of nearest neighbors to return
     * @returns {Candidate<T>[]} K nearest neighbors
     */
    search_by_index(i, K = 5) {
        const elements = this._elements;
        if (i < 0 || i >= elements.length) return [];

        const element = elements[i];
        if (!element) return [];

        return this.search(element, K);
    }
}
