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
 * This is an *approximate* index. Recall above 95% is expected for any `m`; raise `ef` or
 * `ef_construction` if a dataset needs more. Queries overtake the exact {@link KDTree} at around
 * 10 000 points, so below a few thousand prefer {@link spatial_tree}, which is exact and quicker.
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
     * @param {Partial<ParametersHNSW>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(points, parameters = {}) {
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

        // `KNN` keeps the parameter object as handed to it, so defaults are merged here rather
        // than left to a default argument, which would only apply when `parameters` is omitted
        // entirely and drop every unspecified default for `new HNSW(points, { m: 8 })`.
        const merged = Object.assign(
            {
                metric: euclidean,
                heuristic: true,
                m: 16,
                ef_construction: 200,
                m0: null,
                mL: null,
                seed: 1212,
                ef: 50,
            },
            parameters,
        );

        super([firstElement], merged);

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

        /** @private @type {number} */
        this._search_id = 1;

        /** @private @type {Uint32Array} */
        this._visited_stamps = new Uint32Array(1024);

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

            // Phase 1: greedy descent through the layers the element does not join, narrowing the
            // entry point one layer at a time.
            for (let l_c = L; l_c > l; --l_c) {
                const search_result = this._search_layer(element, ep_indices, 1, l_c);
                if (search_result.length > 0) {
                    ep_indices = [search_result[0].index];
                }
            }

            // Phase 2: join every layer from min(L, l) down to 0.
            for (let l_c = Math.min(L, l); l_c >= 0; --l_c) {
                const layer = graph.get(l_c);
                if (!layer) continue;

                layer.point_indices.push(global_index);
                if (!layer.edges.has(global_index)) {
                    layer.edges.set(global_index, []);
                }

                // Candidates are always members of this layer. Deliberately no fallback to a scan
                // over the whole dataset: that would link the element to points not on this layer.
                const W = this._search_layer(element, ep_indices, ef_construction, l_c).filter(
                    (c) => c.index !== global_index,
                );
                if (W.length === 0) continue;

                // The element gets `m` links of its own on every layer, layer 0 included. `m0` is
                // not a second budget for it, but the cap enforced once reverse links accumulate.
                const neighbor_indices = this._select(element, W, m, l_c);
                layer.edges.set(global_index, neighbor_indices.slice());

                const max_conn = l_c === 0 ? m0 : m;
                for (const neighbor_idx of neighbor_indices) {
                    if (neighbor_idx === global_index) continue;

                    // Add the reverse connection from neighbor to element.
                    let neighbor_edges = layer.edges.get(neighbor_idx);
                    if (!neighbor_edges) {
                        neighbor_edges = [];
                        layer.edges.set(neighbor_idx, neighbor_edges);
                    }
                    if (!neighbor_edges.includes(global_index)) {
                        neighbor_edges.push(global_index);
                    }

                    // Prune with the same selector used to build the graph; plain nearest-M would
                    // drop the long-range links that keep the layer navigable.
                    if (neighbor_edges.length > max_conn) {
                        const neighbor_element = elements[neighbor_idx];
                        const neighbor_candidates = neighbor_edges
                            .filter((idx) => idx !== neighbor_idx)
                            .map((idx) => ({
                                element: elements[idx],
                                index: idx,
                                distance: this._metric(neighbor_element, elements[idx]),
                            }));
                        layer.edges.set(
                            neighbor_idx,
                            this._select(neighbor_element, neighbor_candidates, max_conn, l_c),
                        );
                    }
                }

                // Descend from the closest neighbor actually linked on this layer.
                ep_indices = [neighbor_indices.length > 0 ? neighbor_indices[0] : W[0].index];
            }

            // Layers above the previous maximum: the element is their only member and becomes the
            // new entry point. Also covers the first insertion, where `L` is -1.
            if (l > L) {
                for (let i = L + 1; i <= l; ++i) {
                    graph.set(i, {
                        l_c: i,
                        point_indices: [global_index],
                        edges: new Map([[global_index, []]]),
                    });
                }
                this._ep = [global_index];
                this._L = l;
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
     * @param {boolean} [extend_candidates=false] - Whether to extend candidates with their neighbors
     * @param {boolean} [keep_pruned_connections=false] - Whether to add pruned connections back if needed
     * @returns {number[]} Selected neighbor indices, nearest first
     */
    _select_heuristic(q, candidates, M, l_c, extend_candidates = false, keep_pruned_connections = false) {
        const metric = this._metric;
        const elements = this._elements;

        const dist_map = new Map();
        const pool = [];
        for (let i = 0; i < candidates.length; ++i) {
            const c = candidates[i];
            if (dist_map.has(c.index)) continue;
            dist_map.set(c.index, c.distance);
            pool.push(c.index);
        }

        if (extend_candidates) {
            const layer = this._graph.get(l_c);
            if (layer) {
                for (let i = 0; i < candidates.length; ++i) {
                    const edges = layer.edges.get(candidates[i].index) ?? [];
                    for (let j = 0; j < edges.length; ++j) {
                        const idx = edges[j];
                        if (dist_map.has(idx)) continue;
                        const elem = elements[idx];
                        if (!elem || elem.length !== q.length) continue;
                        dist_map.set(idx, metric(elem, q));
                        pool.push(idx);
                    }
                }
            }
        }

        pool.sort((a, b) => dist_map.get(a) - dist_map.get(b));

        // Nothing to choose between: keep everything rather than thinning an already small
        // neighborhood.
        if (pool.length < M) return pool;

        /** @type {number[]} */
        const R = [];
        /** @type {number[]} */
        const discarded = [];

        // Keep `e` only if the query is closer to it than any already-selected neighbor is. This
        // spreads links across directions rather than letting them collapse onto the nearest
        // cluster, which is what makes the layer navigable and not merely a k-nearest-neighbor graph.
        for (let i = 0; i < pool.length; ++i) {
            if (R.length >= M) break;
            const idx = pool[i];
            const dist_to_q = dist_map.get(idx);
            const e = elements[idx];

            let should_add = true;
            for (let j = 0; j < R.length; ++j) {
                if (metric(e, elements[R[j]]) < dist_to_q) {
                    should_add = false;
                    break;
                }
            }

            if (should_add) {
                R.push(idx);
            } else {
                discarded.push(idx);
            }
        }

        // Off by default: refilling to `M` with the nearest rejected candidates undoes the spread
        // above and turns the result back into plain nearest-M.
        if (keep_pruned_connections) {
            for (let i = 0; i < discarded.length && R.length < M; ++i) {
                R.push(discarded[i]);
            }
        }

        return R;
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

        // A layer whose members carry no edges yet is not a failure: the entry points are still the
        // best answer it can give.
        if (!layer || !ep_indices || ep_indices.length === 0) {
            return [];
        }

        const numElements = elements.length;
        if (this._visited_stamps.length < numElements) {
            const nextCap = Math.max(numElements * 2, 1024);
            const old = this._visited_stamps;
            this._visited_stamps = new Uint32Array(nextCap);
            this._visited_stamps.set(old);
        }

        const search_id = ++this._search_id;
        const stamps = this._visited_stamps;

        const init_candidates = [];
        for (let i = 0; i < ep_indices.length; ++i) {
            const idx = ep_indices[i];
            const elem = elements[idx];
            if (elem !== undefined) {
                stamps[idx] = search_id;
                init_candidates.push({
                    element: elem,
                    index: idx,
                    distance: metric(elem, q),
                });
            }
        }

        if (init_candidates.length === 0) {
            return [];
        }

        const C = new Heap(init_candidates, (item) => item.distance, "min");
        const W = new Heap(init_candidates, (item) => item.distance, "max");

        while (!C.empty) {
            const c = C.pop();
            if (!c) break;
            const furthest_dist = W.first?.value ?? Infinity;

            if (c.value > furthest_dist) {
                break;
            }

            const edges = layer.edges.get(c.element.index) ?? [];

            for (let i = 0; i < edges.length; ++i) {
                const neighbor_idx = edges[i];
                if (stamps[neighbor_idx] !== search_id) {
                    stamps[neighbor_idx] = search_id;
                    const neighbor_element = elements[neighbor_idx];
                    if (!neighbor_element || neighbor_element.length !== q.length) continue;
                    if (neighbor_idx === c.element.index) continue;

                    const dist_e = metric(neighbor_element, q);
                    const current_furthest = W.first?.value ?? Infinity;

                    if (dist_e < current_furthest || W.length < ef) {
                        const candidateNode = {
                            element: neighbor_element,
                            index: neighbor_idx,
                            distance: dist_e,
                        };
                        C.push(candidateNode);
                        W.push(candidateNode);

                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }

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
}
