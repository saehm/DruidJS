import { DisjointSet } from "../datastructure/index.js";

/**
 * @typedef {[number, number, number]} WeightedEdge An edge `[u, v, weight]` between the vertices with
 *   index `u` and `v`.
 */

/**
 * Computes a minimum spanning tree of a weighted graph with Kruskal's algorithm.
 *
 * The graph is given as an edge list over vertices `0 … N-1`, so a caller holding a sparse structure
 * (a k-nearest-neighbor graph, say) never has to materialize the `N ⨯ N` matrix. Edges are treated as
 * undirected; passing both `[u, v, w]` and `[v, u, w]` is harmless, the second one is skipped as a
 * cycle.
 *
 * If the graph is disconnected the result is a minimum spanning *forest* — one tree per connected
 * component, and fewer than `N - 1` edges. Callers that need a connected result must guarantee a
 * connected input.
 *
 * @param {WeightedEdge[]} edges - The edges of the graph. Not mutated.
 * @param {number} N - The number of vertices. Vertex indices must be in `[0, N)`.
 * @returns {WeightedEdge[]} The edges of the minimum spanning tree, ascending by weight.
 * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
 */
export function minimum_spanning_tree(edges, N) {
    /** @type {WeightedEdge[]} */
    const F = [];
    if (N <= 1) return F;

    /** @type {DisjointSet<number>} */
    const disjoint_set = new DisjointSet(Array.from({ length: N }, (_, i) => i));
    // Sorting a copy keeps the caller's array intact, which matters because the information graph in
    // MINFOTree is reused after the tree is extracted.
    const sorted = edges.slice().sort((a, b) => a[2] - b[2]);

    for (const edge of sorted) {
        const set_u = disjoint_set.find(edge[0]);
        const set_v = disjoint_set.find(edge[1]);
        if (set_u === null || set_v === null) throw new Error("Edge refers to a vertex outside [0, N)!");
        if (set_u !== set_v) {
            F.push(edge);
            disjoint_set.union(set_u, set_v);
            if (F.length === N - 1) break;
        }
    }

    return F;
}
