import { HierarchicalClustering, KMeans } from "../clustering/index.js";
import { DisjointSet } from "../datastructure/index.js";
import { spatial_tree } from "../knn/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { minimum_spanning_tree } from "../util/index.js";
import { DR } from "./DR.js";
import { KKMDS } from "./KKMDS.js";

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/** @import {ParametersMINFOTree} from "./index.js" */
/** @import {WeightedEdge} from "../util/minimum_spanning_tree.js" */

/** The golden ratio conjugate, the paper's shrinkage factor for intra-cluster edges. */
const ALPHA_GOLDEN = (Math.sqrt(5) - 1) / 2;

/**
 * Minimum Information Trees (MINFO Trees)
 *
 * A visualization method for clustered high-dimensional data. It reads the cluster labels on a
 * k-nearest-neighbor graph as a q-state Potts model, weights the edges by an information-geometric
 * curvature derived from it, and keeps the minimum spanning tree of that weighted graph.
 *
 * {@link projection} returns 2D coordinates like any other method here, but the tree is the actual
 * output — read {@link edges} to draw it.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersMINFOTree>
 * @category Dimensionality Reduction
 * @see {@link https://doi.org/10.1109/ACCESS.2025.3602730}
 * @see {@link TopoMap} for another spanning-tree-based projection.
 */
export class MINFOTree extends DR {
    /**
     * Minimum Information Trees.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersMINFOTree>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1109/ACCESS.2025.3602730}
     */
    constructor(X, parameters = {}) {
        /** @type {ParametersMINFOTree} */
        const defaults = {
            k: -Infinity,
            d: 2,
            metric: euclidean,
            clusters: -Infinity,
            clustering: "hierarchical",
            labels: null,
            alpha: ALPHA_GOLDEN,
            epsilon: 1e-3,
            layout: "kamada_kawai",
            iterations: 300,
            seed: 1212,
            eig_args: {},
        };
        super(X, defaults, parameters);

        if (this._parameters.k === -Infinity) {
            // k = ln(n), as in the paper's experiments and the author's reference implementation.
            this.parameter("k", Math.min(Math.max(2, Math.round(Math.log(this._N))), this._N - 1));
        }

        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * The edges of the Minimum Information Tree, as `[u, v, weight]` over row indices of `X`,
     * ascending by weight. This is the method's real output — {@link projection} is one drawing of
     * it.
     *
     * @returns {WeightedEdge[]}
     */
    get edges() {
        this.check_init();
        return /** @type {WeightedEdge[]} */ (this._edges);
    }

    /**
     * The cluster label per point, whether supplied or computed in step 1.
     *
     * @returns {Int32Array} Labels, remapped to `0 … q-1`.
     */
    get labels() {
        this.check_init();
        return /** @type {Int32Array} */ (this._labels);
    }

    /**
     * The information curvature `S_i` per point, normalised as the edge weighting uses it.
     *
     * @returns {Float64Array}
     */
    get curvature() {
        this.check_init();
        return /** @type {Float64Array} */ (this._curvature);
    }

    /**
     * The maximum pseudo-likelihood estimate of the Potts inverse temperature.
     *
     * @returns {number}
     */
    get beta() {
        this.check_init();
        return /** @type {number} */ (this._beta);
    }

    /**
     * Resolves the labels field: either the caller's, or the outcome of clustering `X`.
     *
     * Labels are remapped to a dense `0 … q-1` so they can index the Potts state vectors directly.
     *
     * @private
     * @returns {Int32Array}
     */
    _make_labels() {
        const N = this._N;
        const supplied = this.parameter("labels");
        /** @type {ArrayLike<any>} */
        let raw;

        if (supplied != null) {
            raw = /** @type {ArrayLike<any>} */ (supplied);
            if (raw.length !== N) {
                throw new Error(`labels has length ${raw.length}, but X has ${N} rows!`);
            }
        } else {
            const clusters = /** @type {number} */ (this.parameter("clusters"));
            if (!(clusters >= 2)) {
                // No default is defensible here. The paper's own limitations section stresses that
                // the tree inherits whatever the clustering got wrong, so silently inventing a
                // cluster count would bury the single most consequential choice in the method.
                throw new Error(
                    "MINFOTree needs a labels field: pass `clusters` (>= 2) to cluster X, or `labels` to supply your own.",
                );
            }
            const clustering = /** @type {"hierarchical" | "kmeans"} */ (this.parameter("clustering"));
            const metric = /** @type {Metric} */ (this.parameter("metric"));
            if (clustering === "kmeans") {
                const km = new KMeans(this.X, { K: clusters, metric, seed: this.parameter("seed") });
                raw = km.get_cluster_list();
            } else {
                // Ward is what the paper's experiments use.
                const hc = new HierarchicalClustering(this.X, { linkage: "ward", metric });
                raw = hc.get_cluster_list(this._dendrogram_cut(hc, clusters));
            }
        }

        /** @type {Map<any, number>} */
        const lookup = new Map();
        const labels = new Int32Array(N);
        for (let i = 0; i < N; ++i) {
            const key = raw[i];
            let id = lookup.get(key);
            if (id === undefined) {
                id = lookup.size;
                lookup.set(key, id);
            }
            labels[i] = id;
        }
        this._q = lookup.size;
        return labels;
    }

    /**
     * Finds the distance threshold that cuts a dendrogram into exactly `k` clusters.
     *
     * {@link HierarchicalClustering.get_cluster_list} cuts at a height, not at a cluster count, so
     * the height has to be recovered: breaking the `k - 1` tallest merges leaves `k` components, and
     * any threshold between the `k`-th and `(k-1)`-th tallest does that. The midpoint is used so
     * ties on either side cannot land the cut on the wrong number.
     *
     * Only `2 … N-1` is reachable this way. A dendrogram has `N-1` merges, and the traversal emits a
     * cluster per subtree it stops at, never a bare leaf — so there is no height that separates all
     * `N` points. Asking for more is rejected rather than approximated: cutting below every merge
     * silently yields *one* cluster, the opposite of what was asked for.
     *
     * @private
     * @param {HierarchicalClustering} hc
     * @param {number} k - At least 2; {@link _make_labels} rejects anything smaller before this runs.
     * @returns {number}
     */
    _dendrogram_cut(hc, k) {
        if (!hc.root) throw new Error("Hierarchical clustering produced no dendrogram!");
        const heights = hc.root
            .descendants()
            .filter((node) => !node.isLeaf)
            .map((node) => node.dist)
            .sort((a, b) => a - b);
        const M = heights.length;
        if (k > M) {
            throw new Error(
                `clusters must be at most ${M} for ${this._N} points with hierarchical clustering; got ${k}. ` +
                    `Use \`clustering: "kmeans"\` or pass \`labels\` for finer partitions.`,
            );
        }
        return (heights[M - k] + heights[M - k + 1]) / 2;
    }

    /**
     * Builds the symmetric k-nearest-neighbor graph, connected.
     *
     * The Potts model is defined over a neighborhood system and the layout needs finite pairwise
     * distances, so a k-NNG that falls into several components would break both. When that happens
     * the components are stitched together with the shortest edge between each pair, which is what
     * the author's reference implementation achieves by unioning in a spanning tree of the complete
     * graph.
     *
     * @private
     * @returns {number[][]} Adjacency list.
     */
    _make_knn_graph() {
        const N = this._N;
        const k = /** @type {number} */ (this.parameter("k"));
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const tree = spatial_tree(this.X.to2dArray(), {
            metric,
            seed: /** @type {number} */ (this.parameter("seed")),
        });

        /** @type {Set<number>[]} */
        const adjacency = Array.from({ length: N }, () => new Set());
        for (let i = 0; i < N; ++i) {
            for (const neighbor of tree.search_by_index(i, k)) {
                const j = neighbor.index;
                if (j === i) continue;
                adjacency[i].add(j);
                adjacency[j].add(i);
            }
        }

        this._connect_components(adjacency);
        return adjacency.map((set) => [...set]);
    }

    /**
     * Adds the fewest edges needed to make `adjacency` connected, preferring short ones.
     *
     * @private
     * @param {Set<number>[]} adjacency - Mutated in place.
     */
    _connect_components(adjacency) {
        const N = this._N;
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        /** @type {DisjointSet<number>} */
        const components = new DisjointSet(Array.from({ length: N }, (_, i) => i));
        for (let i = 0; i < N; ++i) {
            for (const j of adjacency[i]) components.union(i, j);
        }

        const roots = new Set();
        for (let i = 0; i < N; ++i) roots.add(components.find(i));
        if (roots.size <= 1) return;

        // Only now is the O(N²) scan worth it: find the closest pair of points between every two
        // components, then span the components with those edges.
        /** @type {Map<string, WeightedEdge>} */
        const best = new Map();
        for (let i = 0; i < N; ++i) {
            const root_i = /** @type {number} */ (components.find(i));
            const x_i = this.X.row(i);
            for (let j = i + 1; j < N; ++j) {
                const root_j = /** @type {number} */ (components.find(j));
                if (root_i === root_j) continue;
                const key = root_i < root_j ? `${root_i},${root_j}` : `${root_j},${root_i}`;
                const dist = metric(x_i, this.X.row(j));
                const current = best.get(key);
                if (current === undefined || dist < current[2]) best.set(key, [i, j, dist]);
            }
        }

        for (const [u, v] of minimum_spanning_tree([...best.values()], N)) {
            adjacency[u].add(v);
            adjacency[v].add(u);
        }
    }

    /**
     * Counts, for every vertex, how many of its neighbors carry each of the `q` labels.
     *
     * `U[i * q + l]` is `U_i(l)` in the paper.
     *
     * @private
     * @param {number[][]} adjacency
     * @returns {Float64Array}
     */
    _neighbor_counts(adjacency) {
        const N = this._N;
        const q = /** @type {number} */ (this._q);
        const labels = /** @type {Int32Array} */ (this._labels);
        const U = new Float64Array(N * q);
        for (let i = 0; i < N; ++i) {
            const offset = i * q;
            for (const j of adjacency[i]) U[offset + labels[j]] += 1;
        }
        return U;
    }

    /**
     * Gibbs weights `w_i(l) = exp(β ⋅ U_i(l))`, shifted to avoid overflow.
     *
     * Every quantity built from `w` here is a ratio of equally-weighted sums, so scaling the whole
     * vector by `exp(-β ⋅ max_l U_i(l))` leaves φ and ψ untouched while keeping `exp` in range even
     * for a large β on a dense neighborhood.
     *
     * @private
     * @param {Float64Array} U
     * @param {number} offset - Start of vertex `i`'s counts in `U`.
     * @param {number} q
     * @param {number} beta
     * @param {Float64Array} out - Length `q` scratch buffer.
     */
    _gibbs_weights(U, offset, q, beta, out) {
        let max_U = -Infinity;
        for (let l = 0; l < q; ++l) {
            if (U[offset + l] > max_U) max_U = U[offset + l];
        }
        for (let l = 0; l < q; ++l) out[l] = Math.exp(beta * (U[offset + l] - max_U));
    }

    /**
     * Maximum pseudo-likelihood estimate of the inverse temperature β.
     *
     * Solves the pseudo-likelihood equation (3), which sets the observed energy equal to the
     * expected energy:
     *
     * ```
     * f(β) = Σ_i U_i(x_i) - Σ_i E_β[U_i] = 0
     * ```
     *
     * `f` is the negative of a sum of variances, hence non-increasing, so plain bisection is enough
     * and needs no derivatives — the same role Brent's method plays in the paper, without the extra
     * machinery. `f(0) <= 0` means the labels agree with their neighborhoods no better than chance,
     * which leaves β = 0 as the estimate.
     *
     * @private
     * @param {Float64Array} U
     * @returns {number}
     */
    _estimate_beta(U) {
        const N = this._N;
        const q = /** @type {number} */ (this._q);
        const labels = /** @type {Int32Array} */ (this._labels);
        const w = new Float64Array(q);

        let observed = 0;
        for (let i = 0; i < N; ++i) observed += U[i * q + labels[i]];

        /** @param {number} beta */
        const f = (beta) => {
            let expected = 0;
            for (let i = 0; i < N; ++i) {
                const offset = i * q;
                this._gibbs_weights(U, offset, q, beta, w);
                let sum_w = 0;
                let sum_Uw = 0;
                for (let l = 0; l < q; ++l) {
                    sum_w += w[l];
                    sum_Uw += U[offset + l] * w[l];
                }
                expected += sum_Uw / sum_w;
            }
            return observed - expected;
        };

        let lower = 0;
        if (f(lower) <= 0) return lower;

        let upper = 1;
        // β above ~30 is numerically indistinguishable from the zero-temperature limit here.
        while (upper < 32 && f(upper) > 0) upper *= 2;
        if (f(upper) > 0) return upper;

        for (let iter = 0; iter < 64; ++iter) {
            const mid = (lower + upper) / 2;
            if (f(mid) > 0) {
                lower = mid;
            } else {
                upper = mid;
            }
            if (upper - lower < 1e-8) break;
        }
        return (lower + upper) / 2;
    }

    /**
     * Information curvature `S_i = -ψ_i / φ_i` at every vertex, normalised to `(0, 1]`.
     *
     * The paper states φ (17) and ψ (23) as Kronecker products over `q ⨯ q` matrices, but both
     * collapse: since summing every entry of `a ⊗ aᵀ` is just `(Σa)²`, the double sums reduce to
     *
     * ```
     * φ_i = (Σ_l v_l w_l / Σ_l w_l)²                                  — squared observed/expected energy gap
     * ψ_i = [(Σ_l U_l² w_l)(Σ_l w_l) - (Σ_l U_l w_l)²] / (Σ_l w_l)²   — the variance of U under w
     * ```
     *
     * which is O(q) per vertex instead of O(q²), and allocates nothing.
     *
     * `epsilon` floors the denominator, following the author's reference implementation. It is load
     * bearing: both φ and ψ vanish for an interior point as β grows, and by β ≈ 8 ψ has underflowed
     * to exactly 0, so `S` is pinned at the top of the range rather than the bottom. That inverts the
     * interior/boundary reading in §VI-C for anything above β ≈ 2, which well-separated clusters
     * comfortably exceed. The formulas are implemented as published in either regime; see
     * `docs/dimred/minfotree.md` for the measurements.
     *
     * @private
     * @param {Float64Array} U
     * @param {number} beta
     * @returns {Float64Array}
     */
    _information_curvature(U, beta) {
        const N = this._N;
        const q = /** @type {number} */ (this._q);
        const labels = /** @type {Int32Array} */ (this._labels);
        const w = new Float64Array(q);
        const S = new Float64Array(N);
        const EPSILON = /** @type {number} */ (this.parameter("epsilon"));

        for (let i = 0; i < N; ++i) {
            const offset = i * q;
            this._gibbs_weights(U, offset, q, beta, w);

            const U_xi = U[offset + labels[i]];
            let sum_w = 0;
            let sum_Uw = 0;
            let sum_UUw = 0;
            let sum_vw = 0;
            for (let l = 0; l < q; ++l) {
                const U_l = U[offset + l];
                const w_l = w[l];
                sum_w += w_l;
                sum_Uw += U_l * w_l;
                sum_UUw += U_l * U_l * w_l;
                sum_vw += (U_xi - U_l) * w_l;
            }

            const mean_v = sum_vw / sum_w;
            const phi = mean_v * mean_v;
            const psi = (sum_UUw * sum_w - sum_Uw * sum_Uw) / (sum_w * sum_w);
            S[i] = -psi / (phi + EPSILON);
        }

        let min = Infinity;
        let max = -Infinity;
        for (let i = 0; i < N; ++i) {
            if (S[i] < min) min = S[i];
            if (S[i] > max) max = S[i];
        }
        const range = max - min;
        for (let i = 0; i < N; ++i) {
            S[i] = range > 0 ? EPSILON + (S[i] - min) / range : EPSILON;
        }
        return S;
    }

    /**
     * All-pairs shortest paths over the Minimum Information Tree.
     *
     * A tree has exactly one path between any two vertices, so a depth-first walk from each source
     * settles every distance in O(N) — no priority queue and no relaxation, unlike the Dijkstra
     * {@link ISOMAP} needs on its general k-NNG.
     *
     * @private
     * @param {WeightedEdge[]} edges
     * @returns {Matrix}
     */
    _tree_distances(edges) {
        const N = this._N;
        /** @type {[number, number][][]} */
        const adjacency = Array.from({ length: N }, () => []);
        for (const [u, v, weight] of edges) {
            adjacency[u].push([v, weight]);
            adjacency[v].push([u, weight]);
        }

        const D = new Matrix(N, N, 0);
        const stack = new Int32Array(N);
        const visited = new Uint8Array(N);
        for (let source = 0; source < N; ++source) {
            visited.fill(0);
            visited[source] = 1;
            stack[0] = source;
            let top = 1;
            while (top > 0) {
                const u = stack[--top];
                const dist_u = D.entry(source, u);
                for (const [v, weight] of adjacency[u]) {
                    if (visited[v]) continue;
                    visited[v] = 1;
                    D.set_entry(source, v, dist_u + weight);
                    stack[top++] = v;
                }
            }
        }
        return D;
    }

    /** Runs steps 1-5: labels, k-NNG, β, curvature, information graph and its spanning tree. */
    init() {
        const alpha = /** @type {number} */ (this.parameter("alpha"));

        this._labels = this._make_labels();
        const adjacency = this._make_knn_graph();
        const U = this._neighbor_counts(adjacency);
        this._beta = this._estimate_beta(U);
        this._curvature = this._information_curvature(U, this._beta);

        const S = this._curvature;
        const labels = this._labels;
        /** @type {WeightedEdge[]} */
        const information_graph = [];
        for (let u = 0; u < this._N; ++u) {
            for (const v of adjacency[u]) {
                if (v <= u) continue;
                const weight = S[u] + S[v];
                information_graph.push([u, v, labels[u] === labels[v] ? alpha * weight : weight]);
            }
        }

        this._edges = minimum_spanning_tree(information_graph, this._N);
        this._is_initialized = true;
        return this;
    }

    /**
     * Computes the projection.
     *
     * @returns {T}
     */
    transform() {
        this.check_init();
        const d = /** @type {number} */ (this.parameter("d"));
        const layout = /** @type {"kamada_kawai" | "MDS"} */ (this.parameter("layout"));

        const D = this._tree_distances(/** @type {WeightedEdge[]} */ (this._edges));
        const kkmds = new KKMDS(D, {
            d,
            metric: "precomputed",
            init_DR: "MDS",
            // `layout: "MDS"` stops at the classical-MDS warm start, which is what zero gradient
            // steps leaves behind.
            iterations: layout === "MDS" ? 0 : /** @type {number} */ (this.parameter("iterations")),
            seed: /** @type {number} */ (this.parameter("seed")),
        });
        kkmds.transform();
        this.Y = kkmds.Y;

        return this.projection;
    }

    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new MINFOTree(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new MINFOTree(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new MINFOTree(X, parameters);
        return dr.transform_async();
    }
}
