import { Heap } from "../datastructure/index.js";
import { spatial_tree } from "../knn/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { WASM_MIN_PARALLEL_ROWS } from "../wasm/thresholds.js";
import { parallel_available, run_row_range } from "../wasm/worker_pool.js";
import { DR } from "./DR.js";
import { wasmDijkstraAPSP } from "./ISOMAP.wasm.js";
import { SMACOF } from "./SMACOF.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersISOMAP} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/**
 * Isomap (Isometric Mapping)
 *
 * A nonlinear dimensionality reduction algorithm that uses geodesic distances
 * between points on a manifold to perform embedding. It builds a neighborhood
 * graph and uses MDS on the shortest-path distances.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersISOMAP>
 * @category Dimensionality Reduction
 * @see {@link LLE} for another nonlinear alternative
 */
export class ISOMAP extends DR {
    /**
     * Isometric feature mapping (ISOMAP).
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersISOMAP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X, parameters = {}) {
        /** @type {ParametersISOMAP} */
        const defaults = {
            neighbors: -Infinity,
            d: 2,
            metric: euclidean,
            seed: 1212,
            project: "MDS",
            eig_args: {},
        };
        super(X, defaults, parameters);

        this.defaults = defaults;

        if (this._parameters.neighbors === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N - 1));
        }

        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Runs the all-pairs geodesic step across the worker pool.
     *
     * The pool needs a `SharedArrayBuffer` to write into, so the result is copied into `out` rather
     * than computed in place. That copy is far cheaper than the Dijkstra runs it parallelizes.
     *
     * @private
     * @param {Int32Array} neighbors
     * @param {Float64Array} distances
     * @param {Float64Array} out
     * @param {number} rows
     * @param {number} maxK
     * @returns {boolean} True if the pool produced the distances.
     */
    _dijkstra_parallel(neighbors, distances, out, rows, maxK) {
        if (rows < WASM_MIN_PARALLEL_ROWS || !parallel_available()) return false;
        try {
            const shared = new Float64Array(new SharedArrayBuffer(rows * rows * 8));
            const inputs = /** @type {{ data: Float64Array | Int32Array; kind: "f64" | "i32" }[]} */ ([
                { data: neighbors, kind: "i32" },
                { data: distances, kind: "f64" },
            ]);
            if (!run_row_range("dijkstra_apsp_range_f64", inputs, shared, [rows, maxK], rows, rows)) return false;
            out.set(shared);
            return true;
        } catch {
            return false;
        }
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
     * @returns {T}
     */
    transform() {
        this.check_init();
        const X = this.X;
        const rows = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix(rows, rows, 0);
        D.shape = [rows, rows, (i, j) => (i <= j ? metric(X.row(i), X.row(j)) : D.entry(j, i))];

        /** @type {{ index: number; distance: number }[][]} */
        const kNearestNeighbors = [];
        const tree = spatial_tree(X.to2dArray(), {
            metric,
            seed: /** @type {number} */ (this.parameter("seed")),
        });
        for (let i = 0; i < rows; ++i) {
            kNearestNeighbors.push(
                tree.search_by_index(i, neighbors).map((n) => ({
                    index: n.index,
                    distance: n.distance,
                })),
            );
        }

        // ISOMAP requires an undirected/symmetric nearest neighbor graph.
        // If i is a nearest neighbor of j, then j should be connected to i as well.
        for (let i = 0; i < rows; ++i) {
            for (const neighbor of kNearestNeighbors[i]) {
                const j = neighbor.index;
                const d = neighbor.distance;
                const reciprocal_edge = kNearestNeighbors[j].find((n) => n.index === i);
                if (!reciprocal_edge) {
                    kNearestNeighbors[j].push({ index: i, distance: d });
                }
            }
        }

        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths using Dijkstra's algorithm
        // TODO: make extern
        const G = new Matrix(rows, rows, Infinity);

        // Convert kNearestNeighbors list into flattened arrays for WASM Dijkstra execution
        let maxK = 0;
        for (let i = 0; i < rows; ++i) {
            if (kNearestNeighbors[i].length > maxK) maxK = kNearestNeighbors[i].length;
        }

        const flatNeighbors = new Int32Array(rows * maxK).fill(-1);
        const flatDistances = new Float64Array(rows * maxK).fill(Infinity);
        for (let i = 0; i < rows; ++i) {
            const list = kNearestNeighbors[i];
            const i_k = i * maxK;
            for (let idx = 0; idx < list.length; ++idx) {
                flatNeighbors[i_k + idx] = list[idx].index;
                flatDistances[i_k + idx] = list[idx].distance;
            }
        }

        if (
            !this._dijkstra_parallel(flatNeighbors, flatDistances, G.values, rows, maxK) &&
            !wasmDijkstraAPSP(flatNeighbors, flatDistances, G.values, rows, maxK)
        ) {
            // JS fallback using Heap
            for (let i = 0; i < rows; ++i) {
                G.set_entry(i, i, 0);
                const H = new Heap([{ index: i, distance: 0 }], (d) => d.distance, "min");

                while (!H.empty) {
                    const item = H.pop();
                    if (!item) break;

                    const u = item.element.index;
                    const dist_u = item.element.distance;

                    if (dist_u > G.entry(i, u)) continue;

                    for (const neighbor of kNearestNeighbors[u]) {
                        const v = neighbor.index;
                        const alt = dist_u + neighbor.distance;
                        if (alt < G.entry(i, v)) {
                            G.set_entry(i, v, alt);
                            H.push({ index: v, distance: alt });
                        }
                    }
                }
            }
        }

        let max_val = 0;
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < rows; j++) {
                const val = G.entry(i, j);
                if (val !== Infinity && val > max_val) max_val = val;
            }
        }
        const big_val = max_val * 10;

        const project = /** @type {"MDS" | "SMACOF"} */ (this.parameter("project"));

        if (project === "SMACOF") {
            // Apply SMACOF metric MDS to the distance matrix directly
            const D_matrix = new Matrix(rows, rows, (i, j) => {
                const val = G.entry(i, j);
                return val === Infinity ? big_val : val;
            });
            const smacof = new SMACOF(D_matrix, {
                metric: "precomputed",
                d,
                seed: this.parameter("seed"),
            });
            smacof.transform();
            this.Y = smacof.Y;
        } else {
            // "MDS" (Classical MDS) via Eigendecomposition of double-centered squared distance matrix
            const D_sq = new Matrix(rows, rows, (i, j) => {
                let val = G.entry(i, j);
                if (val === Infinity) val = big_val;
                return val * val;
            });

            const ai_ = D_sq.meanCols();
            const a_j = D_sq.meanRows();
            const a__ = D_sq.mean();
            const B = new Matrix(rows, rows, (i, j) => -0.5 * (D_sq.entry(i, j) - ai_[i] - a_j[j] + a__));

            // compute d eigenvectors
            const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
            this.Y = Matrix.from(V).transpose();
        }
        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        return dr.transform_async();
    }
}
