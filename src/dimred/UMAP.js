import { BallTree, NaiveKNN } from "../knn/index.js";
import { linspace, Matrix } from "../matrix/index.js";
import { euclidean, euclidean_squared } from "../metrics/index.js";
import { neumair_sum } from "../numerical/index.js";
import { powell } from "../optimization/index.js";
import { max } from "../util/index.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersUMAP} from "./index.js" */

/**
 * Uniform Manifold Approximation and Projection (UMAP)
 *
 * A novel manifold learning technique for dimensionality reduction. UMAP is constructed
 * from a theoretical framework based on Riemannian geometry and algebraic topology.
 * It is often faster than t-SNE while preserving more of the global structure.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersUMAP>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/1802.03426|UMAP Paper}
 * @see {@link TSNE} for a similar visualization technique
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const umap = new druid.UMAP(X, {
 *     n_neighbors: 15,
 *     min_dist: 0.1,
 *     d: 2,
 *     seed: 42
 * });
 *
 * const Y = umap.transform(500); // 500 iterations
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
export class UMAP extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersUMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            {
                n_neighbors: 15,
                local_connectivity: 1,
                min_dist: 1,
                d: 2,
                metric: euclidean,
                seed: 1212,
                _spread: 1,
                _set_op_mix_ratio: 1,
                _repulsion_strength: 1,
                _negative_sample_rate: 5,
                _n_epochs: 350,
                _initial_alpha: 1,
            },
            parameters,
        );
        [this._N, this._D] = this.X.shape;
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const local_connectivity = /** @type {number} */ (this.parameter("local_connectivity"));
        const d = /** @type {number} */ (this.parameter("d"));
        /* let n_neighbors = Math.min(this._N - 1, parameters.n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(this.parameter("local_connectivity"), n_neighbors - 1)); */
        if (n_neighbors > this._N) {
            throw new Error(
                `Parameter n_neighbors (=${n_neighbors}) needs to be smaller than dataset size (N=${this._N})!`,
            );
        }
        if (local_connectivity > n_neighbors) {
            throw new Error(
                `Parameter local_connectivity (=${local_connectivity}) needs to be smaller than parameter n_neighbors (=${n_neighbors})`,
            );
        }
        this._iter = 0;
        const randomizer = this._randomizer;
        this.Y = new Matrix(this._N, d, () => randomizer.random);
    }

    /**
     * @private
     * @param {number} spread
     * @param {number} min_dist
     * @returns {number[]}
     */
    _find_ab_params(spread, min_dist) {
        /** @type {(x: number, a: number, b: number) => number} */
        const curve = (x, a, b) => 1 / (1 + a * x ** (2 * b));
        const xv = linspace(0, spread * 3, 300);
        const yv = linspace(0, spread * 3, 300);

        for (let i = 0, n = xv.length; i < n; ++i) {
            const xv_i = xv[i];
            yv[i] = xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread);
        }

        /** @type {(p: [number, number]) => number} */
        const err = (p) => {
            const error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map((e) => e * e)));
        };

        return powell(err, [1, 1]);
    }

    /**
     * @private
     * @param {{ element: Float64Array; index: number; distance: number }[][]} distances
     * @param {number[]} sigmas
     * @param {number[]} rhos
     * @returns {{ element: Float64Array; index: number; distance: number }[][]}
     */
    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            const rho = rhos[i];
            const curr_dist = distances[i];
            for (let j = 0, m = curr_dist.length; j < m; ++j) {
                const v = curr_dist[j].distance - rho;
                curr_dist[j].distance = v > 0 ? Math.exp(-v / sigmas[i]) : 1.0;
            }
        }
        return distances;
    }

    /**
     * @private
     * @param {NaiveKNN<Float64Array> | BallTree<Float64Array>} knn
     * @param {number} k
     * @returns {{
     *     distances: { element: Float64Array; index: number; distance: number }[][];
     *     sigmas: number[];
     *     rhos: number[];
     * }}
     */
    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = /** @type {number} */ (this._parameters.local_connectivity);
        const metric = /** @type {Metric | "precomputed"} */ (this._parameters.metric);
        const target = Math.log2(k);
        const rhos = [];
        const sigmas = [];
        const X = this.X;
        const N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());

        /** @type {{ element: Float64Array; index: number; distance: number }[][]} */
        const distances = [];
        if (metric === "precomputed" || knn instanceof NaiveKNN) {
            for (let i = 0; i < N; ++i) {
                distances.push(knn.search_by_index(i, k).reverse());
            }
        } else {
            for (const x_i of X) {
                distances.push(knn.search(x_i, k).reverse());
            }
        }

        const index = Math.floor(local_connectivity);
        const interpolation = local_connectivity - index;
        for (let i = 0; i < N; ++i) {
            let lo = 0;
            let hi = Infinity;
            let mid = 1;
            let rho = 0;

            const search_result = distances[i];
            const non_zero_dist = search_result.filter((d) => d.distance > 0);
            const non_zero_dist_length = non_zero_dist.length;
            if (non_zero_dist_length >= local_connectivity) {
                if (index > 0) {
                    rho = non_zero_dist[index - 1].distance;
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rho += interpolation * (non_zero_dist[index].distance - non_zero_dist[index - 1].distance);
                    }
                } else {
                    rho = interpolation * non_zero_dist[0].distance;
                }
            } else if (non_zero_dist_length > 0) {
                rho = non_zero_dist[non_zero_dist_length - 1].distance;
            }
            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    const d = search_result[j].distance - rho;
                    psum += d > 0 ? Math.exp(-(d / mid)) : 1;
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    [hi, mid] = [mid, (lo + hi) / 2];
                } else {
                    if (hi === Infinity) {
                        [lo, mid] = [mid, mid * 2];
                    } else {
                        [lo, mid] = [mid, (lo + hi) / 2];
                    }
                }
            }

            //let mean_d = null;
            if (rho > 0) {
                const mean_ithd = search_result.reduce((a, b) => a + b.distance, 0) / search_result.length;
                if (mid < MIN_K_DIST_SCALE * mean_ithd) {
                    mid = MIN_K_DIST_SCALE * mean_ithd;
                }
            } else {
                const mean_d = distances.reduce(
                    (acc, res) => acc + res.reduce((a, b) => a + b.distance, 0) / res.length,
                    0,
                );
                if (mid < MIN_K_DIST_SCALE * mean_d) {
                    mid = MIN_K_DIST_SCALE * mean_d;
                }
            }
            rhos[i] = rho;
            sigmas[i] = mid;
        }
        return {
            distances: distances,
            sigmas: sigmas,
            rhos: rhos,
        };
    }

    /**
     * @private
     * @param {Matrix} X
     * @param {number} n_neighbors
     * @returns {Matrix}
     */
    _fuzzy_simplicial_set(X, n_neighbors) {
        const N = X.shape[0];
        const metric = /** @type {Metric | "precomputed"} */ (this._parameters.metric);
        const _set_op_mix_ratio = /** @type {number} */ (this._parameters._set_op_mix_ratio);

        const knn =
            metric === "precomputed"
                ? new NaiveKNN(X.to2dArray(), {
                      metric: "precomputed",
                      seed: /** @type {number} */ (this._parameters.seed),
                  })
                : new BallTree(X.to2dArray(), {
                      metric,
                      seed: /** @type {number} */ (this._parameters.seed),
                  });
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        const result = new Matrix(N, N, "zeros");
        for (let i = 0; i < N; ++i) {
            const distances_i = distances[i];
            for (let j = 0; j < distances_i.length; ++j) {
                result.set_entry(i, distances_i[j].index, distances_i[j].distance);
            }
        }

        const transposed_result = result.T;
        const prod_matrix = result.mult(transposed_result);
        return result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(_set_op_mix_ratio)
            .add(prod_matrix.mult(1 - _set_op_mix_ratio));
    }

    /**
     * @private
     * @param {number} n_epochs
     * @returns {Float32Array}
     */
    _make_epochs_per_sample(n_epochs) {
        if (!this._weights) throw new Error("Call init() first!");
        const weights = this._weights;
        const result = new Float32Array(weights.length).fill(-1);
        const weight_scl = n_epochs / max(weights);
        weights.forEach((w, i) => {
            const sample = w * weight_scl;
            if (sample > 0) result[i] = Math.round(n_epochs / sample);
        });
        return result;
    }

    /**
     * @private
     * @param {Matrix} graph
     * @returns {{ rows: number[]; cols: number[]; data: number[] }}
     */
    _tocoo(graph) {
        const rows = [];
        const cols = [];
        const data = [];
        const [rows_n, cols_n] = graph.shape;
        for (let row = 0; row < rows_n; ++row) {
            for (let col = 0; col < cols_n; ++col) {
                const entry = graph.entry(row, col);
                if (entry !== 0) {
                    rows.push(row);
                    cols.push(col);
                    data.push(entry);
                }
            }
        }
        return {
            rows: rows,
            cols: cols,
            data: data,
        };
    }

    /**
     * Computes all necessary
     *
     * @returns {UMAP<T>}
     */
    init() {
        const _spread = /** @type {number} */ (this._parameters._spread);
        const min_dist = /** @type {number} */ (this._parameters.min_dist);
        const n_neighbors = /** @type {number} */ (this._parameters.n_neighbors);
        const _n_epochs = /** @type {number} */ (this._parameters._n_epochs);
        const _negative_sample_rate = /** @type {number} */ (this._parameters._negative_sample_rate);
        const [a, b] = this._find_ab_params(_spread, min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, n_neighbors);
        const { rows, cols, data: weights } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        this._weights = weights;
        this._epochs_per_sample = this._make_epochs_per_sample(_n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map((d) => d * _negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        return this;
    }

    graph() {
        this.check_init();
        return { cols: this._head, rows: this._tail, weights: this._weights };
    }

    /**
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {T}
     */
    transform(iterations = 350) {
        if (this.parameter("_n_epochs") !== iterations) {
            this.parameter("_n_epochs", iterations);
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {Generator<T, T, void>}
     */
    *generator(iterations = 350) {
        if (this.parameter("_n_epochs") !== iterations) {
            this.parameter("_n_epochs", iterations);
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @private
     * @param {number} x
     * @returns {number}
     */
    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

    /**
     * Performs the optimization step.
     *
     * @private
     * @param {Matrix} head_embedding
     * @param {Matrix} tail_embedding
     * @param {number[]} head
     * @param {number[]} tail
     * @returns {Matrix}
     */
    _optimize_layout(head_embedding, tail_embedding, head, tail) {
        const randomizer = this._randomizer;
        const _repulsion_strength = /** @type {number} */ (this.parameter("_repulsion_strength"));
        const dim = /** @type {number} */ (this.parameter("d"));
        const {
            _alpha: alpha,
            _a: a,
            _b: b,
            _epochs_per_sample: epochs_per_sample,
            _epochs_per_negative_sample: epochs_per_negative_sample,
            _epoch_of_next_negative_sample: epoch_of_next_negative_sample,
            _epoch_of_next_sample: epoch_of_next_sample,
            _clip: clip,
        } = this;
        if (
            alpha === undefined ||
            a === undefined ||
            b === undefined ||
            epochs_per_sample === undefined ||
            epochs_per_negative_sample === undefined ||
            epoch_of_next_negative_sample === undefined ||
            epoch_of_next_sample === undefined ||
            clip === undefined
        ) {
            throw new Error("call init() first!");
        }
        const tail_length = tail.length;

        for (let i = 0, n = epochs_per_sample.length; i < n; ++i) {
            if (epoch_of_next_sample[i] <= this._iter) {
                const j = head[i];
                const k = tail[i];
                const current = head_embedding.row(j);
                const other = tail_embedding.row(k);
                const dist = euclidean_squared(current, other);
                if (dist > 0) {
                    const grad_coeff = (-2 * a * b * dist ** (b - 1)) / (a * dist ** b + 1);
                    for (let d = 0; d < dim; ++d) {
                        const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                        current[d] += grad_d;
                        other[d] -= grad_d;
                    }
                }
                epoch_of_next_sample[i] += epochs_per_sample[i];
                const n_neg_samples = (this._iter - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i];
                for (let p = 0; p < n_neg_samples; ++p) {
                    const k = randomizer.random_int % tail_length;
                    const other = tail_embedding.row(tail[k]);
                    const dist = euclidean_squared(current, other);
                    if (dist > 0) {
                        const grad_coeff = (2 * _repulsion_strength * b) / ((0.01 + dist) * (a * dist ** b + 1));
                        for (let d = 0; d < dim; ++d) {
                            const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                            current[d] += grad_d;
                            other[d] -= grad_d;
                        }
                    } else if (j === k) {
                    }
                }
                epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
            }
        }
        return head_embedding;
    }

    /**
     * @private
     * @returns {Matrix}
     */
    next() {
        if (!this._head || !this._tail) throw new Error("Call init() first!");
        const iter = ++this._iter;
        const Y = this.Y;
        const _initial_alpha = /** @type {number} */ (this._parameters._initial_alpha);
        const _n_epochs = /** @type {number} */ (this._parameters._n_epochs);
        this._alpha = _initial_alpha * (1 - iter / _n_epochs);
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new UMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new UMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new UMAP(X, parameters);
        return dr.transform_async();
    }
}
