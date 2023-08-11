import { Matrix, linspace, norm } from "../matrix/index.js";
import { euclidean, euclidean_squared } from "../metrics/index.js";
import { neumair_sum } from "../numerical/index.js";
import { DR } from "./DR.js";
import { PCA } from "./index.js";

export class SQDMDS extends DR {
    /**
     * SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE and UMAP.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @param {Matrix|number[][]} X
     * @param {object} [parameters]
     * @param {number} [parameters.d=2]
     * @param {function} [parameters.metric = euclidean]
     * @param {number} [parameters.decay_start = 0.1] - Percentage of iterations using exaggeration phase. If random init: it is recommended to start the decay later to give the time for the global config to adjust with big steps.
     * @param {number} [parameters.decay_cte = 0.34] - Controls the decay of the learning parameter.
     * @param {object} [parameters.init_DR]
     * @returns {SQDMDS}
     * @see {@link https://arxiv.org/pdf/2202.12087.pdf}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                d: 2,
                metric: euclidean,
                seed: 1212,
                decay_start: 0.1,
                decay_cte: 0.34, // 0.34
                init_DR: { type: "random" },
            },
            parameters
        );

        return this;
    }

    /**
     * @private
     */
    init() {
        const N = this._N;
        const d = this.parameter("d");

        // initialize helpers.
        this._add = this.__add(d);
        this._sub_div = this.__sub_div(d);
        this._minus = this.__minus(d);
        this._mult = this.__mult(d);
        this._LR_init = Math.max(2, 0.005 * N);
        this._LR = this._LR_init;
        this._offset = -Math.exp(-1 / this.parameter("decay_cte"));
        this._momentums = new Matrix(N, d, 0);
        this._grads = new Matrix(N, d, 0);
        this._indices = linspace(0, N - 1);
        // initialize projection.
        const R = this._randomizer;
        this.Y = new Matrix(N, d, () => R.random - 0.5);

        // preparing metric for optimization.
        const this_metric = this.parameter("metric");
        if (this_metric === "precomputed") {
            this._HD_metric = function (i, j, X) {
                return X.entry(i, j);
            };
            this._HD_metric_exaggeration = function (i, j, X) {
                return Math.pow(X.entry(i, j), 2);
            };
        } else {
            this._HD_metric = function (i, j, X) {
                return this_metric(X.row(i), X.row(j));
            };
            if (this_metric == euclidean) {
                this._HD_metric_exaggeration = function (i, j, X) {
                    return euclidean_squared(X.row(i), X.row(j));
                };
            } else {
                this._HD_metric_exaggeration = function (i, j, X) {
                    return Math.pow(this_metric(X.row(i), X.row(j)), 2);
                };
            }
        }
        return;
    }

    /**
     * Computes the projection.
     * @param {number} [iterations=500] - number of iterations.
     * @returns {Matrix|number[][]} the projection.
     */
    transform(iterations = 500) {
        this.check_init();
        this._decay_start = Math.round(this.parameter("decay_start") * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
        }
        return this.projection;
    }

    /**
     * Computes the projection.
     * @param {number} [iterations=500] - number of iterations.
     * @yields {Matrix|number[][]} the intermediate steps of the projection.
     */
    *generator(iterations = 500) {
        this.check_init();
        this._decay_start = Math.round(this.parameter("decay_start") * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * Performs an optimization step.
     * @private
     * @param {number} i - Acutal iteration.
     * @param {number} iterations - Number of iterations.
     */
    _step(i, iterations) {
        const decay_start = this._decay_start;
        if (i > decay_start) {
            const decay_cte = this.parameter("decay_cte");
            const offset = this._offset;
            const ratio = (i - decay_start) / (iterations - decay_start);
            this._LR = this._LR_init * (Math.exp(-(ratio * ratio) / decay_cte) + offset);
            this._distance_exaggeration = false;
        } else {
            this._distance_exaggeration = true;
        }
        this._nestrov_iteration(this._distance_exaggeration);
    }

    /**
     * Creates quartets of non overlapping indices.
     * @private
     * @returns {number[][]}
     */
    __quartets() {
        const N = this._N;
        const max_N = N - (N % 4);
        const R = this._randomizer;
        const shuffled_indices = R.choice(this._indices, max_N);
        const result = [];
        for (let i = 0; i < max_N; i += 4) {
            result.push(Uint32Array.of(shuffled_indices[i], shuffled_indices[i + 1], shuffled_indices[i + 2], shuffled_indices[i + 3]));
        }
        return result;
    }

    /**
     * Computes and applies gradients, and updates momentum.
     * @private
     * @param {boolean} distance_exaggeration
     */
    _nestrov_iteration(distance_exaggeration) {
        const momentums = this._momentums.mult(0.99, { inline: true });
        const LR = this._LR;
        const grads = this._fill_MDS_grads(this.Y.add(momentums), this._grads, distance_exaggeration);
        const [n, d] = momentums.shape;
        for (let i = 0; i < n; ++i) {
            const g_i = grads.row(i);
            const g_i_norm = norm(g_i);
            if (g_i_norm == 0) continue;
            const mul = LR / g_i_norm;
            const m_i = momentums.row(i);
            for (let j = 0; j < d; ++j) {
                m_i[j] -= mul * g_i[j];
            }
        } // momentums -= (LR / norm) * grads
        this.Y.add(momentums, { inline: true });
    }

    /**
     * Computes the gradients.
     * @param {Matrix} Y - The Projection.
     * @param {Matrix} grads - The gradients.
     * @param {boolean} [exaggeration = false] - Whether or not to use early exaggeration.
     * @param {boolean} [zero_grad = true] - Whether or not to reset the gradient in the beginning.
     * @returns {Matrix} the gradients.
     */
    _fill_MDS_grads(Y, grads, exaggeration = false, zero_grad = true) {
        if (zero_grad) {
            // compute new gradients
            grads.values.fill(0);
        }
        const add = this._add;
        const X = this.X;
        let HD_metric;
        if (exaggeration == true) {
            HD_metric = this._HD_metric_exaggeration;
        } else {
            HD_metric = this._HD_metric;
        }

        const D_quartet = new Float64Array(6);
        const quartets = this.__quartets();
        for (const [i, j, k, l] of quartets) {
            // compute quartet's HD distances.
            D_quartet[0] = HD_metric(i, j, X);
            D_quartet[1] = HD_metric(i, k, X);
            D_quartet[2] = HD_metric(i, l, X);
            D_quartet[3] = HD_metric(j, k, X);
            D_quartet[4] = HD_metric(j, l, X);
            D_quartet[5] = HD_metric(k, l, X);

            const D_quartet_sum = neumair_sum(D_quartet);

            if (D_quartet_sum > 0) {
                for (let i = 0; i < 6; ++i) {
                    D_quartet[i] /= D_quartet_sum;
                    D_quartet[i] += 1e-11;
                }
            }
            const [gi, gj, gk, gl] = this._compute_quartet_grads(Y, [i, j, k, l], D_quartet);

            // add is inline, row acces the matrix
            add(grads.row(i), gi);
            add(grads.row(j), gj);
            add(grads.row(k), gk);
            add(grads.row(l), gl);
        }
        return grads;
    }

    /**
     * Quartet gradients for a projection.
     * @private
     * @param {Matrix} Y - The acutal projection.
     * @param {number[]} quartet - The indices of the quartet.
     * @param {number[]} D_hd - The high-dimensional distances of the quartet.
     * @returns {number[][]} the gradients for the quartet.
     */
    _compute_quartet_grads(Y, quartet, [p_ab, p_ac, p_ad, p_bc, p_bd, p_cd]) {
        const [a, b, c, d] = quartet.map((index) => Y.row(index));
        // LD distances, add a small number just in case
        const d_ab = euclidean(a, b) + 1e-12;
        const d_ac = euclidean(a, c) + 1e-12;
        const d_ad = euclidean(a, d) + 1e-12;
        const d_bc = euclidean(b, c) + 1e-12;
        const d_bd = euclidean(b, d) + 1e-12;
        const d_cd = euclidean(c, d) + 1e-12;
        const sum_LD_dist = neumair_sum([d_ab, d_ac, d_ad, d_bc, d_bd, d_cd]);

        // for each element of the sum: use the same gradient function and just permute the points given in input.
        const [gA1, gB1, gC1, gD1] = this._ABCD_grads(a, b, c, d, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd, p_ab, sum_LD_dist);
        const [gA2, gC2, gB2, gD2] = this._ABCD_grads(a, c, b, d, d_ac, d_ab, d_ad, d_bc, d_cd, d_bd, p_ac, sum_LD_dist);
        const [gA3, gD3, gC3, gB3] = this._ABCD_grads(a, d, c, b, d_ad, d_ac, d_ab, d_cd, d_bd, d_bc, p_ad, sum_LD_dist);
        const [gB4, gC4, gA4, gD4] = this._ABCD_grads(b, c, a, d, d_bc, d_ab, d_bd, d_ac, d_cd, d_ad, p_bc, sum_LD_dist);
        const [gB5, gD5, gA5, gC5] = this._ABCD_grads(b, d, a, c, d_bd, d_ab, d_bc, d_ad, d_cd, d_ac, p_bd, sum_LD_dist);
        const [gC6, gD6, gA6, gB6] = this._ABCD_grads(c, d, a, b, d_cd, d_ac, d_bc, d_ad, d_bd, d_ab, p_cd, sum_LD_dist);

        const add = this._add;
        const gA = add(gA1, gA2, gA3, gA4, gA5, gA6);
        const gB = add(gB1, gB2, gB3, gB4, gB5, gB6);
        const gC = add(gC1, gC2, gC3, gC4, gC5, gC6);
        const gD = add(gD1, gD2, gD3, gD4, gD5, gD6);

        return [gA, gB, gC, gD];
    }

    /**
     * Gradients for one element of the loss function's sum.
     * @private
     */
    _ABCD_grads(a, b, c, d, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd, p_ab, sum_LD_dist) {
        const ratio = d_ab / sum_LD_dist;
        const twice_ratio = 2 * ((p_ab - ratio) / sum_LD_dist);
        const minus = this._minus;
        const add = this._add;
        const mult = this._mult;
        const sub_div = this._sub_div;
        // no side effects because sub_div creates new arrays, and the inline functions work on this new created arrays.
        const gA = mult(minus(mult(add(sub_div(a, b, d_ab), sub_div(a, c, d_ac), sub_div(a, d, d_ad)), ratio), sub_div(a, b, d_ab)), twice_ratio);
        const gB = mult(minus(mult(add(sub_div(b, a, d_ab), sub_div(b, c, d_bc), sub_div(b, d, d_bd)), ratio), sub_div(b, a, d_ab)), twice_ratio);
        const gC = mult(add(sub_div(c, a, d_ac), sub_div(c, b, d_bc), sub_div(c, d, d_cd)), ratio * twice_ratio);
        const gD = mult(add(sub_div(d, a, d_ad), sub_div(d, b, d_bd), sub_div(d, c, d_cd)), ratio * twice_ratio);
        return [gA, gB, gC, gD];
    }

    /**
     * Inline!
     */
    __minus(d) {
        return (a, b) => {
            for (let i = 0; i < d; ++i) {
                a[i] -= b[i];
            }
            return a;
        };
    }

    /**
     * Inline!
     */
    __add(d) {
        return (...summands) => {
            const n = summands.length;
            const s1 = summands[0];
            for (let j = 1; j < n; ++j) {
                const summand = summands[j];
                for (let i = 0; i < d; ++i) {
                    s1[i] += summand[i];
                }
            }
            return s1;
        };
    }

    /**
     * Inline!
     */
    __mult(d) {
        return (a, v) => {
            for (let i = 0; i < d; ++i) {
                a[i] *= v;
            }
            return a;
        };
    }

    /**
     * Creates a new array <code>(x - y) / div</code>
     */
    __sub_div(d) {
        return (x, y, div) => {
            return Float64Array.from({ length: d }, (_, i) => (x[i] - y[i]) / div);
        };
    }
}
