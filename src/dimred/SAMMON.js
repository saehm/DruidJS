import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { MDS, PCA } from "./index.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersPCA, ParametersMDS, ParametersSAMMON} from "./index.js" */
/** @typedef {"PCA" | "MDS" | "random"} AvailableInit */

/** @typedef {{ PCA: ParametersPCA; MDS: ParametersMDS; random: {} }} ChooseDR */

/**
 * Sammon's Mapping
 *
 * A nonlinear dimensionality reduction technique that minimizes a stress
 * function based on the ratio of pairwise distances in high and low dimensional spaces.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersSAMMON<AvailableInit>>
 * @category Dimensionality Reduction
 */
export class SAMMON extends DR {
    /** @type {Matrix | undefined} */
    distance_matrix;

    /**
     * SAMMON's Mapping
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters] - Object containing parameterization of the DR
     *   method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                magic: 0.1,
                d: 2,
                metric: euclidean,
                seed: 1212,
                init_DR: "random",
                init_parameters: {},
            },
            parameters,
        );
    }

    /**
     * Initializes the projection.
     *
     * @param {Matrix | undefined} D
     * @returns {asserts D is Matrix}
     */
    init(D) {
        const N = this.X.shape[0];
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const init_DR = /** @type {AvailableInit} */ (this.parameter("init_DR"));
        const DR_parameters = this.parameter("init_parameters");
        if (init_DR === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (init_DR === "PCA") {
            this.Y = Matrix.from(PCA.transform(this.X, /** @type {ParametersPCA} */ (DR_parameters)));
        } else if (init_DR === "MDS") {
            this.Y = Matrix.from(MDS.transform(this.X, /** @type {ParametersMDS} */ (DR_parameters)));
        } else {
            throw new Error('init_DR needs to be either "random" or a DR method!');
        }
        D = metric === "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        this.distance_matrix = D;
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {T} The projection of `X`.
     */
    transform(max_iter = 200) {
        this.check_init();
        if (!this.distance_matrix) this.init(this.distance_matrix);
        for (let j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimenionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection of
     *   `X`.
     */
    *generator(max_iter = 200) {
        this.check_init();
        if (!this.distance_matrix) this.init(this.distance_matrix);

        for (let j = 0; j < max_iter; ++j) {
            this._step();
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        if (!this.distance_matrix) this.init(this.distance_matrix);
        const MAGIC = /** @type {number} */ (this.parameter("magic"));
        const D = /** @type {Matrix} */ (this.distance_matrix);
        const N = this.X.shape[0];
        const d = /** @type {number} */ (this.parameter("d"));
        const Y = this.Y;

        const G = new Matrix(N, d, 0);

        const sum = new Float64Array(d);
        for (let i = 0; i < N; ++i) {
            const e1 = new Float64Array(d);
            const e2 = new Float64Array(d);
            const Yi = Y.row(i);
            for (let j = 0; j < N; ++j) {
                if (i === j) continue;
                const dX = D.entry(i, j);
                if (dX === 0) continue; // Skip identical points in high-dim

                const Yj = Y.row(j);
                const delta = new Float64Array(d);
                for (let k = 0; k < d; ++k) {
                    delta[k] = Yi[k] - Yj[k];
                }
                const dY = Math.max(euclidean(Yi, Yj), 1e-6);
                const dq = dX - dY;
                const dr = dX * dY;
                for (let k = 0; k < d; ++k) {
                    e1[k] += (delta[k] * dq) / dr;
                    e2[k] += (dq - (delta[k] ** 2 * (1 + dq / dY)) / dY) / dr;
                }
            }
            for (let k = 0; k < d; ++k) {
                const val = Y.entry(i, k) + ((MAGIC * e1[k]) / Math.abs(e2[k]) || 0);
                G.set_entry(i, k, val);
                sum[k] += val;
            }
        }
        for (let k = 0; k < d; ++k) {
            sum[k] /= N;
        }

        for (let i = 0; i < N; ++i) {
            for (let k = 0; k < d; ++k) {
                Y.set_entry(i, k, G.entry(i, k) - sum[k]);
            }
        }
        return Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new SAMMON(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new SAMMON(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new SAMMON(X, parameters);
        return dr.transform_async();
    }
}
