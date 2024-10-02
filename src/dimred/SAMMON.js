import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { PCA, MDS } from "./index.js";
import { distance_matrix } from "../matrix/index.js";

/**
 * @class
 * @alias SAMMON
 * @extends DR
 */
export class SAMMON extends DR {
    /**
     * SAMMON's Mapping
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {"PCA"|"MDS"|"random"} [parameters.init = "random"] - Either "PCA" or "MDS", with which SAMMON initialiates the projection. With "random" a random matrix gets used as starting point.
     * @param {object} [parameters.init_parameters] - Parameters for the {@link init}-DR method.
     * @param {number} [parameters.magic = 0.1] - learning rate for gradient descent.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { magic: 0.1, d: 2, metric: euclidean, seed: 1212, init_DR: "random", init_parameters: {} }, parameters);
        return this;
    }

    /**
     * initializes the projection.
     * @private
     */
    init() {
        if (this._is_initialized) return this;
        const N = this.X.shape[0];
        const { d, metric, init_DR: init_DR, init_parameters: DR_parameters } = this._parameters;
        if (init_DR === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (["PCA", "MDS"].includes(init_DR)) {
            this.Y = Matrix.from(init_DR == "PCA" ? PCA.transform(this.X, DR_parameters) : MDS.transform(this.X, DR_parameters));
        } else {
            throw new Error('init_DR needs to be either "random" or a DR method!');
        }
        this.distance_matrix = metric == "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        this._is_initialized = true;
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Matrix|Array} - The projection of {@link X}.
     */
    transform(max_iter = 200) {
        if (!this._is_initialized) this.init();
        for (let j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Generator} - A generator yielding the intermediate steps of the projection of {@link X}.
     */
    *generator(max_iter = 200) {
        if (!this._is_initialized) this.init();

        for (let j = 0; j < max_iter; ++j) {
            this._step();
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        const MAGIC = this.parameter("magic");
        const D = this.distance_matrix;
        const N = this.X.shape[0];
        const { d, metric } = this._parameters;
        let Y = this.Y;

        let G = new Matrix(N, d, 0);

        let sum = new Float64Array(d);
        for (let i = 0; i < N; ++i) {
            let e1 = new Float64Array(d);
            let e2 = new Float64Array(d);
            const Yi = Y.row(i);
            for (let j = 0; j < N; ++j) {
                if (i === j) continue;
                const Yj = Y.row(j);
                const delta = new Float64Array(d);
                for (let k = 0; k < d; ++k) {
                    delta[k] = Yi[k] - Yj[k];
                }
                const dY = metric(Yi, Yj);
                const dX = D.entry(i, j);
                const dq = dX - dY;
                const dr = Math.max(dX * dY, 1e-2);
                for (let k = 0; k < d; ++k) {
                    e1[k] += (delta[k] * dq) / dr;
                    e2[k] += (dq - (Math.pow(delta[k], 2) * (1 + dq / dY)) / dY) / dr;
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
}
