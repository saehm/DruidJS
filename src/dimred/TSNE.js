import { Matrix } from "../matrix/index.js";
import { euclidean_squared } from "../metrics/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias TSNE
 * @extends DR
 */
export class TSNE extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.perplexity = 50] - perplexity.
     * @param {number} [parameters.epsilon = 10] - learning parameter.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function|"precomputed"} [parameters.metric = euclidean_squared] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TSNE}
     */
    constructor(X, parameters) {
        super(X, { perplexity: 50, epsilon: 10, d: 2, metric: euclidean_squared, seed: 1212 }, parameters);
        [this._N, this._D] = this.X.shape;
        this._iter = 0;
        this.Y = new Matrix(this._N, this.parameter("d"), () => this._randomizer.gauss_random() * 1e-4);
        return this;
    }

    /**
     *
     * @returns {TSNE}
     */
    init() {
        // init
        const Htarget = Math.log(this.parameter("perplexity"));
        const N = this._N;
        const D = this._D;
        const { metric } = this._parameters;
        const X = this.X;
        let Delta;
        if (metric == "precomputed") {
            Delta = druid.Matrix.from(X);
        } else {
            Delta = new Matrix(N, N);
            for (let i = 0; i < N; ++i) {
                const X_i = X.row(i);
                for (let j = i + 1; j < N; ++j) {
                    const distance = metric(X_i, X.row(j));
                    Delta.set_entry(i, j, distance);
                    Delta.set_entry(j, i, distance);
                }
            }
        }

        const P = new Matrix(N, N, 0);

        this._ystep = new Matrix(N, D, 0);
        this._gains = new Matrix(N, D, 1);

        // search for fitting sigma
        const tol = 1e-4;
        const maxtries = 50;
        for (let i = 0; i < N; ++i) {
            const dist_i = Delta.row(i);
            const prow = P.row(i);
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let cnt = maxtries;
            let done = false;
            let psum;

            while (!done && cnt--) {
                // compute entropy and kernel row with beta precision
                psum = 0;
                let dp_sum = 0;
                for (let j = 0; j < N; ++j) {
                    const dist = dist_i[j];
                    const pj = i !== j ? Math.exp(-dist * beta) : 0;
                    dp_sum += dist * pj;
                    prow[j] = pj;
                    psum += pj;
                }
                // compute entropy
                const H = psum > 0 ? Math.log(psum) + (beta * dp_sum) / psum : 0;
                if (H > Htarget) {
                    betamin = beta;
                    beta = betamax === Infinity ? beta * 2 : (beta + betamax) / 2;
                } else {
                    betamax = beta;
                    beta = betamin === -Infinity ? beta / 2 : (beta + betamin) / 2;
                }
                done = Math.abs(H - Htarget) < tol;
            }
            // normalize p
            for (let j = 0; j < N; ++j) {
                prow[j] /= psum;
            }
        }

        // compute probabilities
        const N2 = N * 2;
        for (let i = 0; i < N; ++i) {
            for (let j = i; j < N; ++j) {
                const p = Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100);
                P.set_entry(i, j, p);
                P.set_entry(j, i, p);
            }
        }
        this._P = P;
        return this;
    }

    /**
     *
     * @param {number} [iterations=500] - number of iterations.
     * @returns {Matrix|number[][]} the projection.
     */
    transform(iterations = 500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     *
     * @param {number} [iterations=500] - number of iterations.
     * @yields {Matrix|number[][]} - the projection.
     */
    *generator(iterations = 500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */
    next() {
        const iter = ++this._iter;
        const P = this._P;
        const ystep = this._ystep;
        const gains = this._gains;
        const N = this._N;
        const { d: dim, epsilon } = this._parameters;
        let Y = this.Y;

        //calc cost gradient;
        const pmul = iter < 100 ? 4 : 1;

        // compute Q dist (unnormalized)
        const Qu = new Matrix(N, N, "zeros");
        let qsum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dsum = 0;
                for (let d = 0; d < dim; ++d) {
                    const dhere = Y.entry(i, d) - Y.entry(j, d);
                    dsum += dhere * dhere;
                }
                const qu = 1 / (1 + dsum);
                Qu.set_entry(i, j, qu);
                Qu.set_entry(j, i, qu);
                qsum += 2 * qu;
            }
        }

        // normalize Q dist
        const Q = new Matrix(N, N, 0);
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                const val = Math.max(Qu.entry(i, j) / qsum, 1e-100);
                Q.set_entry(i, j, val);
                Q.set_entry(j, i, val);
            }
        }

        const grad = new Matrix(N, dim, "zeros");
        for (let i = 0; i < N; ++i) {
            for (let j = 0; j < N; ++j) {
                const premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (let d = 0; d < dim; ++d) {
                    grad.add_entry(i, d, premult * (Y.entry(i, d) - Y.entry(j, d)));
                }
            }
        }

        // perform gradient step
        let ymean = new Float64Array(dim);
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const gid = grad.entry(i, d);
                const sid = ystep.entry(i, d);
                const gainid = gains.entry(i, d);

                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * 0.8 : gainid + 0.2;
                if (newgain < 0.01) newgain = 0.01;
                gains.set_entry(i, d, newgain);

                const momval = iter < 250 ? 0.5 : 0.8;
                const newsid = momval * sid - epsilon * newgain * gid;
                ystep.set_entry(i, d, newsid);

                Y.add_entry(i, d, newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                Y.sub_entry(i, d, ymean[d] / N);
            }
        }

        return this.Y;
    }
}
