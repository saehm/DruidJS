import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias TSNE
 */
export class TSNE extends DR {
    static parameter_list = ["perplexity", "epsilon"];

    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [perplexity = 50] - perplexity.
     * @param {Number} [epsilon = 10] - learning parameter.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TSNE}
     */
    
    constructor(X, perplexity=50, epsilon=10, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = TSNE.parameter_list;
        [ this._N, this._D ] = this.X.shape;
        this.parameter("perplexity", Math.min(perplexity, this._N - 1));
        this.parameter("epsilon", epsilon);
        this._iter = 0;
        this.Y = new Matrix(this._N, this._d, () => this._randomizer.random);
        return this;
    }

    init(distance_matrix=null) {
        // init
        const Htarget = Math.log(this._perplexity);
        const N = this._N;
        const D = this._D;
        const metric = this._metric;
        const X = this.X;
        let Delta;
        if (distance_matrix) {
            Delta = distance_matrix;
        } else {
            Delta = new Matrix(N, N);
            for (let i = 0; i < N; ++i) {
                for (let j = i + 1; j < N; ++j) {
                    const distance = metric(X.row(i), X.row(j))
                    Delta.set_entry(i, j, distance);
                    Delta.set_entry(j, i, distance);
                }
            }

        } 
            
        const P = new Matrix(N, N, "zeros");

        this._ystep = new Matrix(N, D, "zeros");
        this._gains = new Matrix(N, D, 1);

        // search for fitting sigma
        let prow = new Array(N).fill(0);
        const tol = 1e-4;
        const maxtries = 50;
        for (let i = 0; i < N; ++i) {
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let done = false;

            let num = 0;
            while(!done) {
                let psum = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = Math.exp(-Delta.entry(i, j) * beta);
                    if (i === j) pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                let Hhere = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = (psum === 0) ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) {
                        Hhere -= pj * Math.log(pj);
                    }
                }
                if (Hhere > Htarget) {
                    betamin = beta;
                    beta = (betamax === Infinity) ? (beta * 2) : ((beta + betamax) / 2);
                } else {
                    betamax = beta;
                    beta = (betamin === -Infinity) ? (beta / 2) : ((beta + betamin) / 2);
                }
                ++num;
                if (Math.abs(Hhere - Htarget) < tol) done = true;
                if (num >= maxtries) done = true;
            }

            for (let j = 0; j < N; ++j) {
                P.set_entry(i, j, prow[j]);
            }
        }

        //compute probabilities
        const Pout = new Matrix(N, N, "zeros")
        const N2 = N * 2;
        for (let i = 0; i < N; ++i) {
            for (let j = i; j < N; ++j) {
                const p = Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100);
                Pout.set_entry(i, j, p);
                Pout.set_entry(j, i, p);
            }
        }
        this._P = Pout;
        return this;
    }

    transform(iterations=500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    * generator() {
        this.check_init();
        while (true) {
            this.next();
            yield {
                "projection": this.projection,
                "iteration": this._iter,
            }
        }
    }

    // perform optimization
    next() {
        const iter = ++this._iter;
        const P = this._P;
        const ystep = this._ystep;
        const gains = this._gains;
        const N = this._N;
        const epsilon = this._epsilon;
        const dim = this._d;
        let Y = this.Y;

        //calc cost gradient;
        const pmul = iter < 100 ? 4 : 1;
        
        // compute Q dist (unnormalized)
        const Qu = new Matrix(N, N, "zeros")
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
        const Q = new Matrix(N, N, 0)
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                const val = Math.max(Qu.entry(i, j) / qsum, 1e-100);
                Q.set_entry(i, j, val);
                Q.set_entry(j, i, val);
            }
        }

        //let cost = 0;
        //let grad = [];
        const grad = new Matrix(N, dim, "zeros");
        for (let i = 0; i < N; ++i) {
            //let gsum = new Float64Array(dim);//Array(dim).fill(0);
            for (let j = 0; j < N; ++j) {
                //cost += -P.entry(i, j) * Math.log(Q.entry(i, j));
                const premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (let d = 0; d < dim; ++d) {
                    //gsum[d] += premult * (Y.entry(i, d) - Y.entry(j, d));
                    grad.set_entry(i, d, grad.entry(i, d) + premult * (Y.entry(i, d) - Y.entry(j, d)));
                }
            }
            //grad.push(gsum);
        }

        // perform gradient step
        let ymean = new Float64Array(dim);
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const gid = grad.entry(i, d);
                const sid = ystep.entry(i, d);
                const gainid = gains.entry(i, d);
                
                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * .8 : gainid + .2;
                if (newgain < .01) newgain = .01;
                gains.set_entry(i, d, newgain);

                const momval = iter < 250 ? .5 : .8;
                const newsid = momval * sid - epsilon * newgain * gid;
                ystep.set_entry(i, d, newsid);

                Y.set_entry(i, d, Y.entry(i, d) + newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < 2; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - ymean[d] / N)
            }
        }

        return this.Y;
    }
} 