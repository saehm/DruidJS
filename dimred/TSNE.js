import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/randomizer";
/*import { simultaneous_poweriteration} from "../linear_algebra/index";
import { k_nearest_neighbors } from "../matrix/index";
import { neumair_sum } from "../numerical/index";
import { norm } from "../matrix/index";
import { linspace } from "../matrix/index";*/

export class TSNE{
    constructor(X, perplexity, epsilon, d=2, metric=euclidean, seed=1212) {
        this._X = X;
        this._d = d;
        [ this._N, this._D ] = X.shape;
        this._perplexity = perplexity;
        this._epsilon = epsilon;
        this._metric = metric;
        this._iter = 0;
        this.randomizer = new Randomizer(seed);
        this._Y = new Matrix(this._N, this._d, () => this.randomizer.random);
    }

    init(distance_matrix=null) {
        // init
        let Htarget = Math.log(this._perplexity);
        let D = distance_matrix || new Matrix(this._N, this._N, (i, j) => this._metric(this._X.row(i), this._X.row(j)));
        let P = new Matrix(this._N, this._N, "zeros");

        this._ystep = new Matrix(this._N, this._D, "zeros").to2dArray;
        this._gains = new Matrix(this._N, this._D, 1).to2dArray;

        // search for fitting sigma
        let prow = new Array(this._N).fill(0);
        for (let i = 0, N = this._N; i < N; ++i) {
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let done = false;
            let maxtries = 50;
            let tol = 1e-4;

            let num = 0;
            while(!done) {
                let psum = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = Math.exp(-D.entry(i, j) * beta);
                    if (i === j) pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                let Hhere = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = (psum === 0) ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) Hhere -= pj * Math.log(pj);
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
        let Pout = new Matrix(this._N, this._N, "zeros")
        let N2 = this._N * 2;
        for (let i = 0, N = this._N; i < N; ++i) {
            for (let j = 0; j < N; ++j) {
                Pout.set_entry(i, j, Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100));
            }
        }
        this._P = Pout;
        return this
    }

    set perplexity(value) {
        this._perplexity = value;
    }

    get perplexity() {
        return this._perplexity;
    }

    set epsilon(value) {
        this._epsilon = value;
    }

    get epsilon() {
        return this._epsilon;
    }

    transform(iterations=1000) {
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this._Y;
    }

    * transform_iter() {
        while (true) {
            this.next();
            yield this._Y;
        }
    }

    // perform optimization
    next() {
        let iter = ++this._iter;
        let P = this._P;
        let ystep = this._ystep;
        let gains = this._gains;
        let Y = this._Y;
        let N = this._N;
        let epsilon = this._epsilon;
        let dim = this._d;

        //calc cost gradient;
        let pmul = iter < 100 ? 4 : 1;
        
        // compute Q dist (unnormalized)
        let Qu = new Matrix(N, N, "zeros")
        let qsum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dsum = 0;
                for (let d = 0; d < dim; ++d) {
                    let dhere = Y.entry(i, d) - Y.entry(j, d);
                    dsum += dhere * dhere;
                }
                let qu = 1 / (1 + dsum);
                Qu.set_entry(i, j, qu);
                Qu.set_entry(j, i, qu);
                qsum += 2 * qu;
            }
        }

        // normalize Q dist
        let Q = new Matrix(N, N, (i, j) => Math.max(Qu.entry(i, j) / qsum, 1e-100));

        let cost = 0;
        let grad = [];
        for (let i = 0; i < N; ++i) {
            let gsum = new Array(dim).fill(0);
            for (let j = 0; j < N; ++j) {
                cost += -P.entry(i, j) * Math.log(Q.entry(i, j));
                let premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (let d = 0; d < dim; ++d) {
                    gsum[d] += premult * (Y.entry(i, d) - Y.entry(j, d));
                }
            }
            grad.push(gsum);
        }

        // perform gradient step
        let ymean = new Array(dim).fill(0);
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                let gid = grad[i][d];
                let sid = ystep[i][d];
                let gainid = gains[i][d];
                
                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * .8 : gainid + .2;
                if (newgain < .01) newgain = .01;
                gains[i][d] = newgain;

                let momval = iter < 250 ? .5 : .8;
                let newsid = momval * sid - epsilon * newgain * grad[i][d];
                ystep[i][d] = newsid;

                Y.set_entry(i, d, Y.entry(i, d) + newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < 2; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - ymean[d] / N)
            }
        }

        return this._Y;
    }

    get projection() {
        return this._Y;
    }
} 