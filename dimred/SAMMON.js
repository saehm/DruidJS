import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { DR as DimRed} from "./DR.js";

/**
 * @class
 * @alias SAMMON
 */
export class SAMMON extends DimRed {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, magic=0.1, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed)
        super.parameter_list = ["magic"];
        this.parameter("magic", magic);
        [ this._N, this._D ] = this.X.shape;
        return this;
    }

    /**
     * initializes SAMMON. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init(DR="random", distance_matrix=null) {
        const N = this._N;
        const d = this._d;

        if (DR === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (DR instanceof DimRed) {
            this.Y = DR.transform(this.X);
        }
        this.distance_matrix = distance_matrix || this.__distance_matrix(this.X);
        return this;
    }

    /**
     * @private
     * @param {Matrix} A
     * @returns {Matrix} 
     */
    __distance_matrix(A) {
        const metric = this._metric;
        const N = A.shape[0];
        const D = new Matrix(N, N);
        for (let i = 0; i < N; ++i) {
            const A_i = A.row(i);
            for (let j = i; j < N; ++j) {
                let distance = (i === j ? 0 : metric(A_i, A.row(j)));
                D.set_entry(i, j, distance);
                D.set_entry(j, i, distance);
            }
        }
        return D;                
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     */
    transform(max_iter=200) {
        if (!this._is_initialized) this.init();
        for (let j = 0; j < max_iter; ++j) {
            this._step()
        }
        return this.projection;
    }

    * generator(max_iter=200) {
        if (!this._is_initialized) this.init();

        for (let j = 0; j < max_iter; ++j) {
            this._step()
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        const MAGIC = this.parameter("magic");
        const D = this.distance_matrix;
        const N = this._N;
        const d = this._d;
        const metric = this._metric;
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
                    delta[k] = Yi[k] - Yj[k]
                }
                const dY = metric(Yi, Yj);
                const dX = D.entry(i, j);
                const dq = dX - dY;
                const dr = Math.max(dX * dY, 1e-2);
                for (let k = 0; k < d; ++k) {
                    e1[k] += delta[k] * dq / dr;
                    e2[k] += (dq - Math.pow(delta[k], 2) * (1 + dq / dY) / dY) / dr;
                }
            }
            for (let k = 0; k < d; ++k) {
                const val = Y.entry(i, k) + (MAGIC * e1[k] / Math.abs(e2[k]) || 0);
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