import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { DR } from "./DR.js";
import { PCA } from "./index.js";

export class SAMMON extends DR {
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
    constructor(X, max_halves=5, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed)
        super.parameter_list = ["max_halves"];
        this.parameter("max_halves", max_halves);
        [ this._N, this._D ] = this.X.shape;
        return this;
    }

    /**
     * initializes SAMMON. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init(DR="random", Distance_matrix=null) {
        const N = this._N;
        const d = this._d;

        if (DR === "random") {
            const randomizer = this._randomizer;
            console.log(randomizer)
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else {
            this.Y = DR.transform(this.X);
        }
        const Y = this.Y;

        if (!Distance_matrix) {
            Distance_matrix = new Matrix(N, N);
        }

        const metric = this._metric;
        let distance_matrix = new Matrix(N, N);
        let distance_inverse_matrix = new Matrix(N, N);
        for (let i = 0; i < N; ++i) {
            const Y_i = Y.row(i);
            for (let j = i; j < N; ++j) {
                let distance = i === j ? 1 : metric(Y_i, Y.row(j));
                let distance_inverse = 1 / distance;
                distance_matrix.set_entry(i, j, distance);
                distance_matrix.set_entry(j, i, distance);
                distance_inverse_matrix.set_entry(i, j, distance_inverse);
                distance_inverse_matrix.set_entry(j, i, distance_inverse);
                if (!Distance_matrix) {
                    let Distance = i === j ? 1 : metric(X.row(i), X.row(j));
                    Distance_matrix.set_entry(i, j, Distance);
                    Distance_matrix.set_entry(j, i, Distance);
                }
            }
        }
        let Distance_inverse_matrix = Distance_matrix._apply(1, (d, v) => v / d);
        let delta = Distance_matrix.sub(distance_matrix);
        let E = delta._apply(2, (d, v) => Math.pow(d, v)).mult(Distance_inverse_matrix);
        console.log(E)
        this._distance_matrix = distance_matrix;
        this._distance_inverse_matrix = distance_inverse_matrix;
        this._Distance_matrix = Distance_matrix;
        this._Distance_inverse_matrix = distance_inverse_matrix;
        this._delta = delta;
        this._ones = new Matrix(N, d, 1);
        this._E = E.sum;

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
                let distance = (i === j ? 1 : metric(A_i, A.row(j)));
                D.set_entry(i, j, distance);
                D.set_entry(j, i, distance);
            }
        }
        return D;                
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     */
    transform(max_iter=20) {
        if (!this._is_initialized) this.init();

        for (let j = 0; j < max_iter; ++j) {
            console.log([...this.Y])
            this._step()
        }

        return this.projection;
    }

    * generator() {
        if (!this._is_initialized) this.init();

        for (let j = 0; j < max_iter; ++j) {
            this._step()
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        const max_halves = this.parameter("max_halves");

        let distance_matrix = this._distance_matrix;
        let distance_inverse_matrix = this._distance_inverse_matrix;
        let Distance_matrix = this._Distance_matrix;
        let Distance_inverse_matrix = this._Distance_inverse_matrix;
        let ones = this._ones;
        let E = this._E;
        let Y = this.Y;

        let delta = distance_inverse_matrix.sub(Distance_inverse_matrix);
        let delta_one = delta.dot(ones);
        let g = delta.dot(Y).sub(Y.mult(delta_one));
        let dinv3 = distance_inverse_matrix._apply(3, (d, v) => Math.pow(d, v));
        let Y2 = Y._apply(2, (d, v) => Math.pow(d, v));
        let H = dinv3.dot(Y2).sub(delta_one).sub(Y.mult(2).mult(dinv3.dot(Y))).add(Y2.mult(dinv3.dot(ones)));
        H = H._apply(null, (d) => Math.abs(d));
        let s = g.divide(H);
        let Y_old = Y.clone();

        for (let j = 0; j < max_halves; ++j) {
            Y = Y_old.add(s);
            distance_matrix = this.__distance_matrix(Y);
            distance_inverse_matrix = distance_matrix._apply(1, (d, v) => v / d);
            delta = Distance_matrix.sub(distance_matrix);
            let E_new = delta._apply(2, (d, v) => Math.pow(d, v)).mult(Distance_inverse_matrix).sum;
            if (E_new < E) {
                break;
            } else {
                s = s.mult(.5);
            }
        }

        this.Y = Y;
    }
} 