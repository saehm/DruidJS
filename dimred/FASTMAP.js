import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/randomizer";
/**
 * @class
 * @alias FASTMAP
 */
export class FASTMAP{
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @returns {FASTMAP}
     */
    constructor(X, d=2, metric=euclidean) {
        this.X = X;
        this.d = d;
        this._metric = metric;
        this._col = -1
        this.randomizer = new Randomizer(1212)
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {function} dist 
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */
    _choose_distant_objects(dist) {
        let X = this.X;
        let N = X.shape[0];
        let a_index = this.randomizer.random_int % N - 1;
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            let d_ai = dist(a_index, i)
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            let d_bi = dist(b_index, i)
            if (d_bi > max_dist) {
                max_dist = d_bi;
                a_index = i;
            }
        }
        return [a_index, b_index, max_dist];
    }

    /**
     * Computes the projection.
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */
    transform() {
        let X = this.X;
        let [ rows, D ] = X.shape;
        let Y = new Matrix(rows, this.d);
        //let PA = [[], []];
        let dist = (a,b) => this._metric(X.row(a), X.row(b));
        let old_dist = dist;

        while(this._col < this.d - 1) {
            this._col += 1;
            let col = this._col;
            // choose pivot objects
            let [a_index, b_index, d_ab] = this._choose_distant_objects(dist);
            // record id of pivot objects
            //PA[0].push(a_index);
            //PA[1].push(b_index);
            if (d_ab === 0) {
                // because all inter-object distances are zeros
                for (let i = 0; i < rows; ++i) {
                    Y.set_entry(i, col, 0);
                }
            } else {
                // project the objects on the line (O_a, O_b)
                for (let i = 0; i < rows; ++i) {
                    let d_ai = dist(a_index, i);
                    let d_bi = dist(b_index, i);
                    let y_i = (d_ai ** 2 + d_ab ** 2 - d_bi ** 2) / (2 * d_ab);
                    Y.set_entry(i, col, y_i);
                }
                // consider the projections of the objects on a
                // hyperplane perpendicluar to the line (a, b);
                // the distance function D'() between two 
                // projections is given by Eq.4
                dist = (a,b) => Math.sqrt((old_dist(a,b) ** 2) - ((Y.entry(a, col) - Y.entry(b, col)) ** 2))
            }
        }
        // return embedding
        this.Y = Y
        return this.Y;
    }

    /**
     * @returns {Matrix}
     */
    get projection() {
        return this.Y
    }

    

}