import { euclidean } from "../metrics/index.js";
import { Matrix } from "../matrix/index.js";
//import { neumair_sum } from "../numerical/index";

/**
 * Computes the norm of a vector, by computing its distance to **0**.
 * @memberof module:matrix
 * @alias norm
 * @param {Matrix|Array<Number>|Float64Array} v - Vector.
 * @param {Function} [metric = euclidean] - Which metric should be used to compute the norm.
 * @returns {Number} - The norm of {@link v}.
 */
export default function (v, metric = euclidean) {
    let vector = null;
    if (v instanceof Matrix) {
        let [rows, cols] = v.shape;
        if (rows === 1) vector = v.row(0);
        else if (cols === 1) vector = v.col(0);
        else throw new Error("Matrix must be 1d!");
    } else {
        vector = v;
    }
    const n = vector.length;
    const zeros = new Float64Array(n);
    return metric(vector, zeros);
}
