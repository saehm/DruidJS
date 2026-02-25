import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * Computes the norm of a vector, by computing its distance to **0**.
 *
 * @category Matrix
 * @param {Matrix | number[] | Float64Array} v - Vector.
 * @param {Metric} [metric=euclidean] - Which metric should be used to compute the norm. Default is `euclidean`
 * @returns {number} - The norm of `v`.
 */
export function norm(v, metric = euclidean) {
    let vector = null;
    if (v instanceof Matrix) {
        const [rows, cols] = v.shape;
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
