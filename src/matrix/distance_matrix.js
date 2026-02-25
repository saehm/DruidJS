import { euclidean } from "../metrics/index.js";
import { Matrix } from "./index.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * @param {Matrix | Float64Array[] | number[][]} A
 * @returns {A is Matrix}
 */
function isMatrix(A) {
    return A instanceof Matrix;
}

/**
 * Computes the distance matrix of datamatrix `A`.
 *
 * @category Matrix
 * @param {Matrix | Float64Array[] | number[][]} A - Matrix.
 * @param {Metric} [metric=euclidean] - The diistance metric. Default is `euclidean`
 * @returns {Matrix} The distance matrix of `A`.
 */
export function distance_matrix(A, metric = euclidean) {
    /** @type {number} */
    const n = isMatrix(A) ? A.shape[0] : A.length;
    const D = new Matrix(n, n);
    for (let i = 0; i < n; ++i) {
        const A_i = isMatrix(A) ? A.row(i) : A[i];
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, isMatrix(A) ? A.row(j) : A[j]);
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}
