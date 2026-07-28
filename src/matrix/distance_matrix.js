import { euclidean } from "../metrics/index.js";
import { wasmDistanceMatrix } from "../wasm/index.js";
import { WASM_MIN_ROWS } from "../wasm/thresholds.js";
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
    const mat = isMatrix(A) ? A : Matrix.from(A);
    const [n, d] = mat.shape;
    const D = new Matrix(n, n);

    if (metric === euclidean && n >= WASM_MIN_ROWS) {
        if (wasmDistanceMatrix(mat.values, n, d, D.values)) {
            return D;
        }
    }

    for (let i = 0; i < n; ++i) {
        const A_i = mat.row(i);
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, mat.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}
