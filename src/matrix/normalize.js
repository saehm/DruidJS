import { euclidean } from "../metrics/index.js";
import { wasmNormalize } from "../wasm/index.js";
import { WASM_MIN_SIMPLE_VECTOR_LENGTH } from "../wasm/thresholds.js";
import { norm } from "./index.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * Normalizes Vector `v`.
 *
 * @category Matrix
 * @param {number[] | Float64Array} v - Vector
 * @param {Metric} metric
 * @returns {number[] | Float64Array} - The normalized vector with length 1.
 */
export function normalize(v, metric = euclidean) {
    if (v instanceof Float64Array && metric === euclidean && v.length >= WASM_MIN_SIMPLE_VECTOR_LENGTH) {
        const copy = new Float64Array(v);
        if (wasmNormalize(copy)) {
            return copy;
        }
    }
    const v_norm = norm(v, metric);
    return v.map((value) => value / v_norm);
}
