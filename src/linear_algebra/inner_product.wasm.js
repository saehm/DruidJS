/**
 * WASM wrapper for the `inner_product` kernel.
 *
 * The copy in, the call and the release are the same for every two-vector reduction, so the wrapper
 * is built by `vector_pair_reduction` rather than written out. The shared runtime -- instance setup,
 * allocation and the persistent buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { vector_pair_reduction } from "../wasm/index.js";

/**
 * Computes Inner Product <A, B> using WASM SIMD.
 *
 * Returns null if WASM is unavailable, so the caller falls back to its JS implementation.
 */
export const wasmInnerProduct = vector_pair_reduction("inner_product_f64");
