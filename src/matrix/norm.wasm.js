/**
 * WASM wrapper for the `norm` kernel.
 *
 * The copy in, the call and the release are the same for every single-vector reduction, so the
 * wrapper is built by `vector_reduction` rather than written out. The shared runtime -- instance
 * setup, allocation and the persistent buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { vector_reduction } from "../wasm/index.js";

/**
 * Computes L2 Norm using WASM SIMD.
 *
 * Returns null if WASM is unavailable, so the caller falls back to its JS implementation.
 */
export const wasmNorm = vector_reduction("norm_f64");
