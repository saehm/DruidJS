/**
 * WASM wrapper for the `canberra` kernel.
 *
 * The copy in, the call and the release are the same for every vector metric, so the wrapper is
 * built by `vector_pair_reduction` rather than written out. The shared runtime -- instance setup,
 * allocation and the persistent buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { vector_pair_reduction } from "../wasm/index.js";

/**
 * Computes Canberra distance between two 1D vectors A and B.
 *
 * Returns null if WASM is unavailable, so the caller falls back to its JS implementation.
 */
export const wasmCanberraDistance = vector_pair_reduction("canberra_distance_f64");
