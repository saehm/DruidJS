/**
 * Problem-size thresholds below which the plain JS implementation is used instead of a WASM kernel.
 *
 * Crossing into WASM is not free: arguments are copied into linear memory and results copied back,
 * and a JIT-warm JS loop over a `Float64Array` is already fast. For small inputs the copy dominates
 * and the JS loop wins outright — a threshold set too low makes the library *slower*.
 *
 * The constants below are break-even points measured with `benchmark/wasm_threshold_calibration.js`
 * (Node 24, x64). Re-run it if the kernels change. Values are rounded up to a power of two so they
 * stay conservative across engines: too high only forgoes a modest speedup, too low costs real time
 * on every call.
 *
 * @module wasm/thresholds
 */

/**
 * Minimum vector length for element-wise kernels that accumulate several running sums per element,
 * such as {@link cosine} (dot product plus both norms).
 *
 * More arithmetic per byte copied means the crossover comes earlier; measured break-even is ~512
 * elements, growing to roughly a 2× speedup by 16k.
 *
 * @type {number}
 */
export const WASM_MIN_VECTOR_LENGTH = 512;

/**
 * Minimum vector length for single-accumulator reductions — `manhattan`, `chebyshev`, `canberra`,
 * `bray_curtis`, `norm`, `normalize`, `inner_product`.
 *
 * These do one multiply-add per element, so the copy is only amortised around ~1024 elements, and
 * even then the gain is modest (~1.3× at 16k).
 *
 * @type {number}
 */
export const WASM_MIN_SIMPLE_VECTOR_LENGTH = 1024;

/**
 * Minimum number of rows for the row-wise kernels that touch every pair of points
 * (`distance_matrix`, k-means assignment, the t-SNE/SMACOF/SAMMON iteration steps).
 *
 * These do O(N²) work for a single O(N·D) copy, so the crossover sits low — measured at ~8 rows.
 *
 * @type {number}
 */
export const WASM_MIN_ROWS = 16;

/**
 * Minimum number of scalar multiply-accumulates (`rows_A * cols_A * cols_B`) for the matrix-product
 * kernels (`dot`, `transDot`, `dotTrans`, `outer`).
 *
 * Measured break-even is ~512 operations (an 8 ⨯ 8 product); beyond it the kernel pulls away fast,
 * reaching ~4× at 32 ⨯ 32.
 *
 * @type {number}
 */
export const WASM_MIN_MATMUL_OPS = 512;

/**
 * Minimum row count before a range kernel is split across the worker pool.
 *
 * Starting the pool costs roughly 100 ms — every worker instantiates its own copy of the module —
 * and a method that dispatches once, such as ISOMAP's geodesic step, cannot amortize that. Measured
 * end to end including start-up, parallel Dijkstra is 0.37× at 200 rows, breaks even near 500, and
 * reaches 1.76× at 1024, 3.99× at 2000 and 5.49× at 4000. Once the pool is warm the speedup is
 * 5-6× from a few hundred rows upward, so this bound is set by start-up, not by the kernel.
 *
 * @type {number}
 */
export const WASM_MIN_PARALLEL_ROWS = 1024;

/**
 * Minimum edge count for the UMAP SGD epoch kernel.
 *
 * Each epoch copies the embedding and four edge-length arrays in and three back out, so tiny graphs
 * are dominated by the transfer.
 *
 * @type {number}
 */
export const WASM_MIN_UMAP_EDGES = 512;
