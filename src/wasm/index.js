import { WASM_BASE64 } from "./wasm_bytes.js";

/** @type {WebAssembly.Instance | null} */
let wasmInstance = null;

/**
 * Set once the module traps. An AssemblyScript `abort` unwinds straight back to the JS caller and
 * leaves the instance in an undefined state, so it must never be used again — every kernel then
 * reports failure and the caller falls back to its JS implementation.
 *
 * @type {boolean}
 */
let wasmAborted = false;

/**
 * Opt-out switch, see {@link setWasmEnabled}. Reads `DRUID_DISABLE_WASM` from the environment so
 * the JS paths can be exercised without changing any code.
 *
 * @type {boolean}
 */
let wasmEnabled = !read_env_flag("DRUID_DISABLE_WASM");

/**
 * Reads a boolean flag from `process.env`, tolerating environments without a `process` global.
 *
 * @param {string} name
 * @returns {boolean}
 */
function read_env_flag(name) {
    try {
        const env = /** @type {any} */ (globalThis).process?.env;
        const value = env?.[name];
        return value !== undefined && value !== "" && value !== "0" && value !== "false";
    } catch {
        return false;
    }
}

/**
 * Enables or disables WASM acceleration at runtime.
 *
 * With WASM disabled every accelerated function takes its JS fallback path. This exists so the two
 * implementations can be compared in tests and benchmarks.
 *
 * For most of the library the two paths agree to within floating point tolerance — the kernels
 * accumulate in a different order than the scalar loops, and `Math.pow`/`Math.exp` are only
 * "implementation-approximated" in ECMAScript, so the last bit may differ. `PCA`, `MDS`, `TSNE` and
 * `TriMap` come out bit-identical in practice.
 *
 * {@link UMAP} and {@link SAMMON} are the exception: their optimisers are chaotic, so a last-bit
 * difference grows into a visibly different (though equally valid) layout. See the note on
 * reproducibility in those classes.
 *
 * @param {boolean} enabled - Whether kernels may run in WASM.
 * @returns {boolean} The new state.
 * @example
 * import { setWasmEnabled } from "@saehrimnir/druidjs";
 * setWasmEnabled(false); // force the pure JS implementations
 */
export function setWasmEnabled(enabled) {
    wasmEnabled = Boolean(enabled);
    return wasmEnabled;
}

function initWasm() {
    if (!wasmEnabled || wasmAborted) return null;
    if (wasmInstance) return wasmInstance;
    try {
        if (typeof WebAssembly === "undefined") return null;
        let binaryString = "";
        const g = /** @type {any} */ (globalThis);
        if (typeof atob === "function") {
            binaryString = atob(WASM_BASE64);
        } else if (g?.Buffer) {
            binaryString = g.Buffer.from(WASM_BASE64, "base64").toString("binary");
        } else {
            return null;
        }
        const bytes = new Uint8Array(binaryString.length);
        for (let i = 0; i < binaryString.length; i++) {
            bytes[i] = binaryString.charCodeAt(i);
        }
        const module = new WebAssembly.Module(bytes);
        wasmInstance = new WebAssembly.Instance(module, {
            env: {
                abort: (
                    /** @type {number} */ _msg,
                    /** @type {number} */ _file,
                    /** @type {number} */ line,
                    /** @type {number} */ column,
                ) => {
                    // The instance cannot be trusted after a trap; retire it and fall back to JS.
                    wasmAborted = true;
                    wasmInstance = null;
                    console.error(`DruidJS: WASM kernel aborted at ${line}:${column}, falling back to JS.`);
                },
            },
        });
        return wasmInstance;
    } catch (_e) {
        return null;
    }
}

/**
 * Returns whether WASM acceleration is available in current environment.
 * @returns {boolean}
 */
export function isWasmAvailable() {
    return initWasm() !== null;
}

/**
 * Returns whether SharedArrayBuffer and WASM multi-threading is supported in current environment (Node.js, Bun, Deno, or Web Browser with COOP/COEP headers).
 * @returns {boolean}
 */
export function isWasmThreadsSupported() {
    try {
        if (typeof SharedArrayBuffer === "undefined") return false;
        const mem = new WebAssembly.Memory({ initial: 1, maximum: 1, shared: true });
        return mem.buffer instanceof SharedArrayBuffer;
    } catch {
        return false;
    }
}

/**
 * Allocates a buffer in the WASM heap and records the pointer in `ptrs`.
 *
 * Every kernel allocates through this and releases through {@link free_all} in a `finally`, so that
 * a failing allocation part-way through a list still frees the ones that already succeeded — WASM
 * linear memory is never reclaimed by the garbage collector, so a leaked pointer is permanent.
 *
 * @private
 * @param {any} exports - The WASM instance exports.
 * @param {number[]} ptrs - Accumulator holding every pointer allocated so far.
 * @param {number} size - Size in bytes.
 * @returns {number} Pointer to the allocated buffer.
 */
function alloc(exports, ptrs, size) {
    const ptr = exports.allocate(size);
    ptrs.push(ptr);
    return ptr;
}

/**
 * Frees every pointer collected by {@link alloc}.
 *
 * @private
 * @param {any} exports - The WASM instance exports.
 * @param {number[]} ptrs - Pointers to release.
 */
function free_all(exports, ptrs) {
    for (let i = 0; i < ptrs.length; ++i) {
        exports.free(ptrs[i]);
    }
    ptrs.length = 0;
}

/**
 * Computes C = A * B using WASM SIMD.
 * A is (rows_A x cols_A), B is (cols_A x cols_B), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmMatMul(A_val, rows_A, cols_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const sizeA = A_val.byteLength;
    const sizeB = B_val.byteLength;
    const sizeC = C_val.byteLength;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, sizeA);
        const ptrB = alloc(exports, ptrs, sizeB);
        const ptrC = alloc(exports, ptrs, sizeC);

        const memA = new Float64Array(memory.buffer, ptrA, A_val.length);
        memA.set(A_val);

        const memB = new Float64Array(memory.buffer, ptrB, B_val.length);
        memB.set(B_val);

        exports.matmul_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based Parallel Matrix Multiplication for worker thread execution.
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @param {number} start_row
 * @param {number} end_row
 * @returns {boolean}
 */
export function wasmMatMulRange(A_val, rows_A, cols_A, B_val, cols_B, C_val, start_row, end_row) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.matmul_range_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B, start_row, end_row);

        // Slice only, so concurrent workers sharing `C_val` do not clobber each other's rows.
        const offset = start_row * cols_B;
        const count = (end_row - start_row) * cols_B;
        C_val.set(new Float64Array(memory.buffer, ptrC + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Pairwise Euclidean Distance Matrix using WASM SIMD.
 *
 * @param {Float64Array} X_val
 * @param {number} n
 * @param {number} d
 * @param {Float64Array} D_val
 * @returns {boolean}
 */
export function wasmDistanceMatrix(X_val, n, d, D_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrD = alloc(exports, ptrs, D_val.byteLength);

        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.euclidean_distance_matrix_f64(ptrX, ptrD, n, d);

        D_val.set(new Float64Array(memory.buffer, ptrD, D_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based Pairwise Euclidean Distance Matrix for worker thread execution.
 *
 * @param {Float64Array} X_val
 * @param {number} n
 * @param {number} d
 * @param {Float64Array} D_val
 * @param {number} start_row
 * @param {number} end_row
 * @returns {boolean}
 */
export function wasmDistanceMatrixRange(X_val, n, d, D_val, start_row, end_row) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrD = alloc(exports, ptrs, D_val.byteLength);

        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.euclidean_distance_matrix_range_f64(ptrX, ptrD, n, d, start_row, end_row);

        // Slice only, so concurrent workers sharing `D_val` do not clobber each other's rows.
        const offset = start_row * n;
        const count = (end_row - start_row) * n;
        D_val.set(new Float64Array(memory.buffer, ptrD + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes cosine distance between two 1D vectors A and B.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null} Distance result, or null if WASM is unavailable.
 */
export function wasmCosineDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        const memA = new Float64Array(memory.buffer, ptrA, len);
        const memB = new Float64Array(memory.buffer, ptrB, len);

        memA.set(A);
        memB.set(B);

        return exports.cosine_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Manhattan (L1) distance between two 1D vectors A and B.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null} Distance result, or null if WASM is unavailable.
 */
export function wasmManhattanDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        const memA = new Float64Array(memory.buffer, ptrA, len);
        const memB = new Float64Array(memory.buffer, ptrB, len);

        memA.set(A);
        memB.set(B);

        return exports.manhattan_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Chebyshev (L_infinity) distance between two 1D vectors A and B.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null} Distance result, or null if WASM is unavailable.
 */
export function wasmChebyshevDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        const memA = new Float64Array(memory.buffer, ptrA, len);
        const memB = new Float64Array(memory.buffer, ptrB, len);

        memA.set(A);
        memB.set(B);

        return exports.chebyshev_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Assigns points to nearest centroids in k-Means.
 *
 * @param {Float64Array} X_val
 * @param {Float64Array} C_val
 * @param {Int32Array} assignments
 * @param {number} n
 * @param {number} k
 * @param {number} d
 * @returns {boolean} True if executed successfully via WASM.
 */
export function wasmKMeansAssign(X_val, C_val, assignments, n, k, d) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);
        const ptrA = alloc(exports, ptrs, assignments.byteLength);

        const memX = new Float64Array(memory.buffer, ptrX, X_val.length);
        const memC = new Float64Array(memory.buffer, ptrC, C_val.length);
        memX.set(X_val);
        memC.set(C_val);

        exports.kmeans_assign_f64(ptrX, ptrC, ptrA, n, k, d);

        const memA = new Int32Array(memory.buffer, ptrA, assignments.length);
        assignments.set(memA);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Performs a single t-SNE iteration gradient and position update step via WASM.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} P_val
 * @param {Float64Array} Qu_val
 * @param {Float64Array} Q_val
 * @param {Float64Array} grad_val
 * @param {Float64Array} ystep_val
 * @param {Float64Array} gains_val
 * @param {number} n
 * @param {number} dim
 * @param {number} pmul
 * @param {number} epsilon
 * @param {number} momval
 * @returns {boolean} True if executed successfully via WASM.
 */
export function wasmTsneStep(
    Y_val,
    P_val,
    Qu_val,
    Q_val,
    grad_val,
    ystep_val,
    gains_val,
    n,
    dim,
    pmul,
    epsilon,
    momval,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrP = alloc(exports, ptrs, P_val.byteLength);
        const ptrQu = alloc(exports, ptrs, Qu_val.byteLength);
        const ptrQ = alloc(exports, ptrs, Q_val.byteLength);
        const ptrGrad = alloc(exports, ptrs, grad_val.byteLength);
        const ptrYstep = alloc(exports, ptrs, ystep_val.byteLength);
        const ptrGains = alloc(exports, ptrs, gains_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrP, P_val.length).set(P_val);
        new Float64Array(memory.buffer, ptrYstep, ystep_val.length).set(ystep_val);
        new Float64Array(memory.buffer, ptrGains, gains_val.length).set(gains_val);

        exports.tsne_step_f64(ptrY, ptrP, ptrQu, ptrQ, ptrGrad, ptrYstep, ptrGains, n, dim, pmul, epsilon, momval);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        ystep_val.set(new Float64Array(memory.buffer, ptrYstep, ystep_val.length));
        gains_val.set(new Float64Array(memory.buffer, ptrGains, gains_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes C = A^T * B using WASM SIMD.
 * A is (cols_A x rows_A), B is (cols_A x cols_B), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} cols_A
 * @param {number} rows_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmTransDot(A_val, cols_A, rows_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.transDot_f64(ptrA, ptrB, ptrC, cols_A, rows_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes C = A * B^T using WASM SIMD.
 * A is (rows_A x cols_A), B is (cols_B x cols_A), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmDotTrans(A_val, rows_A, cols_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.dotTrans_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes outer product C = A x B using WASM SIMD.
 *
 * @param {Float64Array} A_val
 * @param {Float64Array} B_val
 * @param {Float64Array} C_val
 * @param {number} len
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmOuter(A_val, B_val, C_val, len) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.outer_f64(ptrA, ptrB, ptrC, len);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes L2 Norm using WASM SIMD.
 *
 * @param {Float64Array | number[]} V_val
 * @returns {number | null} L2 Norm, or null if WASM is unavailable.
 */
export function wasmNorm(V_val) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = V_val.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrV = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrV, len).set(V_val);
        return exports.norm_f64(ptrV, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Canberra distance.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null}
 */
export function wasmCanberraDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrA, len).set(A);
        new Float64Array(memory.buffer, ptrB, len).set(B);
        return exports.canberra_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Bray-Curtis distance.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null}
 */
export function wasmBrayCurtisDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrA, len).set(A);
        new Float64Array(memory.buffer, ptrB, len).set(B);
        return exports.bray_curtis_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * In-place Vector Normalization using WASM SIMD.
 *
 * @param {Float64Array} V_val
 * @returns {boolean} True if executed via WASM.
 */
export function wasmNormalize(V_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = V_val.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrV = alloc(exports, ptrs, V_val.byteLength);

        new Float64Array(memory.buffer, ptrV, len).set(V_val);
        exports.normalize_f64(ptrV, len);
        V_val.set(new Float64Array(memory.buffer, ptrV, len));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Inner Product <A, B> using WASM SIMD.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null}
 */
export function wasmInnerProduct(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrA, len).set(A);
        new Float64Array(memory.buffer, ptrB, len).set(B);
        return exports.inner_product_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Y = A * X using WASM SIMD.
 *
 * @param {Float64Array} A_val
 * @param {Float64Array} X_val
 * @param {Float64Array} Y_val
 * @param {number} rows
 * @param {number} cols
 * @returns {boolean}
 */
export function wasmMatVecMul(A_val, X_val, Y_val, rows, cols) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.mat_vec_mul_f64(ptrA, ptrX, ptrY, rows, cols);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Neumaier sum using WASM SIMD.
 *
 * @param {Float64Array | number[]} V_val
 * @returns {number | null}
 */
export function wasmNeumaierSum(V_val) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = V_val.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrV = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrV, len).set(V_val);
        return exports.neumair_sum_f64(ptrV, len);
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Sammon mapping single iteration gradient step using WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} D_val
 * @param {Float64Array} grad_val
 * @param {number} n
 * @param {number} dim
 * @param {number} magic
 * @returns {boolean}
 */
export function wasmSammonStep(Y_val, D_val, grad_val, n, dim, magic) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrD = alloc(exports, ptrs, D_val.byteLength);
        const ptrGrad = alloc(exports, ptrs, grad_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrD, D_val.length).set(D_val);

        exports.sammon_step_f64(ptrY, ptrD, ptrGrad, n, dim, magic);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Assigns points to nearest medoids in k-Medoids.
 *
 * @param {Float64Array} D_val
 * @param {Int32Array} medoids
 * @param {Int32Array} assignments
 * @param {number} n
 * @param {number} k
 * @returns {boolean}
 */
export function wasmKMedoidsAssign(D_val, medoids, assignments, n, k) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrD = alloc(exports, ptrs, D_val.byteLength);
        const ptrM = alloc(exports, ptrs, medoids.byteLength);
        const ptrA = alloc(exports, ptrs, assignments.byteLength);

        new Float64Array(memory.buffer, ptrD, D_val.length).set(D_val);
        new Int32Array(memory.buffer, ptrM, medoids.length).set(medoids);

        exports.kmedoids_assign_f64(ptrD, ptrM, ptrA, n, k);

        assignments.set(new Int32Array(memory.buffer, ptrA, assignments.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes All-Pairs Shortest Path geodesic distances for ISOMAP using WASM Min-Heap Repeated Dijkstra.
 *
 * @param {Int32Array} neighbors_val
 * @param {Float64Array} distances_val
 * @param {Float64Array} out_d_val
 * @param {number} n
 * @param {number} k
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmDijkstraAPSP(neighbors_val, distances_val, out_d_val, n, k) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrN = alloc(exports, ptrs, neighbors_val.byteLength);
        const ptrDist = alloc(exports, ptrs, distances_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_d_val.byteLength);

        new Int32Array(memory.buffer, ptrN, neighbors_val.length).set(neighbors_val);
        new Float64Array(memory.buffer, ptrDist, distances_val.length).set(distances_val);

        exports.dijkstra_apsp_f64(ptrN, ptrDist, ptrOut, n, k);

        out_d_val.set(new Float64Array(memory.buffer, ptrOut, out_d_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based All-Pairs Shortest Path calculation for worker thread execution.
 *
 * @param {Int32Array} neighbors_val
 * @param {Float64Array} distances_val
 * @param {Float64Array} out_d_val
 * @param {number} n
 * @param {number} k
 * @param {number} start_src
 * @param {number} end_src
 * @returns {boolean}
 */
export function wasmDijkstraAPSPRange(neighbors_val, distances_val, out_d_val, n, k, start_src, end_src) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrN = alloc(exports, ptrs, neighbors_val.byteLength);
        const ptrDist = alloc(exports, ptrs, distances_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_d_val.byteLength);

        new Int32Array(memory.buffer, ptrN, neighbors_val.length).set(neighbors_val);
        new Float64Array(memory.buffer, ptrDist, distances_val.length).set(distances_val);

        exports.dijkstra_apsp_range_f64(ptrN, ptrDist, ptrOut, n, k, start_src, end_src);

        // Only rows `[start_src, end_src)` were computed. Copying the whole buffer back would
        // overwrite the rows other workers wrote into the same shared output with this instance's
        // untouched (zero) memory.
        const offset = start_src * n;
        const count = (end_src - start_src) * n;
        out_d_val.set(new Float64Array(memory.buffer, ptrOut + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Executes a single SMACOF Guttman Transform iteration step in WASM.
 *
 * @param {Float64Array} Z_val
 * @param {Float64Array} TargetD_val
 * @param {Float64Array} ZNew_val
 * @param {number} n
 * @param {number} d
 * @returns {number | null} Computed stress value, or null if WASM is unavailable.
 */
export function wasmSmacofStep(Z_val, TargetD_val, ZNew_val, n, d) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrZ = alloc(exports, ptrs, Z_val.byteLength);
        const ptrTD = alloc(exports, ptrs, TargetD_val.byteLength);
        const ptrZNew = alloc(exports, ptrs, ZNew_val.byteLength);

        new Float64Array(memory.buffer, ptrZ, Z_val.length).set(Z_val);
        new Float64Array(memory.buffer, ptrTD, TargetD_val.length).set(TargetD_val);

        const stress = exports.smacof_step_f64(ptrZ, ptrTD, ptrZNew, n, d);

        ZNew_val.set(new Float64Array(memory.buffer, ptrZNew, ZNew_val.length));
        return stress;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes SQDMDS quartet gradients in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} X_val
 * @param {Uint32Array} quartets_val
 * @param {Float64Array} grads_val
 * @param {number} n
 * @param {number} d_ld
 * @param {number} d_hd
 * @param {boolean} use_exaggeration
 * @param {boolean} is_precomputed
 * @returns {boolean}
 */
export function wasmSqdmdsFillGrads(
    Y_val,
    X_val,
    quartets_val,
    grads_val,
    n,
    d_ld,
    d_hd,
    use_exaggeration,
    is_precomputed,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const num_quartets = quartets_val.length / 4;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrQ = alloc(exports, ptrs, quartets_val.byteLength);
        const ptrG = alloc(exports, ptrs, grads_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);
        new Int32Array(memory.buffer, ptrQ, quartets_val.length).set(quartets_val);

        exports.sqdmds_fill_grads_f64(
            ptrY,
            ptrX,
            ptrQ,
            num_quartets,
            ptrG,
            n,
            d_ld,
            d_hd,
            use_exaggeration,
            is_precomputed,
        );

        grads_val.set(new Float64Array(memory.buffer, ptrG, grads_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based SQDMDS quartet gradient calculation for multi-threaded worker execution.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} X_val
 * @param {Uint32Array} quartets_val
 * @param {Float64Array} grads_val
 * @param {number} n
 * @param {number} d_ld
 * @param {number} d_hd
 * @param {boolean} use_exaggeration
 * @param {boolean} is_precomputed
 * @param {number} start_q
 * @param {number} end_q
 * @returns {boolean}
 */
export function wasmSqdmdsFillGradsRange(
    Y_val,
    X_val,
    quartets_val,
    grads_val,
    n,
    d_ld,
    d_hd,
    use_exaggeration,
    is_precomputed,
    start_q,
    end_q,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const num_quartets = quartets_val.length / 4;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrQ = alloc(exports, ptrs, quartets_val.byteLength);
        const ptrG = alloc(exports, ptrs, grads_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);
        new Int32Array(memory.buffer, ptrQ, quartets_val.length).set(quartets_val);

        exports.sqdmds_fill_grads_range_f64(
            ptrY,
            ptrX,
            ptrQ,
            num_quartets,
            ptrG,
            n,
            d_ld,
            d_hd,
            use_exaggeration,
            is_precomputed,
            start_q,
            end_q,
        );

        grads_val.set(new Float64Array(memory.buffer, ptrG, grads_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Applies SQDMDS Nesterov momentum update and updates embedding Y in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} M_val
 * @param {Float64Array} Grads_val
 * @param {number} n
 * @param {number} d
 * @param {number} lr
 * @returns {boolean}
 */
export function wasmSqdmdsNestrovStep(Y_val, M_val, Grads_val, n, d, lr) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrM = alloc(exports, ptrs, M_val.byteLength);
        const ptrG = alloc(exports, ptrs, Grads_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrM, M_val.length).set(M_val);
        new Float64Array(memory.buffer, ptrG, Grads_val.length).set(Grads_val);

        exports.sqdmds_nestrov_step_f64(ptrY, ptrM, ptrG, n, d, lr);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        M_val.set(new Float64Array(memory.buffer, ptrM, M_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes TriMap triplet gradients in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Int32Array} triplets_val
 * @param {Float64Array} weights_val
 * @param {Float64Array} grad_val
 * @param {number} n
 * @param {number} dim
 * @param {number} n_inliers
 * @param {number} n_outliers
 * @returns {number | null} Computed loss, or null if WASM is unavailable.
 */
export function wasmTriMapGrad(Y_val, triplets_val, weights_val, grad_val, n, dim, n_inliers, n_outliers) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const num_triplets = triplets_val.length / 3;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrT = alloc(exports, ptrs, triplets_val.byteLength);
        const ptrW = alloc(exports, ptrs, weights_val.byteLength);
        const ptrG = alloc(exports, ptrs, grad_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Int32Array(memory.buffer, ptrT, triplets_val.length).set(triplets_val);
        new Float64Array(memory.buffer, ptrW, weights_val.length).set(weights_val);

        const loss = exports.trimap_grad_f64(ptrY, ptrT, ptrW, ptrG, n, dim, num_triplets, n_inliers, n_outliers);

        grad_val.set(new Float64Array(memory.buffer, ptrG, grad_val.length));
        return loss;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Applies TriMap embedding update in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} grad_val
 * @param {Float64Array} vel_val
 * @param {Float64Array} gain_val
 * @param {number} n
 * @param {number} dim
 * @param {number} gamma
 * @param {number} lr
 * @returns {boolean}
 */
export function wasmTriMapUpdate(Y_val, grad_val, vel_val, gain_val, n, dim, gamma, lr) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrG = alloc(exports, ptrs, grad_val.byteLength);
        const ptrV = alloc(exports, ptrs, vel_val.byteLength);
        const ptrGain = alloc(exports, ptrs, gain_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrG, grad_val.length).set(grad_val);
        new Float64Array(memory.buffer, ptrV, vel_val.length).set(vel_val);
        new Float64Array(memory.buffer, ptrGain, gain_val.length).set(gain_val);

        exports.trimap_update_f64(ptrY, ptrG, ptrV, ptrGain, n, dim, gamma, lr);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        vel_val.set(new Float64Array(memory.buffer, ptrV, vel_val.length));
        gain_val.set(new Float64Array(memory.buffer, ptrGain, gain_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes a single MeanShift iteration step in WASM SIMD.
 *
 * @param {Float64Array} points_val
 * @param {Float64Array} out_points_val
 * @param {number} n
 * @param {number} d
 * @param {number} bandwidth
 * @param {boolean} use_gaussian
 * @returns {number | null} Max shift magnitude, or null if WASM is unavailable.
 */
export function wasmMeanShiftStep(points_val, out_points_val, n, d, bandwidth, use_gaussian) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrIn = alloc(exports, ptrs, points_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_points_val.byteLength);

        new Float64Array(memory.buffer, ptrIn, points_val.length).set(points_val);

        const max_shift = exports.meanshift_step_f64(ptrIn, ptrOut, n, d, bandwidth, use_gaussian);

        out_points_val.set(new Float64Array(memory.buffer, ptrOut, out_points_val.length));
        return max_shift;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based MeanShift iteration step for multi-threaded worker execution.
 *
 * @param {Float64Array} points_val
 * @param {Float64Array} out_points_val
 * @param {number} n
 * @param {number} d
 * @param {number} bandwidth
 * @param {boolean} use_gaussian
 * @param {number} start_i
 * @param {number} end_i
 * @returns {number | null}
 */
export function wasmMeanShiftStepRange(points_val, out_points_val, n, d, bandwidth, use_gaussian, start_i, end_i) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrIn = alloc(exports, ptrs, points_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_points_val.byteLength);

        new Float64Array(memory.buffer, ptrIn, points_val.length).set(points_val);

        const max_shift = exports.meanshift_step_range_f64(
            ptrIn,
            ptrOut,
            n,
            d,
            bandwidth,
            use_gaussian,
            start_i,
            end_i,
        );

        // Slice only, so concurrent workers sharing `out_points_val` do not clobber each other.
        const offset = start_i * d;
        const count = (end_i - start_i) * d;
        out_points_val.set(new Float64Array(memory.buffer, ptrOut + offset * 8, count), offset);
        return max_shift;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Runs one UMAP SGD epoch: attractive forces over the active edges plus repulsive forces over their
 * negative samples.
 *
 * `negative_samples` holds the node indices drawn by the caller's seeded randomizer, in edge order —
 * the RNG stays in JS so both paths consume the identical sequence. `embedding`, `epoch_of_next_sample`
 * and `epoch_of_next_negative_sample` are updated in place.
 *
 * @param {Float64Array} embedding - The `n_points ⨯ dim` embedding, mutated in place.
 * @param {Int32Array} head - Edge source indices.
 * @param {Int32Array} tail - Edge target indices.
 * @param {Float32Array} epochs_per_sample
 * @param {Float32Array} epoch_of_next_sample - Mutated in place.
 * @param {Float32Array} epochs_per_negative_sample
 * @param {Float32Array} epoch_of_next_negative_sample - Mutated in place.
 * @param {Int32Array} negative_samples - Pre-drawn negative sample node indices.
 * @param {number} dim
 * @param {number} iter - Current epoch.
 * @param {number} a
 * @param {number} b
 * @param {number} gamma - Repulsion strength.
 * @param {number} alpha - Current learning rate.
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmUmapOptimizeEpoch(
    embedding,
    head,
    tail,
    epochs_per_sample,
    epoch_of_next_sample,
    epochs_per_negative_sample,
    epoch_of_next_negative_sample,
    negative_samples,
    dim,
    iter,
    a,
    b,
    gamma,
    alpha,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, embedding.byteLength);
        const ptrHead = alloc(exports, ptrs, head.byteLength);
        const ptrTail = alloc(exports, ptrs, tail.byteLength);
        const ptrEps = alloc(exports, ptrs, epochs_per_sample.byteLength);
        const ptrEons = alloc(exports, ptrs, epoch_of_next_sample.byteLength);
        const ptrEpns = alloc(exports, ptrs, epochs_per_negative_sample.byteLength);
        const ptrEonns = alloc(exports, ptrs, epoch_of_next_negative_sample.byteLength);
        // `allocate(0)` is not meaningful, so keep a one-element floor for empty epochs.
        const ptrNeg = alloc(exports, ptrs, Math.max(negative_samples.byteLength, 4));

        new Float64Array(memory.buffer, ptrY, embedding.length).set(embedding);
        new Int32Array(memory.buffer, ptrHead, head.length).set(head);
        new Int32Array(memory.buffer, ptrTail, tail.length).set(tail);
        new Float32Array(memory.buffer, ptrEps, epochs_per_sample.length).set(epochs_per_sample);
        new Float32Array(memory.buffer, ptrEons, epoch_of_next_sample.length).set(epoch_of_next_sample);
        new Float32Array(memory.buffer, ptrEpns, epochs_per_negative_sample.length).set(epochs_per_negative_sample);
        new Float32Array(memory.buffer, ptrEonns, epoch_of_next_negative_sample.length).set(
            epoch_of_next_negative_sample,
        );
        if (negative_samples.length > 0) {
            new Int32Array(memory.buffer, ptrNeg, negative_samples.length).set(negative_samples);
        }

        exports.umap_optimize_epoch_f64(
            ptrY,
            ptrHead,
            ptrTail,
            ptrEps,
            ptrEons,
            ptrEpns,
            ptrEonns,
            ptrNeg,
            head.length,
            negative_samples.length,
            dim,
            iter,
            a,
            b,
            gamma,
            alpha,
        );

        embedding.set(new Float64Array(memory.buffer, ptrY, embedding.length));
        epoch_of_next_sample.set(new Float32Array(memory.buffer, ptrEons, epoch_of_next_sample.length));
        epoch_of_next_negative_sample.set(
            new Float32Array(memory.buffer, ptrEonns, epoch_of_next_negative_sample.length),
        );
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
