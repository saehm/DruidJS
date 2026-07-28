/**
 * Measures where each class of WASM kernel overtakes its JS fallback.
 *
 * The numbers in `src/wasm/thresholds.js` come from this script — re-run it after changing a kernel
 * and update the constants if a crossover has moved:
 *
 *     node benchmark/wasm_threshold_calibration.js
 *
 * Each row reports the mean time per call for the JS loop and for the kernel. The threshold to pick
 * is the first size where the kernel wins by a comfortable margin, rounded up to a power of two.
 */

import { euclidean } from "../src/metrics/index.js";
import { wasmMatMul } from "../src/matrix/Matrix.wasm.js";
import { wasmDistanceMatrix } from "../src/matrix/distance_matrix.wasm.js";
import { wasmCosineDistance } from "../src/metrics/cosine.wasm.js";
import { wasmManhattanDistance } from "../src/metrics/manhattan.wasm.js";
import { wasmNeumaierSum } from "../src/numerical/neumair_sum.wasm.js";
import { isWasmAvailable } from "../src/wasm/index.js";

if (!isWasmAvailable()) {
    console.error("WASM is unavailable in this environment; nothing to calibrate.");
    process.exit(1);
}

/**
 * @param {() => unknown} fn
 * @param {number} reps
 * @returns {number} Mean microseconds per call.
 */
function bench(fn, reps) {
    fn(); // warm up the JIT before timing
    const start = performance.now();
    for (let i = 0; i < reps; ++i) fn();
    return ((performance.now() - start) / reps) * 1e3;
}

/**
 * @param {number} n
 * @returns {Float64Array}
 */
function random_vector(n) {
    const a = new Float64Array(n);
    for (let i = 0; i < n; ++i) a[i] = Math.random();
    return a;
}

/**
 * @param {string} title
 * @param {string} unit
 * @param {number[]} sizes
 * @param {(n: number) => { js: () => unknown; wasm: () => unknown; reps: number }} setup
 */
function table(title, unit, sizes, setup) {
    console.log(`\n${title}`);
    console.log(`${unit.padStart(10)} | ${"js (µs)".padStart(9)} | ${"wasm (µs)".padStart(9)} | speedup`);
    for (const n of sizes) {
        const { js, wasm, reps } = setup(n);
        const t_js = bench(js, reps);
        const t_wasm = bench(wasm, reps);
        const speedup = t_js / t_wasm;
        const marker = speedup > 1 ? "" : "  (JS wins)";
        console.log(
            `${String(n).padStart(10)} | ${t_js.toFixed(2).padStart(9)} | ${t_wasm.toFixed(2).padStart(9)} | ${speedup.toFixed(2)}×${marker}`,
        );
    }
}

table("Multi-accumulator vector kernel — WASM_MIN_VECTOR_LENGTH", "length", [128, 256, 512, 1024, 4096, 16384], (n) => {
    const a = random_vector(n);
    const b = random_vector(n);
    return {
        reps: 5000,
        js: () => {
            let sum = 0;
            let sum_a = 0;
            let sum_b = 0;
            for (let i = 0; i < n; ++i) {
                sum += a[i] * b[i];
                sum_a += a[i] * a[i];
                sum_b += b[i] * b[i];
            }
            return Math.acos(sum / (Math.sqrt(sum_a) * Math.sqrt(sum_b)));
        },
        wasm: () => wasmCosineDistance(a, b),
    };
});

table(
    "Single-accumulator vector kernel — WASM_MIN_SIMPLE_VECTOR_LENGTH",
    "length",
    [128, 256, 512, 1024, 4096, 16384],
    (n) => {
        const a = random_vector(n);
        const b = random_vector(n);
        return {
            reps: 5000,
            js: () => {
                let sum = 0;
                for (let i = 0; i < n; ++i) sum += Math.abs(a[i] - b[i]);
                return sum;
            },
            wasm: () => wasmManhattanDistance(a, b),
        };
    },
);

table("Compensated summation — kept in JS, see neumair_sum", "length", [128, 512, 4096, 16384], (n) => {
    const v = random_vector(n);
    return {
        reps: 5000,
        js: () => {
            let sum = 0;
            let compensation = 0;
            for (let i = 0; i < n; ++i) {
                const summand = v[i];
                const t = sum + summand;
                if (Math.abs(sum) >= Math.abs(summand)) compensation += sum - t + summand;
                else compensation += summand - t + sum;
                sum = t;
            }
            return sum + compensation;
        },
        wasm: () => wasmNeumaierSum(v),
    };
});

table("Matrix product — WASM_MIN_MATMUL_OPS (ops = n³)", "n", [4, 6, 8, 12, 16, 32], (n) => {
    const A = random_vector(n * n);
    const B = random_vector(n * n);
    const C = new Float64Array(n * n);
    return {
        reps: 5000,
        js: () => {
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    let sum = 0;
                    for (let k = 0; k < n; ++k) sum += A[i * n + k] * B[k * n + j];
                    C[i * n + j] = sum;
                }
            }
        },
        wasm: () => wasmMatMul(A, n, n, B, n, C),
    };
});

table("Distance matrix (d = 8) — WASM_MIN_ROWS", "rows", [4, 8, 16, 32, 64, 128], (n) => {
    const d = 8;
    const X = random_vector(n * d);
    const D = new Float64Array(n * n);
    /** @type {Float64Array[]} */
    const rows = [];
    for (let i = 0; i < n; ++i) rows.push(X.subarray(i * d, (i + 1) * d));
    return {
        reps: 3000,
        js: () => {
            for (let i = 0; i < n; ++i) {
                for (let j = i + 1; j < n; ++j) {
                    const dist = euclidean(rows[i], rows[j]);
                    D[i * n + j] = dist;
                    D[j * n + i] = dist;
                }
            }
        },
        wasm: () => wasmDistanceMatrix(X, n, d, D),
    };
});
