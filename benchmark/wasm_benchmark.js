import { inner_product } from "../src/linear_algebra/inner_product.js";
import { distance_matrix } from "../src/matrix/distance_matrix.js";
import { Matrix } from "../src/matrix/Matrix.js";
import { norm } from "../src/matrix/norm.js";
import { bray_curtis } from "../src/metrics/bray_curtis.js";
import { canberra } from "../src/metrics/canberra.js";
import { neumair_sum } from "../src/numerical/neumair_sum.js";
import {
    wasmBrayCurtisDistance,
    wasmCanberraDistance,
    wasmDistanceMatrix,
    wasmInnerProduct,
    wasmMatMul,
    wasmNeumaierSum,
    wasmNorm,
    wasmNormalize,
} from "../src/wasm/index.js";

console.log("==========================================");
console.log("DruidJS WASM SIMD Acceleration Benchmark");
console.log("==========================================\n");

function benchmarkMatMul(size) {
    const A = new Matrix(size, size, () => Math.random());
    const B = new Matrix(size, size, () => Math.random());

    const A_val = A.values;
    const B_val = B.values;
    const C_js = new Matrix(size, size, 0).values;
    const C_wasm = new Matrix(size, size, 0).values;

    const startJS = performance.now();
    const rows_A = size, cols_A = size, cols_B = size;
    for (let i = 0; i < rows_A; ++i) {
        const i_cols_A = i * cols_A;
        const i_cols_B = i * cols_B;
        for (let k = 0; k < cols_A; ++k) {
            const aik = A_val[i_cols_A + k];
            if (aik === 0) continue;
            const k_cols_B = k * cols_B;
            for (let j = 0; j < cols_B; ++j) {
                C_js[i_cols_B + j] += aik * B_val[k_cols_B + j];
            }
        }
    }
    const timeJS = performance.now() - startJS;

    const startWASM = performance.now();
    wasmMatMul(A_val, size, size, B_val, size, C_wasm);
    const timeWASM = performance.now() - startWASM;

    const speedup = (timeJS / timeWASM).toFixed(2);
    console.log(`Matrix.dot (${size} x ${size}):`);
    console.log(`  JS Time:   ${timeJS.toFixed(2)} ms`);
    console.log(`  WASM Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Speedup:   ${speedup}x\n`);
}

function benchmarkDistanceMatrix(n, d) {
    const X = new Matrix(n, d, () => Math.random());
    const D_js = new Matrix(n, n);
    const D_wasm = new Matrix(n, n);

    const startJS = performance.now();
    for (let i = 0; i < n; ++i) {
        const xi = X.row(i);
        for (let j = i + 1; j < n; ++j) {
            const xj = X.row(j);
            let sum = 0;
            for (let k = 0; k < d; ++k) {
                const diff = xi[k] - xj[k];
                sum += diff * diff;
            }
            const dist = Math.sqrt(sum);
            D_js.set_entry(i, j, dist);
            D_js.set_entry(j, i, dist);
        }
    }
    const timeJS = performance.now() - startJS;

    const startWASM = performance.now();
    wasmDistanceMatrix(X.values, n, d, D_wasm.values);
    const timeWASM = performance.now() - startWASM;

    const speedup = (timeJS / timeWASM).toFixed(2);
    console.log(`distance_matrix (${n} samples, ${d} dims):`);
    console.log(`  JS Time:   ${timeJS.toFixed(2)} ms`);
    console.log(`  WASM Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Speedup:   ${speedup}x\n`);
}

function benchmarkInnerProduct(len, iters = 100000) {
    const a = new Float64Array(len).fill(1.5);
    const b = new Float64Array(len).fill(2.5);

    const startJS = performance.now();
    for (let i = 0; i < iters; ++i) {
        let sum = 0;
        for (let k = 0; k < len; ++k) sum += a[k] * b[k];
    }
    const timeJS = performance.now() - startJS;

    const startWASM = performance.now();
    for (let i = 0; i < iters; ++i) {
        wasmInnerProduct(a, b);
    }
    const timeWASM = performance.now() - startWASM;

    console.log(`Inner Product (${len} elements, ${iters} calls):`);
    console.log(`  JS Time:   ${timeJS.toFixed(2)} ms`);
    console.log(`  WASM Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Speedup:   ${(timeJS / timeWASM).toFixed(2)}x\n`);
}

function benchmarkVectorNorm(len, iters = 100000) {
    const v = new Float64Array(len).fill(3.14);

    const startJS = performance.now();
    for (let i = 0; i < iters; ++i) {
        let sum = 0;
        for (let k = 0; k < len; ++k) sum += v[k] * v[k];
        Math.sqrt(sum);
    }
    const timeJS = performance.now() - startJS;

    const startWASM = performance.now();
    for (let i = 0; i < iters; ++i) {
        wasmNorm(v);
    }
    const timeWASM = performance.now() - startWASM;

    console.log(`Vector Norm (${len} elements, ${iters} calls):`);
    console.log(`  JS Time:   ${timeJS.toFixed(2)} ms`);
    console.log(`  WASM Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Speedup:   ${(timeJS / timeWASM).toFixed(2)}x\n`);
}

function benchmarkNeumaierSum(len, iters = 50000) {
    const v = new Float64Array(len).fill(1.0000001);

    const startJS = performance.now();
    for (let i = 0; i < iters; ++i) {
        let sum = 0, c = 0;
        for (let k = 0; k < len; ++k) {
            const val = v[k];
            const t = sum + val;
            if (Math.abs(sum) >= Math.abs(val)) c += (sum - t) + val;
            else c += (val - t) + sum;
            sum = t;
        }
    }
    const timeJS = performance.now() - startJS;

    const startWASM = performance.now();
    for (let i = 0; i < iters; ++i) {
        wasmNeumaierSum(v);
    }
    const timeWASM = performance.now() - startWASM;

    console.log(`Neumaier Sum (${len} elements, ${iters} calls):`);
    console.log(`  JS Time:   ${timeJS.toFixed(2)} ms`);
    console.log(`  WASM Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Speedup:   ${(timeJS / timeWASM).toFixed(2)}x\n`);
}

console.log("--- Matrix Multiplication Benchmarks ---");
benchmarkMatMul(100);
benchmarkMatMul(300);

console.log("--- Pairwise Distance Matrix Benchmarks ---");
benchmarkDistanceMatrix(500, 50);
benchmarkDistanceMatrix(1000, 50);

console.log("--- Vector Operations & Metrics Benchmarks ---");
benchmarkInnerProduct(128, 100000);
benchmarkVectorNorm(128, 100000);
benchmarkNeumaierSum(128, 50000);
