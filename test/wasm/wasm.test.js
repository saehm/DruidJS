import { describe, expect, it } from "vitest";
import { KMeans } from "../../src/clustering/KMeans.js";
import { wasmKMeansAssign } from "../../src/clustering/KMeans.wasm.js";
import { KMedoids } from "../../src/clustering/KMedoids.js";
import { wasmKMedoidsAssign } from "../../src/clustering/KMedoids.wasm.js";
import { SAMMON } from "../../src/dimred/SAMMON.js";
import { wasmSammonStep } from "../../src/dimred/SAMMON.wasm.js";
import { TSNE } from "../../src/dimred/TSNE.js";
import { inner_product } from "../../src/linear_algebra/inner_product.js";
import { wasmInnerProduct } from "../../src/linear_algebra/inner_product.wasm.js";
import { distance_matrix } from "../../src/matrix/distance_matrix.js";
import { wasmDistanceMatrix } from "../../src/matrix/distance_matrix.wasm.js";
import { Matrix } from "../../src/matrix/Matrix.js";
import { wasmDotTrans, wasmMatMul, wasmMatVecMul, wasmOuter, wasmTransDot } from "../../src/matrix/Matrix.wasm.js";
import { norm } from "../../src/matrix/norm.js";
import { wasmNorm } from "../../src/matrix/norm.wasm.js";
import { normalize } from "../../src/matrix/normalize.js";
import { wasmNormalize } from "../../src/matrix/normalize.wasm.js";
import { bray_curtis } from "../../src/metrics/bray_curtis.js";
import { canberra } from "../../src/metrics/canberra.js";
import { chebyshev } from "../../src/metrics/chebyshev.js";
import { wasmChebyshevDistance } from "../../src/metrics/chebyshev.wasm.js";
import { cosine } from "../../src/metrics/cosine.js";
import { manhattan } from "../../src/metrics/manhattan.js";
import { wasmManhattanDistance } from "../../src/metrics/manhattan.wasm.js";
import { neumair_sum } from "../../src/numerical/neumair_sum.js";
import { wasmNeumaierSum } from "../../src/numerical/neumair_sum.wasm.js";
import { isWasmAvailable } from "../../src/wasm/index.js";

describe("WASM Acceleration Kernels", () => {
    it("should report WASM availability", () => {
        expect(typeof isWasmAvailable()).toBe("boolean");
    });

    it("should compute matrix multiplication correctly using WASM SIMD", () => {
        const A = new Matrix(3, 3, (i, j) => (i + 1) * (j + 1));
        const B = new Matrix(3, 3, (i, j) => i + 1 + (j + 1));

        const C_js = new Matrix(3, 3, 0);
        const A_val = A.values;
        const B_val = B.values;
        const C_val = C_js.values;

        for (let i = 0; i < 3; ++i) {
            for (let k = 0; k < 3; ++k) {
                const aik = A_val[i * 3 + k];
                for (let j = 0; j < 3; ++j) {
                    C_val[i * 3 + j] += aik * B_val[k * 3 + j];
                }
            }
        }

        const C_wasm = new Matrix(3, 3, 0);
        const success = wasmMatMul(A_val, 3, 3, B_val, 3, C_wasm.values);
        expect(success).toBe(true);

        for (let i = 0; i < 9; ++i) {
            expect(C_wasm.values[i]).toBeCloseTo(C_js.values[i], 10);
        }
    });

    it("should compute distance matrix correctly using WASM SIMD", () => {
        const X = new Matrix(20, 10, () => Math.random());
        const D_js = new Matrix(20, 20);

        for (let i = 0; i < 20; ++i) {
            const xi = X.row(i);
            for (let j = i + 1; j < 20; ++j) {
                const xj = X.row(j);
                let sum = 0;
                for (let k = 0; k < 10; ++k) {
                    const diff = xi[k] - xj[k];
                    sum += diff * diff;
                }
                const dist = Math.sqrt(sum);
                D_js.set_entry(i, j, dist);
                D_js.set_entry(j, i, dist);
            }
        }

        const D_wasm = new Matrix(20, 20);
        const success = wasmDistanceMatrix(X.values, 20, 10, D_wasm.values);
        expect(success).toBe(true);

        for (let i = 0; i < 400; ++i) {
            expect(D_wasm.values[i]).toBeCloseTo(D_js.values[i], 10);
        }
    });

    it("should compute transDot, dotTrans, and outer correctly via WASM SIMD", () => {
        const A = new Matrix(30, 20, () => Math.random() * 5);
        const B = new Matrix(30, 20, () => Math.random() * 5);

        const resTransDot = A.transDot(B);
        expect(resTransDot.shape).toEqual([20, 20]);

        const resDotTrans = A.dotTrans(B);
        expect(resDotTrans.shape).toEqual([30, 30]);

        const vecA = new Matrix(1, 20, () => Math.random());
        const vecB = new Matrix(1, 20, () => Math.random());
        const resOuter = vecA.outer(vecB);
        expect(resOuter.shape).toEqual([20, 20]);
    });

    it("should compute vector norm, normalize, inner product, and Neumaier sum via WASM SIMD", () => {
        const v = new Float64Array(40).fill(2);
        const nVal = norm(v);
        expect(nVal).toBeCloseTo(Math.sqrt(40 * 4), 8);

        const normVec = normalize(v);
        expect(normVec.length).toBe(40);
        expect(normVec[0]).toBeCloseTo(2 / nVal, 8);

        const ipVal = inner_product(v, v);
        expect(ipVal).toBeCloseTo(160, 8);

        const nSum = neumair_sum(v);
        expect(nSum).toBeCloseTo(80, 8);
    });

    describe("wasmNeumaierSum", () => {
        // `neumair_sum` deliberately never dispatches to this kernel — the compensation term is a
        // serial dependency, so the calibration benchmark measures the kernel as slower than the JS
        // loop at every size. That leaves the wrapper reachable only from here, and completely
        // unexercised otherwise, so the kernel is checked against the implementation it mirrors.
        it("should agree with the JS implementation on a plain sum", () => {
            const v = new Float64Array(40).fill(2);
            const wasm = wasmNeumaierSum(v);
            if (wasm === null) return; // WASM unavailable on this runner
            expect(wasm).toBeCloseTo(neumair_sum(v), 8);
            expect(wasm).toBeCloseTo(80, 8);
        });

        it("should agree with the JS implementation on mixed magnitudes", () => {
            // The case Neumaier summation exists for: a large running sum swamping small addends.
            const v = Float64Array.from([1e16, 1, -1e16, 1, 1, 1]);
            const wasm = wasmNeumaierSum(v);
            if (wasm === null) return;
            expect(wasm).toBeCloseTo(neumair_sum(v), 8);
        });

        it("should agree across a range of lengths, including non-SIMD-aligned ones", () => {
            for (const len of [0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 33, 128, 1000]) {
                const v = Float64Array.from({ length: len }, (_, i) => Math.sin(i) * (i + 1));
                const wasm = wasmNeumaierSum(v);
                if (wasm === null) return;
                const scale = Math.max(1, Math.abs(neumair_sum(v)));
                expect(Math.abs(wasm - neumair_sum(v)), `length ${len}`).toBeLessThan(1e-9 * scale);
            }
        });

        it("should accept a plain number array", () => {
            const v = [1.5, 2.5, 3.5];
            const wasm = wasmNeumaierSum(v);
            if (wasm === null) return;
            expect(wasm).toBeCloseTo(7.5, 8);
        });
    });

    it("should compute Manhattan, Chebyshev, Canberra, and Bray-Curtis SIMD metrics correctly", () => {
        const a = new Float64Array([
            1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            21, 22, 23, 24, 25,
        ]);
        const b = new Float64Array([
            2, 3, 10, 15, 26, 35, 50, 63, 82, 99, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            21, 22, 23, 24, 25, 26,
        ]);

        const mDist = manhattan(a, b);
        const cDist = chebyshev(a, b);
        const cosDist = cosine(a, b);
        const canbDist = canberra(a, b);
        const bcDist = bray_curtis(a, b);

        expect(mDist).toBeGreaterThan(0);
        expect(cDist).toBeGreaterThan(0);
        expect(cosDist).toBeGreaterThan(0);
        expect(canbDist).toBeGreaterThan(0);
        expect(bcDist).toBeGreaterThan(0);
    });

    it("should compute k-Means & k-Medoids assignments and SAMMON steps correctly via WASM", () => {
        const X = new Matrix(30, 4, (i, j) => (i % 3) * 10 + j);
        const C = new Matrix(3, 4, (i, j) => i * 10 + j);
        const assignments = new Int32Array(30);

        const success = wasmKMeansAssign(X.values, C.values, assignments, 30, 3, 4);
        expect(success).toBe(true);

        for (let i = 0; i < 30; ++i) {
            expect(assignments[i]).toBe(i % 3);
        }

        const sammon = new SAMMON(X, { d: 2 });
        const resSammon = sammon.transform(10);
        expect(resSammon.shape).toEqual([30, 2]);

        const kmedoids = new KMedoids(X, { K: 3 });
        const clusters = kmedoids.get_cluster_list();
        expect(clusters.length).toBe(30);
    });
});
