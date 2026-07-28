import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";
import { availableParallelism } from "node:os";
import { fileURLToPath } from "node:url";
import { Matrix } from "../src/matrix/Matrix.js";
import { BallTree } from "../src/knn/BallTree.js";
import { euclidean } from "../src/metrics/euclidean.js";
import { wasmDijkstraAPSP, wasmDijkstraAPSPRange } from "../src/dimred/ISOMAP.wasm.js";
import { wasmMatMul, wasmMatMulRange } from "../src/matrix/Matrix.wasm.js";
import { wasmDistanceMatrix, wasmDistanceMatrixRange } from "../src/matrix/distance_matrix.wasm.js";

const __filename = fileURLToPath(import.meta.url);

if (!isMainThread) {
    const { task, data } = workerData;
    if (task === "dijkstra") {
        const { flatN, flatD, outD, n, k, startSrc, endSrc } = data;
        wasmDijkstraAPSPRange(flatN, flatD, outD, n, k, startSrc, endSrc);
    } else if (task === "distance_matrix") {
        const { X_val, n, d, outD, startRow, endRow } = data;
        wasmDistanceMatrixRange(X_val, n, d, outD, startRow, endRow);
    } else if (task === "matmul") {
        const { A_val, rows_A, cols_A, B_val, cols_B, outC, startRow, endRow } = data;
        wasmMatMulRange(A_val, rows_A, cols_A, B_val, cols_B, outC, startRow, endRow);
    }
    parentPort.postMessage("done");
} else {
    console.log("==========================================");
    console.log("DruidJS Parallel WASM Multi-Threading Benchmark");
    console.log(`Available CPU Threads: ${availableParallelism()}`);
    console.log("==========================================\n");

    async function benchmarkParallelDistanceMatrix(n, d, numThreads = 4) {
        const X = new Matrix(n, d, () => Math.random());
        const D_single = new Matrix(n, n);

        const startSingle = performance.now();
        wasmDistanceMatrix(X.values, n, d, D_single.values);
        const timeSingle = performance.now() - startSingle;

        const sharedBuffer = new SharedArrayBuffer(n * n * 8);
        const outDShared = new Float64Array(sharedBuffer);

        const chunkSize = Math.ceil(n / numThreads);
        const workers = [];
        const startParallel = performance.now();

        for (let t = 0; t < numThreads; ++t) {
            const startRow = t * chunkSize;
            const endRow = Math.min((t + 1) * chunkSize, n);
            if (startRow >= n) break;

            const worker = new Worker(__filename, {
                workerData: {
                    task: "distance_matrix",
                    data: {
                        X_val: X.values,
                        n,
                        d,
                        outD: outDShared,
                        startRow,
                        endRow,
                    },
                },
            });
            workers.push(
                new Promise((resolve) => {
                    worker.on("message", resolve);
                })
            );
        }

        await Promise.all(workers);
        const timeParallel = performance.now() - startParallel;

        console.log(`Pairwise Distance Matrix (${n} points, ${d} dims):`);
        console.log(`  Single-Threaded WASM SIMD: ${timeSingle.toFixed(2)} ms`);
        console.log(`  Multi-Threaded WASM (${numThreads} Workers): ${timeParallel.toFixed(2)} ms`);
        console.log(`  Parallel Speedup:           ${(timeSingle / timeParallel).toFixed(2)}x over Single-Thread WASM\n`);
    }

    async function benchmarkParallelMatMul(size, numThreads = 4) {
        const A = new Matrix(size, size, () => Math.random());
        const B = new Matrix(size, size, () => Math.random());
        const C_single = new Matrix(size, size);

        const startSingle = performance.now();
        wasmMatMul(A.values, size, size, B.values, size, C_single.values);
        const timeSingle = performance.now() - startSingle;

        const sharedBuffer = new SharedArrayBuffer(size * size * 8);
        const outCShared = new Float64Array(sharedBuffer);

        const chunkSize = Math.ceil(size / numThreads);
        const workers = [];
        const startParallel = performance.now();

        for (let t = 0; t < numThreads; ++t) {
            const startRow = t * chunkSize;
            const endRow = Math.min((t + 1) * chunkSize, size);
            if (startRow >= size) break;

            const worker = new Worker(__filename, {
                workerData: {
                    task: "matmul",
                    data: {
                        A_val: A.values,
                        rows_A: size,
                        cols_A: size,
                        B_val: B.values,
                        cols_B: size,
                        outC: outCShared,
                        startRow,
                        endRow,
                    },
                },
            });
            workers.push(
                new Promise((resolve) => {
                    worker.on("message", resolve);
                })
            );
        }

        await Promise.all(workers);
        const timeParallel = performance.now() - startParallel;

        console.log(`Matrix Multiplication (${size} x ${size}):`);
        console.log(`  Single-Threaded WASM SIMD: ${timeSingle.toFixed(2)} ms`);
        console.log(`  Multi-Threaded WASM (${numThreads} Workers): ${timeParallel.toFixed(2)} ms`);
        console.log(`  Parallel Speedup:           ${(timeSingle / timeParallel).toFixed(2)}x over Single-Thread WASM\n`);
    }

    async function runParallelDijkstra(n, k = 10, d = 20, numThreads = 4) {
        const X = new Matrix(n, d, () => Math.random());
        const tree = new BallTree(X.to2dArray(), { metric: euclidean, seed: 1212 });

        const kNN = [];
        for (let i = 0; i < n; ++i) {
            kNN.push(tree.search_by_index(i, k).map((n) => ({ index: n.index, distance: n.distance })));
        }
        for (let i = 0; i < n; ++i) {
            for (const neighbor of kNN[i]) {
                const j = neighbor.index;
                const dist = neighbor.distance;
                if (!kNN[j].find((x) => x.index === i)) kNN[j].push({ index: i, distance: dist });
            }
        }

        let maxK = 0;
        for (let i = 0; i < n; ++i) if (kNN[i].length > maxK) maxK = kNN[i].length;

        const flatN = new Int32Array(n * maxK).fill(-1);
        const flatD = new Float64Array(n * maxK).fill(Infinity);
        for (let i = 0; i < n; ++i) {
            const list = kNN[i];
            const i_k = i * maxK;
            for (let idx = 0; idx < list.length; ++idx) {
                flatN[i_k + idx] = list[idx].index;
                flatD[i_k + idx] = list[idx].distance;
            }
        }

        const G_single = new Matrix(n, n, Infinity);
        const startSingle = performance.now();
        wasmDijkstraAPSP(flatN, flatD, G_single.values, n, maxK);
        const timeSingle = performance.now() - startSingle;

        const sharedBuffer = new SharedArrayBuffer(n * n * 8);
        const outDShared = new Float64Array(sharedBuffer);

        const chunkSize = Math.ceil(n / numThreads);
        const workers = [];
        const startParallel = performance.now();

        for (let t = 0; t < numThreads; ++t) {
            const startSrc = t * chunkSize;
            const endSrc = Math.min((t + 1) * chunkSize, n);
            if (startSrc >= n) break;

            const worker = new Worker(__filename, {
                workerData: {
                    task: "dijkstra",
                    data: {
                        flatN,
                        flatD,
                        outD: outDShared,
                        n,
                        k: maxK,
                        startSrc,
                        endSrc,
                    },
                },
            });
            workers.push(
                new Promise((resolve) => {
                    worker.on("message", resolve);
                })
            );
        }

        await Promise.all(workers);
        const timeParallel = performance.now() - startParallel;

        const speedup = (timeSingle / timeParallel).toFixed(2);
        console.log(`Dijkstra APSP (${n} vertices, ${maxK} neighbors):`);
        console.log(`  Single-Threaded WASM SIMD: ${timeSingle.toFixed(2)} ms`);
        console.log(`  Multi-Threaded WASM (${numThreads} Workers): ${timeParallel.toFixed(2)} ms`);
        console.log(`  Parallel Speedup:           ${speedup}x over Single-Thread WASM\n`);
    }

    async function main() {
        console.log("--- Pairwise Distance Matrix Parallel Benchmarks ---");
        await benchmarkParallelDistanceMatrix(2000, 50, 4);
        await benchmarkParallelDistanceMatrix(2000, 50, 8);

        console.log("--- Matrix Multiplication Parallel Benchmarks ---");
        await benchmarkParallelMatMul(1000, 4);
        await benchmarkParallelMatMul(1000, 8);

        console.log("--- Geodesic Dijkstra APSP Parallel Benchmarks ---");
        await runParallelDijkstra(2000, 10, 20, 4);
        await runParallelDijkstra(2000, 10, 20, 8);
    }

    main();
}
