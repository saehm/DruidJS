import { BallTree } from "../src/knn/BallTree.js";
import { HNSW } from "../src/knn/HNSW.js";
import { NaiveKNN } from "../src/knn/NaiveKNN.js";
import { euclidean } from "../src/metrics/euclidean.js";
import { wasmDistanceMatrix } from "../src/wasm/index.js";

// Pure JS Naive KNN distance matrix calculation (bypassing WASM for benchmarking)
function pureJsNaiveKnn(points, K) {
    const N = points.length;
    const D = points[0].length;
    const neighbors = new Array(N);

    for (let i = 0; i < N; ++i) {
        const p_i = points[i];
        const dists = [];
        for (let j = 0; j < N; ++j) {
            const p_j = points[j];
            let sum = 0;
            for (let d = 0; d < D; ++d) {
                const diff = p_i[d] - p_j[d];
                sum += diff * diff;
            }
            dists.push({ index: j, distance: Math.sqrt(sum) });
        }
        dists.sort((a, b) => a.distance - b.distance);
        neighbors[i] = dists.slice(0, K);
    }
    return neighbors;
}

function runKnnGridBenchmark() {
    console.log("==========================================");
    console.log("DruidJS KNN Grid Runtime Evaluation (N x D Grid)");
    console.log("==========================================\n");

    const N_values = [100, 500, 1000, 2500, 5000];
    const D_values = [2, 10, 64, 784];
    const K = 15;

    console.log("| N | D | WASM NaiveKNN | Pure JS NaiveKNN | Pure JS BallTree | Pure JS HNSW | Best Choice |");
    console.log("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |");

    for (const N of N_values) {
        for (const D of D_values) {
            // Generate uniform dataset N x D
            const points = Array.from({ length: N }, () =>
                Array.from({ length: D }, () => Math.random())
            );

            // 1. WASM SIMD NaiveKNN
            const t0_wasm = performance.now();
            const naiveWasm = new NaiveKNN(points, { metric: euclidean });
            for (let i = 0; i < Math.min(N, 100); ++i) naiveWasm.search(points[i], K);
            const time_wasm = performance.now() - t0_wasm;

            // 2. Pure JS NaiveKNN
            const t0_pure = performance.now();
            for (let i = 0; i < Math.min(N, 100); ++i) {
                const p_i = points[i];
                const dists = [];
                for (let j = 0; j < N; ++j) {
                    const p_j = points[j];
                    let sum = 0;
                    for (let d = 0; d < D; ++d) {
                        const diff = p_i[d] - p_j[d];
                        sum += diff * diff;
                    }
                    dists.push({ index: j, distance: Math.sqrt(sum) });
                }
                dists.sort((a, b) => a.distance - b.distance);
                dists.slice(0, K);
            }
            const time_pure = performance.now() - t0_pure;

            // 3. BallTree (Build + 100 queries)
            const t0_ball = performance.now();
            const ball = new BallTree(points, { metric: euclidean });
            for (let i = 0; i < Math.min(N, 100); ++i) ball.search(points[i], K);
            const time_ball = performance.now() - t0_ball;

            // 4. HNSW (Build + 100 queries)
            const t0_hnsw = performance.now();
            const hnsw = new HNSW(points, { metric: euclidean, m: 16, ef_construction: 200, seed: 42 });
            for (let i = 0; i < Math.min(N, 100); ++i) hnsw.search(points[i], K);
            const time_hnsw = performance.now() - t0_hnsw;

            // Determine winner
            let winner = "WASM NaiveKNN";
            if (time_hnsw < time_wasm && time_hnsw < time_ball) {
                winner = "HNSW";
            } else if (time_ball < time_wasm && time_ball < time_hnsw) {
                winner = "BallTree";
            } else if (time_pure < time_wasm && time_pure < time_ball && time_pure < time_hnsw) {
                winner = "Pure JS NaiveKNN";
            }

            console.log(
                `| ${N.toLocaleString()} | ${D} | ${time_wasm.toFixed(1)} ms | ${time_pure.toFixed(1)} ms | ${time_ball.toFixed(1)} ms | ${time_hnsw.toFixed(1)} ms | **${winner}** |`
            );
        }
    }
}

runKnnGridBenchmark();
