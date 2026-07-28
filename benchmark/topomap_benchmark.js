import { Matrix } from "../src/matrix/Matrix.js";
import { TopoMap } from "../src/dimred/TopoMap.js";
import { DisjointSet } from "../src/datastructure/DisjointSet.js";
import { euclidean } from "../src/metrics/euclidean.js";
import { distance_matrix } from "../src/matrix/distance_matrix.js";

function benchmarkTopoMap() {
    console.log("==========================================");
    console.log("TopoMap Distance Access & Optimization Benchmark");
    console.log("==========================================\n");

    const N = 500;
    const D = 20;
    const X = new Matrix(N, D, () => Math.random());

    // 1. Count actual distance matrix evaluations in lazy mode
    let lazyAccessCount = 0;
    const lazyD = new Matrix(N, N, -1);
    function countLazyDistance(i, j) {
        const D_ij = lazyD.entry(i, j);
        if (D_ij === -1 && i !== j) {
            lazyAccessCount++;
            const dist = euclidean(X.row(i), X.row(j));
            lazyD.set_entry(i, j, dist);
            lazyD.set_entry(j, i, dist);
            return dist;
        }
        return i === j ? 0 : D_ij;
    }

    const E_lazy = [];
    for (let i = 0; i < N; ++i) {
        for (let j = i + 1; j < N; ++j) {
            E_lazy.push([i, j, countLazyDistance(i, j)]);
        }
    }
    const totalEdges = (N * (N - 1)) / 2;

    console.log(`Total Graph Edges N*(N-1)/2: ${totalEdges}`);
    console.log(`Actual Distance Computations in Lazy Mode: ${lazyAccessCount}`);
    console.log(`Percentage of Distances Computed: ${((lazyAccessCount / totalEdges) * 100).toFixed(1)}%\n`);

    // 2. Measure Lazy JS TopoMap Execution Time
    // Temporarily build lazy MST
    const startLazy = performance.now();
    const X_arr = [...X];
    const ds_lazy = new DisjointSet(X_arr);
    const F_lazy = [];
    const E_sorted_lazy = E_lazy.sort((a, b) => a[2] - b[2]);
    for (const [u, v, w] of E_sorted_lazy) {
        const set_u = ds_lazy.find(X_arr[u]);
        const set_v = ds_lazy.find(X_arr[v]);
        if (set_u !== set_v) {
            F_lazy.push([u, v, w]);
            ds_lazy.union(set_u, set_v);
        }
    }
    F_lazy.sort((a, b) => a[2] - b[2]);
    const timeLazy = performance.now() - startLazy;

    // 3. Measure WASM SIMD Precomputed TopoMap Execution Time
    const startWasm = performance.now();
    const D_wasm = distance_matrix(X, euclidean);
    const ds_wasm = new DisjointSet(X_arr);
    const F_wasm = [];
    const E_wasm = [];
    for (let i = 0; i < N; ++i) {
        for (let j = i + 1; j < N; ++j) {
            E_wasm.push([i, j, D_wasm.entry(i, j)]);
        }
    }
    E_wasm.sort((a, b) => a[2] - b[2]);
    for (const [u, v, w] of E_wasm) {
        const set_u = ds_wasm.find(X_arr[u]);
        const set_v = ds_wasm.find(X_arr[v]);
        if (set_u !== set_v) {
            F_wasm.push([u, v, w]);
            ds_wasm.union(set_u, set_v);
        }
    }
    F_wasm.sort((a, b) => a[2] - b[2]);
    const timeWasm = performance.now() - startWasm;

    // 4. Verify Equality of MST Edges
    let maxEdgeDiff = 0;
    let mismatchedEdges = 0;
    if (F_lazy.length !== F_wasm.length) {
        console.error("MST Edge count mismatch!");
    } else {
        for (let i = 0; i < F_lazy.length; ++i) {
            const [u1, v1, w1] = F_lazy[i];
            const [u2, v2, w2] = F_wasm[i];
            if (u1 !== u2 || v1 !== v2) mismatchedEdges++;
            const diff = Math.abs(w1 - w2);
            if (diff > maxEdgeDiff) maxEdgeDiff = diff;
        }
    }

    // 5. Verify Full Projection Equality
    const tm1 = new TopoMap(X);
    const tm2 = new TopoMap(X);
    const Y1 = tm1.transform();
    const Y2 = tm2.transform();

    let maxYDiff = 0;
    for (let i = 0; i < Y1.shape[0] * Y1.shape[1]; ++i) {
        const diff = Math.abs(Y1.values[i] - Y2.values[i]);
        if (diff > maxYDiff) maxYDiff = diff;
    }

    console.log(`Lazy JS MST Execution Time:         ${timeLazy.toFixed(2)} ms`);
    console.log(`WASM SIMD Precomputed MST Time:    ${timeWasm.toFixed(2)} ms`);
    console.log(`Speedup:                            ${(timeLazy / timeWasm).toFixed(2)}x faster\n`);

    console.log(`MST Edge Weight Max Difference:    ${maxEdgeDiff.toExponential(4)}`);
    console.log(`Mismatched MST Edges:               ${mismatchedEdges}`);
    console.log(`Final Projection Y Max Difference: ${maxYDiff.toExponential(4)}`);
}

benchmarkTopoMap();
