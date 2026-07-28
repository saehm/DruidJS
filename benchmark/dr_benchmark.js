import { KMeans } from "../src/clustering/KMeans.js";
import { KMedoids } from "../src/clustering/KMedoids.js";
import { ISOMAP } from "../src/dimred/ISOMAP.js";
import { MDS } from "../src/dimred/MDS.js";
import { PCA } from "../src/dimred/PCA.js";
import { SAMMON } from "../src/dimred/SAMMON.js";
import { SMACOF } from "../src/dimred/SMACOF.js";
import { TSNE } from "../src/dimred/TSNE.js";
import { UMAP } from "../src/dimred/UMAP.js";
import { Matrix } from "../src/matrix/Matrix.js";
import { wasmMatMul } from "../src/matrix/Matrix.wasm.js";
import { wasmDistanceMatrix } from "../src/matrix/distance_matrix.wasm.js";

console.log("==========================================");
console.log("DruidJS DR Algorithms WASM Benchmark");
console.log("==========================================\n");

function benchmarkPCA(n, d) {
    const X = new Matrix(n, d, () => Math.random());

    // Measure PCA with WASM
    const startWASM = performance.now();
    const pcaWASM = new PCA(X, { d: 2 });
    pcaWASM.transform();
    const timeWASM = performance.now() - startWASM;

    // Benchmark pure JS matrix multiplication inside PCA
    // PCA involves X.T.dot(X) and X.dot(V)
    const X_val = X.values;
    const XT = X.T;
    const XT_val = XT.values;
    const Out_js = new Matrix(d, d, 0).values;

    const startJS = performance.now();
    // Simulate JS XT.dot(X)
    for (let i = 0; i < d; ++i) {
        for (let k = 0; k < n; ++k) {
            const aik = XT_val[i * n + k];
            if (aik === 0) continue;
            for (let j = 0; j < d; ++j) {
                Out_js[i * d + j] += aik * X_val[k * d + j];
            }
        }
    }
    const timeJS_dot = performance.now() - startJS;

    const Out_wasm = new Matrix(d, d, 0).values;
    const startWasmDot = performance.now();
    wasmMatMul(XT_val, d, n, X_val, d, Out_wasm);
    const timeWasm_dot = performance.now() - startWasmDot;

    console.log(`PCA Projection (Data: ${n} x ${d}):`);
    console.log(`  PCA Total Execution Time: ${timeWASM.toFixed(2)} ms`);
    console.log(`  Covariance Matrix (X^T * X) JS:   ${timeJS_dot.toFixed(2)} ms`);
    console.log(`  Covariance Matrix (X^T * X) WASM: ${timeWasm_dot.toFixed(2)} ms`);
    console.log(`  MatMul Speedup: ${(timeJS_dot / timeWasm_dot).toFixed(2)}x\n`);
}

function benchmarkMDS(n, d) {
    const X = new Matrix(n, d, () => Math.random());

    // Total MDS execution with WASM distance matrix
    const startMDS = performance.now();
    const mds = new MDS(X, { d: 2 });
    mds.transform();
    const timeMDS = performance.now() - startMDS;

    // Benchmark distance matrix computation inside MDS initialization (JS vs WASM)
    const startJS = performance.now();
    const D_js = new Matrix(n, n);
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
    const timeJS_dist = performance.now() - startJS;

    const startWASM_dist = performance.now();
    const D_wasm = new Matrix(n, n);
    wasmDistanceMatrix(X.values, n, d, D_wasm.values);
    const timeWASM_dist = performance.now() - startWASM_dist;

    console.log(`MDS (Classical Multidimensional Scaling - ${n} points, ${d} dimensions):`);
    console.log(`  Total MDS Execution Time: ${timeMDS.toFixed(2)} ms`);
    console.log(`  Distance Matrix Init (JS):   ${timeJS_dist.toFixed(2)} ms`);
    console.log(`  Distance Matrix Init (WASM): ${timeWASM_dist.toFixed(2)} ms`);
    console.log(`  Distance Matrix Speedup: ${(timeJS_dist / timeWASM_dist).toFixed(2)}x\n`);
}

function benchmarkSMACOF(n, d, iterations = 20) {
    const X = new Matrix(n, d, () => Math.random());

    // SMACOF with WASM distance matrix & WASM B.dot(Z)
    const startSMACOF = performance.now();
    const smacof = new SMACOF(X, { d: 2 });
    smacof.transform(iterations);
    const timeSMACOF = performance.now() - startSMACOF;

    console.log(`SMACOF (${iterations} iterations - ${n} points, ${d} dimensions):`);
    console.log(`  Total SMACOF Execution Time: ${timeSMACOF.toFixed(2)} ms\n`);
}

function benchmarkTSNE(n, d, iterations = 50) {
    const X = new Matrix(n, d, () => Math.random());
    const tsne = new TSNE(X, { d: 2, perplexity: 30 });
    tsne.init();

    const startTSNE = performance.now();
    tsne.transform(iterations);
    const timeTSNE = performance.now() - startTSNE;

    console.log(`t-SNE (${iterations} iterations - ${n} points, ${d} dims):`);
    console.log(`  Total Execution Time: ${timeTSNE.toFixed(2)} ms`);
    console.log(`  Per Iteration Avg:    ${(timeTSNE / iterations).toFixed(2)} ms\n`);
}

function benchmarkKMeans(n, d, k = 5) {
    const X = new Matrix(n, d, () => Math.random());

    const startKMeans = performance.now();
    const kmeans = new KMeans(X, { K: k });
    kmeans.get_cluster_list();
    const timeKMeans = performance.now() - startKMeans;

    console.log(`k-Means (${n} points, ${d} dims, K=${k}):`);
    console.log(`  Total Execution Time: ${timeKMeans.toFixed(2)} ms\n`);
}

function benchmarkISOMAP(n, d) {
    const X = new Matrix(n, d, () => Math.random());
    const isomap = new ISOMAP(X, { d: 2, neighbors: 10 });

    const startISOMAP = performance.now();
    isomap.transform();
    const timeISOMAP = performance.now() - startISOMAP;

    console.log(`ISOMAP (${n} points, ${d} dims, 10 neighbors):`);
    console.log(`  Total Execution Time: ${timeISOMAP.toFixed(2)} ms\n`);
}

function benchmarkUMAP(n, d, epochs = 200) {
    const X = new Matrix(n, d, () => Math.random());
    const umap = new UMAP(X, { d: 2, n_neighbors: 15, n_epochs: epochs });

    const startUMAP = performance.now();
    umap.transform();
    const timeUMAP = performance.now() - startUMAP;

    console.log(`UMAP (${epochs} epochs - ${n} points, ${d} dims):`);
    console.log(`  Total Execution Time: ${timeUMAP.toFixed(2)} ms`);
    console.log(`  Per Epoch Avg:        ${(timeUMAP / epochs).toFixed(2)} ms\n`);
}

function benchmarkSAMMON(n, d, iterations = 30) {
    const X = new Matrix(n, d, () => Math.random());
    const sammon = new SAMMON(X, { d: 2 });

    const startSAMMON = performance.now();
    sammon.transform(iterations);
    const timeSAMMON = performance.now() - startSAMMON;

    console.log(`SAMMON (${iterations} iterations - ${n} points, ${d} dims):`);
    console.log(`  Total Execution Time: ${timeSAMMON.toFixed(2)} ms`);
    console.log(`  Per Iteration Avg:    ${(timeSAMMON / iterations).toFixed(2)} ms\n`);
}

function benchmarkKMedoids(n, d, k = 5) {
    const X = new Matrix(n, d, () => Math.random());

    const startKMedoids = performance.now();
    const kmedoids = new KMedoids(X, { K: k });
    kmedoids.get_cluster_list();
    const timeKMedoids = performance.now() - startKMedoids;

    console.log(`k-Medoids (${n} points, ${d} dims, K=${k}):`);
    console.log(`  Total Execution Time: ${timeKMedoids.toFixed(2)} ms\n`);
}

console.log("--- PCA Benchmarks ---");
benchmarkPCA(2000, 100);

console.log("--- Classical MDS Benchmarks ---");
benchmarkMDS(1000, 50);

console.log("--- SMACOF Iterative MDS Benchmarks ---");
benchmarkSMACOF(1000, 50, 20);

console.log("--- ISOMAP Benchmarks ---");
benchmarkISOMAP(500, 50);

console.log("--- UMAP Benchmarks ---");
benchmarkUMAP(500, 50, 200);

console.log("--- SAMMON Benchmarks ---");
benchmarkSAMMON(300, 50, 30);

console.log("--- t-SNE Benchmarks ---");
benchmarkTSNE(300, 50, 50);

console.log("--- k-Means & k-Medoids Benchmarks ---");
benchmarkKMeans(2000, 50, 5);
benchmarkKMedoids(500, 50, 5);
