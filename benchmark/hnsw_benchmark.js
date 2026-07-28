import { HNSW } from "../src/knn/HNSW.js";
import { euclidean } from "../src/metrics/euclidean.js";

function benchmarkHNSW() {
    console.log("==========================================");
    console.log("DruidJS HNSW WASM & Performance Benchmark");
    console.log("==========================================\n");

    const N = 2000;
    const D = 64; // High dimensional data
    const K = 15;

    console.log(`Generating dataset: N = ${N.toLocaleString()} points, D = ${D} dimensions...`);
    const points = Array.from({ length: N }, () =>
        Array.from({ length: D }, () => Math.random())
    );
    const queries = Array.from({ length: 50 }, () =>
        Array.from({ length: D }, () => Math.random())
    );

    // Build HNSW Index
    const startBuild = performance.now();
    const hnsw = new HNSW(points, {
        metric: euclidean,
        m: 16,
        ef_construction: 200,
        seed: 42,
    });
    const buildTime = performance.now() - startBuild;

    // Search 50 Queries
    const startSearch = performance.now();
    let totalNeighborsFound = 0;
    for (let q = 0; q < queries.length; ++q) {
        const results = hnsw.search(queries[q], K);
        totalNeighborsFound += results.length;
    }
    const searchTime = performance.now() - startSearch;

    console.log(`\nResults:`);
    console.log(`  HNSW Build Time (${N} pts, ${D}-D): ${buildTime.toFixed(2)} ms`);
    console.log(`  HNSW 50 Queries Search Time:      ${searchTime.toFixed(2)} ms (${(searchTime / 50).toFixed(3)} ms / query)`);
    console.log(`  Graph Size:                        ${hnsw.size} nodes, ${hnsw.num_layers} layers`);
}

benchmarkHNSW();
