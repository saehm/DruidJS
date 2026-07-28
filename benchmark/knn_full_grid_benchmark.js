import { Annoy } from "../src/knn/Annoy.js";
import { BallTree } from "../src/knn/BallTree.js";
import { HNSW } from "../src/knn/HNSW.js";
import { KDTree } from "../src/knn/KDTree.js";
import { LSH } from "../src/knn/LSH.js";
import { NaiveKNN } from "../src/knn/NaiveKNN.js";
import { NNDescent } from "../src/knn/NNDescent.js";
import { euclidean } from "../src/metrics/euclidean.js";

/**
 * Runtime and recall comparison across every KNN implementation.
 *
 * Three properties make the output usable as a regression baseline rather than a one-off snapshot:
 *
 * - **The data is seeded and generated here.** Using `Math.random` gave every run a different
 *   dataset, so two runs could not be compared and a regression was invisible. The generator is
 *   deliberately self-contained rather than the library `Randomizer`, so that changing the
 *   library's RNG does not change the benchmark inputs underneath a comparison.
 * - **Build and query are timed separately.** They were previously summed, which for an index built
 *   once and queried many times reports approximately the wrong number: the whole point of an
 *   approximate index is that its build cost amortizes.
 * - **Distance calls are counted.** A structure that quietly degrades to a linear scan still reports
 *   100% recall, because scanning everything is exact. Metric calls per query relative to `N` is
 *   what separates "accurate" from "not actually pruning anything".
 */

/**
 * Self-contained xorshift32, so benchmark inputs never depend on the library under test.
 *
 * @param {number} seed
 * @returns {() => number} Successive floats in [0, 1).
 */
function xorshift32(seed) {
    let x = seed | 0 || 1;
    return () => {
        x ^= x << 13;
        x ^= x >>> 17;
        x ^= x << 5;
        return (x >>> 0) / 4294967296;
    };
}

/**
 * @param {number} n
 * @param {number} d
 * @param {() => number} random
 * @returns {number[][]}
 */
function make_points(n, d, random) {
    return Array.from({ length: n }, () => Array.from({ length: d }, () => random()));
}

/**
 * Exact K nearest neighbors, used as the recall reference.
 *
 * @param {number[][]} points
 * @param {number[][]} queries
 * @param {number} K
 * @returns {Set<number>[]}
 */
function compute_ground_truth(points, queries, K) {
    return queries.map((query) => {
        const distances = points.map((p, index) => {
            let sum = 0;
            for (let d = 0; d < p.length; ++d) {
                const diff = query[d] - p[d];
                sum += diff * diff;
            }
            return { index, distance: sum };
        });
        distances.sort((a, b) => a.distance - b.distance);
        return new Set(distances.slice(0, K).map((d) => d.index));
    });
}

/**
 * @param {{ index: number }[][]} results
 * @param {Set<number>[]} ground_truth
 * @returns {number} Mean fraction of the true neighbors that were returned.
 */
function compute_recall(results, ground_truth) {
    let total = 0;
    const K = ground_truth[0].size;
    for (let q = 0; q < results.length; ++q) {
        let matches = 0;
        for (const candidate of results[q]) {
            if (ground_truth[q].has(candidate.index)) ++matches;
        }
        total += matches / K;
    }
    return total / results.length;
}

/**
 * Times one implementation, then re-runs it with an instrumented metric to count distance
 * computations. The count needs its own pass because wrapping the metric would otherwise distort
 * the timings it is reported alongside.
 *
 * @param {(metric: import("../src/metrics/index.js").Metric) => any} build_index
 * @param {number[][]} queries
 * @param {number} K
 * @param {Set<number>[]} ground_truth
 * @returns {{ build: number; query: number; recall: number; build_calls: number; query_calls: number }}
 */
function measure(build_index, queries, K, ground_truth) {
    const t_build = performance.now();
    const index = build_index(euclidean);
    const build = performance.now() - t_build;

    const t_query = performance.now();
    const results = queries.map((q) => index.search(q, K));
    const query = performance.now() - t_query;

    let calls = 0;
    /** @type {import("../src/metrics/index.js").Metric} */
    const counting = (a, b) => {
        ++calls;
        return euclidean(a, b);
    };
    const counted = build_index(counting);
    const build_calls = calls;
    calls = 0;
    for (const q of queries) counted.search(q, K);

    return { build, query, recall: compute_recall(results, ground_truth), build_calls, query_calls: calls };
}

function run() {
    const N_values = [500, 1000, 2500];
    const D_values = [2, 5, 10, 30, 64, 784];
    const K = 15;
    const NUM_QUERIES = 200;
    const SEED = 42;

    console.log("==========================================");
    console.log("DruidJS KNN Benchmark (runtime, recall, work done)");
    console.log(`K=${K}, ${NUM_QUERIES} held-out queries, seed=${SEED}`);
    console.log("==========================================\n");
    console.log("`dists/query` is metric calls per query as a share of N: ~100% means the structure");
    console.log("is scanning the dataset, so its recall says nothing about how well it prunes.\n");

    console.log("| N | D | Method | Build | Query | µs/query | Recall | dists/query |");
    console.log("| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: |");

    for (const N of N_values) {
        for (const D of D_values) {
            const random = xorshift32(SEED);
            const points = make_points(N, D, random);
            // Held out: querying with points that are in the index makes the first hit trivial.
            const queries = make_points(NUM_QUERIES, D, random);
            const ground_truth = compute_ground_truth(points, queries, K);

            /** @type {[string, (metric: any) => any][]} */
            const implementations = [
                ["NaiveKNN", (metric) => new NaiveKNN(points, { metric })],
                ["BallTree", (metric) => new BallTree(points, { metric })],
                ["KDTree", (metric) => new KDTree(points, { metric })],
                ["Annoy", (metric) => new Annoy(points, { metric, numTrees: 10, maxPointsPerLeaf: 10 })],
                ["LSH", (metric) => new LSH(points, { metric, numHashTables: 5, numHashFunctions: 5 })],
                ["NNDescent", (metric) => new NNDescent(points, { metric, samples: 20 })],
                ["HNSW", (metric) => new HNSW(points, { metric, m: 16, ef_construction: 200, seed: SEED })],
            ];

            for (const [name, build_index] of implementations) {
                let row;
                try {
                    row = measure(build_index, queries, K, ground_truth);
                } catch (error) {
                    // Never swallow this: a failing implementation used to be indistinguishable
                    // from a slow one, both showing up as `N/A`.
                    console.log(`| ${N} | ${D} | ${name} | FAILED: ${error.message} | | | | |`);
                    continue;
                }
                const per_query = (row.query / NUM_QUERIES) * 1000;
                const scanned = row.query_calls / NUM_QUERIES / N;
                console.log(
                    `| ${N} | ${D} | ${name} | ${row.build.toFixed(1)}ms | ${row.query.toFixed(1)}ms | ` +
                        `${per_query.toFixed(1)} | ${(row.recall * 100).toFixed(0)}% | ${(scanned * 100).toFixed(0)}% |`,
                );
            }
        }
    }
}

run();
