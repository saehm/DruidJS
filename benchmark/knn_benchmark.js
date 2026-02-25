#!/usr/bin/env node
/**
 * KNN Algorithms Benchmark
 *
 * This benchmark compares the performance characteristics of all KNN implementations:
 * - Construction/indexing time
 * - Search time for various k values
 * - Recall accuracy for approximate methods
 * - Scalability with dataset size
 *
 * Usage: node benchmark/knn_benchmark.js
 */

import { performance } from "perf_hooks";
import { Annoy, BallTree, HNSW, KDTree, LSH, NaiveKNN, NNDescent } from "../src/knn/index.js";
import { euclidean } from "../src/metrics/index.js";

/**
 * Generate random points in n-dimensional space
 * @param {number} n - Number of points
 * @param {number} dim - Dimensionality
 * @returns {Float64Array[]}
 */
function generateRandomPoints(n, dim) {
  const points = [];
  for (let i = 0; i < n; i++) {
    const point = new Float64Array(dim);
    for (let j = 0; j < dim; j++) {
      point[j] = Math.random() * 100;
    }
    points.push(point);
  }
  return points;
}

/**
 * Generate clustered points for more realistic benchmark
 * @param {number} n - Number of points
 * @param {number} dim - Dimensionality
 * @param {number} numClusters - Number of clusters
 * @returns {Float64Array[]}
 */
function generateClusteredPoints(n, dim, numClusters) {
  const points = [];
  const clusterCenters = [];

  // Generate cluster centers
  for (let c = 0; c < numClusters; c++) {
    const center = new Float64Array(dim);
    for (let j = 0; j < dim; j++) {
      center[j] = Math.random() * 100;
    }
    clusterCenters.push(center);
  }

  // Generate points around cluster centers
  for (let i = 0; i < n; i++) {
    const center = clusterCenters[i % numClusters];
    const point = new Float64Array(dim);
    for (let j = 0; j < dim; j++) {
      point[j] = center[j] + (Math.random() - 0.5) * 10;
    }
    points.push(point);
  }

  return points;
}

/**
 * Measure execution time of a function
 * @template T
 * @param {() => T} fn - Function to measure
 * @returns {{ result: T, time: number }}
 */
function measureTime(fn) {
  const start = performance.now();
  const result = fn();
  const end = performance.now();
  return { result, time: end - start };
}

/**
 * Calculate recall for approximate methods
 * @param {number[]} approximate - Indices from approximate method
 * @param {number[]} exact - Indices from exact method
 * @returns {number}
 */
function calculateRecall(approximate, exact) {
  const exactSet = new Set(exact);
  let found = 0;
  for (const idx of approximate) {
    if (exactSet.has(idx)) {
      found++;
    }
  }
  return found / exact.length;
}

/**
 * Run benchmark for a single KNN algorithm
 * @param {string} name - Algorithm name
 * @param {object} config - Algorithm configuration
 * @param {Float64Array[]} trainPoints - Training points
 * @param {Float64Array[]} testPoints - Test queries
 * @param {number[]} kValues - Different k values to test
 * @param {NaiveKNN} groundTruth - Ground truth for recall calculation
 */
function runBenchmark(name, config, trainPoints, testPoints, kValues, groundTruth) {
  console.log(`\n${"=".repeat(60)}`);
  console.log(`Benchmarking: ${name}`);
  console.log("=".repeat(60));

  // Measure construction time
  console.log("\n📊 Construction/Indexing:");
  const { result: knn, time: buildTime } = measureTime(() => {
    if (name === "HNSW" || name === "LSH" || name === "Annoy") {
      // These support incremental adding
      const instance = new config.class([], config.params);
      instance.add(trainPoints);
      return instance;
    } else {
      return new config.class(trainPoints, config.params);
    }
  });
  console.log(`  Build time: ${buildTime.toFixed(2)}ms`);

  // Measure search time for different k values
  console.log("\n🔍 Search Performance:");
  const results = {};

  for (const k of kValues) {
    const searchTimes = [];
    const recalls = [];

    for (const query of testPoints) {
      const { time: searchTime } = measureTime(() => {
        return knn.search(query, k);
      });
      searchTimes.push(searchTime);

      // Calculate recall for approximate methods
      if (groundTruth && config.approximate) {
        const exactNeighbors = groundTruth.search(query, k);
        const approxNeighbors = knn.search(query, k);
        const exactIndices = exactNeighbors.map((n) => n.index);
        const approxIndices = approxNeighbors.map((n) => n.index);
        recalls.push(calculateRecall(approxIndices, exactIndices));
      }
    }

    const avgTime = searchTimes.reduce((a, b) => a + b, 0) / searchTimes.length;
    const minTime = Math.min(...searchTimes);
    const maxTime = Math.max(...searchTimes);

    console.log(`\n  k=${k}:`);
    console.log(`    Avg search time: ${avgTime.toFixed(3)}ms`);
    console.log(`    Min/Max: ${minTime.toFixed(3)}ms / ${maxTime.toFixed(3)}ms`);
    console.log(`    Queries/sec: ${(1000 / avgTime).toFixed(1)}`);

    if (recalls.length > 0) {
      const avgRecall = recalls.reduce((a, b) => a + b, 0) / recalls.length;
      console.log(`    Recall: ${(avgRecall * 100).toFixed(1)}%`);
      results[`recall_k${k}`] = avgRecall;
    }

    results[`avg_time_k${k}`] = avgTime;
  }

  results.buildTime = buildTime;
  return results;
}

/**
 * Run scalability benchmark
 * @param {string} name - Algorithm name
 * @param {object} config - Algorithm configuration
 * @param {number[]} datasetSizes - Different dataset sizes to test
 * @param {number} dim - Dimensionality
 * @param {number} k - Number of neighbors
 */
function runScalabilityBenchmark(name, config, datasetSizes, dim, k) {
  console.log(`\n${"=".repeat(60)}`);
  console.log(`Scalability Benchmark: ${name}`);
  console.log("=".repeat(60));

  console.log("\n📈 Dataset Size vs Performance:");
  console.log("Size\t\tBuild(ms)\tSearch(ms)");
  console.log("-".repeat(50));

  for (const size of datasetSizes) {
    const points = generateRandomPoints(size, dim);
    const query = points[0];

    // Warm up
    try {
      let knn;
      if (name === "HNSW" || name === "LSH" || name === "Annoy") {
        knn = new config.class([], config.params);
        knn.add(points);
      } else {
        knn = new config.class(points, config.params);
      }
      knn.search(query, k);

      // Measure build time
      const { time: buildTime } = measureTime(() => {
        if (name === "HNSW" || name === "LSH" || name === "Annoy") {
          const instance = new config.class([], config.params);
          instance.add(points);
          return instance;
        } else {
          return new config.class(points, config.params);
        }
      });

      // Measure search time (average of 10 searches)
      const searchTimes = [];
      for (let i = 0; i < 10; i++) {
        const { time } = measureTime(() => knn.search(query, k));
        searchTimes.push(time);
      }
      const avgSearchTime = searchTimes.reduce((a, b) => a + b, 0) / searchTimes.length;

      console.log(`${size}\t\t${buildTime.toFixed(2)}\t\t${avgSearchTime.toFixed(3)}`);
    } catch (error) {
      console.log(`${size}\t\tFAILED: ${error.message}`);
    }
  }
}

/**
 * Main benchmark runner
 */
async function main() {
  console.log("\n" + "=".repeat(60));
  console.log("KNN Algorithms Benchmark Suite");
  console.log("=".repeat(60));

  // Configuration
  const TRAIN_SIZE = 1000;
  const TEST_SIZE = 100;
  const DIMENSIONS = 10;
  const K_VALUES = [1, 5, 10, 20];
  const SEED = 42;

  console.log(`\n📋 Configuration:`);
  console.log(`  Training points: ${TRAIN_SIZE}`);
  console.log(`  Test queries: ${TEST_SIZE}`);
  console.log(`  Dimensions: ${DIMENSIONS}`);
  console.log(`  k values: ${K_VALUES.join(", ")}`);
  console.log(`  Seed: ${SEED}`);

  // Generate dataset
  console.log("\n🎲 Generating dataset...");
  const trainPoints = generateClusteredPoints(TRAIN_SIZE, DIMENSIONS, 10);
  const testPoints = generateRandomPoints(TEST_SIZE, DIMENSIONS);
  console.log(`  Generated ${trainPoints.length} training points`);
  console.log(`  Generated ${testPoints.length} test queries`);

  // Define algorithms to benchmark
  const algorithms = {
    "NaiveKNN (Brute Force)": {
      class: NaiveKNN,
      params: { metric: euclidean },
      approximate: false,
    },
    BallTree: {
      class: BallTree,
      params: { metric: euclidean },
      approximate: false,
    },
    KDTree: {
      class: KDTree,
      params: { metric: euclidean },
      approximate: false,
    },
    HNSW: {
      class: HNSW,
      params: { metric: euclidean, m: 16, ef_construction: 200, ef: 50, seed: SEED },
      approximate: true,
    },
    NNDescent: {
      class: NNDescent,
      params: { metric: euclidean, K: 10, rho: 0.8, delta: 0.001, seed: SEED },
      approximate: true,
    },
    LSH: {
      class: LSH,
      params: { metric: euclidean, numHashTables: 10, numHashFunctions: 10, seed: SEED },
      approximate: true,
    },
    Annoy: {
      class: Annoy,
      params: { metric: euclidean, numTrees: 10, maxPointsPerLeaf: 10, seed: SEED },
      approximate: true,
    },
  };

  // Create ground truth using NaiveKNN
  console.log("\n🎯 Creating ground truth with NaiveKNN...");
  const { result: groundTruth } = measureTime(() => {
    return new NaiveKNN(trainPoints, { metric: euclidean });
  });
  console.log("  Ground truth ready");

  // Run benchmarks
  const allResults = {};
  for (const [name, config] of Object.entries(algorithms)) {
    allResults[name] = runBenchmark(name, config, trainPoints, testPoints, K_VALUES, groundTruth);
  }

  // Summary comparison
  console.log("\n" + "=".repeat(60));
  console.log("📊 Summary Comparison (k=10)");
  console.log("=".repeat(60));
  console.log("\nAlgorithm\t\tBuild(ms)\tSearch(ms)\tRecall(%)");
  console.log("-".repeat(70));

  for (const [name, results] of Object.entries(allResults)) {
    const buildTime = results.buildTime.toFixed(2);
    const searchTime = results.avg_time_k10?.toFixed(3) || "N/A";
    const recall = results.recall_k10 ? (results.recall_k10 * 100).toFixed(1) : "100.0";
    console.log(
      `${name.padEnd(20)}\t${buildTime.padStart(8)}\t${searchTime.padStart(8)}\t${recall.padStart(8)}`,
    );
  }

  // Scalability benchmark for selected algorithms
  console.log("\n" + "=".repeat(60));
  console.log("📈 Scalability Analysis");
  console.log("=".repeat(60));

  const datasetSizes = [100, 500, 1000, 2000];
  const scalabilityAlgorithms = {
    NaiveKNN: algorithms["NaiveKNN (Brute Force)"],
    BallTree: algorithms["BallTree"],
    KDTree: algorithms["KDTree"],
    HNSW: algorithms["HNSW"],
  };

  for (const [name, config] of Object.entries(scalabilityAlgorithms)) {
    runScalabilityBenchmark(name, config, datasetSizes, DIMENSIONS, 10);
  }

  // Recommendations
  console.log("\n" + "=".repeat(60));
  console.log("💡 Recommendations");
  console.log("=".repeat(60));
  console.log(`
Based on the benchmark results, here are some guidelines for choosing a KNN algorithm:

1. Exact Search (Small to Medium Datasets):
   - KDTree: Best for low dimensions (d < 20), balanced performance
   - BallTree: Good for moderate dimensions, flexible metric support
   - NaiveKNN: Use only for very small datasets or when exact results are critical

2. Approximate Search (Large Datasets or High Dimensions):
   - HNSW: Best overall approximate method, excellent recall
   - Annoy: Good balance of speed and accuracy, easy to tune
   - LSH: Best for very high dimensions (d > 50), fast but lower recall
   - NNDescent: Good for graph-based applications

3. Use Cases:
   - d < 10: KDTree or BallTree
   - 10 < d < 50: HNSW or Annoy
   - d > 50: LSH or HNSW with higher ef parameter
   - Real-time applications: HNSW or Annoy
   - Memory constrained: KDTree or LSH
`);

  console.log("\n✅ Benchmark completed!\n");
}

main().catch((error) => {
  console.error("Benchmark failed:", error);
  process.exit(1);
});
