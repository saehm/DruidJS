import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { KNN } from "./KNN.js";

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersLSH } from "./index.js" */

/**
 * Locality Sensitive Hashing (LSH) for approximate nearest neighbor search.
 *
 * LSH uses hash functions that map similar items to the same buckets with high probability.
 * This implementation uses the p-stable scheme of Datar et al. for Euclidean distance: each hash
 * projects onto a Gaussian random vector and quantizes the result into buckets of width w.
 *
 * Key concepts:
 * - Multiple hash tables increase recall probability
 * - Each hash function projects data onto a random Gaussian direction
 * - Points landing in the same quantization bucket are hashed together
 * - Combines results from all tables for better accuracy
 *
 * Best suited for:
 * - High-dimensional data where exact methods fail
 * - Approximate nearest neighbor needs
 * - Large datasets where linear scan is too slow
 * - When some false positives/negatives are acceptable
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersLSH>
 * @see {@link https://en.wikipedia.org/wiki/Locality-sensitive_hashing}
 */
export class LSH extends KNN {
    /**
     * Creates a new LSH index.
     *
     * @param {T[]} elements - Elements to index
     * @param {Partial<ParametersLSH>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(elements, parameters = {}) {
        // Handle empty initialization - use dummy element
        const hasElements = elements && elements.length > 0;
        const firstElement = /** @type {T} */ (hasElements ? elements[0] : new Float64Array([0]));

        // `KNN` keeps the parameter object as handed to it, so defaults are merged here rather
        // than left to a default argument, which would only apply when `parameters` is omitted
        // entirely and drop every unspecified default for `new LSH(elements, { seed: 5 })`.
        super(
            [firstElement],
            Object.assign(
                { metric: euclidean, numHashTables: 10, numHashFunctions: 10, bucketWidth: null, seed: 1212 },
                parameters,
            ),
        );

        this._metric = this._parameters.metric;
        this._numHashTables = this._parameters.numHashTables;
        this._numHashFunctions = this._parameters.numHashFunctions;
        this._seed = this._parameters.seed;
        this._randomizer = new Randomizer(this._seed);

        // Hash tables: array of Maps where key is hash bucket, value is array of element indices
        /** @type {Map<string, number[]>[]} */
        this._hashTables = [];

        // Random projection vectors for each hash table and hash function
        /** @type {Float64Array[][]} */
        this._projections = [];

        // Random offsets for each hash table and hash function (for quantization)
        /** @type {number[][]} */
        this._offsets = [];

        // Store dimensionality for later
        /** @type {number} */
        this._dim = firstElement.length;

        // Bucket width w, estimated from the data when not supplied
        /** @type {number} */
        this._bucketWidth = this._parameters.bucketWidth ?? this._estimateBucketWidth(hasElements ? elements : []);

        // Initialize hash functions
        this._initializeHashFunctions();

        // Reset elements if we were initialized with dummy
        if (!hasElements) {
            /** @type {T[]} */
            this._elements = [];
        } else {
            // Clear and re-add elements properly
            /** @type {T[]} */
            this._elements = [];
            this._hashTables = [];
            this._projections = [];
            this._offsets = [];
            this._initializeHashFunctions();
            this.add(elements);
        }
    }

    /**
     * Estimates a bucket width from the spread of the data.
     *
     * Projected gaps are distributed as `N(0, ||u - v||^2)`, so w must sit near the scale of
     * nearby-point distances: much larger and everything collides, much smaller and nothing does.
     *
     * @private
     * @param {T[]} elements
     * @returns {number} Bucket width, always positive.
     */
    _estimateBucketWidth(elements) {
        const n = elements.length;
        if (n < 2) return 1;

        const neighbor_rank = Math.min(10, n - 1);
        const step = Math.max(1, Math.floor(n / Math.min(n, 32)));

        const scales = [];
        for (let i = 0; i < n; i += step) {
            const distances = new Array(n);
            for (let j = 0; j < n; ++j) distances[j] = this._metric(elements[i], elements[j]);
            distances.sort((a, b) => a - b);
            scales.push(distances[neighbor_rank]);
        }

        scales.sort((a, b) => a - b);
        const scale = scales[scales.length >> 1];
        // Several times the near-neighbor distance: `numHashFunctions` hashes are ANDed per table,
        // so the per-hash collision probability is raised to that power and has to stay high.
        return scale > 0 ? 4 * scale : 1;
    }

    /**
     * Initialize random projection vectors for all hash tables.
     * @private
     */
    _initializeHashFunctions() {
        // From `_dim`, not `_elements`: the constructor clears `_elements` before re-initializing
        const dim = this._dim;

        for (let t = 0; t < this._numHashTables; t++) {
            const tableProjections = [];
            const tableOffsets = [];

            for (let h = 0; h < this._numHashFunctions; h++) {
                // Raw N(0, 1); normalizing to unit length would break p-stability
                const projection = new Float64Array(dim);
                for (let i = 0; i < dim; i++) {
                    // Box-Muller transform for normal distribution
                    const u1 = this._randomizer.random;
                    const u2 = this._randomizer.random;
                    projection[i] = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
                }

                tableProjections.push(projection);
                // Offset b, drawn uniformly from [0, w) as the scheme requires.
                tableOffsets.push(this._randomizer.random * this._bucketWidth);
            }

            this._projections.push(tableProjections);
            this._offsets.push(tableOffsets);
            this._hashTables.push(new Map());
        }
    }

    /**
     * Compute hash signature for an element using random projections.
     * @private
     * @param {T} element
     * @param {number} tableIndex
     * @returns {string} Hash signature
     */
    _computeHash(element, tableIndex) {
        const projections = this._projections[tableIndex];
        const offsets = this._offsets[tableIndex];
        const bits = [];

        for (let i = 0; i < this._numHashFunctions; i++) {
            // Compute dot product
            let dot = 0;
            const proj = projections[i];
            for (let j = 0; j < element.length; j++) {
                dot += element[j] * proj[j];
            }
            // Quantize into buckets of width w
            bits.push(Math.floor((dot + offsets[i]) / this._bucketWidth));
        }

        return bits.join(",");
    }

    /**
     * Add elements to the LSH index.
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements) {
        // Extend elements array
        const startIndex = this._elements.length;
        this._elements = this._elements.concat(elements);

        // Hash each new element and add to tables
        for (let i = 0; i < elements.length; i++) {
            const globalIndex = startIndex + i;
            const element = elements[i];

            for (let t = 0; t < this._numHashTables; t++) {
                const hash = this._computeHash(element, t);
                const table = this._hashTables[t];

                if (!table.has(hash)) {
                    table.set(hash, []);
                }
                const bucket = table.get(hash);
                if (bucket) {
                    bucket.push(globalIndex);
                }
            }
        }

        return this;
    }

    /**
     * Search for k approximate nearest neighbors.
     * @param {T} query
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(query, k = 5) {
        const metric = this._metric;
        const elements = this._elements;

        if (elements.length === 0) return [];

        // Collect candidate indices from all hash tables
        const candidates = new Set();

        for (let t = 0; t < this._numHashTables; t++) {
            const hash = this._computeHash(query, t);
            const table = this._hashTables[t];
            const bucket = table.get(hash);

            if (bucket) {
                for (const idx of bucket) {
                    if (idx !== undefined) {
                        candidates.add(idx);
                    }
                }
            }
        }

        // If insufficient candidates found, fall back to linear search
        if (candidates.size < k) {
            // Add more candidates from all buckets or entire dataset
            //const needed = k - candidates.size;

            // First, try to add from neighboring buckets (different hashes)
            for (let t = 0; t < this._numHashTables && candidates.size < k; t++) {
                const table = this._hashTables[t];
                for (const [, bucket] of table) {
                    for (const idx of bucket) {
                        if (idx !== undefined) {
                            candidates.add(idx);
                            if (candidates.size >= k) break;
                        }
                    }
                    if (candidates.size >= k) break;
                }
            }

            // If still not enough, add from entire dataset
            for (let i = 0; i < elements.length && candidates.size < k; i++) {
                candidates.add(i);
            }
        }

        // Compute exact distances for candidates
        /** @type {Heap<{ index: number; distance: number }>} */
        const best = new Heap(null, (d) => d.distance, "max");

        for (const idx of candidates) {
            const element = elements[idx];
            if (!element || element.length !== query.length) continue;

            const dist = metric(query, element);

            if (best.length < k) {
                best.push({ index: idx, distance: dist });
            } else if (dist < (best.first?.value ?? Infinity)) {
                best.pop();
                best.push({ index: idx, distance: dist });
            }
        }

        // Convert to result format
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (best.length > 0) {
            const item = /** @type {{ element: { index: number; distance: number }; value: number }} */ (best.pop());
            result.push({
                element: elements[item.element.index],
                index: item.element.index,
                distance: item.value,
            });
        }

        return result.reverse();
    }
}
