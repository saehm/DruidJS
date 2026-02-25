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
 * This implementation uses Random Projection hashing (SimHash-style) which works well for
 * cosine similarity and Euclidean distance.
 *
 * Key concepts:
 * - Multiple hash tables increase recall probability
 * - Each hash function projects data onto random hyperplanes
 * - Points on the same side of hyperplanes are hashed together
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
     * @param {ParametersLSH} [parameters={}] - Configuration parameters
     */
    constructor(
        elements,
        parameters = {
            metric: euclidean,
            numHashTables: 10,
            numHashFunctions: 10,
            seed: 1212,
        },
    ) {
        // Handle empty initialization - use dummy element
        const hasElements = elements && elements.length > 0;
        const firstElement = /** @type {T} */ (hasElements ? elements[0] : new Float64Array([0]));

        super([firstElement], parameters);

        this._metric = this._parameters.metric ?? euclidean;
        this._numHashTables = this._parameters.numHashTables ?? 10;
        this._numHashFunctions = this._parameters.numHashFunctions ?? 10;
        this._seed = this._parameters.seed ?? 1212;
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
     * Initialize random projection vectors for all hash tables.
     * @private
     */
    _initializeHashFunctions() {
        const dim = this._elements[0]?.length ?? 0;

        for (let t = 0; t < this._numHashTables; t++) {
            const tableProjections = [];
            const tableOffsets = [];

            for (let h = 0; h < this._numHashFunctions; h++) {
                // Generate random projection vector (normalized)
                const projection = new Float64Array(dim);
                let norm = 0;
                for (let i = 0; i < dim; i++) {
                    // Box-Muller transform for normal distribution
                    const u1 = this._randomizer.random;
                    const u2 = this._randomizer.random;
                    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
                    projection[i] = z;
                    norm += z * z;
                }
                // Normalize
                norm = Math.sqrt(norm);
                for (let i = 0; i < dim; i++) {
                    projection[i] /= norm;
                }

                tableProjections.push(projection);
                // Random offset for quantization buckets
                tableOffsets.push(this._randomizer.random);
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
            // Quantize with offset
            const bucket = Math.floor(dot + offsets[i]);
            bits.push(bucket);
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

    /**
     * @param {number} i
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_by_index(i, k = 5) {
        if (i < 0 || i >= this._elements.length) return [];
        return this.search(this._elements[i], k);
    }
}
