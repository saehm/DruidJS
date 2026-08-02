/**
 * Computes the Bray-Curtis distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Bray-Curtis distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity}
 */
declare function bray_curtis(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the canberra distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The canberra distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
declare function canberra(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The chebyshev distance between `a` and `b`.
 */
declare function chebyshev(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the cosine distance (not similarity) between `a` and `b`.
 *
 * Identical vectors yield exactly `0`, and pairs involving a zero vector yield `π / 2`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The cosine distance between `a` and `b`.
 * @example
 * import { cosine } from "@saehrimnir/druidjs";
 * const a = [1, 2, 3];
 * const b = [4, 5, 6];
 * const distance = cosine(a, b); // 0.9746318461970762
 */
declare function cosine(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the euclidean distance (`l_2`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The euclidean distance between `a` and `b`.
 */
declare function euclidean(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the squared euclidean distance (l_2^2) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The squared euclidean distance between `a` and `b`.

 */
declare function euclidean_squared(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the Goodman-Kruskal gamma coefficient for ordinal association.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First categorical/ordinal variable
 * @param {number[] | Float64Array} b - Second categorical/ordinal variable
 * @returns {number} The Goodman-Kruskal gamma coefficient between `a` and `b` (-1 to 1).
 * @see {@link https://en.wikipedia.org/wiki/Goodman_and_Kruskal%27s_gamma}
 */
declare function goodman_kruskal(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the hamming distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The hamming distance between `a` and `b`.
 */
declare function hamming(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the Haversine distance between two points on a sphere of unit length 1. Multiply the result with the radius of the sphere. (For instance Earth's radius is 6371km)
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - Point [lat1, lon1] in radians
 * @param {number[] | Float64Array} b - Point [lat2, lon2] in radians
 * @returns {number} The Haversine distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Haversine_formula}
 */
declare function haversine(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the jaccard distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The jaccard distance between `a` and `b`.
 */
declare function jaccard(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the manhattan distance (`l_1`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The manhattan distance between `a` and `b`.
 */
declare function manhattan(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the Sokal-Michener distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Sokal-Michener distance between `a` and `b`.

 */
declare function sokal_michener(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the 1D Wasserstein distance (Earth Mover's Distance) between two distributions.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First distribution (histogram or probability mass)
 * @param {number[] | Float64Array} b - Second distribution (histogram or probability mass)
 * @returns {number} The Wasserstein/EMD distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Wasserstein_metric}
 */
declare function wasserstein(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the yule distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The yule distance between `a` and `b`.
 */
declare function yule(a: number[] | Float64Array, b: number[] | Float64Array): number;

type Metric = (a: number[] | Float64Array, b: number[] | Float64Array) => number;

/**
 * Computes the distance matrix of datamatrix `A`.
 *
 * @category Matrix
 * @param {Matrix | Float64Array[] | number[][]} A - Matrix.
 * @param {Metric} [metric=euclidean] - The diistance metric. Default is `euclidean`
 * @returns {Matrix} The distance matrix of `A`.
 */
declare function distance_matrix(A: Matrix | Float64Array[] | number[][], metric?: Metric): Matrix;

/** @import { Metric } from "../metrics/index.js" */
/**
 * Computes the k-nearest neighbors of each row of `A`.
 *
 * @category Matrix
 * @param {Matrix} A - Either the data matrix, or a distance matrix.
 * @param {number} k - The number of neighbors to compute.
 * @param {Metric | "precomputed"} [metric=euclidean] Default is `euclidean`
 * @returns {{ i: number; j: number; distance: number }[][]} The kNN graph.
 */
declare function k_nearest_neighbors(A: Matrix, k: number, metric?: Metric | "precomputed"): {
    i: number;
    j: number;
    distance: number;
}[][];

/**
 * Creates an Array containing `number` numbers from `start` to `end`. If `number = null`.
 *
 * @category Matrix
 * @param {number} start - Start value.
 * @param {number} end - End value.
 * @param {number} [number] - Number of number between `start` and `end`.
 * @returns {number[]} An array with `number` entries, beginning at `start` ending at `end`.
 */
declare function linspace(start: number, end: number, number?: number): number[];

/**
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
declare function max(values: Iterable<number | null>): number;

/**
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
declare function min(values: Iterable<number | null>): number;

/**
 * @typedef {[number, number, number]} WeightedEdge An edge `[u, v, weight]` between the vertices with
 *   index `u` and `v`.
 */
/**
 * Computes a minimum spanning tree of a weighted graph with Kruskal's algorithm.
 *
 * The graph is given as an edge list over vertices `0 … N-1`, so a caller holding a sparse structure
 * (a k-nearest-neighbor graph, say) never has to materialize the `N ⨯ N` matrix. Edges are treated as
 * undirected; passing both `[u, v, w]` and `[v, u, w]` is harmless, the second one is skipped as a
 * cycle.
 *
 * If the graph is disconnected the result is a minimum spanning *forest* — one tree per connected
 * component, and fewer than `N - 1` edges. Callers that need a connected result must guarantee a
 * connected input.
 *
 * @param {WeightedEdge[]} edges - The edges of the graph. Not mutated.
 * @param {number} N - The number of vertices. Vertex indices must be in `[0, N)`.
 * @returns {WeightedEdge[]} The edges of the minimum spanning tree, ascending by weight.
 * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
 */
declare function minimum_spanning_tree(edges: WeightedEdge[], N: number): WeightedEdge[];
/**
 * An edge `[u, v, weight]` between the vertices with
 *   index `u` and `v`.
 */
type WeightedEdge = [number, number, number];

/**
 * Seeded pseudo-random number generator.
 *
 * Implements **sfc32** (Small Fast Counting, Doty-Humphrey), a 128-bit-state counter-based
 * generator that passes TestU01 BigCrush. It is preferred here over the more familiar Mersenne
 * Twister, which needs a 2.5 KB state, is markedly slower to seed, and whose GF(2)-linear output
 * fails BigCrush's matrix-rank and linear-complexity tests.
 *
 * Every operation is 32-bit integer arithmetic (`^ << >>> +` and `Math.imul`), all of which
 * ECMAScript specifies exactly. The stream is therefore identical on every engine for a given seed —
 * unlike floating point transcendentals, which are only "implementation-approximated". See
 * {@link DR} for what that does and does not guarantee about the output of an algorithm.
 *
 * @category Utils
 * @class
 * @see {@link https://pracrand.sourceforge.net/|PractRand, where sfc32 originates}
 * @example
 * const R = new Randomizer(1212);
 * R.random;      // float in [0, 1)
 * R.random_int;  // uint32
 */
declare class Randomizer {
    /**
     * @template T Returns samples from an input Matrix or Array.
     * @param {T[]} A - The input Matrix or Array.
     * @param {number} n - The number of samples.
     * @param {number} seed - The seed for the random number generator.
     * @returns {T[]} - A random selection form `A` of `n` samples.
     */
    static choice<T>(A: T[], n: number, seed?: number): T[];
    /**
     * @param {number} [_seed=new Date().getTime()] - The seed for the random number generator. If `_seed == null` then
     *   the actual time gets used as seed. Default is `new Date().getTime()`
     */
    constructor(_seed?: number);
    /** @type {number} */
    _a: number;
    /** @type {number} */
    _b: number;
    /** @type {number} */
    _c: number;
    /** @type {number} */
    _d: number;
    /** @type {number} */
    _seed: number;
    /** @type {number | null} */
    _val: number | null;
    /** @type {number} seed */
    set seed(_seed: number);
    /**
     * Returns the seed of the random number generator.
     *
     * @returns {number} - The seed.
     */
    get seed(): number;
    /**
     * Returns a float between 0 and 1.
     *
     * @returns {number} - A random number between [0, 1]
     */
    get random(): number;
    /**
     * Returns an integer between 0 and MAX_INTEGER.
     *
     * @returns {number} - A random integer.
     */
    get random_int(): number;
    /**
     * Returns a normally distributed number with mean 0 and standard deviation 1.
     *
     * Uses the Marsaglia polar method, which yields two values per iteration; the spare is cached
     * and returned by the following call.
     *
     * @returns {number} A standard normal variate.
     */
    gauss_random(): number;
    /**
     * Returns `n` samples drawn from `A` without replacement.
     *
     * Uses a partial Fisher-Yates shuffle over a scratch index array, which runs in O(n): `A` itself
     * is never touched, and exactly `n` random values are consumed. Removing each picked index with
     * `Array.prototype.splice` instead would shift the tail on every draw and make the call O(n²) —
     * measurably worse from a few hundred elements upward, and 27× slower at n = 4000.
     *
     * @template T Returns samples from an input Matrix or Array.
     * @param {T[]} A - The input Matrix or Array.
     * @param {number} n - The number of samples.
     * @returns {T[]} A random selection form `A` of `n` samples.
     */
    choice<T>(A: T[], n: number): T[];
}

/**
 * In-place QuickSelect algorithm to partition an array around the k-th smallest element.
 * Runs in O(N) average time complexity (compared to O(N log N) for full Array.prototype.sort).
 *
 * After the call `arr[k]` holds the k-th smallest element, every element left of `k` compares
 * less than or equal to it, and every element right of `k` compares greater than or equal to it.
 *
 * @template T
 * @param {T[]} arr - Array to partition in-place
 * @param {Randomizer} randomizer - Seeded source of randomness for pivot selection.
 * @param {number} k - Target 0-indexed rank to select
 * @param {(a: T, b: T) => number} [compareFn] - Comparison function
 * @param {number} [start_left] - Start index (inclusive)
 * @param {number} [start_right] - End index (inclusive)
 * @returns {T} The k-th smallest element in the array
 */
declare function quickselect<T>(arr: T[], randomizer: Randomizer, k: number, compareFn?: (a: T, b: T) => number, start_left?: number, start_right?: number): T;
/**
 * QuickSelect by specific dimension axis for spatial trees (KDTree, BallTree).
 * Partitions array in-place along `arr[i].element[axis]`.
 *
 * @template {{ element: number[] | Float64Array }} T
 * @param {T[]} arr - Array to partition in-place
 * @param {Randomizer} randomizer - Seeded source of randomness for pivot selection.
 * @param {number} k - Target 0-indexed rank to select
 * @param {number} axis - Dimension coordinate index
 * @param {number} [start_left] - Start index (inclusive)
 * @param {number} [start_right] - End index (inclusive)
 * @returns {T} The k-th element along specified axis
 */
declare function quickselectByAxis<T extends {
    element: number[] | Float64Array;
}>(arr: T[], randomizer: Randomizer, k: number, axis: number, start_left?: number, start_right?: number): T;

/** @typedef {(i: number, j: number) => number} Accessor */
/**
 * @class
 * @category Matrix
 */
declare class Matrix {
    /**
     * Creates a Matrix out of `A`.
     * @param {Matrix | Float64Array[] | number[][]} A - The matrix, array, or number, which should converted to a Matrix.
     * @returns {Matrix}
     * @example
     * let A = Matrix.from([ [1, 0], [0, 1], ]); //creates a two by two identity matrix.
     */
    static from(A: Matrix | Float64Array[] | number[][]): Matrix;
    /**
     * Creates a Matrix with the diagonal being the values of `v`.
     *
     * @example let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     *
     * @param {number[] | Float64Array} v
     * @returns {Matrix}
     */
    static from_diag(v: number[] | Float64Array): Matrix;
    /**
     * Creates a Matrix with the diagonal being the values of `v`.
     *
     * @example let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     *
     * @param {number[] | Float64Array} v
     * @param {"col" | "row"} type
     * @returns {Matrix}
     */
    static from_vector(v: number[] | Float64Array, type: "col" | "row"): Matrix;
    /**
     * Solves the equation `Ax = b` using the conjugate gradient method. Returns the result `x`.
     *
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer | null} [randomizer]
     * @param {number} [tol=1e-3] Default is `1e-3`
     * @returns {Matrix}
     */
    static solve_CG(A: Matrix, b: Matrix, randomizer?: Randomizer | null, tol?: number): Matrix;
    /**
     * Solves the equation `Ax = b`. Returns the result `x`.
     *
     * @param {Matrix | { L: Matrix; U: Matrix }} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A: Matrix | {
        L: Matrix;
        U: Matrix;
    }, b: Matrix): Matrix;
    /**
     * `LU` decomposition of the Matrix `A`. Creates two matrices, so that the dot product `LU` equals `A`.
     *
     * @param {Matrix} A
     * @returns {{ L: Matrix; U: Matrix }} The left triangle matrix `L` and the upper triangle matrix `U`.
     */
    static LU(A: Matrix): {
        L: Matrix;
        U: Matrix;
    };
    /**
     * Computes the determinante of `A`, by using the `LU` decomposition of `A`.
     *
     * @param {Matrix} A
     * @returns {number} The determinate of the Matrix `A`.
     */
    static det(A: Matrix): number;
    /**
     * Computes the `k` components of the SVD decomposition of the matrix `M`.
     *
     * @param {Matrix} M
     * @param {number} [k=2] Default is `2`
     * @returns {{ U: Float64Array[]; Sigma: Float64Array; V: Float64Array[] }}
     */
    static SVD(M: Matrix, k?: number): {
        U: Float64Array[];
        Sigma: Float64Array;
        V: Float64Array[];
    };
    /**
     * @param {unknown} A
     * @returns {A is unknown[]|number[]|Float64Array|Float32Array}
     */
    static isArray(A: unknown): A is unknown[] | number[] | Float64Array | Float32Array;
    /**
     * @param {any[]} A
     * @returns {A is number[][]|Float64Array[]}
     */
    static is2dArray(A: any[]): A is number[][] | Float64Array[];
    /**
     * Creates a new Matrix. Entries are stored in a Float64Array.
     *
     * @example let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix. let B = new
     * Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
     *
     * @param {number} rows - The amount of rows of the matrix.
     * @param {number} cols - The amount of columns of the matrix.
     * @param {Accessor | string | number} value - Can be a function with row and col as parameters, a number, or
     *   "zeros", "identity" or "I", or "center".
     *
     *   - **function**: for each entry the function gets called with the parameters for the actual row and column.
     *   - **string**: allowed are
     *
     *       - "zero", creates a zero matrix.
     *       - "identity" or "I", creates an identity matrix.
     *       - "center", creates an center matrix.
     *   - **number**: create a matrix filled with the given value.
     */
    constructor(rows: number, cols: number, value?: Accessor | string | number);
    /** @type {number} */ _rows: number;
    /** @type {number} */ _cols: number;
    /** @type {Float64Array} */ _data: Float64Array;
    /**
     * Returns the `row`<sup>th</sup> row from the Matrix.
     *
     * @param {number} row
     * @returns {Float64Array}
     */
    row(row: number): Float64Array;
    /**
     * Returns an generator yielding each row of the Matrix.
     *
     * @yields {Float64Array}
     */
    iterate_rows(): Generator<Float64Array<ArrayBufferLike>, void, unknown>;
    /**
     * Sets the entries of `row`<sup>th</sup> row from the Matrix to the entries from `values`.
     *
     * @param {number} row
     * @param {number[]} values
     * @returns {Matrix}
     */
    set_row(row: number, values: number[]): Matrix;
    /**
     * Swaps the rows `row1` and `row2` of the Matrix.
     *
     * @param {number} row1
     * @param {number} row2
     * @returns {Matrix}
     */
    swap_rows(row1: number, row2: number): Matrix;
    /**
     * Returns the col<sup>th</sup> column from the Matrix.
     *
     * @param {number} col
     * @returns {Float64Array}
     */
    col(col: number): Float64Array;
    /**
     * Returns the `col`<sup>th</sup> entry from the `row`<sup>th</sup> row of the Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @returns {number}
     */
    entry(row: number, col: number): number;
    /**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given
     * {@link value}.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    set_entry(row: number, col: number, value: number): Matrix;
    /**
     * Adds a given {@link value} to the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the
     * Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    add_entry(row: number, col: number, value: number): Matrix;
    /**
     * Subtracts a given {@link value} from the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the
     * Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    sub_entry(row: number, col: number, value: number): Matrix;
    /**
     * Returns a new transposed Matrix.
     *
     * @returns {Matrix}
     */
    transpose(): Matrix;
    /**
     * Returns a new transposed Matrix. Short-form of `transpose`.
     *
     * @returns {Matrix}
     */
    get T(): Matrix;
    /**
     * Returns the inverse of the Matrix.
     *
     * @returns {Matrix}
     */
    inverse(): Matrix;
    /**
     * Returns the dot product. If `B` is an Array or Float64Array then an Array gets returned. If `B` is a Matrix then
     * a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    dot(B: Matrix | number[] | Float64Array): Matrix;
    /**
     * Transposes the current matrix and returns the dot product with `B`. If `B` is an Array or Float64Array then an
     * Array gets returned. If `B` is a Matrix then a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    transDot(B: Matrix | number[] | Float64Array): Matrix;
    /**
     * Returns the dot product with the transposed version of `B`. If `B` is an Array or Float64Array then an Array gets
     * returned. If `B` is a Matrix then a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    dotTrans(B: Matrix | number[] | Float64Array): Matrix;
    /**
     * Computes the outer product from `this` and `B`.
     *
     * @param {Matrix} B
     * @returns {Matrix}
     */
    outer(B: Matrix): Matrix;
    /**
     * Appends matrix `B` to the matrix.
     *
     * @example let A = Matrix.from([ [1, 1], [1, 1], ]); // 2 by 2 matrix filled with ones. let B = Matrix.from([ [2,
     * 2], [2, 2], ]); // 2 by 2 matrix filled with twos.
     *
     *     A.concat(B, "horizontal"); // 2 by 4 matrix. [[1, 1, 2, 2], [1, 1, 2, 2]]
     *     A.concat(B, "vertical"); // 4 by 2 matrix. [[1, 1], [1, 1], [2, 2], [2, 2]]
     *     A.concat(B, "diag"); // 4 by 4 matrix. [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]]
     *
     * @param {Matrix} B - Matrix to append.
     * @param {"horizontal" | "vertical" | "diag"} [type="horizontal"] - Type of concatenation. Default is
     *   `"horizontal"`
     * @returns {Matrix}
     */
    concat(B: Matrix, type?: "horizontal" | "vertical" | "diag"): Matrix;
    /**
     * Writes the entries of B in A at an offset position given by `offset_row` and `offset_col`.
     *
     * @param {number} offset_row
     * @param {number} offset_col
     * @param {Matrix} B
     * @returns {Matrix}
     */
    set_block(offset_row: number, offset_col: number, B: Matrix): Matrix;
    /**
     * Extracts the entries from the `start_row`<sup>th</sup> row to the `end_row`<sup>th</sup> row, the
     * `start_col`<sup>th</sup> column to the `end_col`<sup>th</sup> column of the matrix. If `end_row` or `end_col` is
     * empty, the respective value is set to `this.rows` or `this.cols`.
     *
     * @example let A = Matrix.from([ [1, 2, 3], [4, 5, 6], [7, 8, 9], ]); // a 3 by 3 matrix.
     *
     *     A.get_block(1, 1); // [[5, 6], [8, 9]]
     *     A.get_block(0, 0, 1, 1); // [[1]]
     *     A.get_block(1, 1, 2, 2); // [[5]]
     *     A.get_block(0, 0, 2, 2); // [[1, 2], [4, 5]]
     *
     * @param {number} start_row
     * @param {number} start_col
     * @param {number | null} [end_row]
     * @param {number | null} [end_col]
     * @returns {Matrix} Returns a `end_row` - `start_row` times `end_col` - `start_col` matrix, with respective entries
     *   from the matrix.
     */
    get_block(start_row: number, start_col: number, end_row?: number | null, end_col?: number | null): Matrix;
    /**
     * Returns a new array gathering entries defined by the indices given by argument.
     *
     * @param {number[]} row_indices - Array consists of indices of rows for gathering entries of this matrix
     * @param {number[]} col_indices - Array consists of indices of cols for gathering entries of this matrix
     * @returns {Matrix}
     */
    gather(row_indices: number[], col_indices: number[]): Matrix;
    /**
     * Applies a function to each entry of the matrix.
     *
     * @private
     * @param {(d: number, v: number) => number} f Function takes 2 parameters, the value of the actual entry and a
     *   value given by the function `v`. The result of `f` gets writen to the Matrix.
     * @param {Accessor} v Function takes 2 parameters for `row` and `col`, and returns a value witch should be applied
     *   to the `col`<sup>th</sup> entry of the `row`<sup>th</sup> row of the matrix.
     * @returns {Matrix}
     */
    private _apply_array;
    /**
     * @param {number[] | Float64Array} values
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
    _apply_rowwise_array(values: number[] | Float64Array, f: (d: number, v: number) => number): Matrix;
    /**
     * @param {number[] | Float64Array} values
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
    _apply_colwise_array(values: number[] | Float64Array, f: (d: number, v: number) => number): Matrix;
    /**
     * @param {Matrix | number[] | Float64Array | number} value
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
    _apply(value: Matrix | number[] | Float64Array | number, f: (d: number, v: number) => number): Matrix;
    /**
     * Clones the Matrix.
     *
     * @returns {Matrix}
     */
    clone(): Matrix;
    /**
     * Entrywise multiplication with `value`.
     *
     * @example let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;
     *
     *     A.mult(2); // [[2, 4], [6, 8]];
     *     A.mult(B); // [[1, 4], [9, 16]];
     *
     * @param {Matrix | Float64Array | number[] | number} value
     * @param {Object} [options]
     * @param {boolean} [options.inline=false] - If true, applies multiplication to the element, otherwise it creates
     *   first a copy and applies the multiplication on the copy. Default is `false`
     * @returns {Matrix}
     */
    mult(value: Matrix | Float64Array | number[] | number, { inline }?: {
        inline?: boolean | undefined;
    }): Matrix;
    /**
     * Entrywise division with `value`.
     *
     * @example let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;
     *
     *     A.divide(2); // [[0.5, 1], [1.5, 2]];
     *     A.divide(B); // [[1, 1], [1, 1]];
     *
     * @param {Matrix | Float64Array | number[] | number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline=false] - If true, applies division to the element, otherwise it creates first a
     *   copy and applies the division on the copy. Default is `false`
     * @returns {Matrix}
     */
    divide(value: Matrix | Float64Array | number[] | number, { inline }?: {
        inline?: boolean | undefined;
    }): Matrix;
    /**
     * Entrywise addition with `value`.
     *
     * @example let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;
     *
     *     A.add(2); // [[3, 4], [5, 6]];
     *     A.add(B); // [[2, 4], [6, 8]];
     *
     * @param {Matrix | Float64Array | number[] | number} value
     * @param {Object} [options]
     * @param {boolean} [options.inline=false] - If true, applies addition to the element, otherwise it creates first a
     *   copy and applies the addition on the copy. Default is `false`
     * @returns {Matrix}
     */
    add(value: Matrix | Float64Array | number[] | number, { inline }?: {
        inline?: boolean | undefined;
    }): Matrix;
    /**
     * Entrywise subtraction with `value`.
     *
     * @example let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;
     *
     *     A.sub(2); // [[-1, 0], [1, 2]];
     *     A.sub(B); // [[0, 0], [0, 0]];
     *
     * @param {Matrix | Float64Array | number[] | number} value
     * @param {Object} [options]
     * @param {boolean} [options.inline=false] - If true, applies subtraction to the element, otherwise it creates first
     *   a copy and applies the subtraction on the copy. Default is `false`
     * @returns {Matrix}
     */
    sub(value: Matrix | Float64Array | number[] | number, { inline }?: {
        inline?: boolean | undefined;
    }): Matrix;
    /**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     *
     * @param {[number, number, Accessor]} parameter - Takes an Array in the form [rows, cols, value], where rows and
     *   cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row
     *   and col) which has to return a value for the colth entry of the rowth row.
     */
    set shape([rows, cols, value]: [number, number, Accessor]);
    /**
     * Returns the number of rows and columns of the Matrix.
     *
     * @returns {number[]} An Array in the form [rows, columns].
     */
    get shape(): number[];
    /**
     * Returns the Matrix as a Array of Float64Arrays.
     *
     * @returns {Float64Array[]}
     */
    to2dArray(): Float64Array[];
    /**
     * Returns the Matrix as a Array of Arrays.
     *
     * @returns {number[][]}
     */
    asArray(): number[][];
    /**
     * Returns the diagonal of the Matrix.
     *
     * @returns {Float64Array}
     */
    diag(): Float64Array;
    /**
     * Returns the mean of all entries of the Matrix.
     *
     * @returns {number}
     */
    mean(): number;
    /**
     * Returns the sum oof all entries of the Matrix.
     *
     * @returns {number}
     */
    sum(): number;
    /**
     * Returns the entries of the Matrix.
     *
     * @returns {Float64Array}
     */
    get values(): Float64Array;
    /**
     * Returns the mean of each row of the matrix.
     *
     * @returns {Float64Array}
     */
    meanRows(): Float64Array;
    /**
     * Returns the mean of each column of the matrix.
     *
     * @returns {Float64Array}
     */
    meanCols(): Float64Array;
    /**
     * Makes a `Matrix` object an iterable object.
     *
     * @yields {Float64Array}
     */
    [Symbol.iterator](): Generator<Float64Array<ArrayBufferLike>, void, unknown>;
}
type Accessor = (i: number, j: number) => number;

/** @import { Metric } from "../metrics/index.js" */
/**
 * Computes the norm of a vector, by computing its distance to **0**.
 *
 * @category Matrix
 * @param {Matrix | number[] | Float64Array} v - Vector.
 * @param {Metric} [metric=euclidean] - Which metric should be used to compute the norm. Default is `euclidean`
 * @returns {number} - The norm of `v`.
 */
declare function norm(v: Matrix | number[] | Float64Array, metric?: Metric): number;

/** @import { Metric } from "../metrics/index.js" */
/**
 * Normalizes Vector `v`.
 *
 * @category Matrix
 * @param {number[] | Float64Array} v - Vector
 * @param {Metric} metric
 * @returns {number[] | Float64Array} - The normalized vector with length 1.
 */
declare function normalize(v: number[] | Float64Array, metric?: Metric): number[] | Float64Array;

/** @import {InputType} from "../index.js" */
/**
 * Base class for all clustering algorithms.
 * @template Para
 */
declare class Clustering<Para> {
    /**
     * Compute the respective Clustering with given parameters
     * @param {InputType} points
     * @param {Para} parameters
     */
    constructor(points: InputType, parameters: Para);
    /** @type {InputType} */
    _points: InputType;
    /** @type {Para} */
    _parameters: Para;
    /** @type {Matrix} */
    _matrix: Matrix;
    /** @type {number} */
    _N: number;
    /** @type {number} */
    _D: number;
    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[][]} An array with the indices of the clusters.
     */
    get_clusters(...args: unknown[]): number[][];
    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[]} An array with the clusters id's for each point.
     */
    get_cluster_list(...args: unknown[]): number[];
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersCURE } from "./index.js" */
/**
 * CURE (Clustering Using REpresentatives)
 *
 * An efficient clustering algorithm for large databases that is robust to outliers
 * and identifies clusters with non-spherical shapes and wide variances in size.
 *
 * @class
 * @extends Clustering<ParametersCURE>
 * @category Clustering
 */
declare class CURE extends Clustering<ParametersCURE> {
    /**
     * @param {InputType} points
     * @param {Partial<ParametersCURE>} parameters
     */
    constructor(points: InputType, parameters?: Partial<ParametersCURE>);
    /** @type {number} */
    _K: number;
    /** @type {number} */
    _num_representatives: number;
    /** @type {number} */
    _shrink_factor: number;
    /**
     * @private
     * @type {CURECluster[]}
     */
    private _clusters;
    /** @type {number[]} */
    _cluster_ids: number[];
    /**
     * Initialize each point as its own cluster
     * @private
     */
    private _initialize_clusters;
    /**
     * Compute distance between two clusters using representative points
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {number}
     */
    private _cluster_distance;
    /**
     * Find the closest pair of clusters
     * @private
     * @returns {[number, number, number]} [index1, index2, distance]
     */
    private _find_closest_clusters;
    /**
     * Merge two clusters
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {CURECluster}
     */
    private _merge_clusters;
    /**
     * Run CURE clustering algorithm
     * @private
     */
    private _cure;
    /**
     * Build the cluster list (point -> cluster assignment)
     * @private
     */
    private _build_cluster_ids;
    /**
     * @returns {number[][]}
     */
    get_clusters(): number[][];
    /**
     * @returns {number[]}
     */
    get_cluster_list(): number[];
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersHierarchicalClustering } from "./index.js" */
/**
 * Hierarchical Clustering
 *
 * A bottom-up approach (agglomerative) to clustering that builds a tree of clusters (dendrogram).
 * Supports different linkage criteria: single, complete, average, and ward.
 *
 * @class
 * @extends Clustering<ParametersHierarchicalClustering>
 * @category Clustering
 */
declare class HierarchicalClustering extends Clustering<ParametersHierarchicalClustering> {
    /**
     * @param {InputType} points - Data or distance matrix if metric is 'precomputed'
     * @param {Partial<ParametersHierarchicalClustering>} parameters
     */
    constructor(points: InputType, parameters?: Partial<ParametersHierarchicalClustering>);
    /** @type {Cluster | null} */
    root: Cluster | null;
    _id: number;
    _d_min: Float64Array<ArrayBuffer>;
    _distance_matrix: Matrix;
    _clusters: any[];
    _c_size: Uint16Array<ArrayBuffer>;
    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {Cluster[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters_raw(value: number, type?: "distance" | "depth"): Cluster[][];
    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters(value: number, type?: "distance" | "depth"): number[][];
    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[]} - Array of clusters with the indices of the rows in given points.
     */
    get_cluster_list(value: number, type?: "distance" | "depth"): number[];
    /**
     * @private
     * @param {Cluster} node
     * @param {(d: {dist: number, depth: number}) => number} f
     * @param {number} value
     * @param {Cluster[][]} result
     */
    private _traverse;
}

/** @private */
declare class Cluster {
    /**
     *
     * @param {number} id
     * @param {Cluster?} left
     * @param {Cluster?} right
     * @param {number} dist
     * @param {Float64Array?} centroid
     * @param {number} index
     * @param {number} [size]
     * @param {number} [depth]
     */
    constructor(id: number, left: Cluster | null, right: Cluster | null, dist: number, centroid: Float64Array | null, index: number, size?: number, depth?: number);
    /**@type {number} */
    size: number;
    /**@type {number} */
    depth: number;
    /**@type {Cluster | null} */
    parent: Cluster | null;
    id: number;
    left: Cluster | null;
    right: Cluster | null;
    dist: number;
    index: number;
    centroid: Float64Array<ArrayBufferLike>;
    /**
     *
     * @param {Cluster} left
     * @param {Cluster} right
     * @returns {Float64Array}
     */
    _calculate_centroid(left: Cluster, right: Cluster): Float64Array;
    get isLeaf(): boolean;
    /**
     *
     * @returns {Cluster[]}
     */
    leaves(): Cluster[];
    /**
     *
     * @returns {Cluster[]}
     */
    descendants(): Cluster[];
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersKMeans } from "./index.js" */
/**
 * K-Means Clustering
 *
 * A popular clustering algorithm that partitions data into K clusters where each point
 * belongs to the cluster with the nearest mean (centroid).
 *
 * @class
 * @extends Clustering<ParametersKMeans>
 * @category Clustering
 * @see {@link KMedoids} for a more robust alternative
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const points = [[1, 1], [1.5, 1.5], [5, 5], [5.5, 5.5]];
 * const kmeans = new druid.KMeans(points, { K: 2 });
 *
 * const clusters = kmeans.get_cluster_list(); // [0, 0, 1, 1]
 * const centroids = kmeans.centroids; // center points
 */
declare class KMeans extends Clustering<ParametersKMeans> {
    /**
     * @param {InputType} points
     * @param {Partial<ParametersKMeans>} parameters
     */
    constructor(points: InputType, parameters?: Partial<ParametersKMeans>);
    _K: number;
    _randomizer: Randomizer;
    /** @type {number[]} */
    _clusters: number[];
    _cluster_centroids: Float64Array<ArrayBufferLike>[];
    /** @returns {number} The number of clusters */
    get k(): number;
    /** @returns {Float64Array[]} The cluster centroids */
    get centroids(): Float64Array[];
    /** @returns {number[]} The cluster list */
    get_cluster_list(): number[];
    /** @returns {number[][]} An Array of clusters with the indices of the points. */
    get_clusters(): number[][];
    /**
     * @private
     * @param {number[]} point_indices
     * @param {number[]} candidates
     * @returns {number}
     */
    private _furthest_point;
    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    private _get_random_centroids;
    /**
     * @private
     * @param {Float64Array[]} cluster_centroids
     * @returns {{ clusters_changed: boolean; cluster_centroids: Float64Array[] }}
     */
    private _iteration;
    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    private _compute_centroid;
}

/** @import {InputType} from "../index.js" */
/** @import { ParametersKMedoids } from "./index.js" */
/**
 * K-Medoids (PAM - Partitioning Around Medoids)
 *
 * A robust clustering algorithm similar to K-Means, but uses actual data points (medoids)
 * as cluster centers and can work with any distance metric.
 *
 * @class
 * @extends Clustering<ParametersKMedoids>
 * @category Clustering
 * @see {@link KMeans} for a faster but less robust alternative
 */
declare class KMedoids extends Clustering<ParametersKMedoids> {
    /**
     * @param {InputType} points - Data matrix
     * @param {Partial<ParametersKMedoids>} parameters
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(points: InputType, parameters?: Partial<ParametersKMedoids>);
    _A: Float64Array<ArrayBufferLike>[];
    _max_iter: number;
    _distance_matrix: Matrix;
    _randomizer: Randomizer;
    _clusters: any[];
    _cluster_medoids: number[];
    _is_initialized: boolean;
    /** @returns {number[]} The cluster list */
    get_cluster_list(): number[];
    /** @returns {number[][]} - Array of clusters with the indices of the rows in given points. */
    get_clusters(): number[][];
    /** @returns {number} */
    get k(): number;
    /** @returns {number[]} */
    get medoids(): number[];
    /** @returns {number[]} */
    get_medoids(): number[];
    generator(): AsyncGenerator<number[][], void, unknown>;
    /** Algorithm 1. FastPAM1: Improved SWAP algorithm */
    /**
     * FastPAM1: One best swap per iteration
     * @private
     * @returns {boolean}
     */
    private _iteration;
    /**
     * @private
     * Get distance between two points
     * @param {number} i
     * @param {number} j
     * @param {Float64Array?} x_i
     * @param {Float64Array?} x_j
     * @returns {number}
     */
    private _get_distance;
    /**
     * @private
     * @param {Float64Array} x_j
     * @param {number} j
     * @returns
     */
    private _nearest_medoid;
    /**
     * @private
     */
    private _update_clusters;
    /**
     * Computes `K` clusters out of the `matrix`.
     * @param {number} K - Number of clusters.
     * @param {number[]} cluster_medoids
     */
    init(K: number, cluster_medoids: number[]): this;
    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @private
     * @param {number} K - Number of clusters
     * @returns {number[]}
     */
    private _get_random_medoids;
}

/** @import { ParametersMeanShift } from "./index.js" */
/** @import { InputType } from "../index.js" */
/**
 * Mean Shift Clustering
 *
 * A non-parametric clustering technique that does not require prior knowledge of the
 * number of clusters. It identifies centers of density in the data.
 *
 * @class
 * @extends Clustering<ParametersMeanShift>
 * @category Clustering
 */
declare class MeanShift extends Clustering<ParametersMeanShift> {
    /**
     *
     * @param {InputType} points
     * @param {Partial<ParametersMeanShift>} parameters
     */
    constructor(points: InputType, parameters?: Partial<ParametersMeanShift>);
    /**
     * @private
     * @type {number}
     */
    private _bandwidth;
    /**
     * @private
     * @type {number}
     */
    private _max_iter;
    /**
     * @private
     * @type {number}
     */
    private _tolerance;
    /**
     * @private
     * @type {(dist: number) => number}
     */
    private _kernel;
    /**
     * @type {Matrix}
     */
    _points: Matrix;
    /**
     * @private
     * @type {number[] | undefined}
     */
    private _clusters;
    /**
     * @private
     * @type {number[][] | undefined}
     */
    private _cluster_list;
    /**
     * Helper to compute bandwidth if not provided
     * @private
     * @param {Matrix} matrix
     * @returns {number}
     */
    private _compute_bandwidth;
    /**
     * Compute kernel weight
     * @private
     * @param {number} dist
     * @returns {number}
     */
    private _kernel_weight;
    /**
     * Runs one mean shift step across the worker pool, updating `_points` in place.
     *
     * @private
     * @param {number} N
     * @param {number} D
     * @param {boolean} use_gaussian
     * @returns {number | null} The largest shift of the step, or null if the pool did not run it.
     */
    private _mean_shift_parallel;
    /**
     * Perform mean shift iterations
     * @private
     */
    private _mean_shift;
    /**
     * After convergence, assign clusters based on nearest mode
     * @private
     */
    private _assign_clusters;
    /**
     * @returns {number[][]}
     */
    get_clusters(): number[][];
    /**
     *
     * @returns {number[]}
     */
    get_cluster_list(): number[];
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersOptics } from "./index.js" */
/** @typedef {Object} DBEntry
 * @property {Float64Array} element
 * @property {number} index
 * @property {number} [reachability_distance]
 * @property {boolean} processed
 * @property {DBEntry[]} [neighbors]
 */
/**
 * OPTICS (Ordering Points To Identify the Clustering Structure)
 *
 * A density-based clustering algorithm that extends DBSCAN. It handles clusters of varying
 * densities and produces a reachability plot that can be used to extract clusters.
 *
 * @class
 * @extends Clustering<ParametersOptics>
 * @category Clustering
 */
declare class OPTICS extends Clustering<ParametersOptics> {
    /**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     *
     * @param {InputType} points - The data.
     * @param {Partial<ParametersOptics>} [parameters={}]
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
    constructor(points: InputType, parameters?: Partial<ParametersOptics>);
    /**
     * Feeds the quickselect used to find core distances, so that seed alone determines the run.
     *
     * @private
     * @type {Randomizer}
     */
    private _randomizer;
    /**
     * @private
     * @type {DBEntry[]}
     */
    private _ordered_list;
    /** @type {number[][]} */
    _clusters: number[][];
    /**
     * @private
     * @type {DBEntry[]}
     */
    private _DB;
    _cluster_index: number;
    /**
     * @private
     * @param {DBEntry} p - A point of the data.
     * @returns {DBEntry[]} An array consisting of the `epsilon`-neighborhood of `p`.
     */
    private _get_neighbors;
    /**
     * @private
     * @param {DBEntry} p - A point of `matrix`.
     * @returns {number|undefined} The distance to the `min_points`-th nearest point of `p`, or undefined if the
     *   `epsilon`-neighborhood has fewer elements than `min_points`.
     */
    private _core_distance;
    /**
     * Updates the reachability distance of the points.
     *
     * @private
     * @param {DBEntry} p
     * @param {Heap<DBEntry>} seeds
     */
    private _update;
    /**
     * Expands the `cluster` with points in `seeds`.
     *
     * @private
     * @param {Heap<DBEntry>} seeds
     * @param {number[]} cluster
     */
    private _expand_cluster;
    /**
     * Returns an array of clusters.
     *
     * @returns {number[][]} Array of clusters with the indices of the rows in given `matrix`.
     */
    get_clusters(): number[][];
    /**
     * @returns {number[]} Returns an array, where the ith entry defines the cluster affirmation of the ith point of
     *   given data. (-1 stands for outlier)
     */
    get_cluster_list(): number[];
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersXMeans } from "./index.js" */
/**
 * @typedef SplitResult
 * @property {number} index - Index of the cluster being split
 * @property {number} bic_parent - BIC score of the parent cluster
 * @property {number} bic_children - BIC score of the split children
 * @property {number[][]} child_clusters - Clusters after splitting
 * @property {Float64Array[]} child_centroids - Centroids of child clusters
 */
/**
 * @typedef CandidateResult
 * @property {KMeans} kmeans - The KMeans instance for this K
 * @property {number} score - BIC score
 */
/**
 * X-Means Clustering
 *
 * An extension of K-Means that automatically determines the number of clusters (K)
 * using the Bayesian Information Criterion (BIC).
 *
 * @class
 * @extends Clustering<ParametersXMeans>
 * @category Clustering
 */
declare class XMeans extends Clustering<ParametersXMeans> {
    /**
     * XMeans clustering algorithm that automatically determines the optimal number of clusters.
     *
     * X-Means extends K-Means by starting with a minimum number of clusters and iteratively
     * splitting clusters to improve the Bayesian Information Criterion (BIC).
     *
     * Algorithm:
     * 1. Start with K_min clusters using KMeans
     * 2. For each cluster, try splitting it into 2 sub-clusters
     * 3. If BIC improves after splitting, keep the split
     * 4. Run KMeans again with all (old + new) centroids
     * 5. Repeat until K_max is reached or no more improvements
     *
     * @param {InputType} points - The data points to cluster
     * @param {Partial<ParametersXMeans>} [parameters={}] - Configuration parameters
     * @see {@link https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf}
     * @see {@link https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py}
     * @see {@link https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java}
     */
    constructor(points: InputType, parameters?: Partial<ParametersXMeans>);
    _randomizer: Randomizer;
    /** @type {KMeans | null} */
    _best_kmeans: KMeans | null;
    /**
     * Run the XMeans algorithm
     *
     * @private
     */
    private _run;
    /**
     * Select the best candidate based on BIC score
     *
     * @private
     * @param {Map<number, CandidateResult>} candidates
     * @returns {KMeans}
     */
    private _select_best_candidate;
    /**
     * Calculate Bayesian Information Criterion for a set of clusters.
     *
     * Uses Kass's formula for BIC calculation:
     * BIC(θ) = L(D) - 0.5 * p * ln(N)
     *
     * Where:
     * - L(D) is the log-likelihood of the data
     * - p is the number of free parameters: (K-1) + D*K + 1
     * - N is the total number of points
     *
     * @private
     * @param {number[][]} clusters - Array of clusters with point indices
     * @param {Float64Array[]} centroids - Array of centroids
     * @returns {number} BIC score (higher is better)
     */
    private _bic;
    /**
     * Get the computed clusters
     *
     * @returns {number[][]} Array of clusters, each containing indices of points
     */
    get_clusters(): number[][];
    /** @returns {number[]} The cluster list */
    get_cluster_list(): number[];
    /**
     * Get the final centroids
     *
     * @returns {Float64Array[]} Array of centroids
     */
    get centroids(): Float64Array[];
    /**
     * Get the optimal number of clusters found
     *
     * @returns {number} The number of clusters
     */
    get k(): number;
}

type ParametersHierarchicalClustering = {
    linkage: "single" | "complete" | "average" | "ward";
    /**
     * - `"ward"` assumes these are Euclidean distances; its
     * minimum-variance criterion is only meaningful in that geometry.
     */
    metric: Metric | "precomputed";
};
type ParametersKMeans = {
    K: number;
    /**
     * Default is `euclidean`
     */
    metric: Metric;
    /**
     * Default is `1212`
     */
    seed: number;
    /**
     * - Initial centroids. Default is `null`
     */
    initial_centroids?: Float64Array<ArrayBufferLike>[] | number[][] | undefined;
};
type ParametersKMedoids = {
    /**
     * - Number of clusters
     */
    K: number;
    /**
     * - Maximum number of iterations. Default is 10 * Math.log10(N). Default is `null`
     */
    max_iter: number | null;
    /**
     * - Metric defining the dissimilarity. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Seed value for random number generator. Default is `1212`
     */
    seed: number;
};
type ParametersOptics = {
    /**
     * - The minimum distance which defines whether a point is a neighbor or not.
     */
    epsilon: number;
    /**
     * - The minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
     */
    min_points: number;
    /**
     * - The distance metric which defines the distance between two points of the points. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Seed for the random number generator used to select core distances. Default is `1212`
     */
    seed: number;
};
type ParametersXMeans = {
    /**
     * - Minimum number of clusters. Default is `2`
     */
    K_min: number;
    /**
     * - Maximum number of clusters. Default is `10`
     */
    K_max: number;
    /**
     * - Distance metric function. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Random seed. Default is `1212`
     */
    seed: number;
    /**
     * - Minimum points required to consider splitting a cluster. Default is `25`
     */
    min_cluster_size: number;
    /**
     * - Convergence tolerance for KMeans. Default is `0.001`
     */
    tolerance: number;
};
type ParametersMeanShift = {
    /**
     * - bandwidth
     */
    bandwidth: number;
    /**
     * - Metric defining the dissimilarity. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Seed value for random number generator. Default is `1212`
     */
    seed: number;
    /**
     * - Kernel function. Default is `gaussian`
     */
    kernel: "flat" | "gaussian" | ((dist: number) => number);
    /**
     * - Maximum number of iterations. Default is `Math.max(10, Math.floor(10 * Math.log10(N)))`
     */
    max_iter?: number | undefined;
    /**
     * - Convergence tolerance. Default is `1e-3`
     */
    tolerance?: number | undefined;
};
type ParametersCURE = {
    /**
     * - Target number of clusters. Default is `2`
     */
    K: number;
    /**
     * - Number of representative points per cluster. Default is `5`
     */
    num_representatives: number;
    /**
     * - Factor to shrink representatives toward centroid (0-1). Default is `0.5`
     */
    shrink_factor: number;
    /**
     * - Distance metric function. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Random seed. Default is `1212`
     */
    seed: number;
};

/**
 * @template T
 * @typedef {Object} DisjointSetPayload
 * @property {T} parent
 * @property {Set<T>} children
 * @property {number} size
 */
/**
 * @template T
 * @class
 * @category Data Structures
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */
declare class DisjointSet<T> {
    /**
     * @param {T[]?} elements
     */
    constructor(elements?: T[] | null);
    /**
     * @private
     * @type {Map<T, DisjointSetPayload<T>>}
     */
    private _list;
    /**
     * @private
     * @param {T} x
     * @returns {DisjointSet<T>}
     */
    private make_set;
    /**
     * @param {T} x
     * @returns
     */
    find(x: T): T | null;
    /**
     * @param {T} x
     * @param {T} y
     * @returns
     */
    union(x: T, y: T): this;
    /** @param {T} x */
    get_children(x: T): Set<T> | null;
}

/** @import { Comparator } from "./index.js" */
/**
 * @template T
 * @class
 * @category Data Structures
 */
declare class Heap<T> {
    /**
     * Creates a Heap from an Array
     *
     * @template T
     * @param {T[]} elements - Contains the elements for the Heap.
     * @param {(d: T) => number} accessor - Function returns the value of the element.
     * @param {"min" | "max" | Comparator} [comparator="min"] - Function returning true or false
     *   defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
     *   Max_heap). Default is `"min"`
     * @returns {Heap<T>}
     */
    static heapify<T_1>(elements: T_1[], accessor: (d: T_1) => number, comparator?: "min" | "max" | Comparator): Heap<T_1>;
    /**
     * A heap is a datastructure holding its elements in a specific way, so that the top element would be the first
     * entry of an ordered list.
     *
     * @param {T[]?} elements - Contains the elements for the Heap. `elements` can be null.
     * @param {(d: T) => number} accessor - Function returns the value of the element.
     * @param {"min" | "max" | Comparator} [comparator="min"] - Function returning true or false
     *   defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
     *   Max_heap). Default is `"min"`
     * @see {@link https://en.wikipedia.org/wiki/Binary_heap}
     */
    constructor(elements: (T[] | null) | undefined, accessor: (d: T) => number, comparator?: "min" | "max" | Comparator);
    /** @type {{ element: T; value: number }[]} */
    _container: {
        element: T;
        value: number;
    }[];
    /** @type {Comparator} */
    _comparator: Comparator;
    /** @type {(d: T) => number} */
    _accessor: (d: T) => number;
    /**
     * Swaps elements of container array.
     *
     * @private
     * @param {number} index_a
     * @param {number} index_b
     */
    private _swap;
    /** @private */
    private _heapify_up;
    /**
     * Pushes the element to the heap.
     *
     * @param {T} element
     * @returns {Heap<T>}
     */
    push(element: T): Heap<T>;
    /**
     * Non-recursive sift-down implementation.
     *
     * @private
     * @param {number} [start_index=0] Default is `0`
     */
    private _heapify_down;
    /**
     * Removes and returns the top entry of the heap.
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`}).
     */
    pop(): {
        element: T;
        value: number;
    } | null;
    /**
     * Returns the top entry of the heap without removing it.
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`).
     */
    get first(): {
        element: T;
        value: number;
    } | null;
    /**
     * Yields the raw data
     *
     * @yields {T} Object consists of the element and its value (computed by `accessor`}).
     */
    iterate(): Generator<T, void, unknown>;
    /**
     * Returns the heap as ordered array.
     *
     * @returns {T[]} Array consisting the elements ordered by `comparator`.
     */
    toArray(): T[];
    /**
     * Returns elements of container array.
     *
     * @returns {T[]} Array consisting the elements.
     */
    data(): T[];
    /**
     * Returns the container array.
     *
     * @returns {{ element: T; value: number }[]} The container array.
     */
    raw_data(): {
        element: T;
        value: number;
    }[];
    /**
     * The size of the heap.
     *
     * @returns {number}
     */
    get length(): number;
    /**
     * Returns false if the the heap has entries, true if the heap has no entries.
     *
     * @returns {boolean}
     */
    get empty(): boolean;
}

type Comparator = (a: number, b: number) => boolean;

/** @import {InputType} from "../index.js" */
/**
 * @abstract
 * @template {InputType} T
 * @template {{ seed?: number }} Para
 *
 * Base class for all Dimensionality Reduction (DR) algorithms.
 *
 * Provides a common interface for parameters management, data initialization,
 * and transformation (both synchronous and asynchronous).
 *
 * Randomness always comes from a seeded {@link Randomizer}, never from `Math.random`, so the same
 * `seed` reproduces a result exactly within one engine and library build. Across engines it does
 * not: ECMAScript leaves `Math.pow`, `Math.exp` and `Math.log` *implementation-approximated*, so
 * two engines — or a WASM kernel and its JS fallback (see {@link setWasmEnabled}) — may differ in
 * the last bit. That is invisible for most methods, but {@link UMAP} and {@link SAMMON} are chaotic
 * enough to turn it into a visibly different, equally valid layout.
 *
 * @class
 */
declare class DR<T extends InputType, Para extends {
    seed?: number;
}> {
    /**
     * Computes the projection.
     *
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {T} The dimensionality reduced dataset.
     */
    static transform<T_1 extends InputType, Para_1 extends {
        seed?: number;
    }>(X: T_1, parameters: Para_1, ...args: unknown[]): T_1;
    /**
     * Computes the projection.
     *
     * @template {{ seed?: number }} Para
     * @param {InputType} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator<InputType, InputType, void>} A generator yielding the intermediate steps of the dimensionality
     *   reduction method.
     */
    static generator<Para_1 extends {
        seed?: number;
    }>(X: InputType, parameters: Para_1, ...args: unknown[]): Generator<InputType, InputType, void>;
    /**
     * Computes the projection.
     *
     * @template {{ seed?: number }} Para
     * @param {InputType} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise<X>} A promise yielding the dimensionality reduced dataset.
     */
    static transform_async<Para_1 extends {
        seed?: number;
    }>(X: InputType, parameters: Para_1, ...args: unknown[]): Promise<InputType>;
    /**
     * Takes the default parameters and seals them, remembers the type of input `X`, and initializes the random number
     * generator.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Para} default_parameters - Object containing default parameterization of the DR method.
     * @param {Partial<Para>} parameters - Object containing parameterization of the DR method to override defaults.
     */
    constructor(X: T, default_parameters: Para, parameters?: Partial<Para>);
    /** @type {number} */
    _D: number;
    /** @type {number} */
    _N: number;
    /** @type {Randomizer} */
    _randomizer: Randomizer;
    /** @type {boolean} */
    _is_initialized: boolean;
    /** @type {T} */
    __input: T;
    /** @type {Para} */
    _parameters: Para;
    /** @type {"array" | "matrix" | "typed"} */
    _type: "array" | "matrix" | "typed";
    /** @type {Matrix} */
    X: Matrix;
    /** @type {Matrix} */
    Y: Matrix;
    /**
     * Get all Parameters.
     * @overload
     * @returns {Para}
     */
    parameter(): Para;
    /**
     * Get value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @returns {Para[K]}
     */
    parameter<K extends keyof Para>(name: K): Para[K];
    /**
     * Set value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @param {Para[K]} value - Value of the parameter to set.
     * @returns {this}
     */
    parameter<K extends keyof Para>(name: K, value: Para[K]): this;
    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {T} The projection.
     */
    transform(...args: unknown[]): T;
    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {Generator<T, T, void>} The intermediate steps of the projection.
     */
    generator(...args: unknown[]): Generator<T, T, void>;
    /**
     * @abstract
     * @param {...unknown} args
     */
    init(...args: unknown[]): void;
    /**
     * If the respective DR method has an `init` function, call it before `transform`.
     *
     * @returns {DR<T, Para>}
     */
    check_init(): DR<T, Para>;
    /** @returns {T} The projection in the type of input `X`. */
    get projection(): T;
    /**
     * Computes the projection.
     *
     * @param {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<T>} The dimensionality reduced dataset.
     */
    transform_async(...args: unknown[]): Promise<T>;
    /**
     * WASM buffer sessions this method keeps alive between iterations, released when a run ends.
     *
     * Empty for methods with no accelerated iteration step, which is most of them. The accelerated
     * optimisers override it — see {@link TSNE}.
     *
     * @protected
     * @returns {string[]}
     */
    protected get _wasm_session_keys(): string[];
    /**
     * Releases the WASM buffers this run is holding.
     *
     * Called from a `finally` around every iteration loop, so the memory does not outlive the
     * projection. Freeing early is never a correctness problem — the next call reallocates — so
     * this is safe to call at any point, including when WASM never ran.
     *
     * @protected
     * @returns {void}
     */
    protected _release_wasm(): void;
    /**
     * Hands back the WASM buffers this instance is holding.
     *
     * Only needed after driving `generator()` by hand and stopping early — a plain `transform()`,
     * or a `for…of` over `generator()` (including one you `break`), already releases when it ends.
     * It frees only this method's buffers, never another running instance's, and the next run simply
     * reallocates, so it is safe to call at any time, more than once, and while other instances are
     * mid-run.
     *
     * @returns {this}
     * @example
     * const tsne = new TSNE(X, { d: 2 });
     * const steps = tsne.generator(500);
     * steps.next();
     * tsne.release(); // stop early and give the buffers back
     */
    release(): this;
    /**
     * Alias of {@link release} for the `using` declaration, so a hand-driven run frees its buffers
     * when the block exits. Note that {@link transform} and a completed or `break`-ed `generator()`
     * already release on their own — this only matters for a generator abandoned part way.
     *
     * ```js
     * using tsne = new TSNE(X, { d: 2 });
     * const steps = tsne.generator(500);
     * steps.next(); // buffers released when the enclosing block exits
     * ```
     *
     * @returns {void}
     */
    [Symbol.dispose](): void;
}

/** @import { InputType } from "../index.js" */
/** @import { ParametersFASTMAP } from "./index.js"; */
/**
 * FastMap algorithm for dimensionality reduction.
 *
 * A very fast algorithm for projecting high-dimensional data into a lower-dimensional
 * space while preserving pairwise distances. It works similarly to PCA but uses
 * only a subset of the data to find projection axes.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersFASTMAP>
 * @category Dimensionality Reduction
 */
declare class FASTMAP<T extends InputType> extends DR<T, ParametersFASTMAP> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersFASTMAP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersFASTMAP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersFASTMAP>): Promise<T_1>;
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets.
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersFASTMAP>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X: T, parameters: Partial<ParametersFASTMAP>);
    /**
     * Chooses two points which are the most distant in the actual projection.
     *
     * @private
     * @param {(a: number, b: number) => number} dist
     * @returns {[number, number, number]} An array consisting of first index, second index, and distance between the
     *   two points.
     */
    private _choose_distant_objects;
    /**
     * Computes the projection.
     *
     * @returns {T} The `d`-dimensional projection of the data matrix `X`.
     */
    transform(): T;
    generator(): Generator<T, T, unknown>;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersISOMAP} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Isomap (Isometric Mapping)
 *
 * A nonlinear dimensionality reduction algorithm that uses geodesic distances
 * between points on a manifold to perform embedding. It builds a neighborhood
 * graph and uses MDS on the shortest-path distances.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersISOMAP>
 * @category Dimensionality Reduction
 * @see {@link LLE} for another nonlinear alternative
 */
declare class ISOMAP<T extends InputType> extends DR<T, ParametersISOMAP> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersISOMAP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersISOMAP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersISOMAP>): Promise<T_1>;
    /**
     * Isometric feature mapping (ISOMAP).
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersISOMAP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X: T, parameters?: Partial<ParametersISOMAP>);
    defaults: ParametersISOMAP;
    /**
     * Runs the all-pairs geodesic step across the worker pool.
     *
     * The pool needs a `SharedArrayBuffer` to write into, so the result is copied into `out` rather
     * than computed in place. That copy is far cheaper than the Dijkstra runs it parallelizes.
     *
     * @private
     * @param {Int32Array} neighbors
     * @param {Float64Array} distances
     * @param {Float64Array} out
     * @param {number} rows
     * @param {number} maxK
     * @returns {boolean} True if the pool produced the distances.
     */
    private _dijkstra_parallel;
    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * @returns {T}
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersStressMDS, WeightSpec} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/** @typedef {"MDS" | "PCA" | "random"} AvailableInit */
/** Raw stress: every pair counts equally. Recovers the objective {@link SMACOF} minimises. */
declare const WEIGHTS_UNIFORM: 0;
/** Sammon stress. Recovers the objective {@link SAMMON} minimises, up to a constant factor. */
declare const WEIGHTS_SAMMON: -1;
/** Elastic scaling, also known as the Kamada-Kawai energy. See {@link KKMDS}. */
declare const WEIGHTS_ELASTIC: -2;
/**
 * Weighted metric MDS (stress majorization family)
 *
 * Minimises `σ(Y) = Σ_{i<j} w_ij ⋅ (‖y_i - y_j‖ - d_ij)²` for a weighting you choose. An exponent
 * `q` gives `w_ij = d_ij^q`, recovering the classical family: `0` is raw stress ({@link SMACOF}'s
 * objective), `-1` Sammon stress ({@link SAMMON}'s), `-2` elastic scaling ({@link KKMDS}). A matrix
 * or a function works too, where a zero weight drops the pair from the objective — the usual way to
 * express missing dissimilarities.
 *
 * It does not replace {@link SMACOF} or {@link SAMMON}: those minimise the same objectives with
 * their own historical algorithms and reach different local minima.
 *
 * Optimised by Jacobi-preconditioned gradient descent with a backtracking line search, warm-started
 * from classical MDS. O(N²·d) per iteration.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersStressMDS>
 * @category Dimensionality Reduction
 * @see {@link KKMDS} for the `weights: -2` preset.
 * @example
 * // Sammon stress, converged further than SAMMON's own optimizer takes it
 * const Y = new StressMDS(X, { weights: -1 }).transform();
 * @example
 * // Ignore pairs whose dissimilarity was never measured
 * const W = new Matrix(N, N, (i, j) => (observed(i, j) ? 1 : 0));
 * const Y = new StressMDS(D, { metric: "precomputed", weights: W }).transform();
 */
declare class StressMDS<T extends InputType> extends DR<T, ParametersStressMDS> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersStressMDS>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersStressMDS>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersStressMDS>): Promise<T_1>;
    /**
     * Weighted metric MDS.
     *
     * @param {T} X - The high-dimensional data, or a precomputed distance matrix.
     * @param {Partial<ParametersStressMDS>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersStressMDS>);
    /**
     * The target distances. Named apart from the base class's `_D`, which is the *input*
     * dimensionality of `X`, not a matrix.
     *
     * @type {Matrix | undefined}
     */
    _target_distances: Matrix | undefined;
    /** @type {Matrix | undefined} */
    _weights: Matrix | undefined;
    /** @type {number} */
    _energy: number;
    /**
     * The weighted stress `σ(Y)` of the current embedding.
     *
     * @returns {number}
     */
    get energy(): number;
    /**
     * The weight matrix actually in use, whatever form `weights` was given in.
     *
     * @returns {Matrix}
     */
    get weights(): Matrix;
    /**
     * Materialises `weights` into an N ⨯ N matrix.
     *
     * Non-positive and non-finite weights are normalised to exactly 0, which is the kernel's "skip
     * this pair" signal. That is what makes a zero distance harmless under a negative exponent: `0^-2`
     * is `Infinity`, and coincident rows are ordinary in real data, so the pair is dropped rather
     * than given infinite stiffness.
     *
     * @private
     * @param {Matrix} D
     * @returns {Matrix}
     */
    private _build_weights;
    /**
     * Classical MDS on the target distances, scaled by the square roots of the eigenvalues.
     *
     * The scaling matters: the descent starts from here, and an unscaled eigenvector basis would put
     * the configuration at a completely different magnitude from `d_ij`, so the first steps would be
     * spent on scale rather than structure.
     *
     * @private
     * @param {Matrix} D
     * @param {number} d
     * @returns {Matrix}
     */
    private _classical_mds;
    /** Computes the target distances, the weights, and the starting configuration. */
    init(): this;
    /**
     * @private
     * @param {Matrix} Y
     * @param {Matrix} D
     * @param {Matrix} W
     * @returns {number}
     */
    private _compute_energy;
    /**
     * One preconditioned gradient step with a backtracking line search, in JS.
     *
     * Mirrors `stress_mds_step_f64` in `StressMDS.as.ts` loop for loop, so the two paths agree to
     * floating point tolerance. See that file for why the gradient is divided by the weighted degree.
     *
     * @private
     * @param {Matrix} Y - Updated in place on an accepted step.
     * @param {Matrix} D
     * @param {Matrix} W
     * @param {Float64Array} G - Gradient scratch, `N * d`.
     * @param {number} step
     * @param {number} energy_current
     * @returns {{ energy: number; step: number; accepted: boolean }}
     */
    private _step_js;
    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Computes the projection.
     *
     * @returns {T}
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersKKMDS} from "./index.js" */
/**
 * Kamada-Kawai Multidimensional Scaling (KKMDS)
 *
 * {@link StressMDS} fixed at `weights: -2`, weighting each pair by `1 / d_ij²`. That is the
 * Kamada-Kawai energy from graph drawing, known in the MDS literature as *elastic scaling*.
 *
 * Its classic use is laying out a graph from its shortest-path distances — pass those as a
 * `"precomputed"` matrix, as {@link MINFOTree} does.
 *
 * @class
 * @template {InputType} T
 * @extends StressMDS<T>
 * @category Dimensionality Reduction
 * @see {@link https://doi.org/10.1016/0020-0190(89)90102-6}
 * @see {@link StressMDS} for other weightings.
 */
declare class KKMDS<T extends InputType> extends StressMDS<T> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersKKMDS>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersKKMDS>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersKKMDS>): Promise<T_1>;
    /**
     * Kamada-Kawai weighted MDS.
     *
     * @param {T} X - The high-dimensional data, or a precomputed distance matrix.
     * @param {Partial<ParametersKKMDS>} [parameters] - Object containing parameterization of the DR
     *   method. `weights` is not among them; it is fixed at `-2`.
     */
    constructor(X: T, parameters?: Partial<ParametersKKMDS>);
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersLDA} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Linear Discriminant Analysis (LDA)
 *
 * A supervised dimensionality reduction technique that finds the axes that
 * maximize the separation between multiple classes.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLDA>
 * @category Dimensionality Reduction
 */
declare class LDA<T extends InputType> extends DR<T, ParametersLDA> {
    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType, Para extends {
        seed?: number;
    }>(X: T_1, parameters: Para): T_1;
    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType, Para extends {
        seed?: number;
    }>(X: T_1, parameters: Para): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType, Para extends {
        seed?: number;
    }>(X: T_1, parameters: Para): Promise<T_1>;
    /**
     * Linear Discriminant Analysis.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLDA> & { labels: any[] | Float64Array }} parameters - Object containing parameterization of the DR method.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
    constructor(X: T, parameters: Partial<ParametersLDA> & {
        labels: any[] | Float64Array;
    });
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T} - The projected data.
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersLLE} from "./index.js" */
/** @import {KNN} from "../knn/KNN.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Locally Linear Embedding (LLE)
 *
 * A nonlinear dimensionality reduction technique that preserves local
 * linear relationships between points. It represents each point as a linear
 * combination of its neighbors.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLLE>
 * @category Dimensionality Reduction
 * @see {@link ISOMAP} for another nonlinear alternative
 */
declare class LLE<T extends InputType> extends DR<T, ParametersLLE> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLLE>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLLE>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLLE>): Promise<T_1>;
    /**
     * Locally Linear Embedding.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLLE>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X: T, parameters: Partial<ParametersLLE>);
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T}
     */
    transform(): T;
}

/**
 * Pairwise Controlled Manifold Approximation Projection (PaCMAP)
 *
 * A dimensionality reduction technique that uses three types of point pairs — nearest neighbor
 * (NN), mid-near (MN), and further (FP) pairs — with a dynamic three-phase weight schedule and
 * Adam optimization to preserve both local and global structure. Neighbours are scaled by a local
 * density estimate before selection, which is what separates the NN set from a plain KNN graph.
 *
 * @class
 * @template {InputType} T
 * @template {ParametersPaCMAP} [P=ParametersPaCMAP] - Widened by {@link LocalMAP}, which adds a parameter.
 * @extends DR<T, P>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/2012.04456|PaCMAP Paper}
 * @see {@link https://github.com/YingfanWang/PaCMAP|PaCMAP GitHub}
 * @see {@link UMAP} for a related graph-based technique
 * @see {@link LocalMAP} for the local-refinement variant
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const pacmap = new druid.PaCMAP(X, {
 *     n_neighbors: 10,
 *     MN_ratio: 0.5,
 *     FP_ratio: 2.0,
 *     seed: 42
 * });
 *
 * const Y = pacmap.transform(); // 450 iterations (default)
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
declare class PaCMAP<T extends InputType, P extends ParametersPaCMAP = ParametersPaCMAP> extends DR<T, P> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersPaCMAP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersPaCMAP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersPaCMAP>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<P>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<P>);
    _iter: number;
    /**
     * Finds the `n_candidates` nearest neighbors of every point.
     *
     * The reference over-fetches 50 beyond `n_neighbors` so that the density rescaling in
     * {@link _scaled_neighbors} has something to re-rank, which makes this the most expensive part
     * of a run — asking an approximate index for 60 neighbors instead of 10 costs it far more than
     * the extra 50 suggests. The default path is therefore an exact blocked search whose selection
     * happens inside the WASM kernel; it beats a tree at these widths and matches the exact search
     * the reference performs. A `knn` index is used instead when one is supplied.
     *
     * @protected
     * @param {Matrix} X - The matrix to search, already reduced if PCA was applied.
     * @param {number} n_candidates
     * @returns {{ dist: Float64Array; idx: Int32Array }} Row-major `N` ⨯ `n_candidates`, ascending.
     */
    protected _find_candidates(X: Matrix, n_candidates: number): {
        dist: Float64Array;
        idx: Int32Array;
    };
    /**
     * Exact neighbour search in JS, for a non-euclidean metric or when WASM is unavailable.
     *
     * @protected
     * @param {Matrix} X
     * @param {number} n_candidates
     * @returns {{ dist: Float64Array; idx: Int32Array }}
     */
    protected _find_candidates_js(X: Matrix, n_candidates: number): {
        dist: Float64Array;
        idx: Int32Array;
    };
    /**
     * Selects the NN pairs by rescaled distance.
     *
     * Each candidate distance is divided by the local scale of both endpoints — the mean distance
     * to their 4th to 6th neighbors — before the top `n_neighbors` are taken. Without this the NN
     * set is an ordinary KNN graph and dense regions dominate the attractive term.
     *
     * @protected
     * @param {{ dist: Float64Array; idx: Int32Array }} candidates
     * @param {number} n_candidates
     * @param {number} n_neighbors
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    protected _scaled_neighbors(candidates: {
        dist: Float64Array;
        idx: Int32Array;
    }, n_candidates: number, n_neighbors: number): Int32Array;
    /**
     * Draws `n_samples` distinct indices, excluding `self` and everything in `reject`.
     *
     * @protected
     * @param {number} n_samples
     * @param {number} self
     * @param {Int32Array | number[]} reject
     * @param {Int32Array} out - Receives the sample; must hold `n_samples`.
     * @returns {number} How many were drawn, short of `n_samples` only if the pool is exhausted.
     */
    protected _sample_excluding(n_samples: number, self: number, reject: Int32Array | number[], out: Int32Array): number;
    /**
     * Samples mid-near pairs: for each point, six random draws from the whole dataset, of which the
     * second closest is kept. Sampling globally rather than from the neighbour list is what makes
     * these pairs "mid-near" and is the mechanism behind PaCMAP's global structure in phase 1.
     *
     * @protected
     * @param {Matrix} X
     * @param {number} n_MN
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    protected _sample_mn_pairs(X: Matrix, n_MN: number): Int32Array;
    /**
     * Samples further pairs: random points that are not among `i`'s neighbors.
     *
     * @protected
     * @param {Int32Array} nn_pairs
     * @param {number} n_neighbors
     * @param {number} n_FP
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    protected _sample_fp_pairs(nn_pairs: Int32Array, n_neighbors: number, n_FP: number): Int32Array;
    /**
     * Accumulates the gradient of one pair type into `grad_flat`.
     *
     * All three losses share the form `w * f(d_ij)` with `d_ij = 1 + ||y_i - y_j||²`, so one loop
     * covers them: the attractive pairs use `d_ij / (a + d_ij)` and the repulsive ones
     * `1 / (1 + d_ij)`, which is the same denominator with `a = 1` and the sign flipped.
     *
     * `Y` is indexed flat rather than through `Matrix.row`, which would allocate a subarray per
     * endpoint — two per pair, tens of millions per run, and by far the dominant cost of a step.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {Int32Array} pairs - Flat `[i, j, ...]` pair array.
     * @param {number} w - Weight for this pair type.
     * @param {number} a - Denominator constant: 10 for NN, 10000 for MN, 1 for FP.
     * @param {boolean} repulsive
     */
    protected _accumulate_gradients(grad_flat: Float64Array, pairs: Int32Array, w: number, a: number, repulsive: boolean): void;
    /**
     * Returns the weight schedule for the current iteration.
     *
     * @protected
     * @param {number} iter - Current iteration (0-indexed)
     * @returns {{ w_nn: number; w_mn: number; w_fp: number }}
     */
    protected _get_weights(iter: number): {
        w_nn: number;
        w_mn: number;
        w_fp: number;
    };
    /**
     * Applies Adam optimizer update to Y using accumulated gradients.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient
     */
    protected _adam_update(grad_flat: Float64Array): void;
    _adam_t: any;
    /**
     * Initializes PaCMAP: preprocessing, PCA embedding, NN, MN and FP pairs, and Adam state.
     *
     * @returns {this}
     */
    init(): this;
    _X_knn: Matrix | undefined;
    _nn_pairs: Int32Array<ArrayBufferLike> | undefined;
    _mn_pairs: Int32Array<ArrayBufferLike> | undefined;
    _fp_pairs: Int32Array<ArrayBufferLike> | undefined;
    _adam_m: Float64Array<ArrayBuffer> | undefined;
    _adam_v: Float64Array<ArrayBuffer> | undefined;
    _grad: Float64Array<ArrayBuffer> | undefined;
    /**
     * Accumulates the nearest-neighbor gradient for the current step.
     *
     * Split out from {@link next} because it is the only part of a step {@link LocalMAP} changes,
     * and only in phase 3. Everything else — the mid-near and further terms, the Adam update, the
     * iteration bookkeeping — is shared.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {number} w_nn
     */
    protected _accumulate_nn_gradients(grad_flat: Float64Array, w_nn: number): void;
    /**
     * Hook run after the embedding has been updated, with the iteration that just finished.
     *
     * Nothing to do here; {@link LocalMAP} redraws its further pairs from it.
     *
     * @protected
     * @param {number} iter - The 0-indexed iteration that just completed.
     */
    protected _after_step(iter: number): void;
    /**
     * Performs one optimization step.
     *
     * @returns {Matrix}
     */
    next(): Matrix;
    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {T}
     */
    transform(iterations?: number): T;
    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {Generator<T, T, void>}
     */
    generator(iterations?: number): Generator<T, T, void>;
}

/** @import {InputType} from "../index.js" */
/** @import {Matrix} from "../matrix/index.js" */
/** @import {ParametersLocalMAP} from "./index.js" */
/**
 * LocalMAP
 *
 * A variant of {@link PaCMAP} that sharpens local cluster separation in phase 3. Attraction on the
 * nearest-neighbor pairs is scaled by `low_dist_thres / (2 * sqrt(d_ij))`, which strengthens it for
 * pairs already close in the embedding and weakens it for far ones, and the further pairs are
 * periodically redrawn from points that are *near* in the current embedding rather than staying the
 * random set chosen at initialization.
 *
 * @class
 * @template {InputType} T
 * @extends PaCMAP<T, ParametersLocalMAP>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/2012.04456|PaCMAP Paper}
 * @see {@link PaCMAP} for the base algorithm
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const localmap = new druid.LocalMAP(X, {
 *     n_neighbors: 10,
 *     low_dist_thres: 10,
 *     seed: 42
 * });
 *
 * const Y = localmap.transform(); // 450 iterations (default)
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
declare class LocalMAP<T extends InputType> extends PaCMAP<T, ParametersLocalMAP> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLocalMAP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLocalMAP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLocalMAP>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLocalMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersLocalMAP>);
    /**
     * The first iteration of phase 3, past which LocalMAP diverges from PaCMAP.
     *
     * The reference tests `itr > phase_1 + phase_2`, so the first phase 3 step still runs plain
     * PaCMAP; this is that boundary, and "after" means strictly greater.
     *
     * @private
     * @returns {number}
     */
    private get _phase3_start();
    /**
     * Accumulates NN gradients with LocalMAP's local scaling, `nn_scale / sqrt(d_ij)`.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {Int32Array} pairs - Flat `[i, j, ...]` pair array.
     * @param {number} w_nn - NN weight.
     * @param {number} nn_scale - `low_dist_thres / 2`.
     */
    protected _accumulate_gradients_local_nn(grad_flat: Float64Array, pairs: Int32Array, w_nn: number, nn_scale: number): void;
    /**
     * Redraws the further pairs from points within `low_dist_thres` of `i` in the embedding.
     *
     * A point that is far in the input but has drifted close in the embedding is exactly what the
     * repulsive term needs to push apart, so sampling from the current layout rather than the
     * original random set is what sharpens the cluster boundaries. Candidates that are already
     * neighbors are rejected, and a row that finds nothing within the threshold keeps its old
     * partner.
     *
     * @protected
     * @param {number} low_dist_thres
     */
    protected _resample_local_fp_pairs(low_dist_thres: number): void;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersLSP} from "./index.js" */
/**
 * Least Square Projection (LSP)
 *
 * A dimensionality reduction technique that uses a small set of control points
 * (projected with MDS) to define the projection for the rest of the data
 * using a Laplacian-based optimization.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLSP>
 * @category Dimensionality Reduction
 */
declare class LSP<T extends InputType> extends DR<T, ParametersLSP> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLSP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLSP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersLSP>): Promise<T_1>;
    /**
     * Least Squares Projection.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLSP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     */
    constructor(X: T, parameters?: Partial<ParametersLSP>);
    /**
     * @returns {LSP<T>}
     */
    init(): LSP<T>;
    _A: Matrix | undefined;
    _b: Matrix | undefined;
    /**
     * Computes the projection.
     *
     * @returns {T} Returns the projection.
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersLTSA} from "./index.js" */
/** @import {KNN} from "../knn/KNN.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Local Tangent Space Alignment (LTSA)
 *
 * A nonlinear dimensionality reduction algorithm that represents the local
 * geometry of the manifold by tangent spaces and then aligns them to reveal
 * the global structure.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLTSA>
 * @category Dimensionality Reduction
 */
declare class LTSA<T extends InputType> extends DR<T, ParametersLTSA> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLTSA>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLTSA>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersLTSA>): Promise<T_1>;
    /**
     * Local Tangent Space Alignment
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLTSA>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X: T, parameters: Partial<ParametersLTSA>);
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Transforms the inputdata `X` to dimenionality `d`.
     *
     * @returns {T}
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersMDS} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Classical Multidimensional Scaling (MDS)
 *
 * A linear dimensionality reduction technique that seeks to preserve the
 * pairwise distances between points as much as possible in the lower-dimensional
 * space.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersMDS>
 * @category Dimensionality Reduction
 * @see {@link PCA} for another linear alternative
 */
declare class MDS<T extends InputType> extends DR<T, ParametersMDS> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMDS>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMDS>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMDS>): Promise<T_1>;
    /**
     * Classical MDS.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersMDS>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersMDS>);
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T}
     */
    transform(): T;
    _d_X: Matrix | undefined;
    /** @returns {number} - The stress of the projection. */
    stress(): number;
}

/**
 * Minimum Information Trees (MINFO Trees)
 *
 * A visualization method for clustered high-dimensional data. It reads the cluster labels on a
 * k-nearest-neighbor graph as a q-state Potts model, weights the edges by an information-geometric
 * curvature derived from it, and keeps the minimum spanning tree of that weighted graph.
 *
 * {@link projection} returns 2D coordinates like any other method here, but the tree is the actual
 * output — read {@link edges} to draw it.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersMINFOTree>
 * @category Dimensionality Reduction
 * @see {@link https://doi.org/10.1109/ACCESS.2025.3602730}
 * @see {@link TopoMap} for another spanning-tree-based projection.
 */
declare class MINFOTree<T extends InputType> extends DR<T, ParametersMINFOTree> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMINFOTree>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMINFOTree>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMINFOTree>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersMINFOTree>): Promise<T_1>;
    /**
     * Minimum Information Trees.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersMINFOTree>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1109/ACCESS.2025.3602730}
     */
    constructor(X: T, parameters?: Partial<ParametersMINFOTree>);
    /**
     * The edges of the Minimum Information Tree, as `[u, v, weight]` over row indices of `X`,
     * ascending by weight. This is the method's real output — {@link projection} is one drawing of
     * it.
     *
     * @returns {WeightedEdge[]}
     */
    get edges(): WeightedEdge[];
    /**
     * The cluster label per point, whether supplied or computed in step 1.
     *
     * @returns {Int32Array} Labels, remapped to `0 … q-1`.
     */
    get labels(): Int32Array;
    /**
     * The information curvature `S_i` per point, normalised as the edge weighting uses it.
     *
     * @returns {Float64Array}
     */
    get curvature(): Float64Array;
    /**
     * The maximum pseudo-likelihood estimate of the Potts inverse temperature.
     *
     * @returns {number}
     */
    get beta(): number;
    /**
     * Resolves the labels field: either the caller's, or the outcome of clustering `X`.
     *
     * Labels are remapped to a dense `0 … q-1` so they can index the Potts state vectors directly.
     *
     * @private
     * @returns {Int32Array}
     */
    private _make_labels;
    _q: number | undefined;
    /**
     * Finds the distance threshold that cuts a dendrogram into exactly `k` clusters.
     *
     * {@link HierarchicalClustering.get_cluster_list} cuts at a height, not at a cluster count, so
     * the height has to be recovered: breaking the `k - 1` tallest merges leaves `k` components, and
     * any threshold between the `k`-th and `(k-1)`-th tallest does that. The midpoint is used so
     * ties on either side cannot land the cut on the wrong number.
     *
     * Only `2 … N-1` is reachable this way. A dendrogram has `N-1` merges, and the traversal emits a
     * cluster per subtree it stops at, never a bare leaf — so there is no height that separates all
     * `N` points. Asking for more is rejected rather than approximated: cutting below every merge
     * silently yields *one* cluster, the opposite of what was asked for.
     *
     * @private
     * @param {HierarchicalClustering} hc
     * @param {number} k - At least 2; {@link _make_labels} rejects anything smaller before this runs.
     * @returns {number}
     */
    private _dendrogram_cut;
    /**
     * Builds the symmetric k-nearest-neighbor graph, connected.
     *
     * The Potts model is defined over a neighborhood system and the layout needs finite pairwise
     * distances, so a k-NNG that falls into several components would break both. When that happens
     * the components are stitched together with the shortest edge between each pair, which is what
     * the author's reference implementation achieves by unioning in a spanning tree of the complete
     * graph.
     *
     * @private
     * @returns {number[][]} Adjacency list.
     */
    private _make_knn_graph;
    /**
     * Adds the fewest edges needed to make `adjacency` connected, preferring short ones.
     *
     * @private
     * @param {Set<number>[]} adjacency - Mutated in place.
     */
    private _connect_components;
    /**
     * Counts, for every vertex, how many of its neighbors carry each of the `q` labels.
     *
     * `U[i * q + l]` is `U_i(l)` in the paper.
     *
     * @private
     * @param {number[][]} adjacency
     * @returns {Float64Array}
     */
    private _neighbor_counts;
    /**
     * Gibbs weights `w_i(l) = exp(β ⋅ U_i(l))`, shifted to avoid overflow.
     *
     * Every quantity built from `w` here is a ratio of equally-weighted sums, so scaling the whole
     * vector by `exp(-β ⋅ max_l U_i(l))` leaves φ and ψ untouched while keeping `exp` in range even
     * for a large β on a dense neighborhood.
     *
     * @private
     * @param {Float64Array} U
     * @param {number} offset - Start of vertex `i`'s counts in `U`.
     * @param {number} q
     * @param {number} beta
     * @param {Float64Array} out - Length `q` scratch buffer.
     */
    private _gibbs_weights;
    /**
     * Maximum pseudo-likelihood estimate of the inverse temperature β.
     *
     * Solves the pseudo-likelihood equation (3), which sets the observed energy equal to the
     * expected energy:
     *
     * ```
     * f(β) = Σ_i U_i(x_i) - Σ_i E_β[U_i] = 0
     * ```
     *
     * `f` is the negative of a sum of variances, hence non-increasing, so plain bisection is enough
     * and needs no derivatives — the same role Brent's method plays in the paper, without the extra
     * machinery. `f(0) <= 0` means the labels agree with their neighborhoods no better than chance,
     * which leaves β = 0 as the estimate.
     *
     * @private
     * @param {Float64Array} U
     * @returns {number}
     */
    private _estimate_beta;
    /**
     * Information curvature `S_i = -ψ_i / φ_i` at every vertex, normalised to `(0, 1]`.
     *
     * The paper states φ (17) and ψ (23) as Kronecker products over `q ⨯ q` matrices, but both
     * collapse: since summing every entry of `a ⊗ aᵀ` is just `(Σa)²`, the double sums reduce to
     *
     * ```
     * φ_i = (Σ_l v_l w_l / Σ_l w_l)²                                  — squared observed/expected energy gap
     * ψ_i = [(Σ_l U_l² w_l)(Σ_l w_l) - (Σ_l U_l w_l)²] / (Σ_l w_l)²   — the variance of U under w
     * ```
     *
     * which is O(q) per vertex instead of O(q²), and allocates nothing.
     *
     * `epsilon` floors the denominator, following the author's reference implementation. It is load
     * bearing: both φ and ψ vanish for an interior point as β grows, and by β ≈ 8 ψ has underflowed
     * to exactly 0, so `S` is pinned at the top of the range rather than the bottom. That inverts the
     * interior/boundary reading in §VI-C for anything above β ≈ 2, which well-separated clusters
     * comfortably exceed. The formulas are implemented as published in either regime; see
     * `docs/dimred/minfotree.md` for the measurements.
     *
     * @private
     * @param {Float64Array} U
     * @param {number} beta
     * @returns {Float64Array}
     */
    private _information_curvature;
    /**
     * All-pairs shortest paths over the Minimum Information Tree.
     *
     * A tree has exactly one path between any two vertices, so a depth-first walk from each source
     * settles every distance in O(N) — no priority queue and no relaxation, unlike the Dijkstra
     * {@link ISOMAP} needs on its general k-NNG.
     *
     * @private
     * @param {WeightedEdge[]} edges
     * @returns {Matrix}
     */
    private _tree_distances;
    /** Runs steps 1-5: labels, k-NNG, β, curvature, information graph and its spanning tree. */
    init(): this;
    _labels: Int32Array<ArrayBufferLike> | undefined;
    _beta: number | undefined;
    _curvature: Float64Array<ArrayBufferLike> | undefined;
    _edges: WeightedEdge[] | undefined;
    /**
     * Computes the projection.
     *
     * @returns {T}
     */
    transform(): T;
    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersPCA} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */
/**
 * Principal Component Analysis (PCA)
 *
 * A linear dimensionality reduction technique that identifies the axes (principal components)
 * along which the variance of the data is maximized.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersPCA>
 * @category Dimensionality Reduction
 * @see {@link MDS} for another linear alternative
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2], [3, 4], [5, 6]];
 * const pca = new druid.PCA(X, { d: 2 });
 * const Y = pca.transform();
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
declare class PCA<T extends InputType> extends DR<T, ParametersPCA> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersPCA>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {Matrix}
     */
    static principal_components<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersPCA>): Matrix;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersPCA>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersPCA>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersPCA>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersPCA>);
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T} - The projected data.
     */
    transform(): T;
    /**
     * Computes the `d` principal components of Matrix `X`.
     *
     * @returns {Matrix}
     */
    principal_components(): Matrix;
    V: Matrix | undefined;
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersPCA, ParametersMDS, ParametersSAMMON} from "./index.js" */
/** @typedef {"PCA" | "MDS" | "random"} AvailableInit */
/** @typedef {{ PCA: ParametersPCA; MDS: ParametersMDS; random: {} }} ChooseDR */
/**
 * Sammon's Mapping
 *
 * A nonlinear dimensionality reduction technique that minimizes a stress
 * function based on the ratio of pairwise distances in high and low dimensional spaces.
 *
 * A given `seed` reproduces a layout exactly within one engine and library build, but not across
 * browsers or Node versions: the gradient step is chaotic, so a last-bit difference grows into a
 * visibly different — though equally valid — layout. See {@link DR} for why, and {@link SMACOF}
 * for a stable alternative.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersSAMMON<AvailableInit>>
 * @category Dimensionality Reduction
 */
declare class SAMMON<T extends InputType> extends DR<T, ParametersSAMMON<AvailableInit>> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSAMMON<AvailableInit>>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSAMMON<AvailableInit>>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSAMMON<AvailableInit>>): Promise<T_1>;
    /**
     * SAMMON's Mapping
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters] - Object containing parameterization of the DR
     *   method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X: T, parameters?: Partial<ParametersSAMMON<AvailableInit>>);
    /** @type {Matrix | undefined} */
    distance_matrix: Matrix | undefined;
    /**
     * Initializes the projection.
     *
     * @param {Matrix | undefined} D
     * @returns {asserts D is Matrix}
     */
    init(D: Matrix | undefined): asserts D is Matrix;
    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {T} The projection of `X`.
     */
    transform(max_iter?: number): T;
    /**
     * Transforms the inputdata `X` to dimenionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection of
     *   `X`.
     */
    generator(max_iter?: number): Generator<T, T, void>;
    _step(): Matrix;
}
type AvailableInit = "PCA" | "MDS" | "random";
type ChooseDR = {
    PCA: ParametersPCA;
    MDS: ParametersMDS;
    random: {};
};

/** @import {InputType} from "../index.js" */
/** @import {ParametersSMACOF} from "./index.js" */
/**
 * Metric Multidimensional Scaling (MDS) via SMACOF.
 *
 * SMACOF (Scaling by Majorizing a Complicated Function) is an iterative majorization
 * algorithm for solving metric multidimensional scaling problems, which aims to
 * minimize the stress function.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersSMACOF>
 * @category Dimensionality Reduction
 * @see {@link MDS} for the classical approach.
 */
declare class SMACOF<T extends InputType> extends DR<T, ParametersSMACOF> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSMACOF>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSMACOF>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSMACOF>): Promise<T_1>;
    /**
     * SMACOF for MDS.
     *
     * @param {T} X - The high-dimensional data or precomputed distance matrix.
     * @param {Partial<ParametersSMACOF>} [parameters] - Object containing parameterization.
     */
    constructor(X: T, parameters?: Partial<ParametersSMACOF>);
    /**
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    generator(): Generator<T, T, void>;
    /**
     * @returns {T}
     */
    transform(): T;
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersSQDMDS} from "./index.js" */
/**
 * SQuadMDS (Stochastic Quartet MDS)
 *
 * A lean Stochastic Quartet MDS improving global structure preservation in
 * neighbor embedding like t-SNE and UMAP.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersSQDMDS>
 * @category Dimensionality Reduction
 */
declare class SQDMDS<T extends InputType> extends DR<T, ParametersSQDMDS> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSQDMDS>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSQDMDS>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersSQDMDS>): Promise<T_1>;
    /**
     * SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE
     * and UMAP.
     *
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @see {@link https://arxiv.org/pdf/2202.12087.pdf}
     */
    constructor(X: T, parameters?: Partial<ParametersSQDMDS>);
    init(): void;
    _add: ((...summands: Float64Array<ArrayBufferLike>[]) => Float64Array) | undefined;
    _sub_div: ((x: Float64Array<ArrayBufferLike>, y: Float64Array<ArrayBufferLike>, div: number) => Float64Array) | undefined;
    _minus: ((a: Float64Array<ArrayBufferLike>, b: Float64Array<ArrayBufferLike>) => Float64Array) | undefined;
    _mult: ((a: Float64Array<ArrayBufferLike>, v: number) => Float64Array) | undefined;
    _LR_init: number | undefined;
    _LR: number | undefined;
    _offset: number | undefined;
    _momentums: Matrix | undefined;
    _grads: Matrix | undefined;
    _indices: number[] | undefined;
    /** @type {Uint32Array} */
    _flat_quartets: Uint32Array<ArrayBufferLike> | undefined;
    /** @type {(i: number, j: number, X: Matrix) => number} */
    _HD_metric: ((i: number, j: number, X: Matrix) => number) | undefined;
    /** @type {(i: number, j: number, X: Matrix) => number} */
    _HD_metric_exaggeration: ((i: number, j: number, X: Matrix) => number) | undefined;
    /**
     * Computes the projection.
     *
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {T} The projection.
     */
    transform(iterations?: number): T;
    _decay_start: number | undefined;
    /**
     * Computes the projection.
     *
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {Generator<T, T, void>} The intermediate steps of the projection.
     */
    generator(iterations?: number): Generator<T, T, void>;
    /**
     * Performs an optimization step.
     *
     * @private
     * @param {number} i - Acutal iteration.
     * @param {number} iterations - Number of iterations.
     */
    private _step;
    _distance_exaggeration: boolean | undefined;
    /**
     * Creates quartets of non overlapping indices.
     *
     * @private
     * @returns {Uint32Array[]}
     */
    private __quartets;
    /**
     * Computes and applies gradients, and updates momentum.
     *
     * @private
     * @param {boolean} distance_exaggeration
     */
    private _nestrov_iteration;
    /**
     * Computes the gradients.
     *
     * @param {Matrix} Y - The Projection.
     * @param {Matrix} grads - The gradients.
     * @param {boolean} [exaggeration=false] - Whether or not to use early exaggeration. Default is `false`
     * @param {boolean} [zero_grad=true] - Whether or not to reset the gradient in the beginning. Default is `true`
     * @param {Uint32Array[]} [quartets] - Quartets to accumulate over. Pass the ones already drawn for this
     *   iteration; omitting them draws a fresh set, which advances the seeded randomizer.
     * @returns {Matrix} The gradients.
     */
    _fill_MDS_grads(Y: Matrix, grads: Matrix, exaggeration?: boolean, zero_grad?: boolean, quartets?: Uint32Array[]): Matrix;
    /**
     * Quartet gradients for a projection.
     *
     * @private
     * @param {Matrix} Y - The acutal projection.
     * @param {number[]} quartet - The indices of the quartet.
     * @param {Float64Array} D_hd - The high-dimensional distances of the quartet.
     * @returns {Float64Array[]} The gradients for the quartet.
     */
    private _compute_quartet_grads;
    /**
     * Gradients for one element of the loss function's sum.
     *
     * @private
     * @param {Float64Array} a
     * @param {Float64Array} b
     * @param {Float64Array} c
     * @param {Float64Array} d
     * @param {number} d_ab
     * @param {number} d_ac
     * @param {number} d_ad
     * @param {number} d_bc
     * @param {number} d_bd
     * @param {number} d_cd
     * @param {number} p_ab
     * @param {number} sum_LD_dist
     * @returns {Float64Array[]}
     */
    private _ABCD_grads;
    /**
     * Inline!
     *
     * @param {number} d
     */
    __minus(d: number): (a: Float64Array<ArrayBufferLike>, b: Float64Array<ArrayBufferLike>) => Float64Array;
    /**
     * Inline!
     *
     * @param {number} d
     */
    __add(d: number): (...summands: Float64Array<ArrayBufferLike>[]) => Float64Array;
    /**
     * Inline!
     *
     * @param {number} d
     */
    __mult(d: number): (a: Float64Array<ArrayBufferLike>, v: number) => Float64Array;
    /**
     * Creates a new array `(x - y) / div`.
     *
     * @param {number} d
     */
    __sub_div(d: number): (x: Float64Array<ArrayBufferLike>, y: Float64Array<ArrayBufferLike>, div: number) => Float64Array;
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersTopoMap} from "./index.js" */
/**
 * TopoMap
 *
 * A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
 * It aims to preserve the topological structure of the data by maintaining
 * the connectivity of a minimum spanning tree.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersTopoMap>
 * @category Dimensionality Reduction
 */
declare class TopoMap<T extends InputType> extends DR<T, ParametersTopoMap> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersTopoMap>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersTopoMap>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters: Partial<ParametersTopoMap>): Promise<T_1>;
    /**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTopoMap>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X: T, parameters: Partial<ParametersTopoMap>);
    /**
     * Computes the minimum spanning tree, using a given metric
     *
     * @private
     * @param {import("../metrics/index.js").Metric} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    private _make_minimum_spanning_tree;
    /** Initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree. */
    init(): this;
    _Emst: WeightedEdge[] | undefined;
    /**
     * Returns true if Point C is left of line AB.
     *
     * @private
     * @param {Float64Array} PointA - Point A of line AB
     * @param {Float64Array} PointB - Point B of line AB
     * @param {Float64Array} PointC - Point C
     * @returns {boolean}
     */
    private __hull_cross;
    /**
     * Computes the convex hull of the set of Points S
     *
     * @private
     * @param {Float64Array[]} S - Set of Points.
     * @returns {Float64Array[]} Convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     */
    private __hull;
    /**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     *
     * @private
     * @param {Float64Array} PointA
     * @param {Float64Array} PointB
     * @returns {{ sin: number; cos: number }} Object containing the sinus- and cosinus-values for a rotation.
     */
    private __findAngle;
    /**
     * @private
     * @param {Float64Array[]} hull
     * @param {Float64Array} p
     * @param {boolean} topEdge
     * @returns {{ sin: number; cos: number; tx: number; ty: number }}
     */
    private __align_hull;
    /**
     * @private
     * @param {Float64Array} Point - The point which should get transformed.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} Transformation - Contains the values for
     *   translation and rotation.
     */
    private __transform;
    /**
     * Calls `__transform` for each point in Set C
     *
     * @private
     * @param {Float64Array[]} C - Set of points.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} t - Transform object.
     * @param {number} yOffset - Value to offset set C.
     */
    private __transform_component;
    /**
     * @private
     * @param {Float64Array} root_u - Root of component u
     * @param {Float64Array} root_v - Root of component v
     * @param {Float64Array} p_u - Point u
     * @param {Float64Array} p_v - Point v
     * @param {number} w - Edge weight w
     * @param {DisjointSet<Float64Array>} components - The disjoint set containing the components
     */
    private __align_components;
    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {T}
     */
    transform(): T;
    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {Generator<T, T, void>}
     */
    generator(): Generator<T, T, void>;
}

/**
 * Base class for all K-Nearest Neighbors (KNN) search algorithms.
 *
 * Provides a common interface for elements management and search operations.
 *
 * @abstract
 * @category KNN
 * @template {number[] | Float64Array} T - Type of elements
 * @template {Object} Para - Type of parameters
 * @class
 */
declare class KNN<T extends number[] | Float64Array, Para extends Object> {
    /**
     * @param {T[]} elements
     * @param {Para} parameters
     */
    constructor(elements: T[], parameters: Para);
    /** @type {T[]} */
    _elements: T[];
    /** @type {Para} */
    _parameters: Para;
    /** @type {"typed" | "array"} */
    _type: "typed" | "array";
    /**
     * Seeded source of randomness shared by every index. Construction is randomized — the trees
     * pick quickselect pivots from it — so the `seed` parameter is what makes a built index, and
     * therefore its query results, reproducible.
     *
     * @type {Randomizer}
     */
    _randomizer: Randomizer;
    /**
     * @abstract
     * @param {T} t
     * @param {number} k
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(t: T, k: number): {
        element: T;
        index: number;
        distance: number;
    }[];
    /**
     * Searches the `k` nearest neighbors of the element stored at index `i`.
     *
     * The queried element is never part of the result. It is trivially its own closest neighbor at
     * distance 0, which is never what a caller asking "what is this point near?" wants, so every
     * caller used to strip it back out — each in its own, subtly different way. Note the asymmetry
     * with `search`: an arbitrary query point has no "self" to exclude, so there `k` means
     * "k results", while here it means "k neighbors".
     *
     * The self match is removed **by index** — not by position, and not by looking for a zero
     * distance. Position is wrong because an approximate index may order ties differently or miss
     * the element altogether (one extra candidate is requested to cover that), and a zero distance
     * is wrong because genuine duplicate points share it and must survive.
     *
     * @param {number} i - Index of the query element.
     * @param {number} [k=5] - Number of neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} The `k` nearest *other* elements,
     *   closest first. Empty when `i` is out of range.
     */
    search_by_index(i: number, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
}

/**
 * TriMap
 *
 * A dimensionality reduction technique that preserves both local and global
 * structure using triplets. It is designed to be a more robust alternative
 * to t-SNE and UMAP.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersTriMap>
 * @category Dimensionality Reduction
 */
declare class TriMap<T extends InputType> extends DR<T, ParametersTriMap> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTriMap>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTriMap>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTriMap>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTriMap>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X: T, parameters?: Partial<ParametersTriMap>);
    /**
     * @param {Matrix | null} [pca=null] - Initial Embedding (if null then PCA gets used). Default is `null`
     * @param {import("../knn/KNN.js").KNN<number[] | Float64Array, any> | null} [knn=null] - KNN Object (if null then a KDTree or BallTree gets used, depending on the metric). Default is `null`
     */
    init(pca?: Matrix | null, knn?: KNN<number[] | Float64Array, any> | null): this;
    n_inliers: number | undefined;
    n_outliers: number | undefined;
    n_random: number | undefined;
    knn: KNN<number[] | Float64Array<ArrayBufferLike>, any> | undefined;
    triplets: Matrix | undefined;
    weights: Float64Array<ArrayBuffer> | undefined;
    /** @type {Int32Array} */
    _triplets_int32: Int32Array<ArrayBufferLike> | undefined;
    lr: number | undefined;
    C: number | undefined;
    vel: Matrix | undefined;
    gain: Matrix | undefined;
    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     *
     * @param {number} n_inliers
     * @param {number} n_outliers
     * @param {number} n_random
     */
    _generate_triplets(n_inliers: number, n_outliers: number, n_random: number): {
        triplets: Matrix;
        weights: Float64Array<ArrayBuffer>;
    };
    /**
     * Calculates the log-similarity matrix P.
     *
     * Kept in log space: exponentiating here underflows to zero for anything but the closest
     * neighbours, and every consumer either compares entries or takes their difference, both of
     * which are order-preserving under the logarithm.
     *
     * @private
     * @param {Matrix} knn_distances - Matrix of pairwise knn distances
     * @param {Float64Array} sig - Scaling factor for the distances
     * @param {Matrix} nbrs - Nearest neighbors
     * @returns {Matrix} Pairwise log-similarity matrix
     */
    private _find_p;
    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     *
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix
     *   {@link P}. Row i corresponds to the i-th point.
     * @param {number} n_inliers - Number of inlier points.
     * @param {number} n_outliers - Number of outlier points.
     */
    private _sample_knn_triplets;
    /**
     * Should do the same as np.argsort()
     *
     * @private
     * @param {Float64Array | number[]} A
     */
    private __argsort;
    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are
     * in the {@link rejects}.
     *
     * @private
     * @param {number} n_samples
     * @param {number} max_int
     * @param {number[]} rejects
     */
    private _rejection_sample;
    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     *
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - Nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - Scaling factor for the distances.
     */
    private _find_weights;
    /**
     * Sample uniformly ranom triplets
     *
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    private _sample_random_triplets;
    /**
     * Computes the gradient for updating the embedding.
     *
     * @param {Matrix} Y - The embedding
     * @returns {{ grad: Matrix; loss: number }} The gradient and the current loss.
     */
    _grad(Y: Matrix): {
        grad: Matrix;
        loss: number;
    };
    /**
     * @param {number} max_iteration
     * @returns {T}
     */
    transform(max_iteration?: number): T;
    /**
     * @param {number} max_iteration
     * @returns {Generator<T, T, void>}
     */
    generator(max_iteration?: number): Generator<T, T, void>;
    /**
     * Does the iteration step.
     *
     * @private
     * @param {number} iter
     */
    private _next;
    /**
     * Updates the embedding.
     *
     * @private
     * @param {Matrix} Y
     * @param {number} iter
     * @param {Matrix} grad
     */
    private _update_embedding;
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersTSNE} from "./index.js" */
/**
 * t-SNE (t-Distributed Stochastic Neighbor Embedding)
 *
 * A nonlinear dimensionality reduction technique particularly well-suited
 * for visualizing high-dimensional data in 2D or 3D. Preserves local
 * structure while revealing global patterns.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersTSNE>
 * @category Dimensionality Reduction
 * @see {@link https://lvdmaaten.github.io/tsne/|t-SNE Paper}
 * @see {@link UMAP} for faster alternative with similar results
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const tsne = new druid.TSNE(X, {
 *     perplexity: 30,
 *     epsilon: 10,
 *     d: 2,
 *     seed: 42
 * });
 *
 * const Y = tsne.transform(500); // 500 iterations
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
declare class TSNE<T extends InputType> extends DR<T, ParametersTSNE> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTSNE>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTSNE>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersTSNE>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTSNE>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersTSNE>);
    _iter: number;
    init(): this;
    _ystep: Matrix | undefined;
    _gains: Matrix | undefined;
    _P: Matrix | undefined;
    /**
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {T} The projection.
     */
    transform(iterations?: number): T;
    /**
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {Generator<T, T, void>} - The projection.
     */
    generator(iterations?: number): Generator<T, T, void>;
    /**
     * Performs a optimization step
     *
     * @private
     * @returns {Matrix}
     */
    private next;
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersUMAP} from "./index.js" */
/**
 * Uniform Manifold Approximation and Projection (UMAP)
 *
 * A novel manifold learning technique for dimensionality reduction. UMAP is constructed
 * from a theoretical framework based on Riemannian geometry and algebraic topology.
 * It is often faster than t-SNE while preserving more of the global structure.
 *
 * A given `seed` reproduces an embedding exactly within one engine and library build, but not
 * across browsers or Node versions: the gradient descent is chaotic, so a last-bit difference
 * grows into a visibly different — though equally valid — layout. See {@link DR} for why.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersUMAP>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/1802.03426|UMAP Paper}
 * @see {@link TSNE} for a similar visualization technique
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const umap = new druid.UMAP(X, {
 *     n_neighbors: 15,
 *     min_dist: 0.1,
 *     d: 2,
 *     seed: 42
 * });
 *
 * const Y = umap.transform(500); // 500 iterations
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
declare class UMAP<T extends InputType> extends DR<T, ParametersUMAP> {
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {T}
     */
    static transform<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersUMAP>): T_1;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static generator<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersUMAP>): Generator<T_1, T_1, void>;
    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static transform_async<T_1 extends InputType>(X: T_1, parameters?: Partial<ParametersUMAP>): Promise<T_1>;
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersUMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X: T, parameters?: Partial<ParametersUMAP>);
    _iter: number;
    /**
     * @private
     * @param {number} spread
     * @param {number} min_dist
     * @returns {number[]}
     */
    private _find_ab_params;
    /**
     * @private
     * @param {{ element: Float64Array; index: number; distance: number }[][]} distances
     * @param {number[]} sigmas
     * @param {number[]} rhos
     * @returns {{ element: Float64Array; index: number; distance: number }[][]}
     */
    private _compute_membership_strengths;
    /**
     * @private
     * @param {NaiveKNN<Float64Array> | KDTree<Float64Array> | BallTree<Float64Array>} knn
     * @param {number} k
     * @returns {{
     *     distances: { element: Float64Array; index: number; distance: number }[][];
     *     sigmas: number[];
     *     rhos: number[];
     * }}
     */
    private _smooth_knn_dist;
    /**
     * @private
     * @param {Matrix} X
     * @param {number} n_neighbors
     * @returns {Matrix}
     */
    private _fuzzy_simplicial_set;
    /**
     * @private
     * @param {number} n_epochs
     * @returns {Float32Array}
     */
    private _make_epochs_per_sample;
    /**
     * @private
     * @param {Matrix} graph
     * @returns {{ rows: number[]; cols: number[]; data: number[] }}
     */
    private _tocoo;
    /**
     * Computes all necessary
     *
     * @returns {UMAP<T>}
     */
    init(): UMAP<T>;
    _a: number | undefined;
    _b: number | undefined;
    _graph: Matrix | undefined;
    _head: number[] | undefined;
    _tail: number[] | undefined;
    _weights: number[] | undefined;
    _epochs_per_sample: Float32Array<ArrayBufferLike> | undefined;
    _epochs_per_negative_sample: Float32Array<ArrayBuffer> | undefined;
    _epoch_of_next_sample: Float32Array<ArrayBuffer> | undefined;
    _epoch_of_next_negative_sample: Float32Array<ArrayBuffer> | undefined;
    graph(): {
        cols: number[] | undefined;
        rows: number[] | undefined;
        weights: number[] | undefined;
    };
    /**
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {T}
     */
    transform(iterations?: number): T;
    /**
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {Generator<T, T, void>}
     */
    generator(iterations?: number): Generator<T, T, void>;
    /**
     * @private
     * @param {number} x
     * @returns {number}
     */
    private _clip;
    /**
     * Performs the optimization step.
     *
     * @private
     * @param {Matrix} head_embedding
     * @param {Matrix} tail_embedding
     * @param {number[]} head
     * @param {number[]} tail
     * @returns {Matrix}
     */
    private _optimize_layout;
    /**
     * WASM path for {@link _optimize_layout}.
     *
     * Only the forces run in WASM. The negative samples are drawn here first, in the edge order the
     * kernel walks, so both paths consume the seeded randomizer identically — moving the generator
     * into WASM would break that. The kernel updates one embedding in place, which holds because the
     * sole caller passes the same matrix as head and tail.
     *
     * @private
     * @param {Matrix} embedding
     * @param {number[]} head
     * @param {number[]} tail
     * @param {number} dim
     * @param {number} a
     * @param {number} b
     * @param {number} gamma
     * @param {number} alpha
     * @returns {boolean} True if the epoch was executed in WASM.
     */
    private _optimize_layout_wasm;
    _head_int32: Int32Array<ArrayBuffer> | undefined;
    _tail_int32: Int32Array<ArrayBuffer> | undefined;
    /**
     * @private
     * @returns {Matrix}
     */
    private next;
    _alpha: number | undefined;
}

/**
 * Computes the inner product between two arrays of the same length.
 *
 * @category Linear Algebra
 * @param {number[] | Float64Array} a - Array a.
 * @param {number[] | Float64Array} b - Array b.
 * @returns The inner product between `a` and `b`.
 */
declare function inner_product(a: number[] | Float64Array, b: number[] | Float64Array): number;

/**
 * Computes the QR Decomposition of the Matrix `A` using Gram-Schmidt process.
 *
 * @category Linear Algebra
 * @param {Matrix} A
 * @returns {{ R: Matrix; Q: Matrix }}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */
declare function qr(A: Matrix): {
    R: Matrix;
    Q: Matrix;
};

/**
 * Computes the QR Decomposition of the Matrix `A` with householder transformations.
 *
 * @category Linear Algebra
 * @param {Matrix} A
 * @returns {{ R: Matrix; Q: Matrix }}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */
declare function qr_householder(A: Matrix): {
    R: Matrix;
    Q: Matrix;
};

/** @import { EigenArgs } from "./index.js" */
/**
 * Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.
 *
 * @category Linear Algebra
 * @param {Matrix} A - The Matrix
 * @param {number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {EigenArgs} parameters - Object containing parameterization of the simultanious
 *   poweriteration method.
 * @returns {{ eigenvalues: Float64Array; eigenvectors: Float64Array[] }} The `k` biggest eigenvectors and eigenvalues
 *   of Matrix `A`.
 */
declare function simultaneous_poweriteration(A: Matrix, k?: number, { seed, max_iterations, qr, tol }?: EigenArgs): {
    eigenvalues: Float64Array;
    eigenvectors: Float64Array[];
};

type QRDecomposition = (A: Matrix) => {
    R: Matrix;
    Q: Matrix;
};
type EigenArgs = {
    /**
     * - The number of maxiumum iterations the algorithm should run. Default is `100`
     */
    max_iterations?: number | undefined;
    /**
     * - The seed value or a randomizer used in the algorithm. Default is `1212`
     */
    seed?: number | Randomizer | undefined;
    /**
     * - The QR technique to use. Default is `qr_gramschmidt`
     */
    qr?: QRDecomposition | undefined;
    /**
     * - Tolerated error for stopping criteria. Default is `1e-8`
     */
    tol?: number | undefined;
};

type ParametersLSP = {
    /**
     * - number of neighbors to consider.
     */
    neighbors?: number | undefined;
    /**
     * - number of controlpoints
     */
    control_points?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersFASTMAP = {
    /**
     * - The dimensionality of the projection
     */
    d?: number | undefined;
    /**
     * - The metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - The seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersISOMAP = {
    /**
     * - The number of neighbors ISOMAP should use to project the data.
     */
    neighbors?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - Whether to use classical MDS or SMACOF for the final DR.
     */
    project?: "MDS" | "SMACOF" | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
/**
 * How to weight each pair in the {@link StressMDS} objective.
 *
 * A number is an exponent `q`, giving `w_ij = d_ij^q` — `0` for raw stress, `-1` for Sammon stress,
 * `-2` for elastic scaling / Kamada-Kawai. A matrix supplies the weights directly. A function is
 * called per pair with the target distance and the two indices.
 *
 * A weight of zero, or any non-finite value, drops that pair from the objective — which is how
 * missing or untrusted dissimilarities are expressed.
 */
type WeightSpec = number | Matrix | number[][] | ((d_ij: number, i: number, j: number) => number);
type ParametersStressMDS = {
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance
     * between two points. Pass graph shortest-path distances as `"precomputed"` for a graph layout.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - Pair weighting, see {@link WeightSpec}.
     */
    weights?: WeightSpec | undefined;
    /**
     * - maximum number of gradient steps.
     */
    iterations?: number | undefined;
    /**
     * - stop once the relative stress improvement falls below this.
     */
    epsilon?: number | undefined;
    /**
     * - initial step size. Dimensionless: the gradient is
     * preconditioned by the weighted degree, so this needs no rescaling for the data or the weighting.
     * Adapted by the line search, so it only sets where the search starts.
     */
    learning_rate?: number | undefined;
    /**
     * - starting configuration. `"MDS"` runs
     * classical MDS on the same distances, which is what keeps the non-convex descent out of the poor
     * local minima this objective is known for. `"PCA"` needs the original data, not a precomputed
     * matrix.
     */
    init_DR?: "MDS" | "PCA" | "random" | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
/**
 * {@link ParametersStressMDS} without `weights`, which {@link KKMDS} fixes at `-2`.
 */
type ParametersKKMDS = Omit<ParametersStressMDS, "weights">;
type ParametersLDA = {
    /**
     * - The labels / classes for each data point.
     */
    labels: any[] | Float64Array;
    /**
     * - The dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - The seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersLLE = {
    /**
     * - The number of neighbors for LLE.
     */
    neighbors?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Index used to find the
     * neighbors. If `null`, a KDTree or BallTree is built, depending on the metric. Pass an
     * approximate index such as `HNSW`, `Annoy`, or `NNDescent` to avoid the exact O(N^2) search on
     * large datasets. Default is `null`
     */
    knn?: KNN<number[] | Float64Array<ArrayBufferLike>, any> | null | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersLTSA = {
    /**
     * - The number of neighbors for LTSA.
     */
    neighbors?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Index used to find the
     * neighbors. If `null`, a KDTree or BallTree is built, depending on the metric. Pass an
     * approximate index such as `HNSW`, `Annoy`, or `NNDescent` to avoid the exact O(N^2) search on
     * large datasets. Default is `null`
     */
    knn?: KNN<number[] | Float64Array<ArrayBufferLike>, any> | null | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersMDS = {
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersMINFOTree = {
    /**
     * - Neighbors in the k-NN graph. Defaults to `round(ln N)`, as in the paper.
     */
    k?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - Number of clusters to partition `X` into for the labels field.
     * Required unless `labels` is given — the tree inherits whatever the clustering gets wrong, so
     * there is no safe default. With the default hierarchical clustering the usable range is
     * `2 … N-1`, the number of merges a dendrogram cut can separate; `"kmeans"` has no such limit.
     */
    clusters?: number | undefined;
    /**
     * - How to obtain the labels when
     * `labels` is not given. `"hierarchical"` uses Ward linkage, matching the paper's experiments.
     */
    clustering?: "hierarchical" | "kmeans" | undefined;
    /**
     * - Precomputed labels, one per
     * row of `X`. Bypasses the clustering step. Values may be of any type; they are remapped to
     * `0 … q-1`.
     */
    labels?: any[] | Int32Array<ArrayBufferLike> | Float64Array<ArrayBufferLike> | null | undefined;
    /**
     * - Shrinkage applied to intra-cluster edge weights. Defaults to the
     * golden ratio conjugate `(√5-1)/2 ≈ 0.618`, which the paper picks for interpretability rather
     * than by tuning.
     */
    alpha?: number | undefined;
    /**
     * - Floor on the curvature denominator, `S = -ψ / (φ + epsilon)`,
     * needed because φ and ψ both vanish for interior points as β grows. Default matches the author's
     * reference implementation. See `MINFOTree._information_curvature` on how the interior/boundary
     * curvature ordering depends on β.
     */
    epsilon?: number | undefined;
    /**
     * - How to lay the tree out. `"MDS"`
     * stops after the classical-MDS warm start, which is far cheaper and often enough.
     */
    layout?: "MDS" | "kamada_kawai" | undefined;
    /**
     * - Maximum Kamada-Kawai gradient steps.
     */
    iterations?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersPCA = {
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
    /**
     * - Parameters for the eigendecomposition algorithm.
     */
    eig_args?: Partial<EigenArgs> | undefined;
};
type ParametersSAMMON<K extends keyof ChooseDR> = {
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - Either "PCA" or "MDS", with which SAMMON initialiates the projection.
     */
    init_DR?: K | undefined;
    /**
     * - Parameters for the "init"-DR method.
     */
    init_parameters?: ChooseDR[K] | undefined;
    /**
     * - learning rate for gradient descent.
     */
    magic?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersSMACOF = {
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - maximum number of iterations.
     */
    iterations?: number | undefined;
    /**
     * - tolerance for stress difference.
     */
    epsilon?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersSQDMDS = {
    d?: number | undefined;
    metric?: Metric | "precomputed" | undefined;
    /**
     * - Percentage of iterations using exaggeration phase.
     */
    decay_start?: number | undefined;
    /**
     * - Controls the decay of the learning parameter.
     */
    decay_cte?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersTopoMap = {
    /**
     * = euclidean - The metric which defines the distance between
     * two points.
     */
    metric: Metric;
    /**
     * = 1212 - The seed for the random number generator.
     */
    seed: number;
};
type ParametersTriMap = {
    /**
     * - Temperature of the tempered log applied to the triplet
     * weights. `1` recovers the ordinary logarithm; lower values compress large weights harder.
     */
    weight_temp?: number | undefined;
    /**
     * - number of inliers.
     */
    n_inliers?: number | undefined;
    /**
     * - number of outliers.
     */
    n_outliers?: number | undefined;
    /**
     * - number of random points.
     */
    n_random?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - learning rate of the delta-bar-delta optimizer.
     */
    lr?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersPaCMAP = {
    /**
     * - Number of nearest neighbors forming the attractive pairs.
     */
    n_neighbors?: number | undefined;
    /**
     * - Mid-near pairs per point, as a fraction of `n_neighbors`.
     */
    MN_ratio?: number | undefined;
    /**
     * - Further pairs per point, as a multiple of `n_neighbors`.
     */
    FP_ratio?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - learning rate of the Adam optimizer.
     */
    lr?: number | undefined;
    /**
     * - Iterations in each of the three phases.
     */
    num_iters?: number[] | undefined;
    /**
     * - Index used to find the
     * neighbors. If `null`, an exact blocked search runs instead, which is faster than a tree at the
     * `n_neighbors + 50` candidates the density rescaling needs. Pass an approximate index such as
     * `HNSW`, `Annoy`, or `NNDescent` to avoid the O(N^2) search on large datasets. Default is `null`
     */
    knn?: KNN<number[] | Float64Array<ArrayBufferLike>, any> | null | undefined;
    /**
     * - Reduce inputs wider than 100 dimensions to 100 via PCA
     * before the search and the initialization.
     */
    apply_pca?: boolean | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersLocalMAP = {
    /**
     * - Number of nearest neighbors forming the attractive pairs.
     */
    n_neighbors?: number | undefined;
    /**
     * - Mid-near pairs per point, as a fraction of `n_neighbors`.
     */
    MN_ratio?: number | undefined;
    /**
     * - Further pairs per point, as a multiple of `n_neighbors`.
     */
    FP_ratio?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | undefined;
    /**
     * - learning rate of the Adam optimizer.
     */
    lr?: number | undefined;
    /**
     * - Iterations in each of the three phases.
     */
    num_iters?: number[] | undefined;
    /**
     * - Embedding distance below which a point may be redrawn as
     * a further pair in phase 3. Also sets the local attraction scale, `low_dist_thres / 2`.
     */
    low_dist_thres?: number | undefined;
    /**
     * - Index used to find the
     * neighbors. If `null`, an exact blocked search runs instead. Default is `null`
     */
    knn?: KNN<number[] | Float64Array<ArrayBufferLike>, any> | null | undefined;
    /**
     * - Reduce inputs wider than 100 dimensions to 100 via PCA
     * before the search and the initialization.
     */
    apply_pca?: boolean | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersTSNE = {
    /**
     * - perplexity.
     */
    perplexity?: number | undefined;
    /**
     * - learning parameter.
     */
    epsilon?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};
type ParametersUMAP = {
    /**
     * - size of the local neighborhood.
     */
    n_neighbors?: number | undefined;
    /**
     * - number of nearest neighbors connected in the local neighborhood.
     */
    local_connectivity?: number | undefined;
    /**
     * - controls how tightly points get packed together.
     */
    min_dist?: number | undefined;
    /**
     * - the dimensionality of the projection.
     */
    d?: number | undefined;
    /**
     * - the metric which defines the distance between two points in the high-dimensional space.
     */
    metric?: Metric | "precomputed" | undefined;
    /**
     * - The effective scale of embedded points.
     */
    _spread?: number | undefined;
    /**
     * - Interpolate between union and intersection.
     */
    _set_op_mix_ratio?: number | undefined;
    /**
     * - Weighting applied to negative samples.
     */
    _repulsion_strength?: number | undefined;
    /**
     * - The number of negative samples per positive sample.
     */
    _negative_sample_rate?: number | undefined;
    /**
     * - The number of training epochs.
     */
    _n_epochs?: number | undefined;
    /**
     * - The initial learning rate for the optimization.
     */
    _initial_alpha?: number | undefined;
    /**
     * - the seed for the random number generator.
     */
    seed?: number | undefined;
};

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersAnnoy } from "./index.js" */
/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} AnnoyNode
 * @property {boolean} isLeaf - Whether this is a leaf node
 * @property {number[]} indices - Indices of points in this node (leaf) or children (internal)
 * @property {number[]} normal - Hyperplane normal vector (internal nodes only)
 * @property {number} offset - Hyperplane offset (internal nodes only)
 * @property {AnnoyNode<T> | null} left - Left child (internal nodes only)
 * @property {AnnoyNode<T> | null} right - Right child (internal nodes only)
 */
/**
 * Annoy-style (Approximate Nearest Neighbors Oh Yeah) implementation using Random Projection Trees.
 *
 * This implementation builds multiple random projection trees where each tree randomly selects
 * two points and splits the space based on a hyperplane equidistant between them.
 *
 * Key features:
 * - Multiple random projection trees for better recall
 * - Each tree uses random hyperplanes for splitting
 * - Priority queue search for better recall
 * - Combines results from all trees
 *
 * Best suited for:
 * - High-dimensional data
 * - Approximate nearest neighbor search
 * - Large datasets
 * - When high recall is needed with approximate methods
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersAnnoy>
 * @see {@link https://github.com/spotify/annoy}
 * @see {@link https://erikbern.com/2015/09/24/nearest-neighbors-and-vector-models-epilogue-curse-of-dimensionality.html}
 */
declare class Annoy<T extends number[] | Float64Array> extends KNN<T, ParametersAnnoy> {
    /**
     * Creates a new Annoy-style index with random projection trees.
     *
     * @param {T[]} elements - Elements to index
     * @param {Partial<ParametersAnnoy>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(elements: T[], parameters?: Partial<ParametersAnnoy>);
    _metric: Metric;
    _numTrees: number;
    _maxPointsPerLeaf: number;
    _seed: number;
    /**
     * @private
     * @type {AnnoyNode<T>[]}
     */
    private _trees;
    /**
     * Get the number of trees in the index.
     * @returns {number}
     */
    get num_trees(): number;
    /**
     * Get the total number of nodes in all trees.
     * @returns {number}
     */
    get num_nodes(): number;
    /**
     * @private
     * @param {any} node
     * @returns {number}
     */
    private _countNodes;
    /**
     * Add elements to the Annoy index.
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements: T[]): this;
    /**
     * Build all random projection trees.
     * @private
     */
    private _buildTrees;
    /**
     * Recursively build a random projection tree.
     * @private
     * @param {number[]} indices - Indices of elements to include
     * @returns {AnnoyNode<T>}
     */
    private _buildTreeRecursive;
    /**
     * Compute distance from point to hyperplane.
     * @private
     * @param {T} point
     * @param {number[]} normal
     * @param {number} offset
     * @returns {number} Signed distance (positive = right side, negative = left side)
     */
    private _distanceToHyperplane;
    /**
     * Search for k approximate nearest neighbors.
     * @param {T} query
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(query: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
    /**
     * Search tree using priority queue for better recall.
     * Explores nodes in order of distance to hyperplane.
     * @private
     * @param {AnnoyNode<T>} node
     * @param {T} query
     * @param {Set<number>} candidates
     * @param {number} maxCandidates
     */
    private _searchTreePriority;
    /**
     * Alias for search_by_index for backward compatibility.
     *
     * @param {number} i - Index of the query element
     * @param {number} [k=5] - Number of nearest neighbors to return
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_index(i: number, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
}

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersBallTree } from "./index.js" */
/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} ElementWithIndex
 * @property {number} index
 * @property {T} element
 */
/**
 * Ball Tree for efficient nearest neighbor search.
 *
 * A Ball Tree is a metric tree that partitions points into a nested set of
 * hyperspheres (balls). It is particularly effective for high-dimensional
 * data and supports any valid metric.
 *
 * Every ball stores a center and the distance from it to its furthest member, which bounds the
 * subtree from below by `d(t, center) - radius`. That bound holds for any metric by the triangle
 * inequality, so the center does not need to be one of the indexed points — it is the centroid,
 * which keeps radii tighter than an arbitrary member would.
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersBallTree>
 */
declare class BallTree<T extends number[] | Float64Array> extends KNN<T, ParametersBallTree> {
    /**
     * Generates a BallTree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the BallTree
     * @param {Partial<ParametersBallTree>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     */
    constructor(elements: T[], parameters?: Partial<ParametersBallTree>);
    /**
     * @private
     * @type {BallTreeNode<T> | BallTreeLeaf<T>}
     */
    private _root;
    /** @returns {Metric} */
    get _metric(): Metric;
    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @returns {Float64Array} Componentwise mean of `elements`.
     */
    private _centroid;
    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @returns {BallTreeNode<T> | BallTreeLeaf<T>} Root of balltree.
     */
    private _construct;
    /**
     * @private
     * @param {ElementWithIndex<T>[]} B
     * @returns {number}
     */
    private _greatest_spread;
    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
    /**
     * @private
     * @param {T} t - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {Heap<{ point: ElementWithIndex<T>; distance: number }>} Q - Heap consists of the currently found `k`
     *   nearest neighbors.
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} B
     * @param {number} [known_distance] - Distance from `t` to this ball's center, when the caller
     *   already measured it while ordering its children.
     */
    private _search;
}

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersHNSW } from "./index.js" */
/**
 * @typedef {Object} Layer
 * @property {number} l_c - Layer number
 * @property {number[]} point_indices - Global indices of points in this layer
 * @property {Map<number, number[]>} edges - Global index -> array of connected global indices
 */
/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} Candidate
 * @property {T} element - The actual data point
 * @property {number} index - Global index in the dataset
 * @property {number} distance - Distance from query
 */
/**
 * Hierarchical Navigable Small World (HNSW) graph for approximate nearest neighbor search.
 *
 * HNSW builds a multi-layer graph structure where each layer is a navigable small world graph.
 * The top layers serve as "highways" for fast traversal, while lower layers provide accuracy.
 * Each element is assigned to a random level, allowing logarithmic search complexity.
 *
 * Key parameters:
 * - `m`: Controls the number of connections per element (affects accuracy/memory)
 * - `ef_construction`: Controls the quality of the graph during construction (higher = better but slower)
 * - `ef`: Controls the quality of search (higher = better recall but slower)
 *
 * This is an *approximate* index. Recall above 95% is expected for any `m`; raise `ef` or
 * `ef_construction` if a dataset needs more. Queries overtake the exact {@link KDTree} at around
 * 10 000 points, so below a few thousand prefer {@link spatial_tree}, which is exact and quicker.
 *
 * Based on:
 * - "Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs"
 *   by Malkov & Yashunin (2016)
 * - "Approximate Nearest Neighbor Search on High Dimensional Data"
 *   by Li et al. (2019)
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersHNSW>
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const points = [[1, 2], [3, 4], [5, 6], [7, 8]];
 * const hnsw = new druid.HNSW(points, {
 *     metric: druid.euclidean,
 *     m: 16,
 *     ef_construction: 200
 * });
 *
 * const query = [2, 3];
 * const neighbors = hnsw.search(query, 2);
 * // [{ element: [1, 2], index: 0, distance: 1.41 }, ...]
 */
declare class HNSW<T extends number[] | Float64Array> extends KNN<T, ParametersHNSW> {
    /**
     * Creates a new HNSW index.
     *
     * @param {T[]} points - Initial points to add to the index
     * @param {Partial<ParametersHNSW>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(points: T[], parameters?: Partial<ParametersHNSW>);
    /** @type {Metric} */
    _metric: Metric;
    /** @type {Function} */
    _select: Function;
    /**
     * @private
     * @type {Map<number, Layer>}
     */
    private _graph;
    /** @type {number} */
    _next_index: number;
    /** @type {number} */
    _m: number;
    /** @type {number} */
    _ef_construction: number;
    /** @type {number} */
    _ef: number;
    /** @type {number} */
    _m0: number;
    /** @type {number} */
    _mL: number;
    /** @type {number} - Current maximum layer in the graph */
    _L: number;
    /** @type {number[] | null} - Entry point indices for search */
    _ep: number[] | null;
    /** @private @type {number} */
    private _search_id;
    /** @private @type {Uint32Array} */
    private _visited_stamps;
    /**
     * Add a single element to the index.
     *
     * @param {T} element - Element to add
     * @returns {HNSW<T>} This instance for chaining
     */
    addOne(element: T): HNSW<T>;
    /**
     * Add multiple elements to the index.
     *
     * @param {T[]} new_elements - Elements to add
     * @returns {HNSW<T>} This instance for chaining
     */
    add(new_elements: T[]): HNSW<T>;
    /**
     * Select neighbors using the heuristic approach.
     *
     * The heuristic extends candidates with their neighbors and selects
     * points that are closer to the query than to already selected points.
     * This maintains graph connectivity better than simple selection.
     *
     * @private
     * @param {T} q - Query element
     * @param {Candidate<T>[]} candidates - Candidate elements with distances
     * @param {number} M - Maximum number of neighbors to return
     * @param {number} l_c - Layer number
     * @param {boolean} [extend_candidates=false] - Whether to extend candidates with their neighbors
     * @param {boolean} [keep_pruned_connections=false] - Whether to add pruned connections back if needed
     * @returns {number[]} Selected neighbor indices, nearest first
     */
    private _select_heuristic;
    /**
     * Select neighbors using simple distance-based selection.
     *
     * Simply returns the M closest candidates to the query.
     *
     * @private
     * @param {T} q - Query element
     * @param {Candidate<T>[]} C - Candidate elements with distances
     * @param {number} M - Maximum number of neighbors to return
     * @returns {number[]} M nearest candidate indices
     */
    private _select_simple;
    /**
     * Search a single layer for nearest neighbors.
     *
     * Implements the greedy search algorithm: start from entry points,
     * always expand the closest unvisited candidate, maintain a list
     * of the ef closest found neighbors.
     *
     * @private
     * @param {T} q - Query element
     * @param {number[] | null} ep_indices - Entry point indices
     * @param {number} ef - Number of nearest neighbors to find
     * @param {number} l_c - Layer number to search
     * @returns {Candidate<T>[]} ef nearest neighbors found with their distances
     */
    private _search_layer;
    /**
     * Fallback linear search when graph search fails
     * @private
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @returns {Candidate<T>[]}
     */
    private _linear_search;
    /**
     * Iterator for searching the HNSW graph layer by layer.
     *
     * Yields intermediate results at each layer for debugging or visualization.
     *
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @param {number?} [ef] - Size of dynamic candidate list
     * @yields {{layer: number, candidates: Candidate[]}}
     */
    search_iter(q: T, K: number, ef?: number | null): Generator<{
        layer: number;
        candidates: {
            element: T;
            index: number;
            distance: number;
        }[];
    }, void, unknown>;
    /**
     * Get the number of elements in the index.
     *
     * @returns {number} Number of elements
     */
    get size(): number;
    /**
     * Get the number of layers in the graph.
     *
     * @returns {number} Number of layers
     */
    get num_layers(): number;
    /**
     * Get an element by its index.
     *
     * @param {number} index - Element index
     * @returns {T} The element at the given index
     */
    get_element(index: number): T;
}

/** @import { Metric } from "../metrics/index.js" */
/** @import { ParametersKDTree } from "./index.js" */
/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} ElementWithIndex
 * @property {number} index
 * @property {T} element
 */
/**
 * KD-Tree (K-dimensional Tree) for efficient nearest neighbor search.
 *
 * KD-Trees partition k-dimensional space by recursively splitting along coordinate axes.
 * At each level, the tree splits points based on the median of the coordinate with the largest spread.
 * This creates a balanced binary tree structure that enables efficient O(log n) search on average.
 *
 * Best suited for:
 * - Low to moderate dimensional data (d < 20-30)
 * - When exact nearest neighbors are needed
 * - When dimensionality is not too high
 *
 * Performance degrades in high dimensions (curse of dimensionality) where approximate
 * methods like HNSW or LSH become more effective.
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersKDTree>
 * @see {@link https://en.wikipedia.org/wiki/K-d_tree}
 */
declare class KDTree<T extends number[] | Float64Array> extends KNN<T, ParametersKDTree> {
    /**
     * Generates a KD-Tree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KD-Tree
     * @param {Partial<ParametersKDTree>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(elements: T[], parameters?: Partial<ParametersKDTree>);
    /**
     * @private
     * @type {KDTreeNode<T> | KDTreeLeaf<T> | null}
     */
    private _root;
    /** @returns {Metric} */
    get _metric(): Metric;
    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @param {number} depth - Current depth in the tree (determines splitting axis)
     * @returns {KDTreeNode<T> | KDTreeLeaf<T> | null} Root of KD-Tree.
     */
    private _construct;
    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
    /**
     * @private
     * @param {T} target - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} node - Current node.
     * @param {Heap<{ point: ElementWithIndex<T>; distance: number }>} best - Heap of k best found so far.
     */
    private _search_recursive;
}

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
declare class LSH<T extends number[] | Float64Array> extends KNN<T, ParametersLSH> {
    /**
     * Creates a new LSH index.
     *
     * @param {T[]} elements - Elements to index
     * @param {Partial<ParametersLSH>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(elements: T[], parameters?: Partial<ParametersLSH>);
    _metric: Metric;
    _numHashTables: number;
    _numHashFunctions: number;
    _seed: number;
    /** @type {Map<string, number[]>[]} */
    _hashTables: Map<string, number[]>[];
    /** @type {Float64Array[][]} */
    _projections: Float64Array[][];
    /** @type {number[][]} */
    _offsets: number[][];
    /** @type {number} */
    _dim: number;
    /** @type {number} */
    _bucketWidth: number;
    /**
     * Whether the projection geometry still has to be derived from real data.
     *
     * Both things the hash depends on — `_dim` and the bucket width — come from the elements,
     * and an index built empty has none yet. Building hash functions against the placeholder
     * would fix `_dim` at 1 and the width at the `n < 2` fallback, giving one-element projection
     * vectors; `_computeHash` then walks `element.length` components against them, reads past
     * the end, and every real point added later hashes to `"NaN,NaN,…"`. That lands the whole
     * dataset in a single bucket per table and turns every query into a full scan — right
     * answers, no index. So the hash functions wait for `add` to supply real data.
     *
     * @private
     * @type {boolean}
     */
    private _awaiting_data;
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
    private _estimateBucketWidth;
    /**
     * Initialize random projection vectors for all hash tables.
     * @private
     */
    private _initializeHashFunctions;
    /**
     * Compute hash signature for an element using random projections.
     * @private
     * @param {T} element
     * @param {number} tableIndex
     * @returns {string} Hash signature
     */
    private _computeHash;
    /**
     * Add elements to the LSH index.
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements: T[]): this;
    /**
     * Search for k approximate nearest neighbors.
     * @param {T} query
     * @param {number} [k=5]
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(query: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
}

/** @import { ParametersNaiveKNN } from "./index.js" */
/** @import { Metric } from "../metrics/index.js" */
/**
 * Naive KNN implementation performing an exhaustive scan.
 *
 * Every query measures the distance to all `N` elements and then selects the `k` smallest with
 * {@link quickselect}, which is O(N) on average — cheaper than sorting all `N` and far cheaper than
 * the N heaps of N entries this class used to build up front.
 *
 * The N x N distance matrix is only materialized when it actually pays off: a `"precomputed"` index
 * is handed one directly, and {@link NaiveKNN#search_by_index} builds one lazily on first use so
 * that building a whole kNN graph costs a single pass over the pairs instead of one per query.
 * Callers that only ever use {@link NaiveKNN#search} never allocate it.
 *
 * Best suited for small datasets, or when a distance matrix is already available. For larger `N`
 * prefer an approximate index such as {@link HNSW}, {@link Annoy}, or {@link NNDescent}.
 *
 * @template {number[] | Float64Array} T
 * @category KNN
 * @class
 * @extends KNN<T, ParametersNaiveKNN>
 */
declare class NaiveKNN<T extends number[] | Float64Array> extends KNN<T, ParametersNaiveKNN> {
    /**
     * Generates a KNN list with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KNN list
     * @param {Partial<ParametersNaiveKNN>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     */
    constructor(elements: T[], parameters?: Partial<ParametersNaiveKNN>);
    /**
     * Number of indexed elements.
     *
     * @type {number}
     */
    _N: number;
    /**
     * Pairwise distances. Supplied by the caller when `metric` is `"precomputed"`, built on demand
     * by {@link NaiveKNN#search_by_index} otherwise, and `null` until then.
     *
     * @type {Matrix | null}
     */
    _D: Matrix | null;
    /**
     * Reads the element stored at `i`, which may live in a `Matrix` or a plain array.
     *
     * @private
     * @param {number} i
     * @returns {T}
     */
    private _element_at;
    /**
     * Selects the `k` elements closest to a query from its distances to every element.
     *
     * QuickSelect partitions in O(N) average time, after which only the `k` selected entries are
     * sorted. Ties break on the element index so the result never depends on the pivots drawn.
     *
     * @private
     * @param {ArrayLike<number>} distances - Distance from the query to each element, by index.
     * @param {number} k - Number of neighbors to return.
     * @param {number} [exclude=-1] - Index to leave out, or `-1` to keep every element. Default is `-1`
     * @returns {{ element: T; index: number; distance: number }[]} The `k` nearest, closest first.
     */
    private _k_smallest;
    /**
     * Returns the pairwise distance matrix, computing it once on first use.
     *
     * @private
     * @returns {Matrix}
     */
    private _distance_matrix;
    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
}

/** @import {ParametersNNDescent} from "./index.js" */
/**
 *
 * @template {number[] | Float64Array} T
 * @typedef {Object} NNDescentElement
 * @property {T} value
 * @property {number} index
 * @property {boolean} flag
 */
/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} NNDescentNeighbor
 * @property {T} value
 * @property {number} index
 * @property {number} distance
 * @property {boolean} [flag]
 */
/**
 * NN-Descent
 *
 * An efficient graph-based approximate nearest neighbor search algorithm.
 * It works by iteratively improving a neighbor graph using the fact that
 * "neighbors of neighbors are likely to be neighbors".
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersNNDescent>
 * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf|NN-Descent Paper}
 */
declare class NNDescent<T extends number[] | Float64Array> extends KNN<T, ParametersNNDescent> {
    /**
     * @param {T[]} elements - Called V in paper.
     * @param {Partial<ParametersNNDescent>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     */
    constructor(elements: T[], parameters?: Partial<ParametersNNDescent>);
    /**
     * @private
     * @type {KNNHeap<T>[]}
     */
    private _B;
    /**
     * @private
     * @type {NNDescentNeighbor<T>[][]}
     */
    private nn;
    _N: number;
    _sample_size: number;
    _nndescent_elements: {
        value: T;
        index: number;
        flag: boolean;
    }[];
    /**
     * Samples Array A with sample size.
     *
     * @private
     * @template U
     * @param {U[]} A
     * @returns {U[]}
     */
    private _sample;
    /**
     * Offers each of `u1`, `u2` to the other's neighborhood, under the distance between them.
     *
     * The distance has to be computed here. A neighbor object carries the distance to the point
     * whose list it came from, which says nothing about how far apart `u1` and `u2` are; reusing it
     * meant no new distance was ever evaluated and the heaps were ordered by unrelated values.
     *
     * @private
     * @param {KNNHeap<T>[]} B
     * @param {NNDescentNeighbor<T>} u1
     * @param {NNDescentNeighbor<T>} u2
     * @returns {number} How many neighborhoods changed.
     */
    private _join;
    /**
     * @private
     * @param {KNNHeap<T>} B
     * @param {NNDescentNeighbor<T>} u
     * @returns {number}
     */
    private _update;
    /**
     * @private
     * @param {(KNNHeap<T> | null)[]} B
     * @returns {NNDescentNeighbor<T>[][]}
     */
    private _reverse;
    /**
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements: T[]): this;
    /**
     * @param {T} x
     * @param {number} [k=5] Default is `5`
     * @returns {{ element: T, index: number; distance: number }[]}
     */
    search(x: T, k?: number): {
        element: T;
        index: number;
        distance: number;
    }[];
}

/**
 * Picks the fastest spatial index that is still correct for the given metric.
 *
 * Returns a {@link KDTree} for metrics it can prune soundly, and a {@link BallTree} — which only
 * relies on the triangle inequality and therefore works for any metric — otherwise.
 *
 * @template {number[] | Float64Array} T
 * @param {T[]} elements - Elements to index.
 * @param {{ metric: Metric; seed: number }} parameters - Metric and seed for the tree.
 * @returns {KDTree<T> | BallTree<T>} The constructed index.
 */
declare function spatial_tree<T extends number[] | Float64Array>(elements: T[], parameters: {
    metric: Metric;
    seed: number;
}): KDTree<T> | BallTree<T>;

type ParametersAnnoy = {
    /**
     * - Metric to use: (a, b) => distance. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Number of random projection trees to build. Default is `10`
     */
    numTrees: number;
    /**
     * - Maximum points per leaf node. Default is `10`
     */
    maxPointsPerLeaf: number;
    /**
     * - Seed for random number generator. Default is `1212`
     */
    seed: number;
};
type ParametersBallTree = {
    /**
     * - Metric to use: (a, b) => distance. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Members below which a ball stops splitting and is scanned
     * linearly. Default is `32`
     */
    leaf_size: number;
    /**
     * - Seed for random number generator. Default is `1212`
     */
    seed: number;
};
type ParametersHNSW = {
    /**
     * - Metric to use: (a, b) => distance. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Use heuristics or naive selection. Default is `true`
     */
    heuristic: boolean;
    /**
     * - Number of connections a newly inserted element creates on each layer it
     * joins, layer 0 included. Default is `16`
     */
    m: number;
    /**
     * - Size of candidate list during construction. Default is `200`
     */
    ef_construction: number;
    /**
     * - Upper bound on connections at the ground layer (layer 0). Unlike
     * `m` this is not a budget for the inserted element, it is the cap enforced once reverse
     * connections from later insertions have accumulated. Default is `2 * m`
     */
    m0: number | null;
    /**
     * - Normalization factor for level generation. Default is `1 / Math.log(m)`
     */
    mL: number | null;
    /**
     * - Seed for random number generator. Default is `1212`
     */
    seed: number;
    /**
     * - Size of candidate list during search. Default is `50`
     */
    ef: number;
};
type ParametersKDTree = {
    /**
     * - Metric to use: (a, b) => distance. Default is `euclidean`
     */
    metric: Metric;
    seed: number;
};
type ParametersLSH = {
    /**
     * - Metric to use: (a, b) => distance. Default is `euclidean`
     */
    metric: Metric;
    /**
     * - Number of hash tables. Default is `10`
     */
    numHashTables: number;
    /**
     * - Number of hash functions per table. Default is `10`
     */
    numHashFunctions: number;
    /**
     * - Quantization width w. Must match the scale of the data;
     * estimated from the sampled distance distribution when `null`. Default is `null`
     */
    bucketWidth: number | null;
    /**
     * - Seed for random number generator. Default is `1212`
     */
    seed: number;
};
type ParametersNaiveKNN = {
    /**
     * Is either precomputed or a function to use: (a, b) => distance
     */
    metric?: Metric | "precomputed" | undefined;
    seed?: number | undefined;
};
type ParametersNNDescent = {
    /**
     * - Called sigma in paper. Default is `euclidean`
     */
    metric: Metric;
    /**
     * =10 - Number of samples. Default is `10`
     */
    samples: number;
    /**
     * = 1 - Sample rate. Default is `1`
     */
    rho: number;
    /**
     * = 0.001 - Precision parameter. Default is `0.001`
     */
    delta: number;
    /**
     * = 1212 - Seed for the random number generator. Default is `1212`
     */
    seed: number;
};

/**
 * Numerical stable summation with the Kahan summation algorithm.
 *
 * @category Numerical
 * @param {number[] | Float64Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
 */
declare function kahan_sum(summands: number[] | Float64Array): number;

/**
 * Numerical stable summation with the Neumair summation algorithm.
 *
 * Deliberately not WASM accelerated: the compensation term makes each step depend on the previous
 * one, so the kernel cannot vectorise, and it measured slower than this loop at every input size —
 * the argument copy is pure overhead. A `neumair_sum_f64` kernel existed on that basis alone, never
 * called by anything; it and its benchmark row were removed once the measurement had been made.
 *
 * @category Numerical
 * @param {number[] | Float64Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
declare function neumair_sum(summands: number[] | Float64Array): number;

/**
 * @template {Float64Array | number[]} T
 * @category Optimization
 * @param {(d: T) => number} f
 * @param {T} x0
 * @param {number} [max_iter=300] Default is `300`
 * @returns {T}
 * @see http://optimization-js.github.io/optimization-js/optimization.js.html#line438
 */
declare function powell<T extends Float64Array | number[]>(f: (d: T) => number, x0: T, max_iter?: number): T;

/**
 * Enables or disables WASM acceleration at runtime.
 *
 * With WASM disabled every accelerated function takes its JS fallback path. This exists so the two
 * implementations can be compared in tests and benchmarks.
 *
 * For most of the library the two paths agree to within floating point tolerance — the kernels
 * accumulate in a different order than the scalar loops, and `Math.pow`/`Math.exp` are only
 * "implementation-approximated" in ECMAScript, so the last bit may differ. `PCA`, `MDS`, `TSNE` and
 * `TriMap` come out bit-identical in practice.
 *
 * {@link UMAP} and {@link SAMMON} are the exception: their optimisers are chaotic, so a last-bit
 * difference grows into a visibly different (though equally valid) layout. See the note on
 * reproducibility in those classes.
 *
 * @param {boolean} enabled - Whether kernels may run in WASM.
 * @returns {boolean} The new state.
 * @example
 * import { setWasmEnabled } from "@saehrimnir/druidjs";
 * setWasmEnabled(false); // force the pure JS implementations
 */
declare function setWasmEnabled(enabled: boolean): boolean;
/**
 * Returns whether WASM acceleration is available in current environment.
 * @returns {boolean}
 */
declare function isWasmAvailable(): boolean;
/**
 * Returns whether SharedArrayBuffer and WASM multi-threading is supported in current environment (Node.js, Bun, Deno, or Web Browser with COOP/COEP headers).
 * @returns {boolean}
 */
declare function isWasmThreadsSupported(): boolean;

/**
 * Whether row-range kernels can be split across workers in this environment.
 *
 * @returns {boolean} True only on Node, where `Atomics.wait` may block the calling thread.
 * @example
 * import { parallel_available } from "@saehrimnir/druidjs";
 * parallel_available(); // false in a browser
 */
declare function parallel_available(): boolean;
/**
 * Stops the workers. They are unreferenced, so this is only needed to release them early.
 *
 * @returns {void}
 */
declare function terminate_pool(): void;

type InputType = Matrix | Float64Array[] | number[][];
declare const version: string;

export { Annoy, BallTree, CURE, DisjointSet, FASTMAP, HNSW, Heap, HierarchicalClustering, ISOMAP, KDTree, KKMDS, KMeans, KMedoids, LDA, LLE, LSH, LSP, LTSA, LocalMAP, MDS, MINFOTree, Matrix, MeanShift, NNDescent, NaiveKNN, OPTICS, PCA, PaCMAP, Randomizer, SAMMON, SMACOF, SQDMDS, StressMDS, TSNE, TopoMap, TriMap, UMAP, WEIGHTS_ELASTIC, WEIGHTS_SAMMON, WEIGHTS_UNIFORM, XMeans, bray_curtis, canberra, chebyshev, cosine, distance_matrix, euclidean, euclidean_squared, goodman_kruskal, hamming, haversine, inner_product, isWasmAvailable, isWasmThreadsSupported, jaccard, k_nearest_neighbors, kahan_sum, linspace, manhattan, max, min, minimum_spanning_tree, neumair_sum, norm, normalize, parallel_available, powell, qr, qr_householder, quickselect, quickselectByAxis, setWasmEnabled, simultaneous_poweriteration, sokal_michener, spatial_tree, terminate_pool, version, wasserstein, yule };
export type { Comparator, EigenArgs, InputType, Metric, ParametersAnnoy, ParametersBallTree, ParametersCURE, ParametersFASTMAP, ParametersHNSW, ParametersHierarchicalClustering, ParametersISOMAP, ParametersKDTree, ParametersKKMDS, ParametersKMeans, ParametersKMedoids, ParametersLDA, ParametersLLE, ParametersLSH, ParametersLSP, ParametersLTSA, ParametersLocalMAP, ParametersMDS, ParametersMINFOTree, ParametersMeanShift, ParametersNNDescent, ParametersNaiveKNN, ParametersOptics, ParametersPCA, ParametersPaCMAP, ParametersSAMMON, ParametersSMACOF, ParametersSQDMDS, ParametersStressMDS, ParametersTSNE, ParametersTopoMap, ParametersTriMap, ParametersUMAP, ParametersXMeans, QRDecomposition, WeightSpec };
//# sourceMappingURL=druid.d.ts.map
