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
declare function k_nearest_neighbors(
  A: Matrix,
  k: number,
  metric?: Metric | "precomputed",
): {
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
 * @category Utils
 * @class
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
   * Mersenne Twister random number generator.
   *
   * @param {number} [_seed=new Date().getTime()] - The seed for the random number generator. If `_seed == null` then
   *   the actual time gets used as seed. Default is `new Date().getTime()`
   * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
   */
  constructor(_seed?: number);
  _N: number;
  _M: number;
  _MATRIX_A: number;
  _UPPER_MASK: number;
  _LOWER_MASK: number;
  /** @type {number[]} */
  _mt: number[];
  /** @type {number} */
  _mti: number;
  /** @type {number} */
  _seed: number;
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
  gauss_random(): number;
  _val: number | null | undefined;
  /**
   * @template T Returns samples from an input Matrix or Array.
   * @param {T[]} A - The input Matrix or Array.
   * @param {number} n - The number of samples.
   * @returns {T[]} A random selection form `A` of `n` samples.
   */
  choice<T>(A: T[], n: number): T[];
}

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
  static solve(
    A:
      | Matrix
      | {
          L: Matrix;
          U: Matrix;
        },
    b: Matrix,
  ): Matrix;
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
  static SVD(
    M: Matrix,
    k?: number,
  ): {
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
  get_block(
    start_row: number,
    start_col: number,
    end_row?: number | null,
    end_col?: number | null,
  ): Matrix;
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
  _apply_rowwise_array(
    values: number[] | Float64Array,
    f: (d: number, v: number) => number,
  ): Matrix;
  /**
   * @param {number[] | Float64Array} values
   * @param {(d: number, v: number) => number} f
   * @returns {Matrix}
   */
  _apply_colwise_array(
    values: number[] | Float64Array,
    f: (d: number, v: number) => number,
  ): Matrix;
  /**
   * @param {Matrix | number[] | Float64Array | number} value
   * @param {(d: number, v: number) => number} f
   * @returns {Matrix}
   */
  _apply(
    value: Matrix | number[] | Float64Array | number,
    f: (d: number, v: number) => number,
  ): Matrix;
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
  mult(
    value: Matrix | Float64Array | number[] | number,
    {
      inline,
    }?: {
      inline?: boolean | undefined;
    },
  ): Matrix;
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
  divide(
    value: Matrix | Float64Array | number[] | number,
    {
      inline,
    }?: {
      inline?: boolean | undefined;
    },
  ): Matrix;
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
  add(
    value: Matrix | Float64Array | number[] | number,
    {
      inline,
    }?: {
      inline?: boolean | undefined;
    },
  ): Matrix;
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
  sub(
    value: Matrix | Float64Array | number[] | number,
    {
      inline,
    }?: {
      inline?: boolean | undefined;
    },
  ): Matrix;
  /**
   * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
   *
   * @param {[number, number, Accessor]} parameter - Takes an Array in the form [rows, cols, value], where rows and
   *   cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row
   *   and col) which has to return a value for the colth entry of the rowth row.
   * @returns {Matrix}
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
 * Supports different linkage criteria: single, complete, and average.
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
  constructor(
    id: number,
    left: Cluster | null,
    right: Cluster | null,
    dist: number,
    centroid: Float64Array | null,
    index: number,
    size?: number,
    depth?: number,
  );
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
  /** FastPAM1: One best swap per iteration */
  _iteration(): boolean;
  /**
   *
   * @param {number} i
   * @param {number} j
   * @param {Float64Array?} x_i
   * @param {Float64Array?} x_j
   * @returns
   */
  _get_distance(i: number, j: number, x_i?: Float64Array | null, x_j?: Float64Array | null): number;
  /**
   *
   * @param {Float64Array} x_j
   * @param {number} j
   * @returns
   */
  _nearest_medoid(
    x_j: Float64Array,
    j: number,
  ): {
    distance_nearest: number;
    index_nearest: number;
    distance_second: number;
    index_second: number;
  };
  _update_clusters(): void;
  /**
   * Computes `K` clusters out of the `matrix`.
   * @param {number} K - Number of clusters.
   * @param {number[]} cluster_medoids
   */
  init(K: number, cluster_medoids: number[]): this;
  /**
   * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
   *
   * @param {number} K - Number of clusters
   * @returns {number[]}
   */
  _get_random_medoids(K: number): number[];
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
  /** @type {number} */
  _bandwidth: number;
  /** @type {number} */
  _max_iter: number;
  /** @type {number} */
  _tolerance: number;
  /** @type {(dist: number) => number} */
  _kernel: (dist: number) => number;
  /** @type {Matrix} */
  _points: Matrix;
  /** @type {number[] | undefined} */
  _clusters: number[] | undefined;
  /** @type {number[][] | undefined} */
  _cluster_list: number[][] | undefined;
  /**
   * @param {Matrix} matrix
   * @returns {number}
   */
  _compute_bandwidth(matrix: Matrix): number;
  /**
   * @param {number} dist
   * @returns {number}
   */
  _kernel_weight(dist: number): number;
  _mean_shift(): void;
  _assign_clusters(): void;
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
  linkage: "single" | "complete" | "average";
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
  static heapify<T_1>(
    elements: T_1[],
    accessor: (d: T_1) => number,
    comparator?: "min" | "max" | Comparator,
  ): Heap<T_1>;
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
  constructor(
    elements: (T[] | null) | undefined,
    accessor: (d: T) => number,
    comparator?: "min" | "max" | Comparator,
  );
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
   * @private
   * @param {Number} [start_index=0] Default is `0`
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
 * @class
 */
declare class DR<
  T extends InputType,
  Para extends {
    seed?: number;
  },
> {
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
  static transform<
    T_1 extends InputType,
    Para_1 extends {
      seed?: number;
    },
  >(X: T_1, parameters: Para_1, ...args: unknown[]): T_1;
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
  static generator<
    Para_1 extends {
      seed?: number;
    },
  >(X: InputType, parameters: Para_1, ...args: unknown[]): Generator<InputType, InputType, void>;
  /**
   * Computes the projection.
   *
   * @template {{ seed?: number }} Para
   * @param {InputType} X
   * @param {Para} parameters
   * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
   * @returns {Promise<X>} A promise yielding the dimensionality reduced dataset.
   */
  static transform_async<
    Para_1 extends {
      seed?: number;
    },
  >(X: InputType, parameters: Para_1, ...args: unknown[]): Promise<InputType>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersFASTMAP>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersFASTMAP>} parameters
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersFASTMAP>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersISOMAP>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersISOMAP>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersISOMAP>,
  ): Promise<T_1>;
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
  static transform<
    T_1 extends InputType,
    Para extends {
      seed?: number;
    },
  >(X: T_1, parameters: Para): T_1;
  /**
   * @template {InputType} T
   * @template {{ seed?: number }} Para
   * @param {T} X
   * @param {Para} parameters
   * @returns {Generator<T, T, void>}
   */
  static generator<
    T_1 extends InputType,
    Para extends {
      seed?: number;
    },
  >(X: T_1, parameters: Para): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @template {{ seed?: number }} Para
   * @param {T} X
   * @param {Para} parameters
   * @returns {Promise<T>}
   */
  static transform_async<
    T_1 extends InputType,
    Para extends {
      seed?: number;
    },
  >(X: T_1, parameters: Para): Promise<T_1>;
  /**
   * Linear Discriminant Analysis.
   *
   * @param {T} X - The high-dimensional data.
   * @param {Partial<ParametersLDA> & { labels: any[] | Float64Array }} parameters - Object containing parameterization of the DR method.
   * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
   */
  constructor(
    X: T,
    parameters: Partial<ParametersLDA> & {
      labels: any[] | Float64Array;
    },
  );
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersLLE>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersLLE>} parameters
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersLLE>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersLSP>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersLSP>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersLSP>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersLTSA>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersLTSA>} parameters
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersLTSA>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersMDS>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersMDS>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersMDS>,
  ): Promise<T_1>;
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
  static principal_components<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersPCA>,
  ): Matrix;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersPCA>} [parameters]
   * @returns {Generator<T, T, void>}
   */
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersPCA>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersPCA>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersPCA>,
  ): Promise<T_1>;
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
  static transform<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSAMMON<AvailableInit>>,
  ): T_1;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
   * @returns {Generator<T, T, void>}
   */
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSAMMON<AvailableInit>>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSAMMON<AvailableInit>>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSMACOF>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersSMACOF>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSMACOF>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSQDMDS>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersSQDMDS>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersSQDMDS>,
  ): Promise<T_1>;
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
  _sub_div:
    | ((
        x: Float64Array<ArrayBufferLike>,
        y: Float64Array<ArrayBufferLike>,
        div: number,
      ) => Float64Array)
    | undefined;
  _minus:
    | ((a: Float64Array<ArrayBufferLike>, b: Float64Array<ArrayBufferLike>) => Float64Array)
    | undefined;
  _mult: ((a: Float64Array<ArrayBufferLike>, v: number) => Float64Array) | undefined;
  _LR_init: number | undefined;
  _LR: number | undefined;
  _offset: number | undefined;
  _momentums: Matrix | undefined;
  _grads: Matrix | undefined;
  _indices: number[] | undefined;
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
   * @returns {Matrix} The gradients.
   */
  _fill_MDS_grads(Y: Matrix, grads: Matrix, exaggeration?: boolean, zero_grad?: boolean): Matrix;
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
  __minus(
    d: number,
  ): (a: Float64Array<ArrayBufferLike>, b: Float64Array<ArrayBufferLike>) => Float64Array;
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
  __sub_div(
    d: number,
  ): (
    x: Float64Array<ArrayBufferLike>,
    y: Float64Array<ArrayBufferLike>,
    div: number,
  ) => Float64Array;
}

/** @import {InputType} from "../index.js" */
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersTopoMap>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersTopoMap>} parameters
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters: Partial<ParametersTopoMap>,
  ): Promise<T_1>;
  /**
   * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
   *
   * @param {T} X - The high-dimensional data.
   * @param {Partial<ParametersTopoMap>} parameters - Object containing parameterization of the DR method.
   * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
   */
  constructor(X: T, parameters: Partial<ParametersTopoMap>);
  _distance_matrix: Matrix;
  /**
   * @private
   * @param {number} i
   * @param {number} j
   * @param {import("../metrics/index.js").Metric} metric
   * @returns {number}
   */
  private __lazy_distance_matrix;
  /**
   * Computes the minimum spanning tree, using a given metric
   *
   * @private
   * @param {import("../metrics/index.js").Metric} metric
   * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
   */
  private _make_minimum_spanning_tree;
  _disjoint_set: DisjointSet<Float64Array<ArrayBufferLike>> | undefined;
  /** Initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree. */
  init(): this;
  _Emst: number[][] | undefined;
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
   * @abstract
   * @param {T} t
   * @param {number} k
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search(
    t: T,
    k: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @abstract
   * @param {number} i
   * @param {number} k
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search_by_index(
    i: number,
    k: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersTriMap} from "./index.js" */
/** @import {KNN} from "../knn/KNN.js" */
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersTriMap>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersTriMap>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersTriMap>,
  ): Promise<T_1>;
  /**
   * @param {T} X - The high-dimensional data.
   * @param {Partial<ParametersTriMap>} [parameters] - Object containing parameterization of the DR method.
   * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
   * @see {@link https://github.com/eamid/trimap}
   */
  constructor(X: T, parameters?: Partial<ParametersTriMap>);
  /**
   * @param {Matrix | null} [pca=null] - Initial Embedding (if null then PCA gets used). Default is `null`
   * @param {import("../knn/KNN.js").KNN<number[] | Float64Array, any> | null} [knn=null] - KNN Object (if null then BallTree gets used). Default is `null`
   */
  init(pca?: Matrix | null, knn?: KNN<number[] | Float64Array, any> | null): this;
  n_inliers: number | undefined;
  n_outliers: number | undefined;
  n_random: number | undefined;
  knn: KNN<number[] | Float64Array<ArrayBufferLike>, any> | undefined;
  triplets: Matrix | undefined;
  weights: Float64Array<ArrayBuffer> | undefined;
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
  _generate_triplets(
    n_inliers: number,
    n_outliers: number,
    n_random: number,
  ): {
    triplets: Matrix;
    weights: Float64Array<ArrayBuffer>;
  };
  /**
   * Calculates the similarity matrix P
   *
   * @private
   * @param {Matrix} knn_distances - Matrix of pairwise knn distances
   * @param {Float64Array} sig - Scaling factor for the distances
   * @param {Matrix} nbrs - Nearest neighbors
   * @returns {Matrix} Pairwise similarity matrix
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
   */
  _grad(Y: Matrix): {
    grad: Matrix;
    loss: number;
    n_viol: number;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersTSNE>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersTSNE>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersTSNE>,
  ): Promise<T_1>;
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
  static generator<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersUMAP>,
  ): Generator<T_1, T_1, void>;
  /**
   * @template {InputType} T
   * @param {T} X
   * @param {Partial<ParametersUMAP>} [parameters]
   * @returns {Promise<T>}
   */
  static transform_async<T_1 extends InputType>(
    X: T_1,
    parameters?: Partial<ParametersUMAP>,
  ): Promise<T_1>;
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
   * @param {NaiveKNN<Float64Array> | BallTree<Float64Array>} knn
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
declare function simultaneous_poweriteration(
  A: Matrix,
  k?: number,
  { seed, max_iterations, qr, tol }?: EigenArgs,
): {
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
   * - scaling factor.
   */
  weight_adj?: number | undefined;
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
  tol?: number | undefined;
  /**
   * - the metric which defines the distance between two points.
   */
  metric?: Metric | undefined;
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
   * @param {ParametersAnnoy} [parameters={}] - Configuration parameters
   */
  constructor(elements: T[], parameters?: ParametersAnnoy);
  _metric: Metric;
  _numTrees: number;
  _maxPointsPerLeaf: number;
  _seed: number;
  _randomizer: Randomizer;
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
  search(
    query: T,
    k?: number,
  ): {
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
   * @param {number} i
   * @param {number} [k=5]
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * Alias for search_by_index for backward compatibility.
   *
   * @param {number} i - Index of the query element
   * @param {number} [k=5] - Number of nearest neighbors to return
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search_index(
    i: number,
    k?: number,
  ): {
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
   * @param {ParametersBallTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
   * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
   * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
   */
  constructor(elements: T[], parameters?: ParametersBallTree);
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
   * @param {number} i
   * @param {number} k
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @param {T} t - Query element.
   * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
   * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
   */
  search(
    t: T,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @private
   * @param {T} t - Query element.
   * @param {number} k - Number of nearest neighbors to return.
   * @param {Heap<ElementWithIndex<T>>} Q - Heap consists of the currently found `k` nearest neighbors.
   * @param {BallTreeNode<T> | BallTreeLeaf<T>} B
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
   * @param {ParametersHNSW} [parameters={}] - Configuration parameters
   */
  constructor(points: T[], parameters?: ParametersHNSW);
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
  /** @type {Randomizer} */
  _randomizer: Randomizer;
  /** @type {number} - Current maximum layer in the graph */
  _L: number;
  /** @type {number[] | null} - Entry point indices for search */
  _ep: number[] | null;
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
   * @param {boolean} [extend_candidates=true] - Whether to extend candidates with their neighbors
   * @param {boolean} [keep_pruned_connections=true] - Whether to add pruned connections back if needed
   * @returns {number[]} Selected neighbor indices
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
  search_iter(
    q: T,
    K: number,
    ef?: number | null,
  ): Generator<
    {
      layer: number;
      candidates: {
        element: T;
        index: number;
        distance: number;
      }[];
    },
    void,
    unknown
  >;
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
  /**
   * Search for nearest neighbors using an element index as the query.
   *
   * @param {number} i - Index of the query element
   * @param {number} [K=5] - Number of nearest neighbors to return
   * @returns {Candidate<T>[]} K nearest neighbors
   */
  search_by_index(i: number, K?: number): Candidate<T>[];
}
type Candidate<T extends number[] | Float64Array> = {
  /**
   * - The actual data point
   */
  element: T;
  /**
   * - Global index in the dataset
   */
  index: number;
  /**
   * - Distance from query
   */
  distance: number;
};

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
   * @param {ParametersKDTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
   */
  constructor(elements: T[], parameters?: ParametersKDTree);
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
   * @param {number} i
   * @param {number} k
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @param {T} t - Query element.
   * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
   * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
   */
  search(
    t: T,
    k?: number,
  ): {
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
declare class LSH<T extends number[] | Float64Array> extends KNN<T, ParametersLSH> {
  /**
   * Creates a new LSH index.
   *
   * @param {T[]} elements - Elements to index
   * @param {ParametersLSH} [parameters={}] - Configuration parameters
   */
  constructor(elements: T[], parameters?: ParametersLSH);
  _metric: Metric;
  _numHashTables: number;
  _numHashFunctions: number;
  _seed: number;
  _randomizer: Randomizer;
  /** @type {Map<string, number[]>[]} */
  _hashTables: Map<string, number[]>[];
  /** @type {Float64Array[][]} */
  _projections: Float64Array[][];
  /** @type {number[][]} */
  _offsets: number[][];
  /** @type {number} */
  _dim: number;
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
  search(
    query: T,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @param {number} i
   * @param {number} [k=5]
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
}

/** @import { ParametersNaiveKNN } from "./index.js" */
/**
 * Naive KNN implementation using a distance matrix.
 *
 * This implementation pre-computes the entire distance matrix and performs
 * an exhaustive search. Best suited for small datasets or when a distance
 * matrix is already available.
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
   * @param {ParametersNaiveKNN} parameters
   */
  constructor(elements: T[], parameters?: ParametersNaiveKNN);
  _D: Matrix;
  /** @type {Heap<{ value: number; index: number }>[]} */
  KNN: Heap<{
    value: number;
    index: number;
  }>[];
  /**
   * @param {number} i
   * @param {number} k
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @param {T} t - Query element.
   * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
   * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
   */
  search(
    t: T,
    k?: number,
  ): {
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
   * @param {Partial<ParametersNNDescent>} parameters
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
  _randomizer: Randomizer;
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
  search(
    x: T,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
  /**
   * @param {number} i
   * @param {number} [k=5] Default is `5`
   * @returns {{ element: T; index: number; distance: number }[]}
   */
  search_by_index(
    i: number,
    k?: number,
  ): {
    element: T;
    index: number;
    distance: number;
  }[];
}

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
  metric: Metric;
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
   * - Max number of connections per element (excluding ground layer). Default is `16`
   */
  m: number;
  /**
   * - Size of candidate list during construction. Default is `200`
   */
  ef_construction: number;
  /**
   * - Max number of connections for ground layer (layer 0). Default is `2 * m`
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
   * = .8 - Sample rate. Default is `.8`
   */
  rho: number;
  /**
   * = 0.0001 - Precision parameter. Default is `0.0001`
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
declare function powell<T extends Float64Array | number[]>(
  f: (d: T) => number,
  x0: T,
  max_iter?: number,
): T;

type InputType = Matrix | Float64Array[] | number[][];
declare const version: string;

export {
  Annoy,
  BallTree,
  CURE,
  DisjointSet,
  FASTMAP,
  HNSW,
  Heap,
  HierarchicalClustering,
  ISOMAP,
  KDTree,
  KMeans,
  KMedoids,
  LDA,
  LLE,
  LSH,
  LSP,
  LTSA,
  MDS,
  Matrix,
  MeanShift,
  NNDescent,
  NaiveKNN,
  OPTICS,
  PCA,
  Randomizer,
  SAMMON,
  SMACOF,
  SQDMDS,
  TSNE,
  TopoMap,
  TriMap,
  UMAP,
  XMeans,
  bray_curtis,
  canberra,
  chebyshev,
  cosine,
  distance_matrix,
  euclidean,
  euclidean_squared,
  goodman_kruskal,
  hamming,
  haversine,
  inner_product,
  jaccard,
  k_nearest_neighbors,
  kahan_sum,
  linspace,
  manhattan,
  max,
  min,
  neumair_sum,
  norm,
  normalize,
  powell,
  qr,
  qr_householder,
  simultaneous_poweriteration,
  sokal_michener,
  version,
  wasserstein,
  yule,
};
export type {
  Comparator,
  EigenArgs,
  InputType,
  Metric,
  ParametersAnnoy,
  ParametersBallTree,
  ParametersCURE,
  ParametersFASTMAP,
  ParametersHNSW,
  ParametersHierarchicalClustering,
  ParametersISOMAP,
  ParametersKDTree,
  ParametersKMeans,
  ParametersKMedoids,
  ParametersLDA,
  ParametersLLE,
  ParametersLSH,
  ParametersLSP,
  ParametersLTSA,
  ParametersMDS,
  ParametersMeanShift,
  ParametersNNDescent,
  ParametersNaiveKNN,
  ParametersOptics,
  ParametersPCA,
  ParametersSAMMON,
  ParametersSMACOF,
  ParametersSQDMDS,
  ParametersTSNE,
  ParametersTopoMap,
  ParametersTriMap,
  ParametersUMAP,
  ParametersXMeans,
  QRDecomposition,
};
//# sourceMappingURL=druid.d.ts.map
