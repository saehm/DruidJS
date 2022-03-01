/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
export class Matrix {
    /**
     * Creates a Matrix out of {@link A}.
     * @param {(Matrix|Number|Number[]|Number[][])} A - The matrix, array, or number, which should converted to a Matrix.
     * @param {"row"|"col"|"diag"} [type = "row"] - If {@link A} is a Array or Float64Array, then type defines if it is a row- or a column vector.
     * @returns {Matrix}
     *
     * @example
     * let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
     * let S = Matrix.from([1, 2, 3], "diag"); // creates a 3 by 3 matrix with 1, 2, 3 on its diagonal. [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     */
    static from(A: (Matrix | number | number[] | number[][]), type?: "row" | "col" | "diag"): Matrix;
    /**
     * Solves the equation {@link A}x = {@link b} using the conjugate gradient method. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer} [randomizer=null]
     * @param {Number} [tol=1e-3]
     * @returns {Matrix}
     */
    static solve_CG(A: Matrix, b: Matrix, randomizer?: Randomizer, tol?: number): Matrix;
    /**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A: Matrix, b: Matrix): Matrix;
    /**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */
    static LU(A: Matrix): {
        L: Matrix;
        U: Matrix;
    };
    /**
     * Computes the determinante of {@link A}, by using the LU decomposition of {@link A}.
     * @param {Matrix} A
     * @returns {Number} det - Returns the determinate of the Matrix {@link A}.
     */
    static det(A: Matrix): number;
    /**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M
     * @param {Number} [k=2]
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */
    static SVD(M: Matrix, k?: number): {
        U: Matrix;
        Sigma: Matrix;
        V: Matrix;
    };
    /**
     * creates a new Matrix. Entries are stored in a Float64Array.
     * @constructor
     * @memberof module:matrix
     * @alias Matrix
     * @param {Number} rows - The amount of rows of the matrix.
     * @param {Number} cols - The amount of columns of the matrix.
     * @param {(Function|"zero"|"identity"|"center"|"I"|Number)} [value = 0] - Can be a function with row and col as parameters, a number, or "zeros", "identity" or "I", or "center".
     *  - **Function**: for each entry the function gets called with the parameters for the actual row and column.
     *  - **String**: allowed are
     *      - "zero", creates a zero matrix.
     *      - "identity" or "I", creates an identity matrix.
     *      - "center", creates an center matrix.
     *  - **Number**: create a matrix filled with the given value.
     * @example
     *
     * let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix.
     * let B = new Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
     * @returns {Matrix} returns a {@link rows} times {@link cols} Matrix filled with {@link value}.
     */
    constructor(rows?: number, cols?: number, value?: (Function | "zero" | "identity" | "center" | "I" | number));
    _rows: number;
    _cols: number;
    _data: Float64Array;
    /**
     * Returns the {@link row}<sup>th</sup> row from the Matrix.
     * @param {Number} row
     * @returns {Float64Array}
     */
    row(row: number): Float64Array;
    /**
     * Returns an generator yielding each row of the Matrix.
     * @yields {Float64Array}
     */
    iterate_rows(): Generator<Float64Array, void, unknown>;
    /**
     * Sets the entries of {@link row}<sup>th</sup> row from the Matrix to the entries from {@link values}.
     * @param {Number} row
     * @param {Number[]|Matrix} values
     * @returns {Matrix}
     */
    set_row(row: number, values: number[] | Matrix): Matrix;
    /**
     * Returns the {@link col}<sup>th</sup> column from the Matrix.
     * @param {Number} col
     * @returns {Number[]}
     */
    col(col: number): number[];
    /**
     * Returns the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {Number} row
     * @param {Number} col
     * @returns {Number}
     */
    entry(row: number, col: number): number;
    /**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given {@link value}.
     * @param {Number} row
     * @param {Number} col
     * @param {Number} value
     * @returns {Matrix}
     */
    set_entry(row: number, col: number, value: number): Matrix;
    /**
     * Returns a new transposed Matrix.
     * @returns {Matrix}
     */
    transpose(): Matrix;
    /**
     * Returns a new transposed Matrix. Short-form of {@function transpose}.
     * @returns {Matrix}
     */
    get T(): Matrix;
    /**
     * Returns the inverse of the Matrix.
     * @returns {Matrix}
     */
    inverse(): Matrix;
    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {Matrix|Number[]} B the right side
     * @returns {Matrix|Number[]}
     */
    dot(B: Matrix | number[]): Matrix | number[];
    /**
     * Computes the outer product from {@link this} and {@link B}.
     * @param {Matrix} B
     * @returns {Matrix}
     */
    outer(B: Matrix): Matrix;
    /**
     * Appends matrix {@link B} to the matrix.
     * @param {Matrix} B - matrix to append.
     * @param {"horizontal"|"vertical"|"diag"} [type = "horizontal"] - type of concatenation.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 1], [1, 1]]); // 2 by 2 matrix filled with ones.
     * let B = Matrix.from([[2, 2], [2, 2]]); // 2 by 2 matrix filled with twos.
     *
     * A.concat(B, "horizontal"); // 2 by 4 matrix. [[1, 1, 2, 2], [1, 1, 2, 2]]
     * A.concat(B, "vertical"); // 4 by 2 matrix. [[1, 1], [1, 1], [2, 2], [2, 2]]
     * A.concat(B, "diag"); // 4 by 4 matrix. [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]]
     */
    concat(B: Matrix, type?: "horizontal" | "vertical" | "diag"): Matrix;
    /**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {Number} offset_row
     * @param {Number} offset_col
     * @param {Matrix} B
     * @returns {Matrix}
     */
    set_block(offset_row: number, offset_col: number, B: Matrix): Matrix;
    /**
     * Extracts the entries from the {@link start_row}<sup>th</sup> row to the {@link end_row}<sup>th</sup> row, the {@link start_col}<sup>th</sup> column to the {@link end_col}<sup>th</sup> column of the matrix.
     * If {@link end_row} or {@link end_col} is empty, the respective value is set to {@link this.rows} or {@link this.cols}.
     * @param {Number} start_row
     * @param {Number} start_col
     * @param {Number} [end_row = null]
     * @param {Number} [end_col = null]
     * @returns {Matrix} Returns a end_row - start_row times end_col - start_col matrix, with respective entries from the matrix.
     * @example
     *
     * let A = Matrix.from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]); // a 3 by 3 matrix.
     *
     * A.get_block(1, 1); // [[5, 6], [8, 9]]
     * A.get_block(0, 0, 1, 1); // [[1]]
     * A.get_block(1, 1, 2, 2); // [[5]]
     * A.get_block(0, 0, 2, 2); // [[1, 2], [4, 5]]
     */
    get_block(start_row: number, start_col: number, end_row?: number, end_col?: number): Matrix;
    /**
     * Returns a new array gathering entries defined by the indices given by argument.
     * @param {Number[]} row_indices - Array consists of indices of rows for gathering entries of this matrix
     * @param {Number[]} col_indices  - Array consists of indices of cols for gathering entries of this matrix
     * @returns {Matrix}
     */
    gather(row_indices: number[], col_indices: number[]): Matrix;
    /**
     * Applies a function to each entry of the matrix.
     * @private
     * @param {Function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {Function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    private _apply_array;
    _apply_rowwise_array(values: any, f: any): Matrix;
    _apply_colwise_array(values: any, f: any): Matrix;
    _apply(value: any, f: any): Matrix;
    /**
     * Clones the Matrix.
     * @returns {Matrix}
     */
    clone(): Matrix;
    /**
     * Entrywise multiplication with {@link value}.
     * @param {Matrix|Array|Number} value
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.mult(2); // [[2, 4], [6, 8]];
     * A.mult(B); // [[1, 4], [9, 16]];
     */
    mult(value: Matrix | any[] | number): Matrix;
    /**
     * Entrywise division with {@link value}.
     * @param {Matrix|Array|Number} value
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.divide(2); // [[0.5, 1], [1.5, 2]];
     * A.divide(B); // [[1, 1], [1, 1]];
     */
    divide(value: Matrix | any[] | number): Matrix;
    /**
     * Entrywise addition with {@link value}.
     * @param {Matrix|Array|Number} value
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.add(2); // [[3, 4], [5, 6]];
     * A.add(B); // [[2, 4], [6, 8]];
     */
    add(value: Matrix | any[] | number): Matrix;
    /**
     * Entrywise subtraction with {@link value}.
     * @param {Matrix|Array|Number} value
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.sub(2); // [[-1, 0], [1, 2]];
     * A.sub(B); // [[0, 0], [0, 0]];
     */
    sub(value: Matrix | any[] | number): Matrix;
    /**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
     * @returns {Matrix}
     */
    set shape(arg: any[]);
    /**
     * Returns the number of rows and columns of the Matrix.
     * @returns {Array} An Array in the form [rows, columns].
     */
    get shape(): any[];
    /**
     * Returns the Matrix as a Array of Float64Arrays.
     * @returns {Number[][]}
     */
    get to2dArray(): number[][];
    /**
     * Returns the Matrix as a Array of Arrays.
     * @returns {Number[][]}
     */
    get asArray(): number[][];
    /**
     * Returns the diagonal of the Matrix.
     * @returns {Number[]}
     */
    get diag(): number[];
    /**
     * Returns the mean of all entries of the Matrix.
     * @returns {Number}
     */
    get mean(): number;
    /**
     * Returns the sum oof all entries of the Matrix.
     * @returns {Number}
     */
    get sum(): number;
    /**
     * Returns the sum oof all entries of the Matrix.
     * @returns {Float64Array}
     */
    get values(): Float64Array;
    /**
     * Returns the mean of each row of the matrix.
     * @returns {Float64Array}
     */
    get meanRows(): Float64Array;
    /** Returns the mean of each column of the matrix.
     * @returns {Float64Array}
     */
    get meanCols(): Float64Array;
    /**
     * Makes a {@link Matrix} object an iterable object.
     * @yields {Float64Array}
     */
    [Symbol.iterator](): Generator<any, void, unknown>;
}
import { Randomizer } from "../util/index.js";
//# sourceMappingURL=Matrix.d.ts.map