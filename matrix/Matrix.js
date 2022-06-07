import { neumair_sum } from "../numerical/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Randomizer } from "../util/index.js";
/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
export class Matrix {
    /**
     * creates a new Matrix. Entries are stored in a Float64Array.
     * @memberof module:matrix
     * @param {number} rows - The amount of rows of the matrix.
     * @param {number} cols - The amount of columns of the matrix.
     * @param {(function|string|number)} value=0 - Can be a function with row and col as parameters, a number, or "zeros", "identity" or "I", or "center".
     *  - **function**: for each entry the function gets called with the parameters for the actual row and column.
     *  - **string**: allowed are
     *      - "zero", creates a zero matrix.
     *      - "identity" or "I", creates an identity matrix.
     *      - "center", creates an center matrix.
     *  - **number**: create a matrix filled with the given value.
     * @example
     *
     * let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix.
     * let B = new Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
     * @returns {Matrix} returns a {@link rows} times {@link cols} Matrix filled with {@link value}.
     */
    constructor(rows = null, cols = null, value = null) {
        this._rows = rows;
        this._cols = cols;
        this._data = null;
        if (rows && cols) {
            if (!value) {
                this._data = new Float64Array(rows * cols);
                return this;
            }
            if (typeof value === "function") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value(row, col);
                    }
                }
                return this;
            }
            if (typeof value === "string") {
                if (value === "zeros") {
                    return new Matrix(rows, cols, 0);
                }
                if (value === "identity" || value === "I") {
                    this._data = new Float64Array(rows * cols);
                    for (let row = 0; row < rows; ++row) {
                        this._data[row * cols + row] = 1;
                    }
                    return this;
                }
                if (value === "center" && rows == cols) {
                    this._data = new Float64Array(rows * cols);
                    value = (i, j) => (i === j ? 1 : 0) - 1 / rows;
                    for (let row = 0; row < rows; ++row) {
                        for (let col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = value(row, col);
                        }
                    }
                    return this;
                }
            }
            if (typeof value === "number") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value;
                    }
                }
                return this;
            }
        }
        return this;
    }

    /**
     * Creates a Matrix out of {@link A}.
     * @param {(Matrix|Array|Float64Array|number)} A - The matrix, array, or number, which should converted to a Matrix.
     * @param {"row"|"col"|"diag"} [type = "row"] - If {@link A} is a Array or Float64Array, then type defines if it is a row- or a column vector.
     * @returns {Matrix}
     *
     * @example
     * let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
     * let S = Matrix.from([1, 2, 3], "diag"); // creates a 3 by 3 matrix with 1, 2, 3 on its diagonal. [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     */
    static from(A, type = "row") {
        if (A instanceof Matrix) {
            return A.clone();
        } else if (Matrix.isArray(A)) {
            let m = A.length;
            if (m === 0) throw new Error("Array is empty");
            // 1d
            if (!Matrix.isArray(A[0])) {
                if (type === "row") {
                    return new Matrix(1, m, (_, j) => A[j]);
                } else if (type === "col") {
                    return new Matrix(m, 1, (i) => A[i]);
                } else if (type === "diag") {
                    return new Matrix(m, m, (i, j) => (i == j ? A[i] : 0));
                } else {
                    throw new Error("1d array has NaN entries");
                }
                // 2d
            } else {
                let n = A[0].length;
                for (let row = 0; row < m; ++row) {
                    if (A[row].length !== n) {
                        throw new Error("various array lengths");
                    }
                }
                return new Matrix(m, n, (i, j) => A[i][j]);
            }
        } else if (typeof A === "number") {
            return new Matrix(1, 1, A);
        } else {
            throw new Error("error");
        }
    }

    /**
     * Returns the {@link row}<sup>th</sup> row from the Matrix.
     * @param {Number} row
     * @returns {Float64Array}
     */
    row(row) {
        const data = this.values;
        const cols = this._cols;
        return data.subarray(row * cols, (row + 1) * cols);
    }

    /**
     * Returns an generator yielding each row of the Matrix.
     * @yields {Float64Array}
     */
    *iterate_rows() {
        const cols = this._cols;
        const rows = this._rows;
        const data = this.values;
        for (let row = 0; row < rows; ++row) {
            yield data.subarray(row * cols, (row + 1) * cols);
        }
    }

    /**
     * Makes a {@link Matrix} object an iterable object.
     * @yields {Float64Array}
     */
    *[Symbol.iterator]() {
        for (const row of this.iterate_rows()) {
            yield row;
        }
    }

    /**
     * Sets the entries of {@link row}<sup>th</sup> row from the Matrix to the entries from {@link values}.
     * @param {Number} row
     * @param {Array} values
     * @returns {Matrix}
     */
    set_row(row, values) {
        const cols = this._cols;
        if (Matrix.isArray(values) && values.length === cols) {
            const offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this.values[offset + col] = values[col];
            }
        } else if (values instanceof Matrix && values.shape[1] === cols && values.shape[0] === 1) {
            const offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this.values[offset + col] = values._data[col];
            }
        } else {
            throw new Error("Values not valid! Needs to be either an Array, a Float64Array, or a fitting Matrix!")
        }
        return this;
    }

    /**
     * Swaps the rows {@link row1} and {@link row2} of the Matrix.
     * @param {Number} row1
     * @param {Number} row2
     * @returns {Matrix}
     */
    swap_rows(row1, row2) {
        const cols = this._cols;
        const data = this.values;
        for (let i = row1 * cols, j = row2 * cols, col = 0; col < cols; ++col, ++i, ++j) {
            const t = data[i];
            data[i] = data[j];
            data[j] = t;
        }
    }

    /**
     * Returns the {@link col}<sup>th</sup> column from the Matrix.
     * @param {Number} col
     * @returns {Array}
     */
    col(col) {
        const result_col = new Float64Array(this._rows);
        for (let row = 0; row < this._rows; ++row) {
            result_col[row] = this.values[row * this._cols + col];
        }
        return result_col;
    }

    /**
     * Returns the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @returns {float64}
     */
    entry(row, col) {
        return this.values[row * this._cols + col];
    }

    /**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given {@link value}.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */
    set_entry(row, col, value) {
        this.values[row * this._cols + col] = value;
        return this;
    }

    /**
     * Adds a given {@link value} to the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */
    add_entry(row, col, value) {
      this.values[row * this._cols + col] += value;
      return this;
    }

    /**
     * Subtracts a given {@link value} from the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */
    sub_entry(row, col, value) {
      this.values[row * this._cols + col] -= value;
      return this;
    }

    /**
     * Returns a new transposed Matrix.
     * @returns {Matrix}
     */
    transpose() {
        let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
        return B;
    }

    /**
     * Returns a new transposed Matrix. Short-form of {@function transpose}.
     * @returns {Matrix}
     */
    get T() {
        return this.transpose();
    }

    /**
     * Returns the inverse of the Matrix.
     * @returns {Matrix}
     */
    inverse() {
        const rows = this._rows;
        const cols = this._cols;
        const A = this.clone();
        const B = new Matrix(rows, cols, 'I');

        // foreach column
        for (let col = 0; col < cols; ++col) {
            // Search for maximum in this column (pivot)
            let max_idx = col;
            let max_val = Math.abs(A.entry(col, col));
            for (let row = col + 1; row < rows; ++row) {
                const val = Math.abs(A.entry(row, col));
                if (max_val < val) {
                    max_idx = row;
                    max_val = val;
                }
            }
            if (max_val === 0) {
                throw new Error('Cannot compute inverse of Matrix, determinant is zero');
            }
            // Swap maximum row with current row
            if (max_idx !== col) {
                A.swap_rows(col, max_idx);
                B.swap_rows(col, max_idx);
            }

            // eliminate non-zero values on the other rows at column c
            const A_col = A.row(col);
            const B_col = B.row(col);
            for (let row = 0; row < rows; ++row) {
                if (row !== col) {
                    // eliminate value at column c and row r
                    const A_row = A.row(row);
                    const B_row = B.row(row);
                    if (A_row[col] !== 0) {
                        const f = A_row[col] / A_col[col];
                        // sub (f * row c) from row r to eliminate the value at column c
                        for (let s = col; s < cols; ++s) {
                            A_row[s] -= (f * A_col[s]);
                        }
                        for (let s = 0; s < cols; ++s) {
                            B_row[s] -= (f * B_col[s]);
                        }
                    }
                } else {
                    // normalize value at Acc to 1 (diagonal):
                    // divide each value of row r=c by the value at Acc
                    const f = A_col[col];
                    for (let s = col; s < cols; ++s) {
                        A_col[s] /= f;
                    }
                    for (let s = 0; s < cols; ++s) {
                        B_col[s] /= f;
                    }
                }
            }
        }
        return B;
    }

    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    dot(B) {
        if (B instanceof Matrix) {
            let A = this;
            const [rows_A, cols_A] = A.shape;
            const [rows_B, cols_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${A.shape.join(" ⨯ ")}-Matrix, B is a ${B.shape.join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows.
                Must be equal!`);
            }
            const C = new Matrix(rows_A, cols_B, (row, col) => {
                const A_i = A.row(row);
                const B_val = B.values;
                let sum = 0;
                for (let i = 0, j = col; i < cols_A; ++i, j += cols_B) {
                    sum += A_i[i] * B_val[j];
                }
                return sum;
            });
            return C;
        } else if (Matrix.isArray(B)) {
            let rows = this._rows;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            let C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map((e) => e * B[row]));
            }
            return C;
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Transposes the current matrix and returns the dot product with {@link B}.
     * If {@link B} is an Array or Float64Array then an Array gets returned.
     * If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    transDot(B) {
        if (B instanceof Matrix) {
            let A = this;
            const [cols_A, rows_A] = A.shape; // transpose matrix
            const [rows_B, cols_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${[rows_A, cols_A].join(" ⨯ ")}-Matrix, B is a ${B.shape.join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows, which must be equal!`);
            }
            // let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
            // this.values[row * this._cols + col];
            const C = new Matrix(rows_A, cols_B, (row, col) => {
                const A_val = A.values;
                const B_val = B.values;
                let sum = 0;
                for (let i = 0, j = row, k = col; i < cols_A; ++i, j += rows_A, k += cols_B) {
                    sum += A_val[j] * B_val[k];
                }
                return sum;
            });
            return C;
        } else if (Matrix.isArray(B)) {
            let rows = this._cols;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            let C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.col(row).map((e) => e * B[row]));
            }
            return C;
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Returns the dot product with the transposed version of {@link B}.
     * If {@link B} is an Array or Float64Array then an Array gets returned.
     * If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    dotTrans(B) {
        if (B instanceof Matrix) {
            let A = this;
            const [rows_A, cols_A] = A.shape;
            const [cols_B, rows_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${A.shape.join(" ⨯ ")}-Matrix, B is a ${[rows_B, cols_B].join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows, which must be equal!`);
            }
            const C = new Matrix(rows_A, cols_B, (row, col) => {
                const A_i = A.row(row);
                const B_i = B.row(col);
                let sum = 0;
                for (let i = 0; i < cols_A; ++i) {
                    sum += A_i[i] * B_i[i];
                }
                return sum;
            });
            return C;
        } else if (Matrix.isArray(B)) {
            let rows = this._rows;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            let C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map((e) => e * B[row]));
            }
            return C;
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Computes the outer product from {@link this} and {@link B}.
     * @param {Matrix} B
     * @returns {Matrix}
     */
    outer(B) {
        let A = this;
        let l = A._data.length;
        let r = B._data.length;
        if (l != r) return undefined;
        let C = new Matrix();
        C.shape = [
            l,
            l,
            (i, j) => {
                if (i <= j) {
                    return A._data[i] * B._data[j];
                } else {
                    return C.entry(j, i);
                }
            },
        ];
        return C;
    }

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
    concat(B, type = "horizontal") {
        const A = this;
        const [rows_A, cols_A] = A.shape;
        const [rows_B, cols_B] = B.shape;
        if (type == "horizontal") {
            if (rows_A != rows_B) {
                throw new Error(`A.concat(B, "horizontal"): A and B need same number of rows, A has ${rows_A} rows, B has ${rows_B} rows.`);
            }
            const X = new Matrix(rows_A, cols_A + cols_B, "zeros");
            X.set_block(0, 0, A);
            X.set_block(0, cols_A, B);
            return X;
        } else if (type == "vertical") {
            if (cols_A != cols_B) {
                throw new Error(`A.concat(B, "vertical"): A and B need same number of columns, A has ${cols_A} columns, B has ${cols_B} columns.`);
            }
            const X = new Matrix(rows_A + rows_B, cols_A, "zeros");
            X.set_block(0, 0, A);
            X.set_block(rows_A, 0, B);
            return X;
        } else if (type == "diag") {
            const X = new Matrix(rows_A + rows_B, cols_A + cols_B, "zeros");
            X.set_block(0, 0, A);
            X.set_block(rows_A, cols_A, B);
            return X;
        } else {
            throw new Error(`type must be "horizontal" or "vertical", but type is ${type}!`);
        }
    }

    /**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row
     * @param {int} offset_col
     * @param {Matrix} B
     * @returns {Matrix}
     */
    set_block(offset_row, offset_col, B) {
        const rows = Math.min(this._rows - offset_row, B.shape[0]);
        const cols = Math.min(this._cols - offset_col, B.shape[1]);
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < cols; ++col) {
                this.set_entry(row + offset_row, col + offset_col, B.entry(row, col));
            }
        }
        return this;
    }

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
    get_block(start_row, start_col, end_row = null, end_col = null) {
        const [rows, cols] = this.shape;
        end_row = end_row ?? rows;
        end_col = end_col ?? cols;
        if (end_row <= start_row || end_col <= start_col) {
            throw new Error(`
                end_row must be greater than start_row, and
                end_col must be greater than start_col, but
                end_row = ${end_row}, start_row = ${start_row}, end_col = ${end_col}, and start_col = ${start_col}!`);
        }
        const X = new Matrix(end_row - start_row, end_col - start_col, "zeros");
        for (let row = start_row, new_row = 0; row < end_row; ++row, ++new_row) {
            for (let col = start_col, new_col = 0; col < end_col; ++col, ++new_col) {
                X.set_entry(new_row, new_col, this.entry(row, col));
            }
        }
        return X;
        //return new Matrix(end_row - start_row, end_col - start_col, (i, j) => this.entry(i + start_row, j + start_col));
    }

    /**
     * Returns a new array gathering entries defined by the indices given by argument.
     * @param {Array<Number>} row_indices - Array consists of indices of rows for gathering entries of this matrix
     * @param {Array<Number>} col_indices  - Array consists of indices of cols for gathering entries of this matrix
     * @returns {Matrix}
     */
    gather(row_indices, col_indices) {
        const N = row_indices.length;
        const D = col_indices.length;

        const R = new Matrix(N, D);
        for (let i = 0; i < N; ++i) {
            const row_index = row_indices[i];
            for (let j = 0; j < N; ++j) {
                const col_index = col_indices[j];
                R.set_entry(i, j, this.entry(row_index, col_index));
            }
        }

        return R;
    }

    /**
     * Applies a function to each entry of the matrix.
     * @private
     * @param {Function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {Function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    _apply_array(f, v) {
        const data = this.values;
        const [rows, cols] = this.shape;
        for (let i = 0, row = 0; row < rows; ++row) {
            for (let col = 0; col < cols; ++col, ++i) {
                data[i] = f(data[i], v(row, col));
            }
        }
        return this;
    }

    _apply_rowwise_array(values, f) {
        return this._apply_array(f, (_, j) => values[j]);
    }

    _apply_colwise_array(values, f) {
        const data = this.values;
        const [rows, cols] = this.shape;
        for (let i = 0, row = 0; row < rows; ++row) {
            const val = values[row];
            for (let col = 0; col < cols; ++col, ++i) {
                data[i] = f(data[i], val);
            }
        }
        return this;
    }

    _apply(value, f) {
        const data = this.values;
        const [rows, cols] = this.shape;
        if (value instanceof Matrix) {
            const values = value.values;
            const [value_rows, value_cols] = value.shape;
            if (value_rows === 1) {
                if (cols !== value_cols) {
                    throw new Error(`cols !== value_cols`);
                }
                for (let i = 0, row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col, ++i) {
                        data[i] = f(data[i], values[col]);
                    }
                }
            } else if (value_cols === 1) {
                if (rows !== value_rows) {
                    throw new Error(`rows !== value_rows`);
                }
                for (let i = 0, row = 0; row < rows; ++row) {
                    const v =  values[row];
                    for (let col = 0; col < cols; ++col, ++i) {
                        data[i] = f(data[i], v);
                    }
                }
            } else if (rows == value_rows && cols == value_cols) {
                for (let i = 0, n = rows * cols; i < n; ++i) {
                    data[i] = f(data[i], values[i]);
                }
            } else {
                throw new Error(`error`);
            }
        } else if (Matrix.isArray(value)) {
            if (value.length === rows) {
                for (let i = 0, row = 0; row < rows; ++row) {
                    const v = value[row];
                    for (let col = 0; col < cols; ++col, ++i) {
                        data[i] = f(data[i], v);
                    }
                }
            } else if (value.length === cols) {
                for (let i = 0, row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col, ++i) {
                        data[i] = f(data[i], value[col]);
                    }
                }
            } else {
                throw new Error(`error`);
            }
        } else { // scalar value
            for (let i = 0, n = rows * cols; i < n; ++i) {
                data[i] = f(data[i], value);
            }
        }
        return this;
    }

    /**
     * Clones the Matrix.
     * @returns {Matrix}
     */
    clone() {
        let B = new Matrix();
        B._rows = this._rows;
        B._cols = this._cols;
        B._data = this.values.slice(0);
        return B;
    }

    /**
     * Entrywise multiplication with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false]  - If true, applies multiplication to the element, otherwise it creates first a copy and applies the multiplication on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.mult(2); // [[2, 4], [6, 8]];
     * A.mult(B); // [[1, 4], [9, 16]];
     */
    mult(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a * b);
    }

    /**
     * Entrywise division with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false] - If true, applies division to the element, otherwise it creates first a copy and applies the division on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.divide(2); // [[0.5, 1], [1.5, 2]];
     * A.divide(B); // [[1, 1], [1, 1]];
     */
    divide(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a / b);
    }

    /**
     * Entrywise addition with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false]  - If true, applies addition to the element, otherwise it creates first a copy and applies the addition on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.add(2); // [[3, 4], [5, 6]];
     * A.add(B); // [[2, 4], [6, 8]];
     */
    add(value, {inline = false} = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a + b);
    }

    /**
     * Entrywise subtraction with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false] - If true, applies subtraction to the element, otherwise it creates first a copy and applies the subtraction on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.sub(2); // [[-1, 0], [1, 2]];
     * A.sub(B); // [[0, 0], [0, 0]];
     */
    sub(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a - b);
    }

    /**
     * Returns the number of rows and columns of the Matrix.
     * @returns {Array} An Array in the form [rows, columns].
     */
    get shape() {
        return [this._rows, this._cols];
    }

    /**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
     * @returns {Matrix}
     */
    set shape([rows, cols, value = () => 0]) {
        this._rows = rows;
        this._cols = cols;
        this._data = new Float64Array(rows * cols);
        for (let i = 0, row = 0; row < rows; ++row) {
            for (let col = 0; col < cols; ++col, ++i) {
                this._data[i] = value(row, col);
            }
        }
        return this;
    }

    /**
     * Returns the Matrix as a Array of Float64Arrays.
     * @returns {Array<Float64Array>}
     */
    get to2dArray() {
        const result = [];
        for (const row of this.iterate_rows()) {
            result.push(row);
        }
        return result;
    }

    /**
     * Returns the Matrix as a Array of Arrays.
     * @returns {Array<Array>}
     */
    get asArray() {
        const result = [];
        for (const row of this.iterate_rows()) {
            result.push(Array.from(row));
        }
        return result;
    }

    /**
     * Returns the diagonal of the Matrix.
     * @returns {Float64Array}
     */
    get diag() {
        const rows = this._rows;
        const cols = this._cols;
        const min_row_col = Math.min(rows, cols);
        let result = new Float64Array(min_row_col);
        for (let i = 0; i < min_row_col; ++i) {
            result[i] = this.entry(i, i);
        }
        return result;
    }

    /**
     * Returns the mean of all entries of the Matrix.
     * @returns {Number}
     */
    get mean() {
        const sum = this.sum;
        const n = this._rows * this._cols;
        return sum / n;
    }

    /**
     * Returns the sum oof all entries of the Matrix.
     * @returns {Number}
     */
    get sum() {
        const data = this.values;
        return neumair_sum(data);
    }

    /**
     * Returns the entries of the Matrix.
     * @returns {Float64Array}
     */
    get values() {
        const data = this._data;
        return data;
    }

    /**
     * Returns the mean of each row of the matrix.
     * @returns {Float64Array}
     */
    get meanRows() {
        const data = this.values;
        const rows = this._rows;
        const cols = this._cols;
        const result = Float64Array.from({ length: rows });
        for (let i = 0, row = 0; row < rows; ++row) {
            let sum = 0;
            for (let col = 0; col < cols; ++col, ++i) {
                sum += data[i];
            }
            result[row] = sum / cols;
        }
        return result;
    }

    /** Returns the mean of each column of the matrix.
     * @returns {Float64Array}
     */
    get meanCols() {
        const data = this.values;
        const rows = this._rows;
        const cols = this._cols;
        const result = Float64Array.from({ length: cols });
        for (let col = 0; col < cols; ++col) {
            let sum = 0;
            for (let i = col, row = 0; row < rows; ++row, i += cols) {
                sum += data[i];
            }
            result[col] = sum / rows;
        }
        return result;
    }

    /**
     * Solves the equation {@link A}x = {@link b} using the conjugate gradient method. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer} [randomizer=null]
     * @param {Number} [tol=1e-3]
     * @returns {Matrix}
     */
    static solve_CG(A, b, randomizer, tol = 1e-3) {
        if (randomizer === null) {
            randomizer = new Randomizer();
        }
        const rows = A.shape[0];
        const cols = b.shape[1];
        let result = new Matrix(rows, 0);
        for (let i = 0; i < cols; ++i) {
            const b_i = Matrix.from(b.col(i)).T;
            let x = new Matrix(rows, 1, () => randomizer.random);
            let r = b_i.sub(A.dot(x));
            let d = r.clone();
            do {
                const z = A.dot(d);
                const alpha = r.transDot(r).entry(0, 0) / d.transDot(z).entry(0, 0);
                x = x.add(d.mult(alpha));
                const r_next = r.sub(z.mult(alpha));
                const beta = r_next.transDot(r_next).entry(0, 0) / r.transDot(r).entry(0, 0);
                d = r_next.add(d.mult(beta));
                r = r_next;
            } while (Math.abs(r.mean) > tol);
            result = result.concat(x, "horizontal");
        }
        return result;
    }

    /**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A, b) {
        let { L: L, U: U } = "L" in A && "U" in A ? A : Matrix.LU(A);
        let rows = L.shape[0];
        let x = b.clone();

        // forward
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < row - 1; ++col) {
                x.sub_entry(0, row, L.entry(row, col) * x.entry(1, col));
            }
            x.set_entry(0, row, x.entry(0, row) / L.entry(row, row));
        }

        // backward
        for (let row = rows - 1; row >= 0; --row) {
            for (let col = rows - 1; col > row; --col) {
                x.sub_entry(0, row, U.entry(row, col) * x.entry(0, col));
            }
            x.set_entry(0, row, x.entry(0, row) / U.entry(row, row));
        }

        return x;
    }

    /**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */
    static LU(A) {
        const rows = A.shape[0];
        const L = new Matrix(rows, rows, "zeros");
        const U = new Matrix(rows, rows, "identity");

        for (let j = 0; j < rows; ++j) {
            for (let i = j; i < rows; ++i) {
                let sum = 0;
                for (let k = 0; k < j; ++k) {
                    sum += L.entry(i, k) * U.entry(k, j);
                }
                L.set_entry(i, j, A.entry(i, j) - sum);
            }
            for (let i = j; i < rows; ++i) {
                if (L.entry(j, j) === 0) {
                    return undefined;
                }
                let sum = 0;
                for (let k = 0; k < j; ++k) {
                    sum += L.entry(j, k) * U.entry(k, i);
                }
                U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
            }
        }

        return { L: L, U: U };
    }

    /**
     * Computes the determinante of {@link A}, by using the LU decomposition of {@link A}.
     * @param {Matrix} A
     * @returns {Number} det - Returns the determinate of the Matrix {@link A}.
     */
    static det(A) {
        const rows = A.shape[0];
        const { L, U } = Matrix.LU(A);
        const L_diag = L.diag;
        const U_diag = U.diag;
        let det = L_diag[0] * U_diag[0];
        for (let row = 1; row < rows; ++row) {
            det *= L_diag[row] * U_diag[row];
        }
        return det;
    }

    /**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M
     * @param {int} [k=2]
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */
    static SVD(M, k = 2) {
        let MtM = M.transDot(M);
        let MMt = M.dotTrans(M);
        let { eigenvectors: V, eigenvalues: Sigma } = simultaneous_poweriteration(MtM, k);
        let { eigenvectors: U } = simultaneous_poweriteration(MMt, k);
        return { U: U, Sigma: Sigma.map((sigma) => Math.sqrt(sigma)), V: V };

        //Algorithm 1a: Householder reduction to bidiagonal form:
        /* const [m, n] = A.shape;
        let U = new Matrix(m, n, (i, j) => i == j ? 1 : 0);
        console.log(U.to2dArray)
        let V = new Matrix(n, m, (i, j) => i == j ? 1 : 0);
        console.log(V.to2dArray)
        let B = Matrix.bidiagonal(A.clone(), U, V);
        console.log(U,V,B)
        return { U: U, "Sigma": B, V: V }; */
    }

    static isArray(A) {
      return Array.isArray(A) || A instanceof Float64Array || A instanceof Float32Array;
    }
}
