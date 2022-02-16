// https://renecutura.eu v0.4.2 Copyright 2022 Rene Cutura
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.druid = global.druid || {}));
})(this, (function (exports) { 'use strict';

/**
 * Computes the euclidean distance (l<sub>2</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the euclidean distance between {@link a} and {@link b}.
 */
function euclidean (a, b) {
    return Math.sqrt(euclidean_squared(a, b));
}

/**
 * Numerical stable summation with the Kahan summation algorithm.
 * @memberof module:numerical
 * @alias kahan_sum
 * @param {Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
 */
function kahan_sum (summands) {
    let n = summands.length;
    let sum = 0;
    let compensation = 0;
    let y, t;

    for (let i = 0; i < n; ++i) {
        y = summands[i] - compensation;
        t = sum + y;
        compensation = t - sum - y;
        sum = t;
    }
    return sum;
}

/**
 * Numerical stable summation with the Neumair summation algorithm.
 * @memberof module:numerical
 * @alias neumair_sum
 * @param {Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
function neumair_sum (summands) {
    let n = summands.length;
    let sum = 0;
    let compensation = 0;

    for (let i = 0; i < n; ++i) {
        let summand = summands[i];
        let t = sum + summand;
        if (Math.abs(sum) >= Math.abs(summand)) {
            compensation += sum - t + summand;
        } else {
            compensation += summand - t + sum;
        }
        sum = t;
    }
    return sum + compensation;
}

/**
 * Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the squared euclidean distance between {@link a} and {@link b}.
 */
function euclidean_squared (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        let x = a[i];
        let y = b[i];
        s[i] = (x - y) * (x - y);
    }
    return neumair_sum(s);
}

/**
 * Computes the cosine distance (not similarity) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @example
 * druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4
 * @returns {Number} The cosine distance between {@link a} and {@link b}.
 */
function cosine (a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += a[i] * b[i];
        sum_a += a[i] * a[i];
        sum_b += b[i] * b[i];
    }
    return Math.acos(sum / (Math.sqrt(sum_a) * Math.sqrt(sum_b)));
}

/**
 * Computes the manhattan distance (l<sub>1</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the manhattan distance between {@link a} and {@link b}.
 */ 
function manhattan (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum;
}

/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the chebyshev distance between {@link a} and {@link b}.
 */
function chebyshev (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res);
}

/**
 * Computes the canberra distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias canberra
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} The canberra distance between {@link a} and {@link b}.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
function canberra(a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += (Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i])));
    }
    return sum;
}

/**
 * Computes the jaccard distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias jaccard
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the jaccard distance between {@link a} and {@link b}.
 */
function jaccard (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let num_non_zero = 0;
    let num_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_non_zero += x || y;
        num_equal += x && y;
    }
    return (num_non_zero - num_equal) / num_non_zero;
}

/**
 * Computes the hamming distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias hamming
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the hamming distance between {@link a} and {@link b}.
 */
function hamming (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let disagree = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i];
        const y = b[i];
        disagree += x != y;
    }
    return disagree / n;
}

/**
 * Computes the Sokal-Michener distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias sokal_michener
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the Sokal-Michener distance between {@link a} and {@link b}.  
 */
function sokal_michener(a, b) {
    if (a.length != b.length) return undefined
    const n = a.length;
    let num_not_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_not_equal += x != y;
    }
    return (2 * num_not_equal) / (n + num_not_equal);
}

/**
 * Computes the yule distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias yule
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the yule distance between {@link a} and {@link b}.
 */
function yule (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let num_true_true = 0;
    let num_true_false = 0;
    let num_false_true = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_true_true += x && y;
        num_true_false += x && !y;
        num_false_true += !x && x;
    }
    const num_false_false = n - num_true_true - num_true_false - num_false_true;
    return num_true_false == 0 || num_false_true == 0 ? 0 : (2 * num_true_false * num_false_true) / (num_true_true * num_false_false + num_true_false * num_false_true);
}

/**
 *
 * @param {*} A
 * @param {*} k
 * @param {*} distance_matrix
 * @param {*} metric
 */
function k_nearest_neighbors (A, k, distance_matrix$1 = null, metric = euclidean) {
    const rows = A.shape[0];
    let D = distance_matrix$1 ?? distance_matrix(A, metric);
    let nN = new Array(rows);
    for (let row = 0; row < rows; ++row) {
        nN[row] = Array.from(D.row(row))
            .map((distance, col) => {
                return {
                    i: row,
                    j: col,
                    distance: distance,
                };
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1);
    }
    return nN;
}

/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
class Matrix {
    /**
     * creates a new Matrix. Entries are stored in a Float64Array.
     * @constructor
     * @memberof module:matrix
     * @alias Matrix
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
        } else if (Array.isArray(A) || A instanceof Float64Array) {
            let m = A.length;
            if (m === 0) throw new Error("Array is empty");
            // 1d
            if (!Array.isArray(A[0]) && !(A[0] instanceof Float64Array)) {
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
            } else if (Array.isArray(A[0]) || A[0] instanceof Float64Array) {
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
     * @param {int} row
     * @param {Array} values
     * @returns {Matrix}
     */
    set_row(row, values) {
        let cols = this._cols;
        if (Array.isArray(values) && values.length === cols) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this.values[offset + col] = values[col];
            }
        } else if (values instanceof Matrix && values.shape[1] === cols && values.shape[0] === 1) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this.values[offset + col] = values._data[col];
            }
        }
        return this;
    }

    /**
     * Returns the {@link col}<sup>th</sup> column from the Matrix.
     * @param {int} col
     * @returns {Array}
     */
    col(col) {
        let result_col = new Float64Array(this._rows);
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
        let B = new Matrix(rows, 2 * cols, (i, j) => {
            if (j >= cols) {
                return i === j - cols ? 1 : 0;
            } else {
                return this.entry(i, j);
            }
        });
        let h = 0;
        let k = 0;
        while (h < rows && k < cols) {
            var i_max = 0;
            let max_val = -Infinity;
            for (let i = h; i < rows; ++i) {
                let val = Math.abs(B.entry(i, k));
                if (max_val < val) {
                    i_max = i;
                    max_val = val;
                }
            }
            if (B.entry(i_max, k) == 0) {
                k++;
            } else {
                // swap rows
                for (let j = 0; j < 2 * cols; ++j) {
                    let h_val = B.entry(h, j);
                    let i_val = B.entry(i_max, j);
                    B.set_entry(h, j, h_val);
                    B.set_entry(i_max, j, i_val);
                }
                for (let i = h + 1; i < rows; ++i) {
                    let f = B.entry(i, k) / B.entry(h, k);
                    B.set_entry(i, k, 0);
                    for (let j = k + 1; j < 2 * cols; ++j) {
                        B.set_entry(i, j, B.entry(i, j) - B.entry(h, j) * f);
                    }
                }
                h++;
                k++;
            }
        }

        for (let row = 0; row < rows; ++row) {
            let f = B.entry(row, row);
            for (let col = row; col < 2 * cols; ++col) {
                B.set_entry(row, col, B.entry(row, col) / f);
            }
        }

        for (let row = rows - 1; row >= 0; --row) {
            let B_row_row = B.entry(row, row);
            for (let i = 0; i < row; i++) {
                let B_i_row = B.entry(i, row);
                let f = B_i_row / B_row_row;
                for (let j = i; j < 2 * cols; ++j) {
                    let B_i_j = B.entry(i, j);
                    let B_row_j = B.entry(row, j);
                    B_i_j = B_i_j - B_row_j * f;
                    B.set_entry(i, j, B_i_j);
                }
            }
        }

        return new Matrix(rows, cols, (i, j) => B.entry(i, j + cols));
    }

    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    dot(B) {
        if (B instanceof Matrix) {
            let A = this;
            if (A.shape[1] !== B.shape[0]) {
                throw new Error(`A.dot(B): A is a ${A.shape.join(" ⨯ ")}-Matrix, B is a ${B.shape.join(" ⨯ ")}-Matrix: 
                A has ${A.shape[1]} cols and B ${B.shape[0]} rows. 
                Must be equal!`);
            }
            let I = A.shape[1];
            let C = new Matrix(A.shape[0], B.shape[1], (row, col) => {
                const A_i = A.row(row);
                const B_i = B.col(col);
                let sum = 0;
                for (let i = 0; i < I; ++i) {
                    sum += A_i[i] * B_i[i];
                }
                return sum;
            });
            return C;
        } else if (Array.isArray(B) || B instanceof Float64Array) {
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
        let [rows, cols] = B.shape;
        for (let row = 0; row < rows; ++row) {
            if (row > this._rows) {
                continue;
            }
            for (let col = 0; col < cols; ++col) {
                if (col > this._cols) {
                    continue;
                }
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
     * @param {function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    _apply_array(f, v) {
        const data = this.values;
        const [rows, cols] = this.shape;
        for (let row = 0; row < rows; ++row) {
            const offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = offset + col;
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
        for (let row = 0; row < rows; ++row) {
            const offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = offset + col;
                data[i] = f(data[i], values[row]);
            }
        }
        return this;
    }

    _apply(value, f) {
        let data = this.values;
        if (value instanceof Matrix) {
            let [value_rows, value_cols] = value.shape;
            let [rows, cols] = this.shape;
            if (value_rows === 1) {
                if (cols !== value_cols) {
                    throw new Error(`cols !== value_cols`);
                }
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(0, col));
                    }
                }
            } else if (value_cols === 1) {
                if (rows !== value_rows) {
                    throw new Error(`rows !== value_rows`);
                }
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, 0));
                    }
                }
            } else if (rows == value_rows && cols == value_cols) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, col));
                    }
                }
            } else {
                throw new Error(`error`);
            }
        } else if (Array.isArray(value)) {
            let rows = this._rows;
            let cols = this._cols;
            if (value.length === rows) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[row]);
                    }
                }
            } else if (value.length === cols) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[col]);
                    }
                }
            } else {
                throw new Error(`error`);
            }
        } else {
            for (let i = 0, n = this._rows * this._cols; i < n; ++i) {
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
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.mult(2); // [[2, 4], [6, 8]];
     * A.mult(B); // [[1, 4], [9, 16]];
     */
    mult(value) {
        return this.clone()._apply(value, (a, b) => a * b);
    }

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
    divide(value) {
        return this.clone()._apply(value, (a, b) => a / b);
    }

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
    add(value) {
        return this.clone()._apply(value, (a, b) => a + b);
    }

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
    sub(value) {
        return this.clone()._apply(value, (a, b) => a - b);
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
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < cols; ++col) {
                this._data[row * cols + col] = value(row, col);
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
     * Returns the sum oof all entries of the Matrix.
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
        for (let row = 0; row < rows; ++row) {
            result[row] = 0;
            for (let col = 0; col < cols; ++col) {
                result[row] += data[row * cols + col];
            }
            result[row] /= cols;
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
            result[col] = 0;
            for (let row = 0; row < rows; ++row) {
                result[col] += data[row * cols + col];
            }
            result[col] /= rows;
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
                const alpha = r.T.dot(r).entry(0, 0) / d.T.dot(z).entry(0, 0);
                x = x.add(d.mult(alpha));
                const r_next = r.sub(z.mult(alpha));
                const beta = r_next.T.dot(r_next).entry(0, 0) / r.T.dot(r).entry(0, 0);
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
                x.set_entry(0, row, x.entry(0, row) - L.entry(row, col) * x.entry(1, col));
            }
            x.set_entry(0, row, x.entry(0, row) / L.entry(row, row));
        }

        // backward
        for (let row = rows - 1; row >= 0; --row) {
            for (let col = rows - 1; col > row; --col) {
                x.set_entry(0, row, x.entry(0, row) - U.entry(row, col) * x.entry(0, col));
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
        const MT = M.T;
        let MtM = MT.dot(M);
        let MMt = M.dot(MT);
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
}

/**
 * Computes the distance matrix of datamatrix {@link A}.
 * @param {Matrix} A - Matrix
 * @param {Function} [metric=euclidean] - The diistance metric.
 * @returns {Matrix} D - The distance matrix of {@link A}.
 */
function distance_matrix (A, metric = euclidean) {
    let n = A.shape[0];
    const D = new Matrix(n, n);
    for (let i = 0; i < n; ++i) {
        const A_i = A.row(i);
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, A.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}

/**
 * @memberof module:matrix
 * Creates an Array containing {@link number} numbers from {@link start} to {@link end}. If <code>{@link number} = null</null>
 * @param {Number} start
 * @param {Number} end
 * @param {Number} [number = null]
 * @returns {Array}
 */
function linspace (start, end, number = null) {
    if (!number) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return number === 1 ? [start] : [];
    }
    let result = new Array(number);
    number -= 1;
    for (let i = number; i >= 0; --i) {
        result[i] = (i * end + (number - i) * start) / number;
    }
    return result;
}

//import { neumair_sum } from "../numerical/index";

function norm(v, metric = euclidean) {
//export default function(vector, p=2, metric = euclidean) {
    let vector = null;
    if (v instanceof Matrix) {
        let [rows, cols] = v.shape;
        if (rows === 1) vector = v.row(0);
        else if (cols === 1) vector = v.col(0);
        else throw "matrix must be 1d!"
    } else {
        vector = v;
    }
    let n = vector.length;
    let z = new Array(n);
    z.fill(0);
    return metric(vector, z);
    
    
    /*let v;
    if (vector instanceof Matrix) {
        let [ rows, cols ] = v.shape;
        if (rows === 1) {
            v = vector.row(0);
        } else if (cols === 1) {
            v = vector.col(0);
        } else {
            throw "matrix must be 1d"
        }
    } else {
        v = vector;
    }
    return Math.pow(neumair_sum(v.map(e => Math.pow(e, p))), 1 / p)*/
}

/**
 * @class
 * @memberof module:utils
 * @alias Randomizer
 */
class Randomizer {
    /**
     * Mersenne Twister random number generator.
     * @constructor
     * @param {Number} [_seed=new Date().getTime()] - The seed for the random number generator. If <code>_seed == null</code> then the actual time gets used as seed.
     * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
     */
    constructor(_seed) {
        this._N = 624;
        this._M = 397;
        this._MATRIX_A = 0x9908b0df;
        this._UPPER_MASK = 0x80000000;
        this._LOWER_MASK = 0x7fffffff;
        this._mt = new Array(this._N);
        this._mti = this.N + 1;

        this.seed = _seed || new Date().getTime();
        return this;
    }

    set seed(_seed) {
        this._seed = _seed;
        let mt = this._mt;

        mt[0] = _seed >>> 0;
        for (this._mti = 1; this._mti < this._N; this._mti += 1) {
            let mti = this._mti;
            let s = mt[mti - 1] ^ (mt[mti - 1] >>> 30);
            mt[mti] = ((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253 + mti;
            mt[mti] >>>= 0;
        }
    }

    /**
     * Returns the seed of the random number generator.
     * @returns {Number} - The seed.
     */
    get seed() {
        return this._seed;
    }

    /**
     * Returns a float between 0 and 1.
     * @returns {Number} - A random number between [0, 1]
     */
    get random() {
        return this.random_int * (1.0 / 4294967296.0);
    }

    /**
     * Returns an integer between 0 and MAX_INTEGER.
     * @returns {Integer} - A random integer.
     */
    get random_int() {
        let y,
            mag01 = new Array(0x0, this._MATRIX_A);
        if (this._mti >= this._N) {
            let kk;

            /* if (this._mti == this._N + 1) {
                this.seed = 5489;
            } */

            let N_M = this._N - this._M;
            let M_N = this._M - this._N;

            for (kk = 0; kk < N_M; ++kk) {
                y = (this._mt[kk] & this._UPPER_MASK) | (this._mt[kk + 1] & this._LOWER_MASK);
                this._mt[kk] = this._mt[kk + this._M] ^ (y >>> 1) ^ mag01[y & 0x1];
            }
            for (; kk < this._N - 1; ++kk) {
                y = (this._mt[kk] & this._UPPER_MASK) | (this._mt[kk + 1] & this._LOWER_MASK);
                this._mt[kk] = this._mt[kk + M_N] ^ (y >>> 1) ^ mag01[y & 0x1];
            }

            y = (this._mt[this._N - 1] & this._UPPER_MASK) | (this._mt[0] & this._LOWER_MASK);
            this._mt[this._N - 1] = this._mt[this._M - 1] ^ (y >>> 1) ^ mag01[y & 0x1];

            this._mti = 0;
        }

        y = this._mt[(this._mti += 1)];
        y ^= y >>> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >>> 18;

        return y >>> 0;
    }

    /**
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Array|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @returns {Array} - A random selection form {@link A} of {@link n} samples.
     */
    choice(A, n) {
        if (A instanceof Matrix) {
            let rows = A.shape[0];
            if (n > rows) {
                throw new Error("n bigger than A!");
            }
            let sample = new Array(n);
            let index_list = linspace(0, rows - 1);
            for (let i = 0, l = index_list.length; i < n; ++i, --l) {
                let random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map((d) => A.row(d));
        } else if (Array.isArray(A) || A instanceof Float64Array) {
            let rows = A.length;
            if (n > rows) {
                throw new Error("n bigger than A!");
            }
            let sample = new Array(n);
            let index_list = linspace(0, rows - 1);
            for (let i = 0, l = index_list.length; i < n; ++i, --l) {
                let random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map((d) => A[d]);
        }
    }

    /**
     * @static
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Array|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @param {Number} seed - The seed for the random number generator.
     * @returns {Array} - A random selection form {@link A} of {@link n} samples.
     */
    static choice(A, n, seed = 1212) {
        const R = new Randomizer(seed);
        return R.choice(A, n);
        /* let rows = A.shape[0];
        if (n > rows) {
            throw new Error("n bigger than A!");
        }
        let rand = new Randomizer(seed);
        let sample = new Array(n);
        let index_list = linspace(0, rows - 1);
        for (let i = 0, l = index_list.length; i < n; ++i, --l) {
            let random_index = rand.random_int % l;
            sample[i] = index_list.splice(random_index, 1)[0];
        }
        //return result;
        //return new Matrix(n, cols, (row, col) => A.entry(sample[row], col))
        return sample.map((d) => A.row(d)); */
    }
}

/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias max
 * @param {Array} values 
 * @returns {Number}
 */
function max (values) {
    let max;
    for (const value of values) {
        if (value != null && (max < value || (max === undefined && value >= value))) {
            max = value;
        }
    }
    return max;
}

/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias min
 * @param {Array} values
 * @returns {Number}
 */
function min (values) {
    let min;
    for (const value of values) {
        if (value != null && (min > value || (min === undefined && value <= value))) {
            min = value;
        }
    }
    return min;
}

/**
 * @class
 * @alias Heap
 */
class Heap {
    /**
     * A heap is a datastructure holding its elements in a specific way, so that the top element would be the first entry of an ordered list.
     * @constructor
     * @memberof module:datastructure
     * @alias Heap
     * @param {Array=} elements - Contains the elements for the Heap. {@link elements} can be null.
     * @param {Function} [accessor = (d) => d] - Function returns the value of the element.
     * @param {("min"|"max"|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     * @see {@link https://en.wikipedia.org/wiki/Binary_heap}
     */
    constructor(elements = null, accessor = d => d, comparator = "min") {
        if (elements) {
            return Heap.heapify(elements, accessor, comparator);
        } else {
            this._accessor = accessor;
            this._container = [];
            if (comparator == "min") {
                this._comparator = (a, b) => a < b;
            } else if (comparator == "max") {
                this._comparator = (a, b) => a > b;
            } else {
                this._comparator = comparator;
            }
            return this
        }
    }

    /**
     * Creates a Heap from an Array
     * @param {Array|Set} elements - Contains the elements for the Heap.
     * @param {Function=} [accessor = (d) => d] - Function returns the value of the element.
     * @param {(String=|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     */
    static heapify(elements, accessor = d => d, comparator = "min") {
        const heap = new Heap(null, accessor, comparator);
        const container = heap._container;
        for (const e of elements) {
            container.push({
                "element": e,
                "value": accessor(e),
            });
        }
        for (let i = Math.floor((elements.length / 2) - 1); i >= 0; --i) {
            heap._heapify_down(i);
        }
        return heap;
    }

    /**
     * Swaps elements of container array.
     * @private
     * @param {Number} index_a 
     * @param {Number} index_b 
     */
    _swap(index_a, index_b) {
        const container = this._container;
        [container[index_b], container[index_a]] = [container[index_a], container[index_b]];
        return;
    }

    /**
     * @private
     */
    _heapify_up() {
        const container = this._container;
        let index = container.length - 1;
        while (index > 0) {
            let parentIndex = Math.floor((index - 1) / 2);
            if (!this._comparator(container[index].value, container[parentIndex].value)) {
                break;
            } else {
            this._swap(parentIndex, index);
            index = parentIndex;
            }
        }
    }

    /**
     * Pushes the element to the heap.
     * @param {} element
     * @returns {Heap}
     */
    push(element) {
        const value = this._accessor(element);
        //const node = new Node(element, value);
        const node = {"element": element, "value": value};
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    /**
     * @private
     * @param {Number} [start_index = 0] 
     */
    _heapify_down(start_index=0) {
        const container = this._container;
        const comparator = this._comparator;
        const length = container.length;
        let left = 2 * start_index + 1;
        let right = 2 * start_index + 2;
        let index = start_index;
        if (index > length) throw "index higher than length"
        if (left < length && comparator(container[left].value, container[index].value)) {
            index = left;
        }
        if (right < length && comparator(container[right].value, container[index].value)) {
            index = right;
        }
        if (index !== start_index) {
            this._swap(start_index, index);
            this._heapify_down(index);
        }
    }

    /**
     * Removes and returns the top entry of the heap.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    pop() {
        const container = this._container;
        if (container.length === 0) {
            return null;
        } else if (container.length === 1) {
            return container.pop();
        }
        this._swap(0, container.length - 1);
        const item = container.pop();
        this._heapify_down();
        return item;
    }

    /**
     * Returns the top entry of the heap without removing it.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    get first() {
        return this._container.length > 0 ? this._container[0] : null;
    }


    /**
     * Yields the raw data
     * @yields {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    * iterate() {
        for (let i = 0, n = this._container.length; i < n; ++i) {
            yield this._container[i].element;
        }
    }

    /**
     * Returns the heap as ordered array.
     * @returns {Array} Array consisting the elements ordered by {@link comparator}.
     */
    toArray() {
        return this.data()
            .sort((a,b) => this._comparator(a, b) ? -1 : 0)
    }

    /**
     * Returns elements of container array.
     * @returns {Array} Array consisting the elements.
     */
    data() {
        return this._container
            .map(d => d.element)
    }

    /**
     * Returns the container array.
     * @returns {Array} The container array.
     */
    raw_data() {
        return this._container;
    }

    /**
     * The size of the heap.
     * @returns {Number}
     */
    get length() {
        return this._container.length;
    }

    /**
     * Returns false if the the heap has entries, true if the heap has no entries.
     * @returns {Boolean}
     */
    get empty() {
        return this.length === 0;
    }
}

/**
 * @class
 * @alias DisjointSet
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */
class DisjointSet {
    /**
     * @constructor
     * @alias DisjointSet
     * @memberof module:datastructure
     * @param {Array=} elements 
     * @returns {DisjointSet}
     */
    constructor(elements = null) {
        this._list = new Set();
        if (elements) {
            for (const e of elements) {
                this.make_set(e);
            }
        }
        return this;
    }

    make_set(x) {
        const list = this._list;
        if (!list.has(x)) {
            list.add(x);
            x.__disjoint_set = {};
            x.__disjoint_set.parent = x;
            x.__disjoint_set.children = new Set([x]);
            x.__disjoint_set.size = 1;
        }
        return this;
    }

    find(x) {
        const list = this._list;
        if (list.has(x)) {
            if (x.__disjoint_set.parent !== x) {
                x.__disjoint_set.children.add(...x);
                x.__disjoint_set.parent = this.find(x.__disjoint_set.parent);
                return x.__disjoint_set.parent;
            } else {
                return x;
            }
        } else {
            return null;
        }
    }

    union(x, y) {
        let node_x = this.find(x);
        let node_y = this.find(y);

        if (node_x === node_y) return this;
        if (node_x.__disjoint_set.size < node_y.__disjoint_set.size) [node_x, node_y] = [node_y, node_x];

        node_y.__disjoint_set.parent = node_x;
        // keep track of children?
        node_y.__disjoint_set.children.forEach(node_x.__disjoint_set.children.add, node_x.__disjoint_set.children);
        node_x.__disjoint_set.size += node_y.__disjoint_set.size;

        return this;
    }
}

/**
 * @class
 * @alias BallTree
 */
class BallTree {
    /**
     * Generates a BallTree with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias BallTree
     * @param {Array=} elements - Elements which should be added to the BallTree
     * @param {Function} [metric = euclidean] metric to use: (a, b) => distance
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     * @returns {BallTree}
     */
    constructor(elements = null, metric = euclidean) {
        this._Node = class {
            constructor(pivot, child1=null, child2=null, radius=null) {
                this.pivot = pivot;
                this.child1 = child1;
                this.child2 = child2;
                this.radius = radius;
            }
        };
        this._Leaf = class {
            constructor(points) {
                this.points = points;
            }
        };
        this._metric = metric;
        if (elements) {
            this.add(elements);
        }
        return this;
    }

    /**
     * 
     * @param {Array<*>} elements - new elements.
     * @returns {BallTree}
     */
    add(elements) {
        elements = elements.map((element, index) => {
            return {index: index, element: element}
        });
        this._root = this._construct(elements);
        return this;
    }

    /**
     * @private
     * @param {Array<*>} elements 
     * @returns {Node} root of balltree.
     */
    _construct(elements) {
        if (elements.length === 1) {
            return new this._Leaf(elements);
        } else {
            let c = this._greatest_spread(elements);
            let sorted_elements = elements.sort((a, b) => a.element[c] - b.element[c]);
            let n = sorted_elements.length;
            let p_index = Math.floor(n / 2);
            let p = elements[p_index];
            let L = sorted_elements.slice(0, p_index);
            let R = sorted_elements.slice(p_index, n);
            let radius = Math.max(...elements.map(d => this._metric(p.element, d.element)));
            let B;
            if (L.length > 0 && R.length > 0) {         
                B = new this._Node(p, this._construct(L), this._construct(R), radius);
            } else {
                B = new this._Leaf(elements);
            }
            return B;
        }
    }

    /**
     * @private
     * @param {Node} B 
     * @returns {Number}
     */
    _greatest_spread(B) {
        let d = B[0].element.length;
        let start = new Array(d);

        for (let i = 0; i < d; ++i) {
            start[i] = [Infinity, -Infinity];
        }

        let spread = B.reduce((acc, current) => {
            for (let i = 0; i < d; ++i) {
                acc[i][0] = Math.min(acc[i][0], current.element[i]);
                acc[i][1] = Math.max(acc[i][1], current.element[i]);
            }
            return acc;
        }, start);
        spread = spread.map(d => d[1] - d[0]);
        
        let c = 0;
        for (let i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c;
    }

    /**
     * 
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t, k = 5) {
        return this._search(t, k, new Heap(null, d => this._metric(d.element, t), "max"), this._root);
    }

    /**
     * @private
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @param {Heap} Q - Heap consists of the currently found {@link k} nearest neighbors.
     * @param {Node|Leaf} B 
     */
    _search(t, k, Q, B) {
        // B is Node
        if (Q.length >= k && B.pivot && B.radius && this._metric(t, B.pivot.element) - B.radius >= Q.first.value) {
            return Q;
        } 
        if (B.child1) this._search(t, k, Q, B.child1);
        if (B.child2) this._search(t, k, Q, B.child2);
        
        // B is leaf
        if (B.points) {
            for (let i = 0, n = B.points.length; i < n; ++i) {
                let p = B.points[i];
                if (k > Q.length) {
                    Q.push(p);
                } else {
                    Q.push(p);
                    Q.pop();
                }
            }
        }
        return Q;
    }
}

/**
 * @class
 * @alias KNN
 */
class KNN {
    /**
     * Generates a KNN list with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias KNN
     * @param {Array=} elements - Elements which should be added to the KNN list
     * @param {Function|"precomputed"} [metric = euclidean] metric is either precomputed or a function to use: (a, b) => distance
     * @returns {KNN}
     */
    constructor(elements=null, metric=euclidean) {
        this._metric = metric;
        this._elements = elements instanceof Matrix ? elements : Matrix.from(elements);
        const N = this._elements.shape[0];
        if (metric === "precomputed") {
            this._D = this._elements.clone();
        } else {
            this._D = distance_matrix(this._elements, metric);
        }
        this.KNN = [];
        for (let row = 0; row < N; ++row) {
            const distances = this._D.row(row);
            const H = new Heap(null, d => d.value, "min");
            for (let j = 0; j < N; ++j) {
                H.push({
                    value: distances[j],
                    index: j,
                });
            }
            this.KNN.push(H);
        }
    }

    /**
     * 
     * @param {Array|Number} t - query element or index.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t, k = 5) {
        const metric = this._metric;
        const KNN = this.KNN;
        let H;
        if (Array.isArray(t)) {
            if (this._metric == "precomputed") {
                throw "Search by query element is only possible when not using a precomputed distance matrix!"
            } 
            const elements = this._elements;
            const N = KNN.length;
            let nearest_element_index = null;
            let nearest_dist = Infinity;
            for (let i = 0; i < N; ++i) {
                const element = elements.row(i);
                const dist = metric(t, element);
                if (dist < nearest_dist) {
                    nearest_element_index = i;
                    nearest_dist = dist;
                }
            }
            H = KNN[nearest_element_index];
        } else if (Number.isInteger(t)) {
            H = KNN[t];
        }

        let result = [];
        for (let i = 0; i < k; ++i) {
            result.push(H.pop());
        }
        result.forEach(res => H.push(res.element));
        return result
    }    
}

/**
 * Computes the QR Decomposition of the Matrix {@link A} using Gram-Schmidt process.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */
function qr(A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, cols, "identity");
    const R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        let v = A.col(j);
        for (let i = 0; i < j; ++i) {
            const q = Q.col(i);
            const q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            R.set_entry(i,j, q_dot_v);
            v = v.map((v_, k) => v_ - q_dot_v * q[k]);
        }
        const v_norm = norm(v, euclidean);
        for (let k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j,j, v_norm);
    }
    return {"R": R, "Q": Q};
}

/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Number} [max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [seed=1212] - The seed value or a randomizer used in the algorithm.
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */
function simultaneous_poweriteration$1(A, k = 2, max_iterations=100, seed=1212) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q: Q, R: R } = qr(new Matrix(n, k, () => randomizer.random));
    while (max_iterations--) {
        const oldR = R.clone();
        const Z = A.dot(Q);
        const QR = qr(Z); 
        Q = QR.Q;
        R = QR.R;
        if (neumair_sum(R.sub(oldR).diag) / n < 1e-12) {
            max_iterations = 0;
        }        
    }

    const eigenvalues = R.diag;
    const eigenvectors = Q.transpose().to2dArray;
    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
    };
}

/**
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
class DR {
    //static parameter_list = [];
    get parameter_list() {
        return this._parameter_list;
    }

    set parameter_list(list) {
        this._parameter_list = list;
        return this;
    }
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {seed} [seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X, d = 2, metric = euclidean, seed = 1212) {
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("no valid type for X");
        }
        [this._N, this._D] = this.X.shape;
        this._d = d;
        this._metric = metric;
        this._seed = seed;
        this._randomizer = new Randomizer(seed);
        this._is_initialized = false;
        return this;
    }

    /**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {Number} [value = null] - value of the parameter to set, if <code>value == null</code> then return actual parameter value.
     * @memberof DR
     */
    parameter(name, value = null) {
        if (!this.parameter_list.includes(name)) {
            throw new Error(`${name} is not a valid parameter!`);
        }
        if (value) {
            this[`_${name}`] = value;
            this._is_initialized = false;
            return this;
        } else {
            return this[`_${name}`];
        }
    }

    para(name, value = null) {
        return this.parameter(name, value);
    }

    p(name, value = null) {
        return this.parameter(name, value);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        return this.projection;
    }

    *generator() {
        return this.transform();
    }

    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
    }

    /**
     * @returns {Matrix} Returns the projection.
     */
    get projection() {
        return this._type === "matrix" ? this.Y : this.Y.to2dArray;
    }

    async transform_async(...args) {
        const dr = new this(...args);
        return dr.transform();
    }

    static transform(...args) {
        let dr = new this(...args);
        return dr.transform();
    }

    static async transform_async(...args) {
        return this.transform(...args);
    }

    static *generator(...args) {
        const dr = new this(...args);
        const gen = dr.generator();
        for (const res of gen) {
            yield res;
        }
    }
}

/**
 * @class
 * @alias PCA
 * @augments DR
 */
class PCA extends DR {
    /**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @returns {PCA}
     */
    constructor(X, d = 2) {
        super(X, d);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}. If parameter {@link A} is given, then project {@link A} with the principal components of {@link X}.
     * @param {null|Matrix|Array} [A = null] - If given, the data to project.
     * @returns {Matrix|Array} - The projected data.
     */
    transform(A = null) {
        const V = this.principal_components();
        if (A == null) {
            const X = this.X;
            this.Y = X.dot(V);
            return this.projection;
        } else if (Array.isArray(A)) {
            return druid.Matrix.from(A).dot(V).asArray;
        } else if (A instanceof Matrix) {
            return A.dot(V);
        } else {
            throw new Error("No valid type for A!");
        }
    }

    /**
     * Computes the {@link d} principal components of Matrix {@link X}.
     * @returns {Matrix}
     */
    principal_components() {
        if (this.V) {
            return this.V;
        }
        const X = this.X;
        const means = Matrix.from(X.meanCols);
        const X_cent = X.sub(means);
        const C = X_cent.transpose().dot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration$1(C, this._d);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }
}

/**
 * @class
 * @alias MDS
 * @extends DR
 */
class MDS extends DR{
    /**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|Array}
     */
    transform() {
        const X = this.X;
        const rows = X.shape[0];
        const metric = this._metric;
        const A = metric === "precomputed" ? X : distance_matrix(X, metric); 
        const ai_ = A.meanCols;
        const a_j = A.meanRows;
        const a__ = A.mean;

        this._d_X = A;
        const B = new Matrix(rows, rows, (i, j) => (A.entry(i, j) - ai_[i] - a_j[j] + a__));
                
        const { eigenvectors: V } = simultaneous_poweriteration$1(B, this._d);
        this.Y = Matrix.from(V).transpose();
        
        return this.projection;
    }

    /**
     * @returns {Number} - the stress of the projection.
     */
    stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X; /*new Matrix();
        d_X.shape = [N, N, (i, j) => {
            return i < j ? metric(X.row(i), X.row(j)) : d_X.entry(j, i);
        }]*/
        const d_Y = new Matrix();
        d_Y.shape = [N, N, (i, j) => {
            return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
        }];
        let top_sum = 0;
        let bottom_sum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                top_sum += Math.pow(d_X.entry(i, j) - d_Y.entry(i, j), 2);
                bottom_sum += Math.pow(d_X.entry(i, j), 2);
            }
        }
        return Math.sqrt(top_sum / bottom_sum);
    }
}

/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */
class ISOMAP extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [d = 2] - the dimensionality of the projection. 
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points. 
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, neighbors, d = 2, metric = euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k"];
        this.parameter("k", Math.min(neighbors ?? Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N -1));
        return this;
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        const X = this.X;
        const rows = this._N;
        const metric = this._metric;
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix();
        D.shape = [rows, rows, (i,j) => i <= j ? metric(X.row(i), X.row(j)) : D.entry(j,i)];
        const kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            const row = [];
            for (let j = 0; j < rows; ++j) {
                row.push({
                    "index": j,
                    "distance": D.entry(i, j),
                });
            }
            const H = new Heap(row, d => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, this._k + 1));
        }
        
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
        const G = new Matrix(rows, rows, (i,j) => {
            const other = kNearestNeighbors[i].find(n => n.index === j);
            return other ? other.distance : Infinity
        });

        for (let i = 0; i < rows; ++i) {
            for (let j = 0; j < rows; ++j) {
                for (let k = 0; k < rows; ++k) {
                    G.set_entry(i, j, Math.min(G.entry(i, j), G.entry(i, k) + G.entry(k, j)));
                }
            }
        }
        
        let ai_ = new Float64Array(rows);
        let a_j = new Float64Array(rows);
        let a__ = 0;
        let A = new Matrix(rows, rows, (i,j) => {
            let val = G.entry(i, j);
            val = val === Infinity ? 0 : val;
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        });
        
        ai_ = ai_.map(v => v / rows);
        a_j = a_j.map(v => v / rows);
        a__ /= (rows ** 2);
        const B = new Matrix(rows, rows, (i,j) => (A.entry(i,j) - ai_[i] - a_j[j] + a__));
             
        // compute d eigenvectors
        const { eigenvectors: V } = simultaneous_poweriteration$1(B, this._d);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.projection;
    }


}

/**
 * @class
 * @alias FASTMAP
 * @extends DR
 */
class FASTMAP extends DR{
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        return this;
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {function} dist 
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */
    _choose_distant_objects(dist) {
        const X = this.X;
        const N = X.shape[0];
        let a_index = this._randomizer.random_int % N - 1;
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_bi = dist(b_index, i);
            if (d_bi > max_dist) {
                max_dist = d_bi;
                a_index = i;
            }
        }
        return [a_index, b_index, max_dist];
    }

    /**
     * Computes the projection.
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */
    transform() {
        const X = this.X;
        const N = X.shape[0];
        const d = this._d;
        const metric = this._metric;
        const Y = new Matrix(N, d, 0);
        let dist = (a, b) => metric(X.row(a), X.row(b));

        for (let _col = 0; _col < d; ++_col) {
            let old_dist = dist;
            // choose pivot objects
            const [a_index, b_index, d_ab] = this._choose_distant_objects(dist);
            // record id of pivot objects
            //PA[0].push(a_index);
            //PA[1].push(b_index);
            /* if (d_ab === 0) {
                // because all inter-object distances are zeros
                for (let i = 0; i < N; ++i) {
                    Y.set_entry(i, _col, 0);
                }
            } else { */
            if (d_ab !== 0) {
                // project the objects on the line (O_a, O_b)
                for (let i = 0; i < N; ++i) {
                    const d_ai = dist(a_index, i);
                    const d_bi = dist(b_index, i);
                    const y_i = (d_ai ** 2 + d_ab ** 2 - d_bi ** 2) / (2 * d_ab);
                    Y.set_entry(i, _col, y_i);
                }
                // consider the projections of the objects on a
                // hyperplane perpendicluar to the line (a, b);
                // the distance function D'() between two 
                // projections is given by Eq.4
                dist = (a, b) => Math.sqrt(old_dist(a, b) ** 2 - (Y.entry(a, _col) - Y.entry(b, _col)) ** 2);
            }
        }
        // return embedding
        this.Y = Y;
        return this.projection;
    }
}

/**
 * @class
 * @alias LDA
 * @extends DR
 */
class LDA extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Array} labels - the label / class of each data point.
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, labels, d = 2, metric = euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["labels"];
        this.parameter("labels", labels);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        let X = this.X;
        let [ rows, cols ] = X.shape;
        let labels = this._labels;
        let unique_labels = {};
        let label_id = 0;
        labels.forEach((l, i) => {
            if (l in unique_labels) {
                unique_labels[l].count++;
                unique_labels[l].rows.push(X.row(i));
            } else {
                unique_labels[l] = {
                    "id": label_id++,
                    "count": 1,
                    "rows": [X.row(i)]
                };
            }
        });
        
        // create X_mean and vector means;
        let X_mean = X.mean;
        let V_mean = new Matrix(label_id, cols);
        for (let label in unique_labels) {
            let V = Matrix.from(unique_labels[label].rows);
            let v_mean = V.meanCols;
            for (let j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }           
        }
        // scatter_between
        let S_b = new Matrix(cols, cols);
        for (let label in unique_labels) {
            let v = V_mean.row(unique_labels[label].id);
            let m = new Matrix(cols, 1, (j) => v[j] - X_mean);
            let N = unique_labels[label].count;
            S_b = S_b.add(m.dot(m.transpose()).mult(N));
        }

        // scatter_within
        let S_w = new Matrix(cols, cols);
        for (let label in unique_labels) {
            let v = V_mean.row(unique_labels[label].id);
            let m = new Matrix(cols, 1, (j) => v[j]);
            let R = unique_labels[label].rows;
            for (let i = 0, n = unique_labels[label].count; i < n; ++i) {
                let row_v = new Matrix(cols, 1, (j,_) => R[i][j] - m.entry(j, 0));
                S_w = S_w.add(row_v.dot(row_v.transpose()));
            }
        }

        let { eigenvectors: V } = simultaneous_poweriteration$1(S_w.inverse().dot(S_b), this._d);
        V = Matrix.from(V).transpose();
        this.Y = X.dot(V);

        // return embedding
        return this.projection;
    }
}

/**
 * @class
 * @alias LLE
 * @extends DR
 */
class LLE extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, neighbors, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k"];
        this.parameter("k", Math.min(neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const d = this._d;
        const rows = this._N;
        const cols = this._D;
        const k = this.parameter("k");
        const nN = k_nearest_neighbors(X, k, null, this._metric);
        const O = new Matrix(k, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const nN_row = nN[row];
            const Z = new Matrix(k, cols, (i, j) => X.entry(nN_row[i].j, j) - X.entry(row, j));
            const C = Z.dot(Z.T);
            if ( k > cols ) {
                const C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < k; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }
            // reconstruct;
            let w = Matrix.solve_CG(C, O, this._randomizer);
            w = w.divide(w.sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        const I = new Matrix(rows, rows, "identity");
        const IW = I.sub(W);
        const M = IW.T.dot(IW);
        const { eigenvectors: V } = simultaneous_poweriteration$1(M.T.inverse(), d + 1);
        this.Y = Matrix.from(V.slice(1, 1 + d)).T;

        // return embedding
        return this.projection;
    }
}

/**
 * @class
 * @alias LTSA
 * @extends DR
 */
class LTSA extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, neighbors, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k"];
        this.parameter("k", Math.min(neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (this._D <= d) throw `Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${d})!`;
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const d = this._d;
        const [ rows, D ] = X.shape;
        const k = this.parameter("k");
        // 1.1 determine k nearest neighbors
        const nN = k_nearest_neighbors(X, k, null, this._metric);
        // center matrix
        const O = new Matrix(D, D, "center");
        const B = new Matrix(rows, rows, 0);
        
        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            const I_i = [row, ...nN[row].map(n => n.j)];
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            // center X_i
            X_i = X_i.dot(O);
            // correlation matrix
            const C = X_i.dot(X_i.transpose());
            const { eigenvectors: g } = simultaneous_poweriteration$1(C, d);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            const G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            const W_i = G_i_t.transpose().dot(G_i_t).add(1 / Math.sqrt(k + 1));
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < k + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }

        // 3. Aligning global coordinates
        const { eigenvectors: Y } = simultaneous_poweriteration$1(B, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.projection;
    }
}

/**
 * @class
 * @alias TSNE
 * @extends DR
 */
class TSNE extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [perplexity = 50] - perplexity.
     * @param {Number} [epsilon = 10] - learning parameter.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TSNE}
     */
    constructor(X, perplexity=50, epsilon=10, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["perplexity", "epsilon"];
        [ this._N, this._D ] = this.X.shape;
        this.parameter("perplexity", Math.min(perplexity, this._N - 1));
        this.parameter("epsilon", epsilon);
        this._iter = 0;
        this.Y = new Matrix(this._N, this._d, () => this._randomizer.random);
        return this;
    }

    /**
     * 
     * @param {Matrix} distance_matrix - accepts a precomputed distance matrix
     * @returns {TSNE}
     */
    init(distance_matrix=null) {
        // init
        const Htarget = Math.log(this._perplexity);
        const N = this._N;
        const D = this._D;
        const metric = this._metric;
        const X = this.X;
        let Delta;
        if (distance_matrix) {
            Delta = distance_matrix;
        } else {
            Delta = new Matrix(N, N);
            for (let i = 0; i < N; ++i) {
                const X_i = X.row(i);
                for (let j = i + 1; j < N; ++j) {
                    const distance = metric(X_i, X.row(j));
                    Delta.set_entry(i, j, distance);
                    Delta.set_entry(j, i, distance);
                }
            }

        } 
            
        const P = new Matrix(N, N, "zeros");

        this._ystep = new Matrix(N, D, "zeros");
        this._gains = new Matrix(N, D, 1);

        // search for fitting sigma
        let prow = new Array(N).fill(0);
        const tol = 1e-4;
        const maxtries = 50;
        for (let i = 0; i < N; ++i) {
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let done = false;

            let num = 0;
            while(!done) {
                let psum = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = Math.exp(-Delta.entry(i, j) * beta);
                    if (i === j) pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                let Hhere = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = (psum === 0) ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) {
                        Hhere -= pj * Math.log(pj);
                    }
                }
                if (Hhere > Htarget) {
                    betamin = beta;
                    beta = (betamax === Infinity) ? (beta * 2) : ((beta + betamax) / 2);
                } else {
                    betamax = beta;
                    beta = (betamin === -Infinity) ? (beta / 2) : ((beta + betamin) / 2);
                }
                ++num;
                if (Math.abs(Hhere - Htarget) < tol) done = true;
                if (num >= maxtries) done = true;
            }

            for (let j = 0; j < N; ++j) {
                P.set_entry(i, j, prow[j]);
            }
        }

        //compute probabilities
        const Pout = new Matrix(N, N, "zeros");
        const N2 = N * 2;
        for (let i = 0; i < N; ++i) {
            for (let j = i; j < N; ++j) {
                const p = Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100);
                Pout.set_entry(i, j, p);
                Pout.set_entry(j, i, p);
            }
        }
        this._P = Pout;
        return this;
    }

    /**
     * 
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */
    transform(iterations=500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     * 
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */
    * generator(iterations=500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */
    next() {
        const iter = ++this._iter;
        const P = this._P;
        const ystep = this._ystep;
        const gains = this._gains;
        const N = this._N;
        const epsilon = this._epsilon;
        const dim = this._d;
        let Y = this.Y;

        //calc cost gradient;
        const pmul = iter < 100 ? 4 : 1;
        
        // compute Q dist (unnormalized)
        const Qu = new Matrix(N, N, "zeros");
        let qsum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dsum = 0;
                for (let d = 0; d < dim; ++d) {
                    const dhere = Y.entry(i, d) - Y.entry(j, d);
                    dsum += dhere * dhere;
                }
                const qu = 1 / (1 + dsum);
                Qu.set_entry(i, j, qu);
                Qu.set_entry(j, i, qu);
                qsum += 2 * qu;
            }
        }

        // normalize Q dist
        const Q = new Matrix(N, N, 0);
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                const val = Math.max(Qu.entry(i, j) / qsum, 1e-100);
                Q.set_entry(i, j, val);
                Q.set_entry(j, i, val);
            }
        }

        const grad = new Matrix(N, dim, "zeros");
        for (let i = 0; i < N; ++i) {
            for (let j = 0; j < N; ++j) {
                const premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (let d = 0; d < dim; ++d) {
                    grad.set_entry(i, d, grad.entry(i, d) + premult * (Y.entry(i, d) - Y.entry(j, d)));
                }
            }
        }

        // perform gradient step
        let ymean = new Float64Array(dim);
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const gid = grad.entry(i, d);
                const sid = ystep.entry(i, d);
                const gainid = gains.entry(i, d);
                
                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * .8 : gainid + .2;
                if (newgain < .01) newgain = .01;
                gains.set_entry(i, d, newgain);

                const momval = iter < 250 ? .5 : .8;
                const newsid = momval * sid - epsilon * newgain * gid;
                ystep.set_entry(i, d, newsid);

                Y.set_entry(i, d, Y.entry(i, d) + newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < 2; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - ymean[d] / N);
            }
        }

        return this.Y;
    }
}

/**
 *
 * @memberof module:optimization
 * @alias powell
 * @param {Function} f
 * @param {Array} x0
 * @param {Number} [max_iter = 300]
 * @returns {Array}
 * @see http://optimization-js.github.io/optimization-js/optimization.js.html#line438
 */
function powell (f, x0, max_iter = 300) {
    const epsilon = 1e-2;
    const n = x0.length;
    let alpha = 1e-3;
    let pfx = 10000;
    let x = x0.slice();
    let fx = f(x);
    let convergence = false;

    while (max_iter-- >= 0 && !convergence) {
        convergence = true;
        for (let i = 0; i < n; ++i) {
            x[i] += 1e-6;
            let fxi = f(x);
            x[i] -= 1e-6;
            let dx = (fxi - fx) / 1e-6;
            if (Math.abs(dx) > epsilon) {
                convergence = false;
            }
            x[i] -= alpha * dx;
            fx = f(x);
        }
        alpha *= pfx >= fx ? 1.05 : 0.4;
        pfx = fx;
    }
    return x;
}

/**
 * @class
 * @alias UMAP
 * @extends DR
 */
class UMAP extends DR {

    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias UMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {UMAP}
     */
    constructor(X, n_neighbors=15, local_connectivity=1, min_dist=1, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["n_neighbors", "local_connectivity", "min_dist"];
        [ this._N, this._D ] = this.X.shape;
        n_neighbors = Math.min(this._N - 1, n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(local_connectivity, n_neighbors - 1));
        this.parameter("min_dist", min_dist);
        this._iter = 0;
        this._spread = 1;
        this._set_op_mix_ratio = 1;
        this._repulsion_strength = 1;
        this._negative_sample_rate = 5;
        this._n_epochs = 350;
        this._initial_alpha = 1;
        this.Y = new Matrix(this._N, this._d, () => this._randomizer.random);
        return this;
    }

    /**
     * @private
     * @param {Number} spread 
     * @param {Number} min_dist 
     * @returns {Array}
     */
    _find_ab_params(spread, min_dist) {
        const curve = (x, a, b) => 1 / (1 + a * Math.pow(x, 2 * b));
        const xv = linspace(0, spread * 3, 300);
        const yv = linspace(0, spread * 3, 300);
        
        for (let i = 0, n = xv.length; i < n; ++i) {
            const xv_i = xv[i];
            yv[i] = (xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread));
        }
      
        const err = (p) => {
            const error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map(e => e * e)));
        };
      
        return powell(err, [1, 1]);
    }

    /**
     * @private
     * @param {Array<Array>} distances 
     * @param {Array<Number>} sigmas 
     * @param {Array<Number>} rhos 
     * @returns {Array}
     */
    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            for (let j = 0, m = distances[i].length; j < m; ++j) {
                const v = distances[i][j].value - rhos[i];
                distances[i][j].value = v > 0 ? Math.exp(-v / sigmas[i]) : 1;
            }
        }
        return distances;
    }

    /**
     * @private
     * @param {KNN|BallTree} knn 
     * @param {Number} k 
     * @returns {Object}
     */
    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = this._local_connectivity;
        const target = Math.log2(k);
        const rhos = [];
        const sigmas = [];
        const X = this.X;
        const N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());

        const distances = [];
        if (this._metric === "precomputed") {
            for (let i = 0; i < N; ++i) {
                distances.push(knn.search(i, k).reverse());
            }
        } else {
           for (const x_i of X) {
                distances.push(knn.search(x_i, k).raw_data().reverse());
            }
        }

        for (let i = 0; i < N; ++i) {
            let lo = 0;
            let hi = Infinity;
            let mid = 1;

            const search_result = distances[i];
            const non_zero_dist = search_result.filter(d => d.value > 0);
            const non_zero_dist_length = non_zero_dist.length;
            if (non_zero_dist_length >= local_connectivity) {
                const index = Math.floor(local_connectivity);
                const interpolation = local_connectivity - index;
                if (index > 0) {
                    rhos.push(non_zero_dist[index - 1]);
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rhos[i].value += interpolation * (non_zero_dist[index].value - non_zero_dist[index - 1]);
                    }
                } else {
                    rhos[i].value = interpolation * non_zero_dist[0].value;
                }
            } else if (non_zero_dist_length > 0) {
                rhos[i] = non_zero_dist[non_zero_dist_length - 1].value;
            }
            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    const d = search_result[j].value - rhos[i];
                    psum += (d > 0 ? Math.exp(-(d / mid)) : 1);
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    [hi, mid] = [mid, (lo + hi) / 2];
                } else {
                    if (hi === Infinity) {
                        [lo, mid] = [mid, mid * 2];
                    } else {
                        [lo, mid] = [mid, (lo + hi) / 2];
                    }
                }
            }
            sigmas[i] = mid;

            const mean_ithd = search_result.reduce((a, b) => a + b.value, 0) / search_result.length;
            //let mean_d = null;
            if (rhos[i] > 0) {
                if (sigmas[i] < MIN_K_DIST_SCALE * mean_ithd) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_ithd;
                }
            } else {
                const mean_d = distances.reduce((acc, res) => acc + res.reduce((a, b) => a + b.value, 0) / res.length);
                if (sigmas[i] > MIN_K_DIST_SCALE * mean_d) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_d;
                }
                
            }
        }
        return {
            "distances": distances, 
            "sigmas": sigmas, 
            "rhos": rhos
        }
    }

    /**
     * @private
     * @param {Matrix} X 
     * @param {Number} n_neighbors 
     * @returns {Matrix}
     */
    _fuzzy_simplicial_set(X, n_neighbors) {
        const N = X.shape[0];
        const metric = this._metric;
        const knn = metric === "precomputed" ? new KNN(X, "precomputed") : new BallTree(X.to2dArray, metric);
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        const result = new Matrix(N, N, "zeros");
        for (let i = 0; i < N; ++i) {
            const distances_i = distances[i];
            for (let j = 0; j < distances_i.length; ++j) {
                result.set_entry(i, distances_i[j].element.index, distances_i[j].value);
            }
        }
        const transposed_result = result.T;
        const prod_matrix = result.mult(transposed_result);
        return result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(this._set_op_mix_ratio)
            .add(prod_matrix.mult(1 - this._set_op_mix_ratio));
    }

    /**
     * @private
     * @param {Number} n_epochs 
     * @returns {Array}
     */
    _make_epochs_per_sample(n_epochs) {
        const weights = this._weights;
        const result = new Float32Array(weights.length).fill(-1);
        const weights_max = max(weights);
        const n_samples = weights.map(w => n_epochs * (w / weights_max));
        for (let i = 0; i < result.length; ++i) 
          if (n_samples[i] > 0) result[i] = Math.round(n_epochs / n_samples[i]);
        return result;
    }

    /**
     * @private
     * @param {Matrix} graph 
     * @returns {Object}
     */
    _tocoo(graph) {
        const rows = [];
        const cols = [];
        const data = [];
        const [ rows_n, cols_n ] = graph.shape;
        for (let row = 0; row < rows_n; ++row) {
            for (let col = 0; col < cols_n; ++col) {
                const entry = graph.entry(row, col);
                if (entry !== 0) {
                    rows.push(row);
                    cols.push(col);
                    data.push(entry);
                }
            }
        }
        return {
            "rows": rows, 
            "cols": cols, 
            "data": data
        };
    }

    /**
     * Computes all necessary 
     * @returns {UMAP}
     */
    init() {
        const [ a, b ] = this._find_ab_params(this._spread, this._min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, this._n_neighbors);
        const { rows, cols, data: weights } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        this._weights = weights;
        this._epochs_per_sample = this._make_epochs_per_sample(this._n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map(d => d * this._negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        return this;
    }

    set local_connectivity(value) {
        this._local_connectivity = value;
    }

    get local_connectivity() {
        return this._local_connectivity;
    }

    set min_dist(value) {
        this._min_dist = value;
    }

    get min_dist() {
        return this._min_dist;
    }

    graph() {
        this.check_init();
        return { cols: this._head, rows: this._tail, weights: this._weights };
    }

    /**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    transform(iterations=350) {
        if (this._n_epochs != iterations) {
            this._n_epochs = iterations;
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }


    /**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    * generator(iterations=350) {
        if (this._n_epochs != iterations) {
            this._n_epochs = iterations;
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @private
     * @param {Number} x 
     * @returns {Number}
     */
    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

    /**
     * performs the optimization step.
     * @private
     * @param {Matrix} head_embedding 
     * @param {Matrix} tail_embedding 
     * @param {Matrix} head 
     * @param {Matrix} tail 
     * @returns {Matrix}
     */
    _optimize_layout(head_embedding, tail_embedding, head, tail) {
        const { 
            _d: dim, 
            _alpha: alpha, 
            _repulsion_strength: repulsion_strength, 
            _a: a, 
            _b: b,
            _epochs_per_sample: epochs_per_sample,
            _epochs_per_negative_sample: epochs_per_negative_sample,
            _epoch_of_next_negative_sample: epoch_of_next_negative_sample,
            _epoch_of_next_sample: epoch_of_next_sample,
            _clip: clip
        } = this;
        const tail_length = tail.length;

        for (let i = 0, n = epochs_per_sample.length; i < n; ++i) {
            if (epoch_of_next_sample[i] <= this._iter) {
                const j = head[i];
                const k = tail[i];
                const current = head_embedding.row(j);
                const other = tail_embedding.row(k);
                const dist = euclidean_squared(current, other);
                let grad_coeff = 0;
                if (dist > 0) {
                    grad_coeff = (-2 * a * b * Math.pow(dist, b - 1)) / (a * Math.pow(dist, b) + 1);
                }
                for (let d = 0; d < dim; ++d) {
                    const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                    const c = current[d] + grad_d;
                    const o = other[d] - grad_d;
                    current[d] = c;
                    other[d] = o;
                    head_embedding.set_entry(j, d, c);
                    tail_embedding.set_entry(k, d, o);
                }
                epoch_of_next_sample[i] += epochs_per_sample[i];
                const n_neg_samples = (this._iter - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i];
                for (let p = 0; p < n_neg_samples; ++p) {
                    const k = Math.floor(this._randomizer.random * tail_length);
                    const other = tail_embedding.row(tail[k]);
                    const dist = euclidean_squared(current, other);
                    let grad_coeff = 0;
                    if (dist > 0) {
                        grad_coeff = (2 * repulsion_strength * b) / ((.01 + dist) * (a * Math.pow(dist, b) + 1));
                    } else if (j === k) {
                        continue;
                    }
                    for (let d = 0; d < dim; ++d) {
                        const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                        const c = current[d] + grad_d;
                        const o = other[d] - grad_d;
                        current[d] = c;
                        other[d] = o;
                        head_embedding.set_entry(j, d, c);
                        tail_embedding.set_entry(tail[k], d, o);
                    }
                }
                epoch_of_next_negative_sample[i] += (n_neg_samples * epochs_per_negative_sample[i]);
            }
        }
        return head_embedding;
    }

    /**
     * @private
     * @returns {Matrix}
     */
    next() {
        let iter = ++this._iter;
        let Y = this.Y;

        this._alpha = (this._initial_alpha * (1 - iter / this._n_epochs));
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }
}

/**
 * @class
 * @alias TriMap
 * @extends DR
 */
class TriMap extends DR{
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [weight_adj = 500] - scaling factor.
     * @param {Number} [c = 5] - number of triplets multiplier.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X, weight_adj = 500, c = 5, d = 2, metric = euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["weight_adj", "c"];
        this.parameter("weight_adj", weight_adj);
        this.parameter("c", c);
        return this;
    }

    /**
     * 
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used). 
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used). 
     */
    init(pca = null, knn = null) {
        const X = this.X;
        const N = X.shape[0];
        const d = this._d;
        const metric = this._metric;
        const c = this._c;
        this.n_inliers = 2 * c;
        this.n_outliers = 1 * c;
        this.n_random = 1 * c;
        this.Y = pca || new PCA(X, d).transform();//.mult(.01);
        this.knn = knn || new BallTree(X.to2dArray, metric);
        const {triplets, weights} = this._generate_triplets(this.n_inliers, this.n_outliers, this.n_random);
        this.triplets = triplets;
        this.weights = weights;
        this.lr = 1000 * N / triplets.shape[0];
        this.C = Infinity;
        this.tol = 1e-7;
        this.vel = new Matrix(N, d, 0);
        this.gain = new Matrix(N, d, 1);
        return this;
    }

    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {Number} n_inliers 
     * @param {Number} n_outliers 
     * @param {Number} n_random 
     */
    _generate_triplets(n_inliers, n_outliers, n_random) {
        const metric = this._metric;
        const weight_adj = this._weight_adj;
        const X = this.X;
        const N = X.shape[0];
        const knn = this.knn;
        const n_extra = Math.min(n_inliers + 20, N);
        const nbrs = new Matrix(N, n_extra);
        const knn_distances = new Matrix(N, n_extra);
        for (let i = 0; i < N; ++i) {
            knn.search(X.row(i), n_extra + 1)
                .raw_data()
                .filter(d => d.value != 0)
                .sort((a, b) => a.value - b.value)
                .forEach((d, j) => {
                    nbrs.set_entry(i, j, d.element.index);
                    knn_distances.set_entry(i, j, d.value);
                });
        }
        // scale parameter
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            sig[i] = Math.max(
                   (knn_distances.entry(i, 3) +
                    knn_distances.entry(i, 4) +
                    knn_distances.entry(i, 5) +
                    knn_distances.entry(i, 6)) / 4,
                    1e-10);
        }
        
        const P = this._find_p(knn_distances, sig, nbrs);
        
        let triplets = this._sample_knn_triplets(P, nbrs, n_inliers, n_outliers);
        let n_triplets = triplets.shape[0];
        const outlier_distances = new Float64Array(n_triplets);
        for (let i = 0; i < n_triplets; ++i) {
            const j = triplets.entry(i, 0);
            const k = triplets.entry(i, 2);
            outlier_distances[i] = metric(X.row(j), X.row(k));
        }
        let weights = this._find_weights(triplets, P, nbrs, outlier_distances, sig);
        
        if (n_random > 0) {
            const {random_triplets, random_weights} = this._sample_random_triplets(X, n_random, sig);
            triplets = triplets.concat(random_triplets, "vertical");
            weights = Float64Array.from([...weights, ...random_weights]);
        }
        n_triplets = triplets.shape[0];
        let max_weight = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            if (isNaN(weights[i])) {weights[i] = 0;}
            if (max_weight < weights[i]) max_weight = weights[i];
        }
        let max_weight_2 = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight;
            weights[i] += .0001;
            weights[i] = Math.log(1 + weight_adj * weights[i]);
            if (max_weight_2 < weights[i]) max_weight_2 = weights[i];
        }
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight_2;
        }
        return {
            "triplets": triplets,
            "weights": weights,
        }
    }

    /**
     * Calculates the similarity matrix P
     * @private
     * @param {Matrix} knn_distances - matrix of pairwise knn distances
     * @param {Float64Array} sig - scaling factor for the distances
     * @param {Matrix} nbrs - nearest neighbors
     * @returns {Matrix} pairwise similarity matrix
     */
    _find_p(knn_distances, sig, nbrs) {
        const [N, n_neighbors] = knn_distances.shape;
        return new Matrix(N, n_neighbors, (i, j) => {
            return Math.exp(-((knn_distances.entry(i, j) ** 2) / sig[i] / sig[nbrs.entry(i, j)]));
        });
    }

    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix {@link P}. Row i corresponds to the i-th point.
     * @param {Number} n_inliers - Number of inlier points.
     * @param {Number} n_outliers - Number of outlier points.
     * 
     */
    _sample_knn_triplets(P, nbrs, n_inliers, n_outliers) {
        const N = nbrs.shape[0];
        const triplets = new Matrix(N * n_inliers * n_outliers, 3);
        for (let i = 0; i < N; ++i) {
            let n_i = i * n_inliers * n_outliers;
            const sort_indices = this.__argsort(P.row(i).map(d => -d));
            for (let j = 0; j < n_inliers; ++j) {
                let n_j = j * n_outliers;
                const sim = nbrs.entry(i, sort_indices[j]);
                const samples = this._rejection_sample(n_outliers, N, sort_indices.slice(0, j + 1));
                for (let k = 0; k < n_outliers; ++k) {
                    const index = n_i + n_j + k;
                    const out = samples[k];
                    triplets.set_entry(index, 0, i);
                    triplets.set_entry(index, 1, sim);
                    triplets.set_entry(index, 2, out);
                }
            }
        }
        return triplets;
    }

    /**
     * Should do the same as np.argsort()
     * @private
     * @param {Array} A 
     */
    __argsort(A) {
        return A
            .map((d, i) => {return {d: d, i: i};})
            .sort((a, b) => a.d - b.d)
            .map((d) => d.i);
    }

    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are in the {@link rejects}.
     * @private
     * @param {*} n_samples 
     * @param {*} max_int 
     * @param {*} rejects 
     */
    _rejection_sample(n_samples, max_int, rejects) {
        const randomizer = this._randomizer;
        const interval = linspace(0, max_int - 1).filter(d => rejects.indexOf(d) < 0);
        return randomizer.choice(interval, Math.min(n_samples, interval.length - 2));
    }

    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - scaling factor for the distances.
     */
    _find_weights(triplets, P, nbrs, outlier_distances, sig) {
        const n_triplets = triplets.shape[0];
        const weights = new Float64Array(n_triplets);
        for (let t = 0; t < n_triplets; ++t) {
            const i = triplets.entry(t, 0);
            const sim = nbrs.row(i).indexOf(triplets.entry(t, 1));
            const p_sim = P.entry(i, sim);
            let p_out = Math.exp(-(outlier_distances[t] ** 2 / (sig[i] * sig[triplets.entry(t, 2)])));
            if (p_out < 1e-20) p_out = 1e-20;
            weights[t] = p_sim / p_out;
        }
        return weights;
    }

    /**
     * Sample uniformly ranom triplets
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {Number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    _sample_random_triplets(X, n_random, sig) {
        const metric = this._metric;
        const randomizer = this._randomizer;
        const N = X.shape[0];
        const random_triplets = new Matrix(N * n_random, 3);
        const random_weights = new Float64Array(N * n_random);
        for (let i = 0; i < N; ++i) {
            const n_i = i * n_random;
            const indices = [...linspace(0, i - 1), ...linspace(i + 1, N - 1)];
            for (let j = 0; j < n_random; ++j) {
                let [sim, out] = randomizer.choice(indices, 2);
                let p_sim = Math.exp(-((metric(X.row(i), X.row(sim)) ** 2) / (sig[i] * sig[sim])));
                if (p_sim < 1e-20) p_sim = 1e-20;
                let p_out = Math.exp(-((metric(X.row(i), X.row(out)) ** 2) / (sig[i] * sig[out]))); 
                if (p_out < 1e-20) p_out = 1e-20;

                if (p_sim < p_out) {
                    [sim, out] = [out, sim];
                    [p_sim, p_out] = [p_out, p_sim];
                }
                const index = n_i + j;
                random_triplets.set_entry(index, 0, i);
                random_triplets.set_entry(index, 1, sim);
                random_triplets.set_entry(index, 2, out);
                random_weights[index] = p_sim / p_out;
            }
        }
        return {
            "random_triplets": random_triplets,
            "random_weights": random_weights,
        }
    }

    /**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */
    _grad(Y) {
        const n_inliers = this.n_inliers;
        const n_outliers = this.n_outliers;
        const triplets = this.triplets;
        const weights = this.weights;
        const [N, dim] = Y.shape;
        const n_triplets = triplets.shape[0];
        const grad = new Matrix(N, dim, 0);
        let y_ij = new Array(dim).fill(0);
        let y_ik = new Array(dim).fill(0);
        let d_ij = 1;
        let d_ik = 1;
        let n_viol = 0;
        let loss = 0;
        const n_knn_triplets = N * n_inliers * n_outliers;

        for (let t = 0; t < n_triplets; ++t) {
            const [i, j, k] = triplets.row(t);
            // update y_ij, y_ik, d_ij, d_ik
            if (t % n_outliers == 0 || t >= n_knn_triplets) {
                d_ij = 1;
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_jd = Y.entry(j, d);
                    const Y_kd = Y.entry(k, d);
                    y_ij[d] = Y_id - Y_jd;
                    y_ik[d] = Y_id - Y_kd;
                    d_ij += (y_ij[d] ** 2);
                    d_ik += (y_ik[d] ** 2);
                }
            // update y_ik and d_ik only
            } else {
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_kd = Y.entry(k, d);
                    y_ik[d] = Y_id - Y_kd;
                    d_ik += (y_ik[d] ** 2);
                }
            }

            if (d_ij > d_ik) ++n_viol;
            loss += weights[t] / (1 + d_ik / d_ij);
            const w = (weights[t] / (d_ij + d_ik)) ** 2;
            for (let d = 0; d < dim; ++d) {
                const gs = y_ij[d] * d_ik * w;
                const go = y_ik[d] * d_ij * w;
                grad.set_entry(i, d, grad.entry(i, d) + gs - go);
                grad.set_entry(j, d, grad.entry(j, d) - gs);
                grad.set_entry(k, d, grad.entry(k, d) + go);
            }
        }
        return {
            "grad": grad,
            "loss": loss,
            "n_viol": n_viol,
        };
    }

    /**
     * 
     * @param {Number} max_iteration 
     */
    transform(max_iteration = 400) {
        this.check_init();
        for (let iter = 0; iter < max_iteration; ++iter) {
            this._next(iter);
        }
        return this.projection;
    }

    /**
     * @yields {Matrix}
     * @returns {Matrix}
     */
    * generator() {
        this.check_init();
        for (let iter = 0; iter < 800; ++iter) {
            this._next(iter);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * Does the iteration step.
     * @private
     * @param {Number} iter 
     */
    _next(iter) {
        const gamma = iter > 150 ? .5 : .3;
        const old_C = this.C;
        const vel = this.vel;
        const Y = this.Y.add(vel.mult(gamma));
        const {grad, loss, n_viol} = this._grad(Y);
        this.C = loss;
        this.Y = this._update_embedding(Y, iter, grad);
        this.lr *= (old_C > loss + this.tol)  ? 1.01 : .9;
        return this.Y;
    }

    /**
     * Updates the embedding.
     * @private
     * @param {Matrix} Y 
     * @param {Number} iter 
     * @param {Matrix} grad 
     */
    _update_embedding(Y, iter, grad) {
        const [N, dim] = Y.shape;
        const gamma = iter > 150 ? .9 : .5; // moment parameter
        const min_gain = .01;
        const gain = this.gain;
        const vel = this.vel;
        const lr = this.lr;
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const new_gain = (Math.sign(vel.entry(i, d)) != Math.sign(grad.entry(i, d))) ? gain.entry(i, d) + .2 : Math.max(gain.entry(i, d) * .8, min_gain);
                gain.set_entry(i, d, new_gain);
                vel.set_entry(i, d, gamma * vel.entry(i, d) - lr * gain.entry(i, d) * grad.entry(i, d));
                Y.set_entry(i, d, Y.entry(i, d) + vel.entry(i, d));
            }
        }
        return Y;
    }
}

/**
 * @class
 * @alias Hierarchical_Clustering
 */
class Hierarchical_Clustering {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias Hierarchical_Clustering
     * @todo needs restructuring.
     * @param {Matrix} - Data or distance matrix if metric is 'precomputed'
     * @param {("single"|"complete"|"average")} [linkage = "complete"]
     * @param {Function|"precomputed"} [metric = euclidean]
     * @returns {Hierarchical_Clustering}
     */
    constructor(matrix, linkage = "complete", metric = euclidean) {
        this._id = 0;
        this._matrix = matrix instanceof Matrix ? matrix : Matrix.from(matrix);
        this._metric = metric;
        this._linkage = linkage;
        if (metric === "precomputed" && this._matrix.shape[0] !== this._matrix.shape[1]) {
            throw new Error("If metric is 'precomputed', then matrix has to be square!");
        }
        this.init();
        this.root = this.do();
        return this;
    }

    /**
     *
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(value, type = "distance") {
        let clusters = [];
        let accessor;
        switch (type) {
            case "distance":
                accessor = (d) => d.dist;
                break;
            case "depth":
                accessor = (d) => d.depth;
                break;
            default:
                throw new Error("invalid type");
        }
        this._traverse(this.root, accessor, value, clusters);
        return clusters;
    }

    /**
     * @private
     * @param {} node
     * @param {*} f
     * @param {*} value
     * @param {*} result
     */
    _traverse(node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        } else {
            this._traverse(node.left, f, value, result);
            this._traverse(node.right, f, value, result);
        }
    }

    /**
     * computes the tree.
     */
    init() {
        const metric = this._metric;
        const A = this._matrix;
        const n = (this._n = A.shape[0]);
        const d_min = (this._d_min = new Float64Array(n));
        let distance_matrix;
        if (metric !== "precomputed") {
            distance_matrix = new Matrix(n, n, 0); //new Array(n);
            for (let i = 0; i < n; ++i) {
                d_min[i] = 0;
                //distance_matrix[i] = new Float64Array(n);
                for (let j = 0; j < n; ++j) {
                    distance_matrix.set_entry(i, j, i === j ? Infinity : metric(A.row(i), A.row(j)));
                    if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        } else {
            distance_matrix = this._matrix.clone();
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    if (i === j) {
                        distance_matrix.set_entry(i, j, Infinity);
                    } else if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        this._distance_matrix = distance_matrix;
        const clusters = (this._clusters = new Array(n));
        const c_size = (this._c_size = new Uint16Array(n));
        for (let i = 0; i < n; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        return this;
    }

    /**
     * computes the tree.
     */
    do() {
        const n = this._n;
        const d_min = this._d_min;
        const D = this._distance_matrix;
        const clusters = this._clusters;
        const c_size = this._c_size;
        const linkage = this._linkage;
        let root = null;
        for (let p = 0, p_max = n - 1; p < p_max; ++p) {
            let c1 = 0;
            for (let i = 0; i < n; ++i) {
                let D_i_min = D.entry(i, d_min[i]);
                for (let j = i + 1; j < n; ++j) {
                    if (D_i_min > D.entry(i, j)) {
                        d_min[i] = j;
                        D_i_min = D.entry(i, d_min[i]);
                    }
                }
            }
            for (let i = 0; i < n; ++i) {
                if (D.entry(i, d_min[i]) < D.entry(c1, d_min[c1])) {
                    c1 = i;
                }
            }
            let c2 = d_min[c1];
            let c1_cluster = clusters[c1][0];
            let c2_cluster = clusters[c2][0];
            let c1_cluster_indices = c1_cluster.isLeaf ? [c1_cluster.index] : c1_cluster.index;
            let c2_cluster_indices = c2_cluster.isLeaf ? [c2_cluster.index] : c2_cluster.index;
            let indices = c1_cluster_indices.concat(c2_cluster_indices);
            let new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D.entry(c1, c2), null, indices);
            c1_cluster.parent = new_cluster;
            c2_cluster.parent = new_cluster;
            clusters[c1].unshift(new_cluster);
            c_size[c1] += c_size[c2];
            for (let j = 0; j < n; ++j) {
                const D_c1_j = D.entry(c1, j);
                const D_c2_j = D.entry(c2, j);
                let value;
                switch (linkage) {
                    case "single":
                        value = Math.min(D_c1_j, D_c2_j);
                        break;
                    case "complete":
                        value = Math.max(D_c1_j, D_c2_j);
                        break;
                    case "average":
                        value = (c_size[c1] * D_c1_j + c_size[c2] * D_c2_j) / (c_size[c1] + c_size[j]);
                        break;
                }
                D.set_entry(j, c1, value);
                D.set_entry(c1, j, value);
            }

            D.set_entry(c1, c1, Infinity);
            for (let i = 0; i < n; ++i) {
                D.set_entry(i, c2, Infinity);
                D.set_entry(c2, i, Infinity);
            }

            /* for (let j = 0; j < n; ++j) {
                if (d_min[j] === c2) {
                    d_min[j] = c1;
                }
                if (D.entry(c1, j) < D.entry(c1, d_min[c1])) {
                    d_min[c1] = j;
                }
            } */
            root = new_cluster;
        }
        return root;
    }
}

class Cluster {
    constructor(id, left, right, dist, centroid, index, size, depth) {
        this.id = id;
        this.left = left;
        this.right = right;
        this.dist = dist;
        this.index = index;
        this.size = size ?? left.size + right.size;
        this.depth = depth ?? 1 + Math.max(left.depth, right.depth);
        this.centroid = centroid ?? this._calculate_centroid(left, right);
        this.parent = null;
        return this;
    }

    _calculate_centroid(left, right) {
        const l_size = left.size;
        const r_size = right.size;
        const l_centroid = left.centroid;
        const r_centroid = right.centroid;
        const size = this.size;
        const n = left.centroid.length;
        const new_centroid = new Float64Array(n);
        for (let i = 0; i < n; ++i) {
            new_centroid[i] = (l_size * l_centroid[i] + r_size * r_centroid[i]) / size;
        }
        return new_centroid;
    }

    get isLeaf() {
        return this.depth === 0;
    }

    leaves() {
        if (this.isLeaf) return [this];
        const left = this.left;
        const right = this.right;
        return (left.isLeaf ? [left] : left.leaves()).concat(right.isLeaf ? [right] : right.leaves());
    }

    descendants() {
        if (this.isLeaf) return [this];
        const left_descendants = this.left.descendants();
        const right_descendants = this.right.descendants();
        return left_descendants.concat(right_descendants).concat([this]);
    }
}

/**
 * @class
 * @alias KMeans
 */
class KMeans {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias KMeans
     * @todo needs restructuring. 
     * @param {Matrix} matrix 
     * @param {Numbers} K 
     * @param {Function} [metric = euclidean] 
     * @param {Number} [seed = 1987]
     * @param {Boolean} [init = true]
     * @returns {KMeans}
     */
    constructor(matrix, K, metric = euclidean, seed=1987, init = true) {
        this._metric = metric;
        this._matrix = matrix;
        this._K = K;
        const [N, D] = matrix.shape;
        this._N = N;
        this._D = D;
        if (K > N) K = N;
        this._randomizer = new Randomizer(seed);
        this._clusters = new Array(N).fill(undefined);
        this._cluster_centroids = this._get_random_centroids(K);
        if (init) this.init(K, this._cluster_centroids);
        return this;
    }

    /**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */
    get_clusters() {
        const K = this._K;
        const clusters = this._clusters;
        const result = new Array(K).fill().map(() => new Array());
        clusters.forEach((c, i) => result[c].push(i));
        return result;
    }

    /**
     * @private
     * @param {Array} points 
     * @param {Array} candidates 
     */
    _furthest_point(points, candidates) {
        const A = this._matrix;
        const metric = this._metric;
        let i = points.length;
        let H = Heap.heapify(
            candidates, 
            (d) => {
                const Ad = A.row(d);
                let sum = 0;
                for (let j = 0; j < i; ++j) {
                    sum += metric(Ad, points[j]);
                }
                return sum;
            }, 
            "max"
        );
        return H.pop().element;
    }

    _get_random_centroids(K) {
        const N = this._N;
        const randomizer = this._randomizer;
        const A = this._matrix;
        const cluster_centroids = new Array(K).fill();
        const indices = linspace(0, N - 1);
        const random_point = randomizer.random_int % (N - 1);
        cluster_centroids[0] = A.row(random_point);
        const init_points = [random_point];
        const sample_size = Math.floor((N - K) / K);// / K
        for (let i = 1; i < K; ++i) {
            // sampling + kmeans++ improvement?
            const sample = randomizer.choice(indices.filter(d => init_points.indexOf(d) == -1), sample_size);
            const furthest_point = this._furthest_point(cluster_centroids.slice(0, i), sample);
            init_points.push(furthest_point);
            cluster_centroids[i] = A.row(furthest_point);
        }
        return cluster_centroids;
    }

    _iteration(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const metric = this._metric;
        const clusters = this._clusters;
        let clusters_changed = false;
        // find nearest cluster centroid.
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            let min_dist = Infinity;
            let min_cluster = null;
            for (let j = 0; j < K; ++j) {
                let d = metric(cluster_centroids[j], Ai);
                if (d < min_dist) {
                    min_dist = d;
                    min_cluster = j; 
                }
            }
            if (clusters[i] !== min_cluster) {
                clusters_changed = true;
            }
            clusters[i] = min_cluster;
        }
        // update cluster centroid
        // reset cluster centroids to 0
        for (let i = 0; i < K; ++i) {
            const centroid = cluster_centroids[i];
            for (let j = 0; j < D; ++j) {
                centroid[j] = 0;
            }
        }
        // compute centroid
        this._compute_centroid(cluster_centroids);

        return {   
            "clusters_changed": clusters_changed,
            "cluster_centroids": cluster_centroids
        };
    }

    _compute_centroid(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const clusters = this._clusters;
        const cluster_counter = new Array(K).fill(0);

        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            const ci = clusters[i];
            cluster_counter[ci]++;
            const centroid = cluster_centroids[ci];
            for (let j = 0; j < D; ++j) {
                centroid[j] += Ai[j];
            }
        }
        for (let i = 0; i < K; ++i) {
            const n = cluster_counter[i];
            cluster_centroids[i] = cluster_centroids[i].map(c => c / n);
        }
        
    }

    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K, cluster_centroids) {
        if (!K) K = this._K;
        if (!cluster_centroids) cluster_centroids = this._get_random_centroids(K);
        let clusters_changed = false;
        do {
            const iteration_result = this._iteration(cluster_centroids);
            cluster_centroids = iteration_result.cluster_centroids;
            clusters_changed = iteration_result.clusters_changed;
        } while (clusters_changed)
    }
    
}

/**
 * @class
 * @alias KMedoids
 */
class KMedoids {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias KMedoids
     * @todo needs restructuring. 
     * @param {Matrix} matrix - data matrix
     * @param {Numbers} K - number of clusters
     * @param {number} [max_iter=null] - maximum number of iterations. Default is 10 * Math.log10(N)
     * @param {Function} [metric = euclidean] - metric defining the dissimilarity 
     * @param {Number} [seed = 1212] - seed value for random number generator
     * @returns {KMedoids}
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(matrix, K, max_iter=null, metric = euclidean, seed=1212) {
        this._metric = metric;
        this._matrix = matrix;
        this._A = this._matrix.to2dArray;
        this._K = K;
        const [N, D] = matrix.shape;
        this._N = N;
        this._D = D;
        this._max_iter = max_iter || 10 * Math.log10(N); 
        this._distance_matrix = new Matrix(N, N, "zeros");
        /* for (let i = 1; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dist = metric(this._A[i], this._A[j]);
                this._distance_matrix.set_entry(i, j, dist);
                this._distance_matrix.set_entry(j, i, dist)
            }
        } */
        if (K > N) K = N;
        this._randomizer = new Randomizer(seed);
        this._clusters = new Array(N).fill(undefined);
        this._cluster_medoids = this._get_random_medoids(K);
        //if (init) this.init(K, this._cluster_medoids);
        this._is_initialized = false;
        return this;
    }

    /**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */
    get_clusters() {
        const K = this._K;
        const A = this._A;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        const result = new Array(K).fill().map(() => new Array());
        A.forEach((x_j, j) => {
            result[this._nearest_medoid(x_j, j).index_nearest].push(j);
        });
        result.medoids = this._cluster_medoids;
        return result;
    }

    async* generator() {
        const max_iter = this._max_iter;
        yield this.get_clusters();
        let finish = false;
        let i = 0;
        do {
            finish = this._iteration();
            yield this.get_clusters();
        } while (!finish && ++i < max_iter)
    }

    /**
     * Algorithm 1. FastPAM1: Improved SWAP algorithm
     */
    /* _iteration_1() {
        const A = this._A;
        const N = this._N;
        const K = this._K;
        const medoids = this._cluster_medoids;
        let DeltaTD = 0;
        let m0 = null;
        let x0 = null;
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const nearest_medoid = this._nearest_medoid(x_j, j);
                const d_j = nearest_medoid.distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    // disance to new medoid
                    const d_oj = this._get_distance(o, j, x_o, x_j);
                    const {
                        "index_nearest": n,
                        "distance_nearest": d_n,
                        "distance_second": d_s,
                    } = this._nearest_medoid(x_o, o); 
                    this._clusters[o] = n; // cached values
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change
                    if (d_oj < d_n) { // reassignment check
                        deltaTD.forEach((d_i, i) => {
                            if (n !== i) {
                                deltaTD[i] = d_i + d_oj - d_n; // update loss change
                            }
                        });
                    }
                });
                // choose best medoid i;
                const i = deltaTD
                    .map((d, i) => [d, i])
                    .sort((d1, d2) => d1[0] - d2[0])[0][1];
                const deltaTD_i = deltaTD[i];
                // store
                if (deltaTD_i < DeltaTD) {
                    DeltaTD = deltaTD_i;
                    m0 = i;
                    x0 = j;
                }
            }
        });

        if (DeltaTD >= 0) {
            return true // break loop if DeltaTD >= 0
        }
        // swap roles of medoid m and non-medoid x;
        medoids[m0] = x0;
        this._cluster_medoids = medoids;
        return false
    } */

    /** Algorithm 2. FastPAM2: SWAP with multiple candidates
     * 
     */
    _iteration() {
        const A = this._A;
        const K = this._K;
        const medoids = this._cluster_medoids;
        const cache = A.map((x_o, o) => this._nearest_medoid(x_o, o));
        // empty best candidates array
        const DeltaTD = new Array(K).fill(0);
        const xs = new Array(K).fill(null);
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const d_j = cache[j].distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    if (j === o) return;
                    const d_oj = this._get_distance(o, j, x_o, x_j); // distance to new medoid
                    const {"index_nearest": n, "distance_nearest": d_n, "distance_second": d_s} = cache[o]; // cached
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change for x_o
                    // Reassignment check
                    if (d_oj < d_n) { 
                        // update loss change
                        for (let i = 0; i < K; ++i) {
                            if (i !== n) deltaTD[i] += d_oj - d_n;
                        }
                    }
                });
                // remember best swap for i;
                deltaTD
                    .map((d, i) => [d, i])
                    .filter(([d, i]) => d < DeltaTD[i])
                    .forEach(([d, i]) => {
                        if (d < DeltaTD[i]) {
                            DeltaTD[i] = d;
                            xs[i] = j;
                        }
                    });
            }
        });
        // stop if no improvements were found
        if (min(DeltaTD) >= 0) return true; 

        // execute all improvements
        while (min(DeltaTD) < 0) {
            // swap roles of medoid m_i and non_medoid xs_i
            const i = DeltaTD
                .map((d, i) => [d, i])
                .sort(([a], [b]) => a - b)[0][1];
            if (medoids.filter(m => m == xs[i]).length == 0) {
                medoids[i] = xs[i];
            }
            // disable the swap just performed
            DeltaTD[i] = 0; 
            // recompute TD for remaining swap candidates
            DeltaTD
                .map((d_j, j) => [d_j, j])
                .filter(([d_j]) => d_j < 0)
                .forEach(([_, j]) => {
                    const x_j = A[j];
                    let sum = 0;
                    A.forEach((x_o, o) => {
                        if (medoids.findIndex(m => m != j && m == o) >= 0) return;
                        if (i == j) return;
                        if (cache[o].index_nearest === medoids[j])
                            sum += (Math.min(this._get_distance(o, j, x_o, x_j), cache[o].distance_second) - cache[o].distance_nearest); 
                        else {
                            sum += (Math.min(this._get_distance(o, j, x_o, x_j) - cache[o].distance_nearest, 0));
                        }
                    });
                    DeltaTD[j] = sum;
                });
        }
        this._cluster_medoids = medoids;
        return false;
    }

    _get_distance(i, j, x_i=null, x_j=null) {
        if (i === j) return 0;
        const D = this._distance_matrix;
        const A = this._A;
        const metric = this._metric;
        let d_ij = D.entry(i, j);
        if (d_ij === 0) {
            d_ij = metric(x_i || A[i], x_j || A[j]);
            D.set_entry(i, j, d_ij);
            D.set_entry(j, i, d_ij);
        }
        return d_ij;
    }

    _nearest_medoid(x_j, j) {
        const medoids = this._cluster_medoids;
        const A = this._A;
        const [nearest, second] = medoids
            .map((m, i) => {
                const x_m = A[m]; 
                return [this._get_distance(j, m, x_j, x_m), i];
            })
            .sort((m1, m2) => m1[0] - m2[0]);
        
        return { 
            "distance_nearest": nearest[0], 
            "index_nearest": nearest[1],
            "distance_second": second[0],
            "index_second": second[1],
        };
    }

    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K, cluster_medoids) {
        if (!K) K = this._K;
        if (!cluster_medoids) cluster_medoids = this._get_random_medoids(K);
        const max_iter = this._max_iter;
        let finish = false;
        let i = 0;
        do {
            finish = this._iteration();
        } while (!finish && ++i < max_iter)
        return this;
    }

    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @param {number} K - number of clusters
     * 
     */
    _get_random_medoids(K) {
        const N = this._N;
        const A = this._A;
        const indices = linspace(0, N - 1);
        const randomizer = this._randomizer;
        const n = Math.min(N, 10 + Math.ceil(Math.sqrt(N)));
        const TD = new Array(n).fill(Infinity);
        const medoids = [];
        // first medoid
        let TD0 = Infinity;
        let S = randomizer.choice(indices, n);
        for (let j = 0; j < n; ++j) {
            const S_j = S[j];
            const x_j = A[S_j];
            for (let o = 0; o < n; ++o) {
                if (o === j) continue;
                const x_o = A[S[o]];
                TD[j] += this._get_distance(j, o, x_j, x_o);
            }
            if (TD[j] < TD0) {
                TD0 = TD[j]; // smallest distance sum
                medoids.push(S_j);
            }
        }
        // other medoids
        for (let i = 1; i < K; ++i) {
            let DeltaTD = Infinity;
            S = randomizer.choice(indices.filter(index => medoids.findIndex(d => d === index) < 0), n);
            for (let j = 0; j < n; ++j) {
                let deltaTD = 0;
                const S_j = S[j];
                const x_j = A[S_j];
                for (let o = 0; o < n; ++o) {
                    if (o === j) continue;
                    const S_o = S[o];
                    const x_o = A[S_o];
                    let delta = this._get_distance(S_j, S_o, x_j, x_o) - min(medoids.map(m => this._get_distance(S_o, m, x_o)));
                    if (delta < 0) {
                        deltaTD = deltaTD + delta;
                    }
                }
                // best reduction
                if (deltaTD < DeltaTD) {
                    DeltaTD = deltaTD;
                    medoids.push(S_j);
                }
            }
            TD0 += DeltaTD;
        }
        return medoids.slice(0, K);
    }
    
}

/**
 * @class
 * @alias OPTICS
 */
class OPTICS {
    /**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     * @constructor
     * @memberof module:clustering
     * @alias OPTICS
     * @todo needs restructuring. 
     * @param {Matrix} matrix - the data.
     * @param {Number} epsilon - the minimum distance which defines whether a point is a neighbor or not.
     * @param {Number} min_points - the minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
     * @param {Function} [metric = euclidean] - the distance metric which defines the distance between two points of the {@link matrix}.
     * @returns {OPTICS}
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
    constructor(matrix, epsilon, min_points, metric = euclidean) {
        this._matrix = matrix;
        this._epsilon = epsilon;
        this._min_points = min_points;
        this._metric = metric;

        this._ordered_list = [];
        this._clusters = [];
        this._DB = new Array(matrix.shape[0]).fill();
        this.init();
        return this;
    }

    /**
     * Computes the clustering.
     */
    init() {
        const ordered_list = this._ordered_list;
        const matrix = this._matrix;
        const N = matrix.shape[0];
        const DB = this._DB;
        const clusters = this._clusters;
        let cluster_index = this._cluster_index = 0;

        for (let i = 0; i < N; ++i) {
            DB[i] = {
                "element": matrix.row(i),
                "index": i,
                "reachability_distance": undefined,
                "processed": false,
            };
        }
        for (const p of DB) {
            if (p.processed) continue;
            p.neighbors = this._get_neighbors(p);
            p.processed = true;
            clusters.push([p.index]);
            cluster_index = clusters.length - 1;
            ordered_list.push(p);
            if (this._core_distance(p) != undefined) {
                const seeds = new Heap(null, d => d.reachability_distance, "min");
                this._update(p, seeds);
                this._expand_cluster(seeds, clusters[cluster_index]);
            }
        }
        return this;
    }

    /**
     * 
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Array} An array consisting of the {@link epsilon}-neighborhood of {@link p}.
     */
    _get_neighbors(p) {
        if ("neighbors" in p) return p.neighbors;
        const DB = this._DB;
        const metric = this._metric;
        const epsilon = this._epsilon;
        const neighbors = [];
        for (const q of DB) {
            if (q.index == p.index) continue;
            if (metric(p.element, q.element) < epsilon) {
                neighbors.push(q);
            }
        }
        return neighbors;
    }

    /**
     * 
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Number} The distance to the {@link min_points}-th nearest point of {@link p}, or undefined if the {@link epsilon}-neighborhood has fewer elements than {@link min_points}.
     */
    _core_distance(p) {
        const min_points = this._min_points;
        const metric = this._metric;
        if (p.neighbors && p.neighbors.length <= min_points) {
            return undefined;
        }
        return metric(p.element, p.neighbors[min_points].element);
    }

    /**
     * Updates the reachability distance of the points.
     * @private
     * @param {Object} p 
     * @param {Heap} seeds 
     */
    _update(p, seeds) {
        const metric = this._metric;
        const core_distance = this._core_distance(p);
        const neighbors = this._get_neighbors(p);//p.neighbors;
        for (const q of neighbors) {
            if (q.processed) continue;
            const new_reachability_distance = Math.max(core_distance, metric(p.element, q.element));
            //if (q.reachability_distance == undefined) { // q is not in seeds
            if (seeds.raw_data().findIndex(d => d.element == q) < 0) {
                q.reachability_distance = new_reachability_distance;
                seeds.push(q);
            } else { // q is in seeds
                if (new_reachability_distance < q.reachability_distance) {
                    q.reachability_distance = new_reachability_distance;
                    seeds = Heap.heapify(seeds.data(), d => d.reachability_distance, "min"); // seeds change key =/
                }
            }
        }
    }

    /**
     * Expands the {@link cluster} with points in {@link seeds}.
     * @private
     * @param {Heap} seeds 
     * @param {Array} cluster 
     */
    _expand_cluster(seeds, cluster) {
        const ordered_list = this._ordered_list;
        while (!seeds.empty) {
            const q = seeds.pop().element;
            q.neighbors = this._get_neighbors(q);
            q.processed = true;
            cluster.push(q.index);
            ordered_list.push(q);
            if (this._core_distance(q) != undefined) {
                this._update(q, seeds);
                this._expand_cluster(seeds, cluster);
            }
        }
    }

    /**
     * Returns an array of clusters.
     * @returns {Array<Array>} Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters() {
        const clusters = [];
        const outliers = [];
        const min_points = this._min_points;
        for (const cluster of this._clusters) {
            if (cluster.length < min_points) {
                outliers.push(...cluster);
            } else {
                clusters.push(cluster);
            }
        }
        clusters.push(outliers);
        return clusters;
    }

    /**
     * @returns {Array} Returns an array, where the ith entry defines the cluster affirmation of the ith point of {@link matrix}. (-1 stands for outlier)
     */
    get_cluster_affirmation() {
        const N = this._matrix.shape[0];
        const result = new Array(N).fill();
        const clusters = this.get_clusters();
        for (let i = 0, n = clusters.length; i < n; ++i) {
            const cluster = clusters[i];
            for (const index of cluster) {
                result[index] = (i < n - 1) ? i : -1;
            }
        }
        return result;
    }
}

/**
 * @class
 * @alias LSP
 * @extends DR
 */
class LSP extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {number} [k = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {number} [control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     */
    constructor(X, k, control_points, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k", "control_points"];
        this.parameter("k", Math.min(k ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        this.parameter("control_points", Math.min(control_points ?? Math.ceil(Math.sqrt(this._N)), this._N - 1));
        this._is_initialized = false;
        return this;
    }

    /**
     * 
     * @param {DR} DR - method used for position control points.
     * @param {DR_parameters} DR_parameters - array containing parameters for the DR method which projects the control points
     * @returns {LSP} 
     */
    init(DR=MDS, DR_parameters=[], KNN=BallTree) {
        if (this._is_initialized) return this;
        const X = this.X;
        const N = this._N;
        const K = this.parameter("k");
        const d = this._d;
        const metric = this._metric;
        const nc = this.parameter("control_points");
        const control_points = new KMedoids(X, nc, null, metric).get_clusters().medoids;
        const C = new Matrix(nc, N, "zeros");
        control_points.forEach((c_i, i) => {
            C.set_entry(i, c_i, 1);
        });
        const Y_C = new DR(Matrix.from(control_points.map(c_i => X.row(c_i))), ...DR_parameters, d).transform();
        
        const XA = X.to2dArray;
        const knn = new KNN(XA, metric);
        const L = new Matrix(N, N, "I");
        const alpha = -1/K;
        XA.forEach((x_i, i) => {
            for (const {"index": j} of knn.search(x_i, K).iterate()) {
                if (i === j) continue;
                L.set_entry(i, j, alpha);
            }
        });
        const A = L.concat(C, "vertical");

        const z = new Matrix(N, d, "zeros");
        const b = z.concat(Y_C, "vertical");
        
        this._A = A;
        this._b = b;
        this._is_initialized = true;
        return this;
    }


    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        const A = this._A;
        const AT = A.T;
        const b = this._b;
        const ATA = AT.dot(A);
        const ATb = AT.dot(b);
        this.Y = Matrix.solve_CG(ATA, ATb, this._randomizer);
        return this.projection;
    }
}

/**
 * @class
 * @alias TopoMap
 * @memberof module:dimensionality_reduction
 * @extends DR
 */
class TopoMap extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TopoMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, d = 2, metric = euclidean, seed = 1212) {
        super(X, d, metric, seed);
        super.parameter_list = [];
        [this._N, this._D] = this.X.shape;
        this._distance_matrix = new Matrix(this._N, this._N, 0);
        return this;
    }

    /**
     * @private
     */
    __lazy_distance_matrix(i, j, metric) {
        const D = this._distance_matrix;
        const X = this.X;
        const D_ij = D.entry(i, j);
        if (D_ij === 0) {
            let dist = metric(X.row(i), X.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
            return dist;
        }
        return D_ij;
    }

    /**
     * Computes the minimum spanning tree, using a given metric
     * @private
     * @param {Function} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    _make_minimum_spanning_tree(metric = euclidean) {
        const N = this._N;
        const X = [...this.X];

        let disjoint_set = new DisjointSet(X);
        const F = [];
        let E = [];
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                E.push([i, j, this.__lazy_distance_matrix(i, j, metric)]);
            }
        }
        E = E.sort((a, b) => a[2] - b[2]);

        for (const [u, v, w] of E) {
            const set_u = disjoint_set.find(X[u]);
            const set_v = disjoint_set.find(X[v]);
            if (set_u !== set_v) {
                F.push([u, v, w]);
                disjoint_set.union(set_u, set_v);
            }
        }

        return F.sort((a, b) => a[2] - b[2]);
    }

    /**
     * initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init() {
        this.Y = new Matrix(this._N, this._d, 0);
        this._Emst = this._make_minimum_spanning_tree(this._metric);
        this._is_initialized = true;
        return this;
    }

    /**
     * Returns true if Point C is left of line AB.
     * @private
     * @param {Array} PointA - Point A of line AB
     * @param {Array} PointB - Point B of line AB
     * @param {Array} PointC - Point C
     * @returns {Boolean}
     */
    __hull_cross([ax, ay], [bx, by], [sx, sy]) {
        return (bx - ax) * (sy - ay) - (by - ay) * (sx - ax) <= 0;
    }

    /**
     * Computes the convex hull of the set of Points S
     * @private
     * @param {Array} S - Set of Points.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     * @returns {Array} convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     */
    __hull(S) {
        const points = S.sort(([x1, y1], [x2, y2]) => y1 - y2 || x1 - x2);
        const N = points.length;
        if (N <= 2) return points;

        const lower = [];
        for (let i = 0; i < N; ++i) {
            while (lower.length >= 2 && this.__hull_cross(lower[lower.length - 2], lower[lower.length - 1], points[i])) {
                lower.pop();
            }
            lower.push(points[i]);
        }
        const upper = [];
        for (let i = N - 1; i >= 0; --i) {
            while (upper.length >= 2 && this.__hull_cross(upper[upper.length - 2], upper[upper.length - 1], points[i])) {
                upper.pop();
            }
            upper.push(points[i]);
        }
        upper.pop();
        lower.pop();
        return lower.concat(upper);
    }

    /**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     * @private
     * @param {Array} PointA
     * @param {Array} PointB
     * @return {Object} Object containing the sinus- and cosinus-values for a rotation.
     */
    __findAngle([p1x, p1y], [p2x, p2y]) {
        const n = euclidean([p1x, p1y], [p2x, p2y]);
        if (n === 0)
            return {
                sin: 0,
                cos: 1,
            };
        const vec = [(p2x - p1x) / n, (p2y - p1y) / n];
        const cos = vec[0];
        let sin = Math.sqrt(1 - cos * cos);
        sin = vec[1] >= 0 ? -sin : sin;
        return {
            sin: sin,
            cos: cos,
        };
    }

    /**
     * @private
     * @param {Array} hull
     * @param {Array} p
     * @param {Bool} topEdge
     */
    __align_hull(hull, p, topEdge) {
        let v = -1;
        let d2;
        for (let i = 0; i < hull.length; ++i) {
            const d = euclidean(hull[i], p);
            if (v === -1) {
                d2 = d;
                v = i;
            } else {
                if (d2 > d) {
                    d2 = d;
                    v = i;
                }
            }
        }

        let v1;
        let v2;
        if (topEdge) {
            v1 = hull[v];
            v2 = hull[(v + 1) % hull.length];
        } else {
            if (v == 0) v = hull.length - 1;
            v1 = hull[v];
            v2 = hull[(v - 1) % hull.length];
        }

        const transformation = {
            tx: -hull[v][0],
            ty: -hull[v][1],
        };

        if (hull.length >= 2) {
            const { sin, cos } = this.__findAngle(v1, v2);
            transformation.sin = sin;
            transformation.cos = cos;
        } else {
            transformation.sin = 0;
            transformation.cos = 1;
        }

        return transformation;
    }

    /**
     * @private
     * @param {Array} Point - The point which should get transformed.
     * @param {Object} Transformation - contains the values for translation and rotation.
     */
    __transform([px, py], { tx, ty, sin, cos }) {
        let x = px + tx;
        let y = py + ty;
        let xx = x * cos - y * sin;
        let yy = x * sin + y * cos;
        return [xx, yy];
    }

    /**
     * Calls {@link __transform} for each point in Set C
     * @private
     * @param {Array} C - Set of points.
     * @param {Object} t - Transform object.
     * @param {Number} yOffset - value to offset set C.
     */
    __transform_component(C, t, yOffset) {
        const N = C.length;
        for (let i = 0; i < N; ++i) {
            const c = C[i];
            const [cx, cy] = this.__transform(c, t);
            c[0] = cx;
            c[1] = cy + yOffset;
        }
    }

    /**
     * @private
     * @param {Array} u - point u
     * @param {Array} v - point v
     * @param {Number} w - edge weight w
     */
    __align_components(u, v, w) {
        const points_u = [...u.__disjoint_set.children];
        const points_v = [...v.__disjoint_set.children];

        const hull_u = this.__hull(points_u);
        const hull_v = this.__hull(points_v);

        const t_u = this.__align_hull(hull_u, u, false);
        const t_v = this.__align_hull(hull_v, v, true);

        this.__transform_component(points_u, t_u, 0);
        this.__transform_component(points_v, t_v, w);
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality 2.
     */
    transform() {
        if (!this._is_initialized) this.init();
        const Emst = this._Emst;
        const Y = [...this.Y];
        const components = new DisjointSet(
            Y.map((y, i) => {
                y.i = i;
                return y;
            })
        );

        for (const [u, v, w] of Emst) {
            const component_u = components.find(Y[u]);
            const component_v = components.find(Y[v]);
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, w);
            components.union(component_u, component_v);
        }
        return this.projection;
    }

    *generator() {
        if (!this._is_initialized) this.init();
        const Emst = this._Emst;
        const Y = [...this.Y];
        const components = new DisjointSet(
            Y.map((y, i) => {
                y.i = i;
                return y;
            })
        );

        for (const [u, v, w] of Emst) {
            const component_u = components.find(Y[u]);
            const component_v = components.find(Y[v]);
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, w);
            components.union(component_u, component_v);
            /* let ok = true
            Y.forEach(([x, y]) => ok = ok && !isNaN(x) && !isNaN(y))
            if (!ok) {
                console.log(...Y) 
                throw "error" 
            } */
            yield this.projection;
        }
        return this.projection;
    }
}

/**
 * @class
 * @alias SAMMON
 * @extends DR
 */
class SAMMON extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, magic=0.1, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["magic"];
        this.parameter("magic", magic);
        [ this._N, this._D ] = this.X.shape;
        return this;
    }

    /**
     * initializes SAMMON. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init(DR$1="random", distance_matrix=null) {
        const N = this._N;
        const d = this._d;

        if (DR$1 === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (DR$1 instanceof DR) {
            this.Y = DR$1.transform(this.X);
        }
        this.distance_matrix = distance_matrix || this.__distance_matrix(this.X);
        return this;
    }

    /**
     * @private
     * @param {Matrix} A
     * @returns {Matrix} 
     */
    __distance_matrix(A) {
        const metric = this._metric;
        const N = A.shape[0];
        const D = new Matrix(N, N);
        for (let i = 0; i < N; ++i) {
            const A_i = A.row(i);
            for (let j = i; j < N; ++j) {
                let distance = (i === j ? 0 : metric(A_i, A.row(j)));
                D.set_entry(i, j, distance);
                D.set_entry(j, i, distance);
            }
        }
        return D;                
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     */
    transform(max_iter=200) {
        if (!this._is_initialized) this.init();
        for (let j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    }

    * generator(max_iter=200) {
        if (!this._is_initialized) this.init();

        for (let j = 0; j < max_iter; ++j) {
            this._step();
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        const MAGIC = this.parameter("magic");
        const D = this.distance_matrix;
        const N = this._N;
        const d = this._d;
        const metric = this._metric;
        let Y = this.Y;
        
        let G = new Matrix(N, d, 0);

        let sum = new Float64Array(d);
        for (let i = 0; i < N; ++i) {
            let e1 = new Float64Array(d);
            let e2 = new Float64Array(d);
            const Yi = Y.row(i);
            for (let j = 0; j < N; ++j) {
                if (i === j) continue;
                const Yj = Y.row(j);
                const delta = new Float64Array(d);
                for (let k = 0; k < d; ++k) {
                    delta[k] = Yi[k] - Yj[k];
                }
                const dY = metric(Yi, Yj);
                const dX = D.entry(i, j);
                const dq = dX - dY;
                const dr = Math.max(dX * dY, 1e-2);
                for (let k = 0; k < d; ++k) {
                    e1[k] += delta[k] * dq / dr;
                    e2[k] += (dq - Math.pow(delta[k], 2) * (1 + dq / dY) / dY) / dr;
                }
            }
            for (let k = 0; k < d; ++k) {
                const val = Y.entry(i, k) + (MAGIC * e1[k] / Math.abs(e2[k]) || 0);
                G.set_entry(i, k, val);
                sum[k] += val;
            }
        }
        for (let k = 0; k < d; ++k) {
            sum[k] /= N;
        }

        for (let i = 0; i < N; ++i) {
            for (let k = 0; k < d; ++k) {
                Y.set_entry(i, k, G.entry(i, k) - sum[k]);
            }
        }
        return Y;
    }
}

var version="0.4.2";

exports.BallTree = BallTree;
exports.DisjointSet = DisjointSet;
exports.FASTMAP = FASTMAP;
exports.Heap = Heap;
exports.Hierarchical_Clustering = Hierarchical_Clustering;
exports.ISOMAP = ISOMAP;
exports.KMeans = KMeans;
exports.KMedoids = KMedoids;
exports.KNN = KNN;
exports.LDA = LDA;
exports.LLE = LLE;
exports.LSP = LSP;
exports.LTSA = LTSA;
exports.MDS = MDS;
exports.Matrix = Matrix;
exports.OPTICS = OPTICS;
exports.PCA = PCA;
exports.Randomizer = Randomizer;
exports.SAMMON = SAMMON;
exports.TSNE = TSNE;
exports.TopoMap = TopoMap;
exports.TriMap = TriMap;
exports.UMAP = UMAP;
exports.canberra = canberra;
exports.chebyshev = chebyshev;
exports.cosine = cosine;
exports.distance_matrix = distance_matrix;
exports.euclidean = euclidean;
exports.euclidean_squared = euclidean_squared;
exports.hamming = hamming;
exports.jaccard = jaccard;
exports.k_nearest_neighbors = k_nearest_neighbors;
exports.kahan_sum = kahan_sum;
exports.linspace = linspace;
exports.manhattan = manhattan;
exports.max = max;
exports.min = min;
exports.neumair_sum = neumair_sum;
exports.norm = norm;
exports.powell = powell;
exports.qr = qr;
exports.simultaneous_poweriteration = simultaneous_poweriteration$1;
exports.sokal_michener = sokal_michener;
exports.version = version;
exports.yule = yule;

Object.defineProperty(exports, '__esModule', { value: true });

}));
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiZHJ1aWQuanMiLCJzb3VyY2VzIjpbIi4uL21ldHJpY3MvZXVjbGlkZWFuLmpzIiwiLi4vbnVtZXJpY2FsL2thaGFuX3N1bS5qcyIsIi4uL251bWVyaWNhbC9uZXVtYWlyX3N1bS5qcyIsIi4uL21ldHJpY3MvZXVjbGlkZWFuX3NxdWFyZWQuanMiLCIuLi9tZXRyaWNzL2Nvc2luZS5qcyIsIi4uL21ldHJpY3MvbWFuaGF0dGFuLmpzIiwiLi4vbWV0cmljcy9jaGVieXNoZXYuanMiLCIuLi9tZXRyaWNzL2NhbmJlcnJhLmpzIiwiLi4vbWV0cmljcy9qYWNjYXJkLmpzIiwiLi4vbWV0cmljcy9oYW1taW5nLmpzIiwiLi4vbWV0cmljcy9zb2thbF9taWNoZW5lci5qcyIsIi4uL21ldHJpY3MveXVsZS5qcyIsIi4uL21hdHJpeC9rX25lYXJlc3RfbmVpZ2hib3JzLmpzIiwiLi4vbWF0cml4L01hdHJpeC5qcyIsIi4uL21hdHJpeC9kaXN0YW5jZV9tYXRyaXguanMiLCIuLi9tYXRyaXgvbGluc3BhY2UuanMiLCIuLi9tYXRyaXgvbm9ybS5qcyIsIi4uL3V0aWwvcmFuZG9taXplci5qcyIsIi4uL3V0aWwvbWF4LmpzIiwiLi4vdXRpbC9taW4uanMiLCIuLi9kYXRhc3RydWN0dXJlL0hlYXAuanMiLCIuLi9kYXRhc3RydWN0dXJlL0Rpc2pvaW50U2V0LmpzIiwiLi4va25uL0JhbGxUcmVlLmpzIiwiLi4va25uL0tOTi5qcyIsIi4uL2xpbmVhcl9hbGdlYnJhL3FyLmpzIiwiLi4vbGluZWFyX2FsZ2VicmEvc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uLmpzIiwiLi4vZGltcmVkL0RSLmpzIiwiLi4vZGltcmVkL1BDQS5qcyIsIi4uL2RpbXJlZC9NRFMuanMiLCIuLi9kaW1yZWQvSVNPTUFQLmpzIiwiLi4vZGltcmVkL0ZBU1RNQVAuanMiLCIuLi9kaW1yZWQvTERBLmpzIiwiLi4vZGltcmVkL0xMRS5qcyIsIi4uL2RpbXJlZC9MVFNBLmpzIiwiLi4vZGltcmVkL1RTTkUuanMiLCIuLi9vcHRpbWl6YXRpb24vcG93ZWxsLmpzIiwiLi4vZGltcmVkL1VNQVAuanMiLCIuLi9kaW1yZWQvVHJpTWFwLmpzIiwiLi4vY2x1c3RlcmluZy9IaWVyYXJjaGljYWxfQ2x1c3RlcmluZy5qcyIsIi4uL2NsdXN0ZXJpbmcvS01lYW5zLmpzIiwiLi4vY2x1c3RlcmluZy9LTWVkb2lkcy5qcyIsIi4uL2NsdXN0ZXJpbmcvT1BUSUNTLmpzIiwiLi4vZGltcmVkL0xTUC5qcyIsIi4uL2RpbXJlZC9Ub3BvTWFwLmpzIiwiLi4vZGltcmVkL1NBTU1PTi5qcyJdLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgeyBldWNsaWRlYW5fc3F1YXJlZCB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG4vKipcbiAqIENvbXB1dGVzIHRoZSBldWNsaWRlYW4gZGlzdGFuY2UgKGw8c3ViPjI8L3N1Yj4pIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm1ldHJpY3NcbiAqIEBhbGlhcyBldWNsaWRlYW5cbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYVxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBiXG4gKiBAcmV0dXJucyB7TnVtYmVyfSB0aGUgZXVjbGlkZWFuIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChhLCBiKSB7XG4gICAgcmV0dXJuIE1hdGguc3FydChldWNsaWRlYW5fc3F1YXJlZChhLCBiKSk7XG59XG4iLCIvKipcbiAqIE51bWVyaWNhbCBzdGFibGUgc3VtbWF0aW9uIHdpdGggdGhlIEthaGFuIHN1bW1hdGlvbiBhbGdvcml0aG0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm51bWVyaWNhbFxuICogQGFsaWFzIGthaGFuX3N1bVxuICogQHBhcmFtIHtBcnJheX0gc3VtbWFuZHMgLSBBcnJheSBvZiB2YWx1ZXMgdG8gc3VtIHVwLlxuICogQHJldHVybnMge251bWJlcn0gVGhlIHN1bS5cbiAqIEBzZWUge0BsaW5rIGh0dHBzOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL0thaGFuX3N1bW1hdGlvbl9hbGdvcml0aG19XG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChzdW1tYW5kcykge1xuICAgIGxldCBuID0gc3VtbWFuZHMubGVuZ3RoO1xuICAgIGxldCBzdW0gPSAwO1xuICAgIGxldCBjb21wZW5zYXRpb24gPSAwO1xuICAgIGxldCB5LCB0O1xuXG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgeSA9IHN1bW1hbmRzW2ldIC0gY29tcGVuc2F0aW9uO1xuICAgICAgICB0ID0gc3VtICsgeTtcbiAgICAgICAgY29tcGVuc2F0aW9uID0gdCAtIHN1bSAtIHk7XG4gICAgICAgIHN1bSA9IHQ7XG4gICAgfVxuICAgIHJldHVybiBzdW07XG59XG4iLCIvKipcbiAqIE51bWVyaWNhbCBzdGFibGUgc3VtbWF0aW9uIHdpdGggdGhlIE5ldW1haXIgc3VtbWF0aW9uIGFsZ29yaXRobS5cbiAqIEBtZW1iZXJvZiBtb2R1bGU6bnVtZXJpY2FsXG4gKiBAYWxpYXMgbmV1bWFpcl9zdW1cbiAqIEBwYXJhbSB7QXJyYXl9IHN1bW1hbmRzIC0gQXJyYXkgb2YgdmFsdWVzIHRvIHN1bSB1cC5cbiAqIEByZXR1cm5zIHtudW1iZXJ9IFRoZSBzdW0uXG4gKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9LYWhhbl9zdW1tYXRpb25fYWxnb3JpdGhtI0Z1cnRoZXJfZW5oYW5jZW1lbnRzfVxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAoc3VtbWFuZHMpIHtcbiAgICBsZXQgbiA9IHN1bW1hbmRzLmxlbmd0aDtcbiAgICBsZXQgc3VtID0gMDtcbiAgICBsZXQgY29tcGVuc2F0aW9uID0gMDtcblxuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbjsgKytpKSB7XG4gICAgICAgIGxldCBzdW1tYW5kID0gc3VtbWFuZHNbaV07XG4gICAgICAgIGxldCB0ID0gc3VtICsgc3VtbWFuZDtcbiAgICAgICAgaWYgKE1hdGguYWJzKHN1bSkgPj0gTWF0aC5hYnMoc3VtbWFuZCkpIHtcbiAgICAgICAgICAgIGNvbXBlbnNhdGlvbiArPSBzdW0gLSB0ICsgc3VtbWFuZDtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGNvbXBlbnNhdGlvbiArPSBzdW1tYW5kIC0gdCArIHN1bTtcbiAgICAgICAgfVxuICAgICAgICBzdW0gPSB0O1xuICAgIH1cbiAgICByZXR1cm4gc3VtICsgY29tcGVuc2F0aW9uO1xufVxuIiwiaW1wb3J0IHsgbmV1bWFpcl9zdW0gfSBmcm9tIFwiLi4vbnVtZXJpY2FsL2luZGV4LmpzXCI7XG4vKipcbiAqIENvbXB1dGVzIHRoZSBzcXVhcmVkIGV1Y2xpZGVhbiBkaXN0YW5jZSAobDxzdWI+Mjwvc3ViPjxzdXA+Mjwvc3VwPikgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqIEBtZW1iZXJvZiBtb2R1bGU6bWV0cmljc1xuICogQGFsaWFzIGV1Y2xpZGVhbl9zcXVhcmVkXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGFcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYlxuICogQHJldHVybnMge051bWJlcn0gdGhlIHNxdWFyZWQgZXVjbGlkZWFuIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChhLCBiKSB7XG4gICAgaWYgKGEubGVuZ3RoICE9IGIubGVuZ3RoKSByZXR1cm4gdW5kZWZpbmVkO1xuICAgIGxldCBuID0gYS5sZW5ndGg7XG4gICAgbGV0IHMgPSBuZXcgQXJyYXkobik7XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgbGV0IHggPSBhW2ldO1xuICAgICAgICBsZXQgeSA9IGJbaV07XG4gICAgICAgIHNbaV0gPSAoeCAtIHkpICogKHggLSB5KTtcbiAgICB9XG4gICAgcmV0dXJuIG5ldW1haXJfc3VtKHMpO1xufVxuIiwiLyoqXG4gKiBDb21wdXRlcyB0aGUgY29zaW5lIGRpc3RhbmNlIChub3Qgc2ltaWxhcml0eSkgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqIEBtZW1iZXJvZiBtb2R1bGU6bWV0cmljc1xuICogQGFsaWFzIGNvc2luZVxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBhXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGJcbiAqIEBleGFtcGxlXG4gKiBkcnVpZC5jb3NpbmUoWzEsMF0sWzEsMV0pID09IDAuNzg1Mzk4MTYzMzk3NDQ4NCA9PSDPgC80XG4gKiBAcmV0dXJucyB7TnVtYmVyfSBUaGUgY29zaW5lIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChhLCBiKSB7XG4gICAgaWYgKGEubGVuZ3RoICE9PSBiLmxlbmd0aCkgcmV0dXJuIHVuZGVmaW5lZDtcbiAgICBsZXQgbiA9IGEubGVuZ3RoO1xuICAgIGxldCBzdW0gPSAwO1xuICAgIGxldCBzdW1fYSA9IDA7XG4gICAgbGV0IHN1bV9iID0gMDtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICBzdW0gKz0gYVtpXSAqIGJbaV07XG4gICAgICAgIHN1bV9hICs9IGFbaV0gKiBhW2ldO1xuICAgICAgICBzdW1fYiArPSBiW2ldICogYltpXTtcbiAgICB9XG4gICAgcmV0dXJuIE1hdGguYWNvcyhzdW0gLyAoTWF0aC5zcXJ0KHN1bV9hKSAqIE1hdGguc3FydChzdW1fYikpKTtcbn1cbiIsIi8qKlxuICogQ29tcHV0ZXMgdGhlIG1hbmhhdHRhbiBkaXN0YW5jZSAobDxzdWI+MTwvc3ViPikgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqIEBtZW1iZXJvZiBtb2R1bGU6bWV0cmljc1xuICogQGFsaWFzIG1hbmhhdHRhblxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBhXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGJcbiAqIEByZXR1cm5zIHtOdW1iZXJ9IHRoZSBtYW5oYXR0YW4gZGlzdGFuY2UgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqLyBcbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChhLCBiKSB7XG4gICAgaWYgKGEubGVuZ3RoICE9IGIubGVuZ3RoKSByZXR1cm4gdW5kZWZpbmVkO1xuICAgIGxldCBuID0gYS5sZW5ndGg7XG4gICAgbGV0IHN1bSA9IDA7XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgc3VtICs9IE1hdGguYWJzKGFbaV0gLSBiW2ldKTtcbiAgICB9XG4gICAgcmV0dXJuIHN1bTtcbn1cbiIsIi8qKlxuICogQ29tcHV0ZXMgdGhlIGNoZWJ5c2hldiBkaXN0YW5jZSAoTDxzdWI+4oiePC9zdWI+KSBiZXR3ZWVuIHtAbGluayBhfSBhbmQge0BsaW5rIGJ9LlxuICogQG1lbWJlcm9mIG1vZHVsZTptZXRyaWNzXG4gKiBAYWxpYXMgY2hlYnlzaGV2XG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGFcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYlxuICogQHJldHVybnMge051bWJlcn0gdGhlIGNoZWJ5c2hldiBkaXN0YW5jZSBiZXR3ZWVuIHtAbGluayBhfSBhbmQge0BsaW5rIGJ9LlxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAoYSwgYikge1xuICAgIGlmIChhLmxlbmd0aCAhPSBiLmxlbmd0aCkgcmV0dXJuIHVuZGVmaW5lZDtcbiAgICBsZXQgbiA9IGEubGVuZ3RoO1xuICAgIGxldCByZXMgPSBbXTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICByZXMucHVzaChNYXRoLmFicyhhW2ldIC0gYltpXSkpO1xuICAgIH1cbiAgICByZXR1cm4gTWF0aC5tYXgoLi4ucmVzKTtcbn1cbiIsIi8qKlxuICogQ29tcHV0ZXMgdGhlIGNhbmJlcnJhIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm1ldHJpY3NcbiAqIEBhbGlhcyBjYW5iZXJyYVxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBhIFxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBiIFxuICogQHJldHVybnMge051bWJlcn0gVGhlIGNhbmJlcnJhIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9DYW5iZXJyYV9kaXN0YW5jZX1cbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24oYSwgYikge1xuICAgIGlmIChhLmxlbmd0aCAhPT0gYi5sZW5ndGgpIHJldHVybiB1bmRlZmluZWQ7XG4gICAgbGV0IG4gPSBhLmxlbmd0aDtcbiAgICBsZXQgc3VtID0gMDtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICBzdW0gKz0gKE1hdGguYWJzKGFbaV0gLSBiW2ldKSAvIChNYXRoLmFicyhhW2ldKSArIE1hdGguYWJzKGJbaV0pKSlcbiAgICB9XG4gICAgcmV0dXJuIHN1bTtcbn0iLCIvKipcbiAqIENvbXB1dGVzIHRoZSBqYWNjYXJkIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm1ldHJpY3NcbiAqIEBhbGlhcyBqYWNjYXJkXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGFcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYlxuICogQHJldHVybnMge051bWJlcn0gdGhlIGphY2NhcmQgZGlzdGFuY2UgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24gKGEsIGIpIHtcbiAgICBpZiAoYS5sZW5ndGggIT0gYi5sZW5ndGgpIHJldHVybiB1bmRlZmluZWQ7XG4gICAgY29uc3QgbiA9IGEubGVuZ3RoO1xuICAgIGxldCBudW1fbm9uX3plcm8gPSAwO1xuICAgIGxldCBudW1fZXF1YWwgPSAwO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbjsgKytpKSB7XG4gICAgICAgIGNvbnN0IHggPSBhW2ldICE9IDA7XG4gICAgICAgIGNvbnN0IHkgPSBiW2ldICE9IDA7XG4gICAgICAgIG51bV9ub25femVybyArPSB4IHx8IHk7XG4gICAgICAgIG51bV9lcXVhbCArPSB4ICYmIHk7XG4gICAgfVxuICAgIHJldHVybiAobnVtX25vbl96ZXJvIC0gbnVtX2VxdWFsKSAvIG51bV9ub25femVybztcbn1cbiIsIi8qKlxuICogQ29tcHV0ZXMgdGhlIGhhbW1pbmcgZGlzdGFuY2UgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqIEBtZW1iZXJvZiBtb2R1bGU6bWV0cmljc1xuICogQGFsaWFzIGhhbW1pbmdcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYVxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBiXG4gKiBAcmV0dXJucyB7TnVtYmVyfSB0aGUgaGFtbWluZyBkaXN0YW5jZSBiZXR3ZWVuIHtAbGluayBhfSBhbmQge0BsaW5rIGJ9LlxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAoYSwgYikge1xuICAgIGlmIChhLmxlbmd0aCAhPSBiLmxlbmd0aCkgcmV0dXJuIHVuZGVmaW5lZDtcbiAgICBjb25zdCBuID0gYS5sZW5ndGg7XG4gICAgbGV0IGRpc2FncmVlID0gMDtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICBjb25zdCB4ID0gYVtpXTtcbiAgICAgICAgY29uc3QgeSA9IGJbaV07XG4gICAgICAgIGRpc2FncmVlICs9IHggIT0geTtcbiAgICB9XG4gICAgcmV0dXJuIGRpc2FncmVlIC8gbjtcbn1cbiIsIi8qKlxuICogQ29tcHV0ZXMgdGhlIFNva2FsLU1pY2hlbmVyIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm1ldHJpY3NcbiAqIEBhbGlhcyBzb2thbF9taWNoZW5lclxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBhIFxuICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBiIFxuICogQHJldHVybnMge051bWJlcn0gdGhlIFNva2FsLU1pY2hlbmVyIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uICBcbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24oYSwgYikge1xuICAgIGlmIChhLmxlbmd0aCAhPSBiLmxlbmd0aCkgcmV0dXJuIHVuZGVmaW5lZFxuICAgIGNvbnN0IG4gPSBhLmxlbmd0aDtcbiAgICBsZXQgbnVtX25vdF9lcXVhbCA9IDA7XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgY29uc3QgeCA9IGFbaV0gIT0gMDtcbiAgICAgICAgY29uc3QgeSA9IGJbaV0gIT0gMDtcbiAgICAgICAgbnVtX25vdF9lcXVhbCArPSB4ICE9IHk7XG4gICAgfVxuICAgIHJldHVybiAoMiAqIG51bV9ub3RfZXF1YWwpIC8gKG4gKyBudW1fbm90X2VxdWFsKTtcbn0iLCIvKipcbiAqIENvbXB1dGVzIHRoZSB5dWxlIGRpc3RhbmNlIGJldHdlZW4ge0BsaW5rIGF9IGFuZCB7QGxpbmsgYn0uXG4gKiBAbWVtYmVyb2YgbW9kdWxlOm1ldHJpY3NcbiAqIEBhbGlhcyB5dWxlXG4gKiBAcGFyYW0ge0FycmF5PE51bWJlcj59IGFcbiAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gYlxuICogQHJldHVybnMge051bWJlcn0gdGhlIHl1bGUgZGlzdGFuY2UgYmV0d2VlbiB7QGxpbmsgYX0gYW5kIHtAbGluayBifS5cbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24gKGEsIGIpIHtcbiAgICBpZiAoYS5sZW5ndGggIT0gYi5sZW5ndGgpIHJldHVybiB1bmRlZmluZWQ7XG4gICAgY29uc3QgbiA9IGEubGVuZ3RoO1xuICAgIGxldCBudW1fdHJ1ZV90cnVlID0gMDtcbiAgICBsZXQgbnVtX3RydWVfZmFsc2UgPSAwO1xuICAgIGxldCBudW1fZmFsc2VfdHJ1ZSA9IDA7XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgY29uc3QgeCA9IGFbaV0gIT0gMDtcbiAgICAgICAgY29uc3QgeSA9IGJbaV0gIT0gMDtcbiAgICAgICAgbnVtX3RydWVfdHJ1ZSArPSB4ICYmIHk7XG4gICAgICAgIG51bV90cnVlX2ZhbHNlICs9IHggJiYgIXk7XG4gICAgICAgIG51bV9mYWxzZV90cnVlICs9ICF4ICYmIHg7XG4gICAgfVxuICAgIGNvbnN0IG51bV9mYWxzZV9mYWxzZSA9IG4gLSBudW1fdHJ1ZV90cnVlIC0gbnVtX3RydWVfZmFsc2UgLSBudW1fZmFsc2VfdHJ1ZTtcbiAgICByZXR1cm4gbnVtX3RydWVfZmFsc2UgPT0gMCB8fCBudW1fZmFsc2VfdHJ1ZSA9PSAwID8gMCA6ICgyICogbnVtX3RydWVfZmFsc2UgKiBudW1fZmFsc2VfdHJ1ZSkgLyAobnVtX3RydWVfdHJ1ZSAqIG51bV9mYWxzZV9mYWxzZSArIG51bV90cnVlX2ZhbHNlICogbnVtX2ZhbHNlX3RydWUpO1xufVxuIiwiaW1wb3J0IHsgZGlzdGFuY2VfbWF0cml4IGFzIGRtYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBldWNsaWRlYW4gfSBmcm9tIFwiLi4vbWV0cmljcy9pbmRleC5qc1wiO1xuXG4vKipcbiAqXG4gKiBAcGFyYW0geyp9IEFcbiAqIEBwYXJhbSB7Kn0ga1xuICogQHBhcmFtIHsqfSBkaXN0YW5jZV9tYXRyaXhcbiAqIEBwYXJhbSB7Kn0gbWV0cmljXG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChBLCBrLCBkaXN0YW5jZV9tYXRyaXggPSBudWxsLCBtZXRyaWMgPSBldWNsaWRlYW4pIHtcbiAgICBjb25zdCByb3dzID0gQS5zaGFwZVswXTtcbiAgICBsZXQgRCA9IGRpc3RhbmNlX21hdHJpeCA/PyBkbWF0cml4KEEsIG1ldHJpYyk7XG4gICAgbGV0IG5OID0gbmV3IEFycmF5KHJvd3MpO1xuICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3M7ICsrcm93KSB7XG4gICAgICAgIG5OW3Jvd10gPSBBcnJheS5mcm9tKEQucm93KHJvdykpXG4gICAgICAgICAgICAubWFwKChkaXN0YW5jZSwgY29sKSA9PiB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHtcbiAgICAgICAgICAgICAgICAgICAgaTogcm93LFxuICAgICAgICAgICAgICAgICAgICBqOiBjb2wsXG4gICAgICAgICAgICAgICAgICAgIGRpc3RhbmNlOiBkaXN0YW5jZSxcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgfSlcbiAgICAgICAgICAgIC5zb3J0KChhLCBiKSA9PiBhLmRpc3RhbmNlIC0gYi5kaXN0YW5jZSlcbiAgICAgICAgICAgIC5zbGljZSgxLCBrICsgMSk7XG4gICAgfVxuICAgIHJldHVybiBuTjtcbn1cbiIsImltcG9ydCB7IG5ldW1haXJfc3VtIH0gZnJvbSBcIi4uL251bWVyaWNhbC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgUmFuZG9taXplciB9IGZyb20gXCIuLi91dGlsL2luZGV4LmpzXCI7XG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIE1hdHJpeFxuICogQHJlcXVpcmVzIG1vZHVsZTpudW1lcmljYWwvbmV1bWFpcl9zdW1cbiAqL1xuZXhwb3J0IGNsYXNzIE1hdHJpeCB7XG4gICAgLyoqXG4gICAgICogY3JlYXRlcyBhIG5ldyBNYXRyaXguIEVudHJpZXMgYXJlIHN0b3JlZCBpbiBhIEZsb2F0NjRBcnJheS5cbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOm1hdHJpeFxuICAgICAqIEBhbGlhcyBNYXRyaXhcbiAgICAgKiBAcGFyYW0ge251bWJlcn0gcm93cyAtIFRoZSBhbW91bnQgb2Ygcm93cyBvZiB0aGUgbWF0cml4LlxuICAgICAqIEBwYXJhbSB7bnVtYmVyfSBjb2xzIC0gVGhlIGFtb3VudCBvZiBjb2x1bW5zIG9mIHRoZSBtYXRyaXguXG4gICAgICogQHBhcmFtIHsoZnVuY3Rpb258c3RyaW5nfG51bWJlcil9IHZhbHVlPTAgLSBDYW4gYmUgYSBmdW5jdGlvbiB3aXRoIHJvdyBhbmQgY29sIGFzIHBhcmFtZXRlcnMsIGEgbnVtYmVyLCBvciBcInplcm9zXCIsIFwiaWRlbnRpdHlcIiBvciBcIklcIiwgb3IgXCJjZW50ZXJcIi5cbiAgICAgKiAgLSAqKmZ1bmN0aW9uKio6IGZvciBlYWNoIGVudHJ5IHRoZSBmdW5jdGlvbiBnZXRzIGNhbGxlZCB3aXRoIHRoZSBwYXJhbWV0ZXJzIGZvciB0aGUgYWN0dWFsIHJvdyBhbmQgY29sdW1uLlxuICAgICAqICAtICoqc3RyaW5nKio6IGFsbG93ZWQgYXJlXG4gICAgICogICAgICAtIFwiemVyb1wiLCBjcmVhdGVzIGEgemVybyBtYXRyaXguXG4gICAgICogICAgICAtIFwiaWRlbnRpdHlcIiBvciBcIklcIiwgY3JlYXRlcyBhbiBpZGVudGl0eSBtYXRyaXguXG4gICAgICogICAgICAtIFwiY2VudGVyXCIsIGNyZWF0ZXMgYW4gY2VudGVyIG1hdHJpeC5cbiAgICAgKiAgLSAqKm51bWJlcioqOiBjcmVhdGUgYSBtYXRyaXggZmlsbGVkIHdpdGggdGhlIGdpdmVuIHZhbHVlLlxuICAgICAqIEBleGFtcGxlXG4gICAgICpcbiAgICAgKiBsZXQgQSA9IG5ldyBNYXRyaXgoMTAsIDEwLCAoKSA9PiBNYXRoLnJhbmRvbSgpKTsgLy9jcmVhdGVzIGEgMTAgdGltZXMgMTAgcmFuZG9tIG1hdHJpeC5cbiAgICAgKiBsZXQgQiA9IG5ldyBNYXRyaXgoMywgMywgXCJJXCIpOyAvLyBjcmVhdGVzIGEgMyB0aW1lcyAzIGlkZW50aXR5IG1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fSByZXR1cm5zIGEge0BsaW5rIHJvd3N9IHRpbWVzIHtAbGluayBjb2xzfSBNYXRyaXggZmlsbGVkIHdpdGgge0BsaW5rIHZhbHVlfS5cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3Rvcihyb3dzID0gbnVsbCwgY29scyA9IG51bGwsIHZhbHVlID0gbnVsbCkge1xuICAgICAgICB0aGlzLl9yb3dzID0gcm93cztcbiAgICAgICAgdGhpcy5fY29scyA9IGNvbHM7XG4gICAgICAgIHRoaXMuX2RhdGEgPSBudWxsO1xuICAgICAgICBpZiAocm93cyAmJiBjb2xzKSB7XG4gICAgICAgICAgICBpZiAoIXZhbHVlKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5fZGF0YSA9IG5ldyBGbG9hdDY0QXJyYXkocm93cyAqIGNvbHMpO1xuICAgICAgICAgICAgICAgIHJldHVybiB0aGlzO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKHR5cGVvZiB2YWx1ZSA9PT0gXCJmdW5jdGlvblwiKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5fZGF0YSA9IG5ldyBGbG9hdDY0QXJyYXkocm93cyAqIGNvbHMpO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3M7ICsrcm93KSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLl9kYXRhW3JvdyAqIGNvbHMgKyBjb2xdID0gdmFsdWUocm93LCBjb2wpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHJldHVybiB0aGlzO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKHR5cGVvZiB2YWx1ZSA9PT0gXCJzdHJpbmdcIikge1xuICAgICAgICAgICAgICAgIGlmICh2YWx1ZSA9PT0gXCJ6ZXJvc1wiKSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiBuZXcgTWF0cml4KHJvd3MsIGNvbHMsIDApO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBpZiAodmFsdWUgPT09IFwiaWRlbnRpdHlcIiB8fCB2YWx1ZSA9PT0gXCJJXCIpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5fZGF0YSA9IG5ldyBGbG9hdDY0QXJyYXkocm93cyAqIGNvbHMpO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5fZGF0YVtyb3cgKiBjb2xzICsgcm93XSA9IDE7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGlmICh2YWx1ZSA9PT0gXCJjZW50ZXJcIiAmJiByb3dzID09IGNvbHMpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5fZGF0YSA9IG5ldyBGbG9hdDY0QXJyYXkocm93cyAqIGNvbHMpO1xuICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IChpLCBqKSA9PiAoaSA9PT0gaiA/IDEgOiAwKSAtIDEgLyByb3dzO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgICAgICAgICAgICAgZm9yIChsZXQgY29sID0gMDsgY29sIDwgY29sczsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB0aGlzLl9kYXRhW3JvdyAqIGNvbHMgKyBjb2xdID0gdmFsdWUocm93LCBjb2wpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiB0aGlzO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmICh0eXBlb2YgdmFsdWUgPT09IFwibnVtYmVyXCIpIHtcbiAgICAgICAgICAgICAgICB0aGlzLl9kYXRhID0gbmV3IEZsb2F0NjRBcnJheShyb3dzICogY29scyk7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgICAgICAgICAgZm9yIChsZXQgY29sID0gMDsgY29sIDwgY29sczsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMuX2RhdGFbcm93ICogY29scyArIGNvbF0gPSB2YWx1ZTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICByZXR1cm4gdGhpcztcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDcmVhdGVzIGEgTWF0cml4IG91dCBvZiB7QGxpbmsgQX0uXG4gICAgICogQHBhcmFtIHsoTWF0cml4fEFycmF5fEZsb2F0NjRBcnJheXxudW1iZXIpfSBBIC0gVGhlIG1hdHJpeCwgYXJyYXksIG9yIG51bWJlciwgd2hpY2ggc2hvdWxkIGNvbnZlcnRlZCB0byBhIE1hdHJpeC5cbiAgICAgKiBAcGFyYW0ge1wicm93XCJ8XCJjb2xcInxcImRpYWdcIn0gW3R5cGUgPSBcInJvd1wiXSAtIElmIHtAbGluayBBfSBpcyBhIEFycmF5IG9yIEZsb2F0NjRBcnJheSwgdGhlbiB0eXBlIGRlZmluZXMgaWYgaXQgaXMgYSByb3ctIG9yIGEgY29sdW1uIHZlY3Rvci5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqXG4gICAgICogQGV4YW1wbGVcbiAgICAgKiBsZXQgQSA9IE1hdHJpeC5mcm9tKFtbMSwgMF0sIFswLCAxXV0pOyAvL2NyZWF0ZXMgYSB0d28gYnkgdHdvIGlkZW50aXR5IG1hdHJpeC5cbiAgICAgKiBsZXQgUyA9IE1hdHJpeC5mcm9tKFsxLCAyLCAzXSwgXCJkaWFnXCIpOyAvLyBjcmVhdGVzIGEgMyBieSAzIG1hdHJpeCB3aXRoIDEsIDIsIDMgb24gaXRzIGRpYWdvbmFsLiBbWzEsIDAsIDBdLCBbMCwgMiwgMF0sIFswLCAwLCAzXV1cbiAgICAgKi9cbiAgICBzdGF0aWMgZnJvbShBLCB0eXBlID0gXCJyb3dcIikge1xuICAgICAgICBpZiAoQSBpbnN0YW5jZW9mIE1hdHJpeCkge1xuICAgICAgICAgICAgcmV0dXJuIEEuY2xvbmUoKTtcbiAgICAgICAgfSBlbHNlIGlmIChBcnJheS5pc0FycmF5KEEpIHx8IEEgaW5zdGFuY2VvZiBGbG9hdDY0QXJyYXkpIHtcbiAgICAgICAgICAgIGxldCBtID0gQS5sZW5ndGg7XG4gICAgICAgICAgICBpZiAobSA9PT0gMCkgdGhyb3cgbmV3IEVycm9yKFwiQXJyYXkgaXMgZW1wdHlcIik7XG4gICAgICAgICAgICAvLyAxZFxuICAgICAgICAgICAgaWYgKCFBcnJheS5pc0FycmF5KEFbMF0pICYmICEoQVswXSBpbnN0YW5jZW9mIEZsb2F0NjRBcnJheSkpIHtcbiAgICAgICAgICAgICAgICBpZiAodHlwZSA9PT0gXCJyb3dcIikge1xuICAgICAgICAgICAgICAgICAgICByZXR1cm4gbmV3IE1hdHJpeCgxLCBtLCAoXywgaikgPT4gQVtqXSk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmICh0eXBlID09PSBcImNvbFwiKSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiBuZXcgTWF0cml4KG0sIDEsIChpKSA9PiBBW2ldKTtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT09IFwiZGlhZ1wiKSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiBuZXcgTWF0cml4KG0sIG0sIChpLCBqKSA9PiAoaSA9PSBqID8gQVtpXSA6IDApKTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoXCIxZCBhcnJheSBoYXMgTmFOIGVudHJpZXNcIik7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIC8vIDJkXG4gICAgICAgICAgICB9IGVsc2UgaWYgKEFycmF5LmlzQXJyYXkoQVswXSkgfHwgQVswXSBpbnN0YW5jZW9mIEZsb2F0NjRBcnJheSkge1xuICAgICAgICAgICAgICAgIGxldCBuID0gQVswXS5sZW5ndGg7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgbTsgKytyb3cpIHtcbiAgICAgICAgICAgICAgICAgICAgaWYgKEFbcm93XS5sZW5ndGggIT09IG4pIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcInZhcmlvdXMgYXJyYXkgbGVuZ3Roc1wiKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICByZXR1cm4gbmV3IE1hdHJpeChtLCBuLCAoaSwgaikgPT4gQVtpXVtqXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSBpZiAodHlwZW9mIEEgPT09IFwibnVtYmVyXCIpIHtcbiAgICAgICAgICAgIHJldHVybiBuZXcgTWF0cml4KDEsIDEsIEEpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKFwiZXJyb3JcIik7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRoZSB7QGxpbmsgcm93fTxzdXA+dGg8L3N1cD4gcm93IGZyb20gdGhlIE1hdHJpeC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gcm93XG4gICAgICogQHJldHVybnMge0Zsb2F0NjRBcnJheX1cbiAgICAgKi9cbiAgICByb3cocm93KSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLnZhbHVlcztcbiAgICAgICAgY29uc3QgY29scyA9IHRoaXMuX2NvbHM7XG4gICAgICAgIHJldHVybiBkYXRhLnN1YmFycmF5KHJvdyAqIGNvbHMsIChyb3cgKyAxKSAqIGNvbHMpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgYW4gZ2VuZXJhdG9yIHlpZWxkaW5nIGVhY2ggcm93IG9mIHRoZSBNYXRyaXguXG4gICAgICogQHlpZWxkcyB7RmxvYXQ2NEFycmF5fVxuICAgICAqL1xuICAgICppdGVyYXRlX3Jvd3MoKSB7XG4gICAgICAgIGNvbnN0IGNvbHMgPSB0aGlzLl9jb2xzO1xuICAgICAgICBjb25zdCByb3dzID0gdGhpcy5fcm93cztcbiAgICAgICAgY29uc3QgZGF0YSA9IHRoaXMudmFsdWVzO1xuICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgeWllbGQgZGF0YS5zdWJhcnJheShyb3cgKiBjb2xzLCAocm93ICsgMSkgKiBjb2xzKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIE1ha2VzIGEge0BsaW5rIE1hdHJpeH0gb2JqZWN0IGFuIGl0ZXJhYmxlIG9iamVjdC5cbiAgICAgKiBAeWllbGRzIHtGbG9hdDY0QXJyYXl9XG4gICAgICovXG4gICAgKltTeW1ib2wuaXRlcmF0b3JdKCkge1xuICAgICAgICBmb3IgKGNvbnN0IHJvdyBvZiB0aGlzLml0ZXJhdGVfcm93cygpKSB7XG4gICAgICAgICAgICB5aWVsZCByb3c7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBTZXRzIHRoZSBlbnRyaWVzIG9mIHtAbGluayByb3d9PHN1cD50aDwvc3VwPiByb3cgZnJvbSB0aGUgTWF0cml4IHRvIHRoZSBlbnRyaWVzIGZyb20ge0BsaW5rIHZhbHVlc30uXG4gICAgICogQHBhcmFtIHtpbnR9IHJvd1xuICAgICAqIEBwYXJhbSB7QXJyYXl9IHZhbHVlc1xuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgc2V0X3Jvdyhyb3csIHZhbHVlcykge1xuICAgICAgICBsZXQgY29scyA9IHRoaXMuX2NvbHM7XG4gICAgICAgIGlmIChBcnJheS5pc0FycmF5KHZhbHVlcykgJiYgdmFsdWVzLmxlbmd0aCA9PT0gY29scykge1xuICAgICAgICAgICAgbGV0IG9mZnNldCA9IHJvdyAqIGNvbHM7XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIHRoaXMudmFsdWVzW29mZnNldCArIGNvbF0gPSB2YWx1ZXNbY29sXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIGlmICh2YWx1ZXMgaW5zdGFuY2VvZiBNYXRyaXggJiYgdmFsdWVzLnNoYXBlWzFdID09PSBjb2xzICYmIHZhbHVlcy5zaGFwZVswXSA9PT0gMSkge1xuICAgICAgICAgICAgbGV0IG9mZnNldCA9IHJvdyAqIGNvbHM7XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIHRoaXMudmFsdWVzW29mZnNldCArIGNvbF0gPSB2YWx1ZXMuX2RhdGFbY29sXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRoZSB7QGxpbmsgY29sfTxzdXA+dGg8L3N1cD4gY29sdW1uIGZyb20gdGhlIE1hdHJpeC5cbiAgICAgKiBAcGFyYW0ge2ludH0gY29sXG4gICAgICogQHJldHVybnMge0FycmF5fVxuICAgICAqL1xuICAgIGNvbChjb2wpIHtcbiAgICAgICAgbGV0IHJlc3VsdF9jb2wgPSBuZXcgRmxvYXQ2NEFycmF5KHRoaXMuX3Jvd3MpO1xuICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCB0aGlzLl9yb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgcmVzdWx0X2NvbFtyb3ddID0gdGhpcy52YWx1ZXNbcm93ICogdGhpcy5fY29scyArIGNvbF07XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHJlc3VsdF9jb2w7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUge0BsaW5rIGNvbH08c3VwPnRoPC9zdXA+IGVudHJ5IGZyb20gdGhlIHtAbGluayByb3d9PHN1cD50aDwvc3VwPiByb3cgb2YgdGhlIE1hdHJpeC5cbiAgICAgKiBAcGFyYW0ge2ludH0gcm93XG4gICAgICogQHBhcmFtIHtpbnR9IGNvbFxuICAgICAqIEByZXR1cm5zIHtmbG9hdDY0fVxuICAgICAqL1xuICAgIGVudHJ5KHJvdywgY29sKSB7XG4gICAgICAgIHJldHVybiB0aGlzLnZhbHVlc1tyb3cgKiB0aGlzLl9jb2xzICsgY29sXTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBTZXRzIHRoZSB7QGxpbmsgY29sfTxzdXA+dGg8L3N1cD4gZW50cnkgZnJvbSB0aGUge0BsaW5rIHJvd308c3VwPnRoPC9zdXA+IHJvdyBvZiB0aGUgTWF0cml4IHRvIHRoZSBnaXZlbiB7QGxpbmsgdmFsdWV9LlxuICAgICAqIEBwYXJhbSB7aW50fSByb3dcbiAgICAgKiBAcGFyYW0ge2ludH0gY29sXG4gICAgICogQHBhcmFtIHtmbG9hdDY0fSB2YWx1ZVxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgc2V0X2VudHJ5KHJvdywgY29sLCB2YWx1ZSkge1xuICAgICAgICB0aGlzLnZhbHVlc1tyb3cgKiB0aGlzLl9jb2xzICsgY29sXSA9IHZhbHVlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIGEgbmV3IHRyYW5zcG9zZWQgTWF0cml4LlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgdHJhbnNwb3NlKCkge1xuICAgICAgICBsZXQgQiA9IG5ldyBNYXRyaXgodGhpcy5fY29scywgdGhpcy5fcm93cywgKHJvdywgY29sKSA9PiB0aGlzLmVudHJ5KGNvbCwgcm93KSk7XG4gICAgICAgIHJldHVybiBCO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgYSBuZXcgdHJhbnNwb3NlZCBNYXRyaXguIFNob3J0LWZvcm0gb2Yge0BmdW5jdGlvbiB0cmFuc3Bvc2V9LlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgZ2V0IFQoKSB7XG4gICAgICAgIHJldHVybiB0aGlzLnRyYW5zcG9zZSgpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIGludmVyc2Ugb2YgdGhlIE1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqL1xuICAgIGludmVyc2UoKSB7XG4gICAgICAgIGNvbnN0IHJvd3MgPSB0aGlzLl9yb3dzO1xuICAgICAgICBjb25zdCBjb2xzID0gdGhpcy5fY29scztcbiAgICAgICAgbGV0IEIgPSBuZXcgTWF0cml4KHJvd3MsIDIgKiBjb2xzLCAoaSwgaikgPT4ge1xuICAgICAgICAgICAgaWYgKGogPj0gY29scykge1xuICAgICAgICAgICAgICAgIHJldHVybiBpID09PSBqIC0gY29scyA/IDEgOiAwO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gdGhpcy5lbnRyeShpLCBqKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgICAgIGxldCBoID0gMDtcbiAgICAgICAgbGV0IGsgPSAwO1xuICAgICAgICB3aGlsZSAoaCA8IHJvd3MgJiYgayA8IGNvbHMpIHtcbiAgICAgICAgICAgIHZhciBpX21heCA9IDA7XG4gICAgICAgICAgICBsZXQgbWF4X3ZhbCA9IC1JbmZpbml0eTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSBoOyBpIDwgcm93czsgKytpKSB7XG4gICAgICAgICAgICAgICAgbGV0IHZhbCA9IE1hdGguYWJzKEIuZW50cnkoaSwgaykpO1xuICAgICAgICAgICAgICAgIGlmIChtYXhfdmFsIDwgdmFsKSB7XG4gICAgICAgICAgICAgICAgICAgIGlfbWF4ID0gaTtcbiAgICAgICAgICAgICAgICAgICAgbWF4X3ZhbCA9IHZhbDtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAoQi5lbnRyeShpX21heCwgaykgPT0gMCkge1xuICAgICAgICAgICAgICAgIGsrKztcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgLy8gc3dhcCByb3dzXG4gICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCAyICogY29sczsgKytqKSB7XG4gICAgICAgICAgICAgICAgICAgIGxldCBoX3ZhbCA9IEIuZW50cnkoaCwgaik7XG4gICAgICAgICAgICAgICAgICAgIGxldCBpX3ZhbCA9IEIuZW50cnkoaV9tYXgsIGopO1xuICAgICAgICAgICAgICAgICAgICBCLnNldF9lbnRyeShoLCBqLCBoX3ZhbCk7XG4gICAgICAgICAgICAgICAgICAgIEIuc2V0X2VudHJ5KGlfbWF4LCBqLCBpX3ZhbCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGZvciAobGV0IGkgPSBoICsgMTsgaSA8IHJvd3M7ICsraSkge1xuICAgICAgICAgICAgICAgICAgICBsZXQgZiA9IEIuZW50cnkoaSwgaykgLyBCLmVudHJ5KGgsIGspO1xuICAgICAgICAgICAgICAgICAgICBCLnNldF9lbnRyeShpLCBrLCAwKTtcbiAgICAgICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IGsgKyAxOyBqIDwgMiAqIGNvbHM7ICsraikge1xuICAgICAgICAgICAgICAgICAgICAgICAgQi5zZXRfZW50cnkoaSwgaiwgQi5lbnRyeShpLCBqKSAtIEIuZW50cnkoaCwgaikgKiBmKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBoKys7XG4gICAgICAgICAgICAgICAgaysrO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGxldCBmID0gQi5lbnRyeShyb3csIHJvdyk7XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSByb3c7IGNvbCA8IDIgKiBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIEIuc2V0X2VudHJ5KHJvdywgY29sLCBCLmVudHJ5KHJvdywgY29sKSAvIGYpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgZm9yIChsZXQgcm93ID0gcm93cyAtIDE7IHJvdyA+PSAwOyAtLXJvdykge1xuICAgICAgICAgICAgbGV0IEJfcm93X3JvdyA9IEIuZW50cnkocm93LCByb3cpO1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCByb3c7IGkrKykge1xuICAgICAgICAgICAgICAgIGxldCBCX2lfcm93ID0gQi5lbnRyeShpLCByb3cpO1xuICAgICAgICAgICAgICAgIGxldCBmID0gQl9pX3JvdyAvIEJfcm93X3JvdztcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gaTsgaiA8IDIgKiBjb2xzOyArK2opIHtcbiAgICAgICAgICAgICAgICAgICAgbGV0IEJfaV9qID0gQi5lbnRyeShpLCBqKTtcbiAgICAgICAgICAgICAgICAgICAgbGV0IEJfcm93X2ogPSBCLmVudHJ5KHJvdywgaik7XG4gICAgICAgICAgICAgICAgICAgIEJfaV9qID0gQl9pX2ogLSBCX3Jvd19qICogZjtcbiAgICAgICAgICAgICAgICAgICAgQi5zZXRfZW50cnkoaSwgaiwgQl9pX2opO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIHJldHVybiBuZXcgTWF0cml4KHJvd3MsIGNvbHMsIChpLCBqKSA9PiBCLmVudHJ5KGksIGogKyBjb2xzKSk7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgZG90IHByb2R1Y3QuIElmIHtAbGluayBCfSBpcyBhbiBBcnJheSBvciBGbG9hdDY0QXJyYXkgdGhlbiBhbiBBcnJheSBnZXRzIHJldHVybmVkLiBJZiB7QGxpbmsgQn0gaXMgYSBNYXRyaXggdGhlbiBhIE1hdHJpeCBnZXRzIHJldHVybmVkLlxuICAgICAqIEBwYXJhbSB7KE1hdHJpeHxBcnJheXxGbG9hdDY0QXJyYXkpfSBCIHRoZSByaWdodCBzaWRlXG4gICAgICogQHJldHVybnMgeyhNYXRyaXh8QXJyYXkpfVxuICAgICAqL1xuICAgIGRvdChCKSB7XG4gICAgICAgIGlmIChCIGluc3RhbmNlb2YgTWF0cml4KSB7XG4gICAgICAgICAgICBsZXQgQSA9IHRoaXM7XG4gICAgICAgICAgICBpZiAoQS5zaGFwZVsxXSAhPT0gQi5zaGFwZVswXSkge1xuICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgQS5kb3QoQik6IEEgaXMgYSAke0Euc2hhcGUuam9pbihcIiDiqK8gXCIpfS1NYXRyaXgsIEIgaXMgYSAke0Iuc2hhcGUuam9pbihcIiDiqK8gXCIpfS1NYXRyaXg6IFxuICAgICAgICAgICAgICAgIEEgaGFzICR7QS5zaGFwZVsxXX0gY29scyBhbmQgQiAke0Iuc2hhcGVbMF19IHJvd3MuIFxuICAgICAgICAgICAgICAgIE11c3QgYmUgZXF1YWwhYCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBsZXQgSSA9IEEuc2hhcGVbMV07XG4gICAgICAgICAgICBsZXQgQyA9IG5ldyBNYXRyaXgoQS5zaGFwZVswXSwgQi5zaGFwZVsxXSwgKHJvdywgY29sKSA9PiB7XG4gICAgICAgICAgICAgICAgY29uc3QgQV9pID0gQS5yb3cocm93KTtcbiAgICAgICAgICAgICAgICBjb25zdCBCX2kgPSBCLmNvbChjb2wpO1xuICAgICAgICAgICAgICAgIGxldCBzdW0gPSAwO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgSTsgKytpKSB7XG4gICAgICAgICAgICAgICAgICAgIHN1bSArPSBBX2lbaV0gKiBCX2lbaV07XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHJldHVybiBzdW07XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgICAgIHJldHVybiBDO1xuICAgICAgICB9IGVsc2UgaWYgKEFycmF5LmlzQXJyYXkoQikgfHwgQiBpbnN0YW5jZW9mIEZsb2F0NjRBcnJheSkge1xuICAgICAgICAgICAgbGV0IHJvd3MgPSB0aGlzLl9yb3dzO1xuICAgICAgICAgICAgaWYgKEIubGVuZ3RoICE9PSByb3dzKSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGBBLmRvdChCKTogQSBoYXMgJHtyb3dzfSBjb2xzIGFuZCBCIGhhcyAke0IubGVuZ3RofSByb3dzLiBNdXN0IGJlIGVxdWFsIWApO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgbGV0IEMgPSBuZXcgQXJyYXkocm93cyk7XG4gICAgICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgICAgIENbcm93XSA9IG5ldW1haXJfc3VtKHRoaXMucm93KHJvdykubWFwKChlKSA9PiBlICogQltyb3ddKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gQztcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgQiBtdXN0IGJlIE1hdHJpeCBvciBBcnJheWApO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQ29tcHV0ZXMgdGhlIG91dGVyIHByb2R1Y3QgZnJvbSB7QGxpbmsgdGhpc30gYW5kIHtAbGluayBCfS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gQlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgb3V0ZXIoQikge1xuICAgICAgICBsZXQgQSA9IHRoaXM7XG4gICAgICAgIGxldCBsID0gQS5fZGF0YS5sZW5ndGg7XG4gICAgICAgIGxldCByID0gQi5fZGF0YS5sZW5ndGg7XG4gICAgICAgIGlmIChsICE9IHIpIHJldHVybiB1bmRlZmluZWQ7XG4gICAgICAgIGxldCBDID0gbmV3IE1hdHJpeCgpO1xuICAgICAgICBDLnNoYXBlID0gW1xuICAgICAgICAgICAgbCxcbiAgICAgICAgICAgIGwsXG4gICAgICAgICAgICAoaSwgaikgPT4ge1xuICAgICAgICAgICAgICAgIGlmIChpIDw9IGopIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIEEuX2RhdGFbaV0gKiBCLl9kYXRhW2pdO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiBDLmVudHJ5KGosIGkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0sXG4gICAgICAgIF07XG4gICAgICAgIHJldHVybiBDO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEFwcGVuZHMgbWF0cml4IHtAbGluayBCfSB0byB0aGUgbWF0cml4LlxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBCIC0gbWF0cml4IHRvIGFwcGVuZC5cbiAgICAgKiBAcGFyYW0ge1wiaG9yaXpvbnRhbFwifFwidmVydGljYWxcInxcImRpYWdcIn0gW3R5cGUgPSBcImhvcml6b250YWxcIl0gLSB0eXBlIG9mIGNvbmNhdGVuYXRpb24uXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKiBAZXhhbXBsZVxuICAgICAqXG4gICAgICogbGV0IEEgPSBNYXRyaXguZnJvbShbWzEsIDFdLCBbMSwgMV1dKTsgLy8gMiBieSAyIG1hdHJpeCBmaWxsZWQgd2l0aCBvbmVzLlxuICAgICAqIGxldCBCID0gTWF0cml4LmZyb20oW1syLCAyXSwgWzIsIDJdXSk7IC8vIDIgYnkgMiBtYXRyaXggZmlsbGVkIHdpdGggdHdvcy5cbiAgICAgKlxuICAgICAqIEEuY29uY2F0KEIsIFwiaG9yaXpvbnRhbFwiKTsgLy8gMiBieSA0IG1hdHJpeC4gW1sxLCAxLCAyLCAyXSwgWzEsIDEsIDIsIDJdXVxuICAgICAqIEEuY29uY2F0KEIsIFwidmVydGljYWxcIik7IC8vIDQgYnkgMiBtYXRyaXguIFtbMSwgMV0sIFsxLCAxXSwgWzIsIDJdLCBbMiwgMl1dXG4gICAgICogQS5jb25jYXQoQiwgXCJkaWFnXCIpOyAvLyA0IGJ5IDQgbWF0cml4LiBbWzEsIDEsIDAsIDBdLCBbMSwgMSwgMCwgMF0sIFswLCAwLCAyLCAyXSwgWzAsIDAsIDIsIDJdXVxuICAgICAqL1xuICAgIGNvbmNhdChCLCB0eXBlID0gXCJob3Jpem9udGFsXCIpIHtcbiAgICAgICAgY29uc3QgQSA9IHRoaXM7XG4gICAgICAgIGNvbnN0IFtyb3dzX0EsIGNvbHNfQV0gPSBBLnNoYXBlO1xuICAgICAgICBjb25zdCBbcm93c19CLCBjb2xzX0JdID0gQi5zaGFwZTtcbiAgICAgICAgaWYgKHR5cGUgPT0gXCJob3Jpem9udGFsXCIpIHtcbiAgICAgICAgICAgIGlmIChyb3dzX0EgIT0gcm93c19CKSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGBBLmNvbmNhdChCLCBcImhvcml6b250YWxcIik6IEEgYW5kIEIgbmVlZCBzYW1lIG51bWJlciBvZiByb3dzLCBBIGhhcyAke3Jvd3NfQX0gcm93cywgQiBoYXMgJHtyb3dzX0J9IHJvd3MuYCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjb25zdCBYID0gbmV3IE1hdHJpeChyb3dzX0EsIGNvbHNfQSArIGNvbHNfQiwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgICAgIFguc2V0X2Jsb2NrKDAsIDAsIEEpO1xuICAgICAgICAgICAgWC5zZXRfYmxvY2soMCwgY29sc19BLCBCKTtcbiAgICAgICAgICAgIHJldHVybiBYO1xuICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gXCJ2ZXJ0aWNhbFwiKSB7XG4gICAgICAgICAgICBpZiAoY29sc19BICE9IGNvbHNfQikge1xuICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgQS5jb25jYXQoQiwgXCJ2ZXJ0aWNhbFwiKTogQSBhbmQgQiBuZWVkIHNhbWUgbnVtYmVyIG9mIGNvbHVtbnMsIEEgaGFzICR7Y29sc19BfSBjb2x1bW5zLCBCIGhhcyAke2NvbHNfQn0gY29sdW1ucy5gKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGNvbnN0IFggPSBuZXcgTWF0cml4KHJvd3NfQSArIHJvd3NfQiwgY29sc19BLCBcInplcm9zXCIpO1xuICAgICAgICAgICAgWC5zZXRfYmxvY2soMCwgMCwgQSk7XG4gICAgICAgICAgICBYLnNldF9ibG9jayhyb3dzX0EsIDAsIEIpO1xuICAgICAgICAgICAgcmV0dXJuIFg7XG4gICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PSBcImRpYWdcIikge1xuICAgICAgICAgICAgY29uc3QgWCA9IG5ldyBNYXRyaXgocm93c19BICsgcm93c19CLCBjb2xzX0EgKyBjb2xzX0IsIFwiemVyb3NcIik7XG4gICAgICAgICAgICBYLnNldF9ibG9jaygwLCAwLCBBKTtcbiAgICAgICAgICAgIFguc2V0X2Jsb2NrKHJvd3NfQSwgY29sc19BLCBCKTtcbiAgICAgICAgICAgIHJldHVybiBYO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGB0eXBlIG11c3QgYmUgXCJob3Jpem9udGFsXCIgb3IgXCJ2ZXJ0aWNhbFwiLCBidXQgdHlwZSBpcyAke3R5cGV9IWApO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogV3JpdGVzIHRoZSBlbnRyaWVzIG9mIEIgaW4gQSBhdCBhbiBvZmZzZXQgcG9zaXRpb24gZ2l2ZW4gYnkge0BsaW5rIG9mZnNldF9yb3d9IGFuZCB7QGxpbmsgb2Zmc2V0X2NvbH0uXG4gICAgICogQHBhcmFtIHtpbnR9IG9mZnNldF9yb3dcbiAgICAgKiBAcGFyYW0ge2ludH0gb2Zmc2V0X2NvbFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBCXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBzZXRfYmxvY2sob2Zmc2V0X3Jvdywgb2Zmc2V0X2NvbCwgQikge1xuICAgICAgICBsZXQgW3Jvd3MsIGNvbHNdID0gQi5zaGFwZTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGlmIChyb3cgPiB0aGlzLl9yb3dzKSB7XG4gICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIGlmIChjb2wgPiB0aGlzLl9jb2xzKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB0aGlzLnNldF9lbnRyeShyb3cgKyBvZmZzZXRfcm93LCBjb2wgKyBvZmZzZXRfY29sLCBCLmVudHJ5KHJvdywgY29sKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogRXh0cmFjdHMgdGhlIGVudHJpZXMgZnJvbSB0aGUge0BsaW5rIHN0YXJ0X3Jvd308c3VwPnRoPC9zdXA+IHJvdyB0byB0aGUge0BsaW5rIGVuZF9yb3d9PHN1cD50aDwvc3VwPiByb3csIHRoZSB7QGxpbmsgc3RhcnRfY29sfTxzdXA+dGg8L3N1cD4gY29sdW1uIHRvIHRoZSB7QGxpbmsgZW5kX2NvbH08c3VwPnRoPC9zdXA+IGNvbHVtbiBvZiB0aGUgbWF0cml4LlxuICAgICAqIElmIHtAbGluayBlbmRfcm93fSBvciB7QGxpbmsgZW5kX2NvbH0gaXMgZW1wdHksIHRoZSByZXNwZWN0aXZlIHZhbHVlIGlzIHNldCB0byB7QGxpbmsgdGhpcy5yb3dzfSBvciB7QGxpbmsgdGhpcy5jb2xzfS5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gc3RhcnRfcm93XG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IHN0YXJ0X2NvbFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbZW5kX3JvdyA9IG51bGxdXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtlbmRfY29sID0gbnVsbF1cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fSBSZXR1cm5zIGEgZW5kX3JvdyAtIHN0YXJ0X3JvdyB0aW1lcyBlbmRfY29sIC0gc3RhcnRfY29sIG1hdHJpeCwgd2l0aCByZXNwZWN0aXZlIGVudHJpZXMgZnJvbSB0aGUgbWF0cml4LlxuICAgICAqIEBleGFtcGxlXG4gICAgICpcbiAgICAgKiBsZXQgQSA9IE1hdHJpeC5mcm9tKFtbMSwgMiwgM10sIFs0LCA1LCA2XSwgWzcsIDgsIDldXSk7IC8vIGEgMyBieSAzIG1hdHJpeC5cbiAgICAgKlxuICAgICAqIEEuZ2V0X2Jsb2NrKDEsIDEpOyAvLyBbWzUsIDZdLCBbOCwgOV1dXG4gICAgICogQS5nZXRfYmxvY2soMCwgMCwgMSwgMSk7IC8vIFtbMV1dXG4gICAgICogQS5nZXRfYmxvY2soMSwgMSwgMiwgMik7IC8vIFtbNV1dXG4gICAgICogQS5nZXRfYmxvY2soMCwgMCwgMiwgMik7IC8vIFtbMSwgMl0sIFs0LCA1XV1cbiAgICAgKi9cbiAgICBnZXRfYmxvY2soc3RhcnRfcm93LCBzdGFydF9jb2wsIGVuZF9yb3cgPSBudWxsLCBlbmRfY29sID0gbnVsbCkge1xuICAgICAgICBjb25zdCBbcm93cywgY29sc10gPSB0aGlzLnNoYXBlO1xuICAgICAgICBlbmRfcm93ID0gZW5kX3JvdyA/PyByb3dzO1xuICAgICAgICBlbmRfY29sID0gZW5kX2NvbCA/PyBjb2xzO1xuICAgICAgICBpZiAoZW5kX3JvdyA8PSBzdGFydF9yb3cgfHwgZW5kX2NvbCA8PSBzdGFydF9jb2wpIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgXG4gICAgICAgICAgICAgICAgZW5kX3JvdyBtdXN0IGJlIGdyZWF0ZXIgdGhhbiBzdGFydF9yb3csIGFuZCBcbiAgICAgICAgICAgICAgICBlbmRfY29sIG11c3QgYmUgZ3JlYXRlciB0aGFuIHN0YXJ0X2NvbCwgYnV0XG4gICAgICAgICAgICAgICAgZW5kX3JvdyA9ICR7ZW5kX3Jvd30sIHN0YXJ0X3JvdyA9ICR7c3RhcnRfcm93fSwgZW5kX2NvbCA9ICR7ZW5kX2NvbH0sIGFuZCBzdGFydF9jb2wgPSAke3N0YXJ0X2NvbH0hYCk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3QgWCA9IG5ldyBNYXRyaXgoZW5kX3JvdyAtIHN0YXJ0X3JvdywgZW5kX2NvbCAtIHN0YXJ0X2NvbCwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gc3RhcnRfcm93LCBuZXdfcm93ID0gMDsgcm93IDwgZW5kX3JvdzsgKytyb3csICsrbmV3X3Jvdykge1xuICAgICAgICAgICAgZm9yIChsZXQgY29sID0gc3RhcnRfY29sLCBuZXdfY29sID0gMDsgY29sIDwgZW5kX2NvbDsgKytjb2wsICsrbmV3X2NvbCkge1xuICAgICAgICAgICAgICAgIFguc2V0X2VudHJ5KG5ld19yb3csIG5ld19jb2wsIHRoaXMuZW50cnkocm93LCBjb2wpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gWDtcbiAgICAgICAgLy9yZXR1cm4gbmV3IE1hdHJpeChlbmRfcm93IC0gc3RhcnRfcm93LCBlbmRfY29sIC0gc3RhcnRfY29sLCAoaSwgaikgPT4gdGhpcy5lbnRyeShpICsgc3RhcnRfcm93LCBqICsgc3RhcnRfY29sKSk7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyBhIG5ldyBhcnJheSBnYXRoZXJpbmcgZW50cmllcyBkZWZpbmVkIGJ5IHRoZSBpbmRpY2VzIGdpdmVuIGJ5IGFyZ3VtZW50LlxuICAgICAqIEBwYXJhbSB7QXJyYXk8TnVtYmVyPn0gcm93X2luZGljZXMgLSBBcnJheSBjb25zaXN0cyBvZiBpbmRpY2VzIG9mIHJvd3MgZm9yIGdhdGhlcmluZyBlbnRyaWVzIG9mIHRoaXMgbWF0cml4XG4gICAgICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBjb2xfaW5kaWNlcyAgLSBBcnJheSBjb25zaXN0cyBvZiBpbmRpY2VzIG9mIGNvbHMgZm9yIGdhdGhlcmluZyBlbnRyaWVzIG9mIHRoaXMgbWF0cml4XG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBnYXRoZXIocm93X2luZGljZXMsIGNvbF9pbmRpY2VzKSB7XG4gICAgICAgIGNvbnN0IE4gPSByb3dfaW5kaWNlcy5sZW5ndGg7XG4gICAgICAgIGNvbnN0IEQgPSBjb2xfaW5kaWNlcy5sZW5ndGg7XG5cbiAgICAgICAgY29uc3QgUiA9IG5ldyBNYXRyaXgoTiwgRCk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCByb3dfaW5kZXggPSByb3dfaW5kaWNlc1tpXTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgY29sX2luZGV4ID0gY29sX2luZGljZXNbal07XG4gICAgICAgICAgICAgICAgUi5zZXRfZW50cnkoaSwgaiwgdGhpcy5lbnRyeShyb3dfaW5kZXgsIGNvbF9pbmRleCkpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIFI7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQXBwbGllcyBhIGZ1bmN0aW9uIHRvIGVhY2ggZW50cnkgb2YgdGhlIG1hdHJpeC5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7ZnVuY3Rpb259IGYgZnVuY3Rpb24gdGFrZXMgMiBwYXJhbWV0ZXJzLCB0aGUgdmFsdWUgb2YgdGhlIGFjdHVhbCBlbnRyeSBhbmQgYSB2YWx1ZSBnaXZlbiBieSB0aGUgZnVuY3Rpb24ge0BsaW5rIHZ9LiBUaGUgcmVzdWx0IG9mIHtAbGluayBmfSBnZXRzIHdyaXRlbiB0byB0aGUgTWF0cml4LlxuICAgICAqIEBwYXJhbSB7ZnVuY3Rpb259IHYgZnVuY3Rpb24gdGFrZXMgMiBwYXJhbWV0ZXJzIGZvciByb3cgYW5kIGNvbCwgYW5kIHJldHVybnMgYSB2YWx1ZSB3aXRjaCBzaG91bGQgYmUgYXBwbGllZCB0byB0aGUgY29sdGggZW50cnkgb2YgdGhlIHJvd3RoIHJvdyBvZiB0aGUgbWF0cml4LlxuICAgICAqL1xuICAgIF9hcHBseV9hcnJheShmLCB2KSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLnZhbHVlcztcbiAgICAgICAgY29uc3QgW3Jvd3MsIGNvbHNdID0gdGhpcy5zaGFwZTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGNvbnN0IG9mZnNldCA9IHJvdyAqIGNvbHM7XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IGkgPSBvZmZzZXQgKyBjb2w7XG4gICAgICAgICAgICAgICAgZGF0YVtpXSA9IGYoZGF0YVtpXSwgdihyb3csIGNvbCkpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIF9hcHBseV9yb3d3aXNlX2FycmF5KHZhbHVlcywgZikge1xuICAgICAgICByZXR1cm4gdGhpcy5fYXBwbHlfYXJyYXkoZiwgKF8sIGopID0+IHZhbHVlc1tqXSk7XG4gICAgfVxuXG4gICAgX2FwcGx5X2NvbHdpc2VfYXJyYXkodmFsdWVzLCBmKSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLnZhbHVlcztcbiAgICAgICAgY29uc3QgW3Jvd3MsIGNvbHNdID0gdGhpcy5zaGFwZTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGNvbnN0IG9mZnNldCA9IHJvdyAqIGNvbHM7XG4gICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IGkgPSBvZmZzZXQgKyBjb2w7XG4gICAgICAgICAgICAgICAgZGF0YVtpXSA9IGYoZGF0YVtpXSwgdmFsdWVzW3Jvd10pO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIF9hcHBseSh2YWx1ZSwgZikge1xuICAgICAgICBsZXQgZGF0YSA9IHRoaXMudmFsdWVzO1xuICAgICAgICBpZiAodmFsdWUgaW5zdGFuY2VvZiBNYXRyaXgpIHtcbiAgICAgICAgICAgIGxldCBbdmFsdWVfcm93cywgdmFsdWVfY29sc10gPSB2YWx1ZS5zaGFwZTtcbiAgICAgICAgICAgIGxldCBbcm93cywgY29sc10gPSB0aGlzLnNoYXBlO1xuICAgICAgICAgICAgaWYgKHZhbHVlX3Jvd3MgPT09IDEpIHtcbiAgICAgICAgICAgICAgICBpZiAoY29scyAhPT0gdmFsdWVfY29scykge1xuICAgICAgICAgICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYGNvbHMgIT09IHZhbHVlX2NvbHNgKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgICAgICAgICAgZm9yIChsZXQgY29sID0gMDsgY29sIDwgY29sczsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRhdGFbcm93ICogY29scyArIGNvbF0gPSBmKGRhdGFbcm93ICogY29scyArIGNvbF0sIHZhbHVlLmVudHJ5KDAsIGNvbCkpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfSBlbHNlIGlmICh2YWx1ZV9jb2xzID09PSAxKSB7XG4gICAgICAgICAgICAgICAgaWYgKHJvd3MgIT09IHZhbHVlX3Jvd3MpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGByb3dzICE9PSB2YWx1ZV9yb3dzYCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3M7ICsrcm93KSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXRhW3JvdyAqIGNvbHMgKyBjb2xdID0gZihkYXRhW3JvdyAqIGNvbHMgKyBjb2xdLCB2YWx1ZS5lbnRyeShyb3csIDApKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0gZWxzZSBpZiAocm93cyA9PSB2YWx1ZV9yb3dzICYmIGNvbHMgPT0gdmFsdWVfY29scykge1xuICAgICAgICAgICAgICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3M7ICsrcm93KSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXRhW3JvdyAqIGNvbHMgKyBjb2xdID0gZihkYXRhW3JvdyAqIGNvbHMgKyBjb2xdLCB2YWx1ZS5lbnRyeShyb3csIGNvbCkpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYGVycm9yYCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSBpZiAoQXJyYXkuaXNBcnJheSh2YWx1ZSkpIHtcbiAgICAgICAgICAgIGxldCByb3dzID0gdGhpcy5fcm93cztcbiAgICAgICAgICAgIGxldCBjb2xzID0gdGhpcy5fY29scztcbiAgICAgICAgICAgIGlmICh2YWx1ZS5sZW5ndGggPT09IHJvd3MpIHtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgICAgICAgICBmb3IgKGxldCBjb2wgPSAwOyBjb2wgPCBjb2xzOyArK2NvbCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgZGF0YVtyb3cgKiBjb2xzICsgY29sXSA9IGYoZGF0YVtyb3cgKiBjb2xzICsgY29sXSwgdmFsdWVbcm93XSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2UgaWYgKHZhbHVlLmxlbmd0aCA9PT0gY29scykge1xuICAgICAgICAgICAgICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3M7ICsrcm93KSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXRhW3JvdyAqIGNvbHMgKyBjb2xdID0gZihkYXRhW3JvdyAqIGNvbHMgKyBjb2xdLCB2YWx1ZVtjb2xdKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGBlcnJvcmApO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDAsIG4gPSB0aGlzLl9yb3dzICogdGhpcy5fY29sczsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgICAgIGRhdGFbaV0gPSBmKGRhdGFbaV0sIHZhbHVlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDbG9uZXMgdGhlIE1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqL1xuICAgIGNsb25lKCkge1xuICAgICAgICBsZXQgQiA9IG5ldyBNYXRyaXgoKTtcbiAgICAgICAgQi5fcm93cyA9IHRoaXMuX3Jvd3M7XG4gICAgICAgIEIuX2NvbHMgPSB0aGlzLl9jb2xzO1xuICAgICAgICBCLl9kYXRhID0gdGhpcy52YWx1ZXMuc2xpY2UoMCk7XG4gICAgICAgIHJldHVybiBCO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEVudHJ5d2lzZSBtdWx0aXBsaWNhdGlvbiB3aXRoIHtAbGluayB2YWx1ZX0uXG4gICAgICogQHBhcmFtIHtNYXRyaXh8QXJyYXl8TnVtYmVyfSB2YWx1ZVxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICogQGV4YW1wbGVcbiAgICAgKlxuICAgICAqIGxldCBBID0gTWF0cml4LmZyb20oW1sxLCAyXSwgWzMsIDRdXSk7IC8vIGEgMiBieSAyIG1hdHJpeC5cbiAgICAgKiBsZXQgQiA9IEEuY2xvbmUoKTsgLy8gQiA9PSBBO1xuICAgICAqXG4gICAgICogQS5tdWx0KDIpOyAvLyBbWzIsIDRdLCBbNiwgOF1dO1xuICAgICAqIEEubXVsdChCKTsgLy8gW1sxLCA0XSwgWzksIDE2XV07XG4gICAgICovXG4gICAgbXVsdCh2YWx1ZSkge1xuICAgICAgICByZXR1cm4gdGhpcy5jbG9uZSgpLl9hcHBseSh2YWx1ZSwgKGEsIGIpID0+IGEgKiBiKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBFbnRyeXdpc2UgZGl2aXNpb24gd2l0aCB7QGxpbmsgdmFsdWV9LlxuICAgICAqIEBwYXJhbSB7TWF0cml4fEFycmF5fE51bWJlcn0gdmFsdWVcbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqIEBleGFtcGxlXG4gICAgICpcbiAgICAgKiBsZXQgQSA9IE1hdHJpeC5mcm9tKFtbMSwgMl0sIFszLCA0XV0pOyAvLyBhIDIgYnkgMiBtYXRyaXguXG4gICAgICogbGV0IEIgPSBBLmNsb25lKCk7IC8vIEIgPT0gQTtcbiAgICAgKlxuICAgICAqIEEuZGl2aWRlKDIpOyAvLyBbWzAuNSwgMV0sIFsxLjUsIDJdXTtcbiAgICAgKiBBLmRpdmlkZShCKTsgLy8gW1sxLCAxXSwgWzEsIDFdXTtcbiAgICAgKi9cbiAgICBkaXZpZGUodmFsdWUpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuY2xvbmUoKS5fYXBwbHkodmFsdWUsIChhLCBiKSA9PiBhIC8gYik7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogRW50cnl3aXNlIGFkZGl0aW9uIHdpdGgge0BsaW5rIHZhbHVlfS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeHxBcnJheXxOdW1iZXJ9IHZhbHVlXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKiBAZXhhbXBsZVxuICAgICAqXG4gICAgICogbGV0IEEgPSBNYXRyaXguZnJvbShbWzEsIDJdLCBbMywgNF1dKTsgLy8gYSAyIGJ5IDIgbWF0cml4LlxuICAgICAqIGxldCBCID0gQS5jbG9uZSgpOyAvLyBCID09IEE7XG4gICAgICpcbiAgICAgKiBBLmFkZCgyKTsgLy8gW1szLCA0XSwgWzUsIDZdXTtcbiAgICAgKiBBLmFkZChCKTsgLy8gW1syLCA0XSwgWzYsIDhdXTtcbiAgICAgKi9cbiAgICBhZGQodmFsdWUpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuY2xvbmUoKS5fYXBwbHkodmFsdWUsIChhLCBiKSA9PiBhICsgYik7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogRW50cnl3aXNlIHN1YnRyYWN0aW9uIHdpdGgge0BsaW5rIHZhbHVlfS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeHxBcnJheXxOdW1iZXJ9IHZhbHVlXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKiBAZXhhbXBsZVxuICAgICAqXG4gICAgICogbGV0IEEgPSBNYXRyaXguZnJvbShbWzEsIDJdLCBbMywgNF1dKTsgLy8gYSAyIGJ5IDIgbWF0cml4LlxuICAgICAqIGxldCBCID0gQS5jbG9uZSgpOyAvLyBCID09IEE7XG4gICAgICpcbiAgICAgKiBBLnN1YigyKTsgLy8gW1stMSwgMF0sIFsxLCAyXV07XG4gICAgICogQS5zdWIoQik7IC8vIFtbMCwgMF0sIFswLCAwXV07XG4gICAgICovXG4gICAgc3ViKHZhbHVlKSB7XG4gICAgICAgIHJldHVybiB0aGlzLmNsb25lKCkuX2FwcGx5KHZhbHVlLCAoYSwgYikgPT4gYSAtIGIpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIG51bWJlciBvZiByb3dzIGFuZCBjb2x1bW5zIG9mIHRoZSBNYXRyaXguXG4gICAgICogQHJldHVybnMge0FycmF5fSBBbiBBcnJheSBpbiB0aGUgZm9ybSBbcm93cywgY29sdW1uc10uXG4gICAgICovXG4gICAgZ2V0IHNoYXBlKCkge1xuICAgICAgICByZXR1cm4gW3RoaXMuX3Jvd3MsIHRoaXMuX2NvbHNdO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIG1hdHJpeCBpbiB0aGUgZ2l2ZW4gc2hhcGUgd2l0aCB0aGUgZ2l2ZW4gZnVuY3Rpb24gd2hpY2ggcmV0dXJucyB2YWx1ZXMgZm9yIHRoZSBlbnRyaWVzIG9mIHRoZSBtYXRyaXguXG4gICAgICogQHBhcmFtIHtBcnJheX0gcGFyYW1ldGVyIC0gdGFrZXMgYW4gQXJyYXkgaW4gdGhlIGZvcm0gW3Jvd3MsIGNvbHMsIHZhbHVlXSwgd2hlcmUgcm93cyBhbmQgY29scyBhcmUgdGhlIG51bWJlciBvZiByb3dzIGFuZCBjb2x1bW5zIG9mIHRoZSBtYXRyaXgsIGFuZCB2YWx1ZSBpcyBhIGZ1bmN0aW9uIHdoaWNoIHRha2VzIHR3byBwYXJhbWV0ZXJzIChyb3cgYW5kIGNvbCkgd2hpY2ggaGFzIHRvIHJldHVybiBhIHZhbHVlIGZvciB0aGUgY29sdGggZW50cnkgb2YgdGhlIHJvd3RoIHJvdy5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqL1xuICAgIHNldCBzaGFwZShbcm93cywgY29scywgdmFsdWUgPSAoKSA9PiAwXSkge1xuICAgICAgICB0aGlzLl9yb3dzID0gcm93cztcbiAgICAgICAgdGhpcy5fY29scyA9IGNvbHM7XG4gICAgICAgIHRoaXMuX2RhdGEgPSBuZXcgRmxvYXQ2NEFycmF5KHJvd3MgKiBjb2xzKTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5fZGF0YVtyb3cgKiBjb2xzICsgY29sXSA9IHZhbHVlKHJvdywgY29sKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRoZSBNYXRyaXggYXMgYSBBcnJheSBvZiBGbG9hdDY0QXJyYXlzLlxuICAgICAqIEByZXR1cm5zIHtBcnJheTxGbG9hdDY0QXJyYXk+fVxuICAgICAqL1xuICAgIGdldCB0bzJkQXJyYXkoKSB7XG4gICAgICAgIGNvbnN0IHJlc3VsdCA9IFtdO1xuICAgICAgICBmb3IgKGNvbnN0IHJvdyBvZiB0aGlzLml0ZXJhdGVfcm93cygpKSB7XG4gICAgICAgICAgICByZXN1bHQucHVzaChyb3cpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgTWF0cml4IGFzIGEgQXJyYXkgb2YgQXJyYXlzLlxuICAgICAqIEByZXR1cm5zIHtBcnJheTxBcnJheT59XG4gICAgICovXG4gICAgZ2V0IGFzQXJyYXkoKSB7XG4gICAgICAgIGNvbnN0IHJlc3VsdCA9IFtdO1xuICAgICAgICBmb3IgKGNvbnN0IHJvdyBvZiB0aGlzLml0ZXJhdGVfcm93cygpKSB7XG4gICAgICAgICAgICByZXN1bHQucHVzaChBcnJheS5mcm9tKHJvdykpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgZGlhZ29uYWwgb2YgdGhlIE1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7RmxvYXQ2NEFycmF5fVxuICAgICAqL1xuICAgIGdldCBkaWFnKCkge1xuICAgICAgICBjb25zdCByb3dzID0gdGhpcy5fcm93cztcbiAgICAgICAgY29uc3QgY29scyA9IHRoaXMuX2NvbHM7XG4gICAgICAgIGNvbnN0IG1pbl9yb3dfY29sID0gTWF0aC5taW4ocm93cywgY29scyk7XG4gICAgICAgIGxldCByZXN1bHQgPSBuZXcgRmxvYXQ2NEFycmF5KG1pbl9yb3dfY29sKTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBtaW5fcm93X2NvbDsgKytpKSB7XG4gICAgICAgICAgICByZXN1bHRbaV0gPSB0aGlzLmVudHJ5KGksIGkpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgbWVhbiBvZiBhbGwgZW50cmllcyBvZiB0aGUgTWF0cml4LlxuICAgICAqIEByZXR1cm5zIHtOdW1iZXJ9XG4gICAgICovXG4gICAgZ2V0IG1lYW4oKSB7XG4gICAgICAgIGNvbnN0IHN1bSA9IHRoaXMuc3VtO1xuICAgICAgICBjb25zdCBuID0gdGhpcy5fcm93cyAqIHRoaXMuX2NvbHM7XG4gICAgICAgIHJldHVybiBzdW0gLyBuO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIHN1bSBvb2YgYWxsIGVudHJpZXMgb2YgdGhlIE1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7TnVtYmVyfVxuICAgICAqL1xuICAgIGdldCBzdW0oKSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLnZhbHVlcztcbiAgICAgICAgcmV0dXJuIG5ldW1haXJfc3VtKGRhdGEpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIHN1bSBvb2YgYWxsIGVudHJpZXMgb2YgdGhlIE1hdHJpeC5cbiAgICAgKiBAcmV0dXJucyB7RmxvYXQ2NEFycmF5fVxuICAgICAqL1xuICAgIGdldCB2YWx1ZXMoKSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLl9kYXRhO1xuICAgICAgICByZXR1cm4gZGF0YTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRoZSBtZWFuIG9mIGVhY2ggcm93IG9mIHRoZSBtYXRyaXguXG4gICAgICogQHJldHVybnMge0Zsb2F0NjRBcnJheX1cbiAgICAgKi9cbiAgICBnZXQgbWVhblJvd3MoKSB7XG4gICAgICAgIGNvbnN0IGRhdGEgPSB0aGlzLnZhbHVlcztcbiAgICAgICAgY29uc3Qgcm93cyA9IHRoaXMuX3Jvd3M7XG4gICAgICAgIGNvbnN0IGNvbHMgPSB0aGlzLl9jb2xzO1xuICAgICAgICBjb25zdCByZXN1bHQgPSBGbG9hdDY0QXJyYXkuZnJvbSh7IGxlbmd0aDogcm93cyB9KTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIHJlc3VsdFtyb3ddID0gMDtcbiAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICAgICAgcmVzdWx0W3Jvd10gKz0gZGF0YVtyb3cgKiBjb2xzICsgY29sXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJlc3VsdFtyb3ddIC89IGNvbHM7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG5cbiAgICAvKiogUmV0dXJucyB0aGUgbWVhbiBvZiBlYWNoIGNvbHVtbiBvZiB0aGUgbWF0cml4LlxuICAgICAqIEByZXR1cm5zIHtGbG9hdDY0QXJyYXl9XG4gICAgICovXG4gICAgZ2V0IG1lYW5Db2xzKCkge1xuICAgICAgICBjb25zdCBkYXRhID0gdGhpcy52YWx1ZXM7XG4gICAgICAgIGNvbnN0IHJvd3MgPSB0aGlzLl9yb3dzO1xuICAgICAgICBjb25zdCBjb2xzID0gdGhpcy5fY29scztcbiAgICAgICAgY29uc3QgcmVzdWx0ID0gRmxvYXQ2NEFycmF5LmZyb20oeyBsZW5ndGg6IGNvbHMgfSk7XG4gICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHM7ICsrY29sKSB7XG4gICAgICAgICAgICByZXN1bHRbY29sXSA9IDA7XG4gICAgICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgICAgIHJlc3VsdFtjb2xdICs9IGRhdGFbcm93ICogY29scyArIGNvbF07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXN1bHRbY29sXSAvPSByb3dzO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogU29sdmVzIHRoZSBlcXVhdGlvbiB7QGxpbmsgQX14ID0ge0BsaW5rIGJ9IHVzaW5nIHRoZSBjb25qdWdhdGUgZ3JhZGllbnQgbWV0aG9kLiBSZXR1cm5zIHRoZSByZXN1bHQgeC5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gQSAtIE1hdHJpeFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBiIC0gTWF0cml4XG4gICAgICogQHBhcmFtIHtSYW5kb21pemVyfSBbcmFuZG9taXplcj1udWxsXVxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbdG9sPTFlLTNdXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBzdGF0aWMgc29sdmVfQ0coQSwgYiwgcmFuZG9taXplciwgdG9sID0gMWUtMykge1xuICAgICAgICBpZiAocmFuZG9taXplciA9PT0gbnVsbCkge1xuICAgICAgICAgICAgcmFuZG9taXplciA9IG5ldyBSYW5kb21pemVyKCk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3Qgcm93cyA9IEEuc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IGNvbHMgPSBiLnNoYXBlWzFdO1xuICAgICAgICBsZXQgcmVzdWx0ID0gbmV3IE1hdHJpeChyb3dzLCAwKTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBjb2xzOyArK2kpIHtcbiAgICAgICAgICAgIGNvbnN0IGJfaSA9IE1hdHJpeC5mcm9tKGIuY29sKGkpKS5UO1xuICAgICAgICAgICAgbGV0IHggPSBuZXcgTWF0cml4KHJvd3MsIDEsICgpID0+IHJhbmRvbWl6ZXIucmFuZG9tKTtcbiAgICAgICAgICAgIGxldCByID0gYl9pLnN1YihBLmRvdCh4KSk7XG4gICAgICAgICAgICBsZXQgZCA9IHIuY2xvbmUoKTtcbiAgICAgICAgICAgIGRvIHtcbiAgICAgICAgICAgICAgICBjb25zdCB6ID0gQS5kb3QoZCk7XG4gICAgICAgICAgICAgICAgY29uc3QgYWxwaGEgPSByLlQuZG90KHIpLmVudHJ5KDAsIDApIC8gZC5ULmRvdCh6KS5lbnRyeSgwLCAwKTtcbiAgICAgICAgICAgICAgICB4ID0geC5hZGQoZC5tdWx0KGFscGhhKSk7XG4gICAgICAgICAgICAgICAgY29uc3Qgcl9uZXh0ID0gci5zdWIoei5tdWx0KGFscGhhKSk7XG4gICAgICAgICAgICAgICAgY29uc3QgYmV0YSA9IHJfbmV4dC5ULmRvdChyX25leHQpLmVudHJ5KDAsIDApIC8gci5ULmRvdChyKS5lbnRyeSgwLCAwKTtcbiAgICAgICAgICAgICAgICBkID0gcl9uZXh0LmFkZChkLm11bHQoYmV0YSkpO1xuICAgICAgICAgICAgICAgIHIgPSByX25leHQ7XG4gICAgICAgICAgICB9IHdoaWxlIChNYXRoLmFicyhyLm1lYW4pID4gdG9sKTtcbiAgICAgICAgICAgIHJlc3VsdCA9IHJlc3VsdC5jb25jYXQoeCwgXCJob3Jpem9udGFsXCIpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogU29sdmVzIHRoZSBlcXVhdGlvbiB7QGxpbmsgQX14ID0ge0BsaW5rIGJ9LiBSZXR1cm5zIHRoZSByZXN1bHQgeC5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gQSAtIE1hdHJpeCBvciBMVSBEZWNvbXBvc2l0aW9uXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IGIgLSBNYXRyaXhcbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fVxuICAgICAqL1xuICAgIHN0YXRpYyBzb2x2ZShBLCBiKSB7XG4gICAgICAgIGxldCB7IEw6IEwsIFU6IFUgfSA9IFwiTFwiIGluIEEgJiYgXCJVXCIgaW4gQSA/IEEgOiBNYXRyaXguTFUoQSk7XG4gICAgICAgIGxldCByb3dzID0gTC5zaGFwZVswXTtcbiAgICAgICAgbGV0IHggPSBiLmNsb25lKCk7XG5cbiAgICAgICAgLy8gZm9yd2FyZFxuICAgICAgICBmb3IgKGxldCByb3cgPSAwOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgZm9yIChsZXQgY29sID0gMDsgY29sIDwgcm93IC0gMTsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICB4LnNldF9lbnRyeSgwLCByb3csIHguZW50cnkoMCwgcm93KSAtIEwuZW50cnkocm93LCBjb2wpICogeC5lbnRyeSgxLCBjb2wpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHguc2V0X2VudHJ5KDAsIHJvdywgeC5lbnRyeSgwLCByb3cpIC8gTC5lbnRyeShyb3csIHJvdykpO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gYmFja3dhcmRcbiAgICAgICAgZm9yIChsZXQgcm93ID0gcm93cyAtIDE7IHJvdyA+PSAwOyAtLXJvdykge1xuICAgICAgICAgICAgZm9yIChsZXQgY29sID0gcm93cyAtIDE7IGNvbCA+IHJvdzsgLS1jb2wpIHtcbiAgICAgICAgICAgICAgICB4LnNldF9lbnRyeSgwLCByb3csIHguZW50cnkoMCwgcm93KSAtIFUuZW50cnkocm93LCBjb2wpICogeC5lbnRyeSgwLCBjb2wpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHguc2V0X2VudHJ5KDAsIHJvdywgeC5lbnRyeSgwLCByb3cpIC8gVS5lbnRyeShyb3csIHJvdykpO1xuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIHg7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICoge0BsaW5rIEx9e0BsaW5rIFV9IGRlY29tcG9zaXRpb24gb2YgdGhlIE1hdHJpeCB7QGxpbmsgQX0uIENyZWF0ZXMgdHdvIG1hdHJpY2VzLCBzbyB0aGF0IHRoZSBkb3QgcHJvZHVjdCBMVSBlcXVhbHMgQS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gQVxuICAgICAqIEByZXR1cm5zIHt7TDogTWF0cml4LCBVOiBNYXRyaXh9fSByZXN1bHQgLSBSZXR1cm5zIHRoZSBsZWZ0IHRyaWFuZ2xlIG1hdHJpeCB7QGxpbmsgTH0gYW5kIHRoZSB1cHBlciB0cmlhbmdsZSBtYXRyaXgge0BsaW5rIFV9LlxuICAgICAqL1xuICAgIHN0YXRpYyBMVShBKSB7XG4gICAgICAgIGNvbnN0IHJvd3MgPSBBLnNoYXBlWzBdO1xuICAgICAgICBjb25zdCBMID0gbmV3IE1hdHJpeChyb3dzLCByb3dzLCBcInplcm9zXCIpO1xuICAgICAgICBjb25zdCBVID0gbmV3IE1hdHJpeChyb3dzLCByb3dzLCBcImlkZW50aXR5XCIpO1xuXG4gICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgcm93czsgKytqKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gajsgaSA8IHJvd3M7ICsraSkge1xuICAgICAgICAgICAgICAgIGxldCBzdW0gPSAwO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGsgPSAwOyBrIDwgajsgKytrKSB7XG4gICAgICAgICAgICAgICAgICAgIHN1bSArPSBMLmVudHJ5KGksIGspICogVS5lbnRyeShrLCBqKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgTC5zZXRfZW50cnkoaSwgaiwgQS5lbnRyeShpLCBqKSAtIHN1bSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gajsgaSA8IHJvd3M7ICsraSkge1xuICAgICAgICAgICAgICAgIGlmIChMLmVudHJ5KGosIGopID09PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGxldCBzdW0gPSAwO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGsgPSAwOyBrIDwgajsgKytrKSB7XG4gICAgICAgICAgICAgICAgICAgIHN1bSArPSBMLmVudHJ5KGosIGspICogVS5lbnRyeShrLCBpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgVS5zZXRfZW50cnkoaiwgaSwgKEEuZW50cnkoaiwgaSkgLSBzdW0pIC8gTC5lbnRyeShqLCBqKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4geyBMOiBMLCBVOiBVIH07XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQ29tcHV0ZXMgdGhlIGRldGVybWluYW50ZSBvZiB7QGxpbmsgQX0sIGJ5IHVzaW5nIHRoZSBMVSBkZWNvbXBvc2l0aW9uIG9mIHtAbGluayBBfS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gQVxuICAgICAqIEByZXR1cm5zIHtOdW1iZXJ9IGRldCAtIFJldHVybnMgdGhlIGRldGVybWluYXRlIG9mIHRoZSBNYXRyaXgge0BsaW5rIEF9LlxuICAgICAqL1xuICAgIHN0YXRpYyBkZXQoQSkge1xuICAgICAgICBjb25zdCByb3dzID0gQS5zaGFwZVswXTtcbiAgICAgICAgY29uc3QgeyBMLCBVIH0gPSBNYXRyaXguTFUoQSk7XG4gICAgICAgIGNvbnN0IExfZGlhZyA9IEwuZGlhZztcbiAgICAgICAgY29uc3QgVV9kaWFnID0gVS5kaWFnO1xuICAgICAgICBsZXQgZGV0ID0gTF9kaWFnWzBdICogVV9kaWFnWzBdO1xuICAgICAgICBmb3IgKGxldCByb3cgPSAxOyByb3cgPCByb3dzOyArK3Jvdykge1xuICAgICAgICAgICAgZGV0ICo9IExfZGlhZ1tyb3ddICogVV9kaWFnW3Jvd107XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIGRldDtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyB0aGUge0BsaW5rIGt9IGNvbXBvbmVudHMgb2YgdGhlIFNWRCBkZWNvbXBvc2l0aW9uIG9mIHRoZSBtYXRyaXgge0BsaW5rIE19XG4gICAgICogQHBhcmFtIHtNYXRyaXh9IE1cbiAgICAgKiBAcGFyYW0ge2ludH0gW2s9Ml1cbiAgICAgKiBAcmV0dXJucyB7e1U6IE1hdHJpeCwgU2lnbWE6IE1hdHJpeCwgVjogTWF0cml4fX1cbiAgICAgKi9cbiAgICBzdGF0aWMgU1ZEKE0sIGsgPSAyKSB7XG4gICAgICAgIGNvbnN0IE1UID0gTS5UO1xuICAgICAgICBsZXQgTXRNID0gTVQuZG90KE0pO1xuICAgICAgICBsZXQgTU10ID0gTS5kb3QoTVQpO1xuICAgICAgICBsZXQgeyBlaWdlbnZlY3RvcnM6IFYsIGVpZ2VudmFsdWVzOiBTaWdtYSB9ID0gc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uKE10TSwgayk7XG4gICAgICAgIGxldCB7IGVpZ2VudmVjdG9yczogVSB9ID0gc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uKE1NdCwgayk7XG4gICAgICAgIHJldHVybiB7IFU6IFUsIFNpZ21hOiBTaWdtYS5tYXAoKHNpZ21hKSA9PiBNYXRoLnNxcnQoc2lnbWEpKSwgVjogViB9O1xuXG4gICAgICAgIC8vQWxnb3JpdGhtIDFhOiBIb3VzZWhvbGRlciByZWR1Y3Rpb24gdG8gYmlkaWFnb25hbCBmb3JtOlxuICAgICAgICAvKiBjb25zdCBbbSwgbl0gPSBBLnNoYXBlO1xuICAgICAgICBsZXQgVSA9IG5ldyBNYXRyaXgobSwgbiwgKGksIGopID0+IGkgPT0gaiA/IDEgOiAwKTtcbiAgICAgICAgY29uc29sZS5sb2coVS50bzJkQXJyYXkpXG4gICAgICAgIGxldCBWID0gbmV3IE1hdHJpeChuLCBtLCAoaSwgaikgPT4gaSA9PSBqID8gMSA6IDApO1xuICAgICAgICBjb25zb2xlLmxvZyhWLnRvMmRBcnJheSlcbiAgICAgICAgbGV0IEIgPSBNYXRyaXguYmlkaWFnb25hbChBLmNsb25lKCksIFUsIFYpO1xuICAgICAgICBjb25zb2xlLmxvZyhVLFYsQilcbiAgICAgICAgcmV0dXJuIHsgVTogVSwgXCJTaWdtYVwiOiBCLCBWOiBWIH07ICovXG4gICAgfVxufVxuIiwiaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IE1hdHJpeCB9IGZyb20gXCIuL01hdHJpeC5qc1wiO1xuXG4vKipcbiAqIENvbXB1dGVzIHRoZSBkaXN0YW5jZSBtYXRyaXggb2YgZGF0YW1hdHJpeCB7QGxpbmsgQX0uXG4gKiBAcGFyYW0ge01hdHJpeH0gQSAtIE1hdHJpeFxuICogQHBhcmFtIHtGdW5jdGlvbn0gW21ldHJpYz1ldWNsaWRlYW5dIC0gVGhlIGRpaXN0YW5jZSBtZXRyaWMuXG4gKiBAcmV0dXJucyB7TWF0cml4fSBEIC0gVGhlIGRpc3RhbmNlIG1hdHJpeCBvZiB7QGxpbmsgQX0uXG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uIChBLCBtZXRyaWMgPSBldWNsaWRlYW4pIHtcbiAgICBsZXQgbiA9IEEuc2hhcGVbMF07XG4gICAgY29uc3QgRCA9IG5ldyBNYXRyaXgobiwgbik7XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgY29uc3QgQV9pID0gQS5yb3coaSk7XG4gICAgICAgIGZvciAobGV0IGogPSBpICsgMTsgaiA8IG47ICsraikge1xuICAgICAgICAgICAgY29uc3QgZGlzdCA9IG1ldHJpYyhBX2ksIEEucm93KGopKTtcbiAgICAgICAgICAgIEQuc2V0X2VudHJ5KGksIGosIGRpc3QpO1xuICAgICAgICAgICAgRC5zZXRfZW50cnkoaiwgaSwgZGlzdCk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIEQ7XG59XG4iLCIvKipcbiAqIEBtZW1iZXJvZiBtb2R1bGU6bWF0cml4XG4gKiBDcmVhdGVzIGFuIEFycmF5IGNvbnRhaW5pbmcge0BsaW5rIG51bWJlcn0gbnVtYmVycyBmcm9tIHtAbGluayBzdGFydH0gdG8ge0BsaW5rIGVuZH0uIElmIDxjb2RlPntAbGluayBudW1iZXJ9ID0gbnVsbDwvbnVsbD5cbiAqIEBwYXJhbSB7TnVtYmVyfSBzdGFydFxuICogQHBhcmFtIHtOdW1iZXJ9IGVuZFxuICogQHBhcmFtIHtOdW1iZXJ9IFtudW1iZXIgPSBudWxsXVxuICogQHJldHVybnMge0FycmF5fVxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAoc3RhcnQsIGVuZCwgbnVtYmVyID0gbnVsbCkge1xuICAgIGlmICghbnVtYmVyKSB7XG4gICAgICAgIG51bWJlciA9IE1hdGgubWF4KE1hdGgucm91bmQoZW5kIC0gc3RhcnQpICsgMSwgMSk7XG4gICAgfVxuICAgIGlmIChudW1iZXIgPCAyKSB7XG4gICAgICAgIHJldHVybiBudW1iZXIgPT09IDEgPyBbc3RhcnRdIDogW107XG4gICAgfVxuICAgIGxldCByZXN1bHQgPSBuZXcgQXJyYXkobnVtYmVyKTtcbiAgICBudW1iZXIgLT0gMTtcbiAgICBmb3IgKGxldCBpID0gbnVtYmVyOyBpID49IDA7IC0taSkge1xuICAgICAgICByZXN1bHRbaV0gPSAoaSAqIGVuZCArIChudW1iZXIgLSBpKSAqIHN0YXJ0KSAvIG51bWJlcjtcbiAgICB9XG4gICAgcmV0dXJuIHJlc3VsdDtcbn1cbiIsImltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBNYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG4vL2ltcG9ydCB7IG5ldW1haXJfc3VtIH0gZnJvbSBcIi4uL251bWVyaWNhbC9pbmRleFwiO1xuXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbih2LCBtZXRyaWMgPSBldWNsaWRlYW4pIHtcbi8vZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24odmVjdG9yLCBwPTIsIG1ldHJpYyA9IGV1Y2xpZGVhbikge1xuICAgIGxldCB2ZWN0b3IgPSBudWxsO1xuICAgIGlmICh2IGluc3RhbmNlb2YgTWF0cml4KSB7XG4gICAgICAgIGxldCBbcm93cywgY29sc10gPSB2LnNoYXBlO1xuICAgICAgICBpZiAocm93cyA9PT0gMSkgdmVjdG9yID0gdi5yb3coMCk7XG4gICAgICAgIGVsc2UgaWYgKGNvbHMgPT09IDEpIHZlY3RvciA9IHYuY29sKDApO1xuICAgICAgICBlbHNlIHRocm93IFwibWF0cml4IG11c3QgYmUgMWQhXCJcbiAgICB9IGVsc2Uge1xuICAgICAgICB2ZWN0b3IgPSB2O1xuICAgIH1cbiAgICBsZXQgbiA9IHZlY3Rvci5sZW5ndGg7XG4gICAgbGV0IHogPSBuZXcgQXJyYXkobilcbiAgICB6LmZpbGwoMCk7XG4gICAgcmV0dXJuIG1ldHJpYyh2ZWN0b3IsIHopO1xuICAgIFxuICAgIFxuICAgIC8qbGV0IHY7XG4gICAgaWYgKHZlY3RvciBpbnN0YW5jZW9mIE1hdHJpeCkge1xuICAgICAgICBsZXQgWyByb3dzLCBjb2xzIF0gPSB2LnNoYXBlO1xuICAgICAgICBpZiAocm93cyA9PT0gMSkge1xuICAgICAgICAgICAgdiA9IHZlY3Rvci5yb3coMCk7XG4gICAgICAgIH0gZWxzZSBpZiAoY29scyA9PT0gMSkge1xuICAgICAgICAgICAgdiA9IHZlY3Rvci5jb2woMCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0aHJvdyBcIm1hdHJpeCBtdXN0IGJlIDFkXCJcbiAgICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICAgIHYgPSB2ZWN0b3I7XG4gICAgfVxuICAgIHJldHVybiBNYXRoLnBvdyhuZXVtYWlyX3N1bSh2Lm1hcChlID0+IE1hdGgucG93KGUsIHApKSksIDEgLyBwKSovXG59IiwiaW1wb3J0IHsgbGluc3BhY2UsIE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcblxuLyoqXG4gKiBAY2xhc3NcbiAqIEBtZW1iZXJvZiBtb2R1bGU6dXRpbHNcbiAqIEBhbGlhcyBSYW5kb21pemVyXG4gKi9cbmV4cG9ydCBjbGFzcyBSYW5kb21pemVyIHtcbiAgICAvKipcbiAgICAgKiBNZXJzZW5uZSBUd2lzdGVyIHJhbmRvbSBudW1iZXIgZ2VuZXJhdG9yLlxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbX3NlZWQ9bmV3IERhdGUoKS5nZXRUaW1lKCldIC0gVGhlIHNlZWQgZm9yIHRoZSByYW5kb20gbnVtYmVyIGdlbmVyYXRvci4gSWYgPGNvZGU+X3NlZWQgPT0gbnVsbDwvY29kZT4gdGhlbiB0aGUgYWN0dWFsIHRpbWUgZ2V0cyB1c2VkIGFzIHNlZWQuXG4gICAgICogQHNlZSBodHRwczovL2dpdGh1Yi5jb20vYm11cnJheTcvbWVyc2VubmUtdHdpc3Rlci1leGFtcGxlcy9ibG9iL21hc3Rlci9qYXZhc2NyaXB0LW1lcnNlbm5lLXR3aXN0ZXIuanNcbiAgICAgKi9cbiAgICBjb25zdHJ1Y3Rvcihfc2VlZCkge1xuICAgICAgICB0aGlzLl9OID0gNjI0O1xuICAgICAgICB0aGlzLl9NID0gMzk3O1xuICAgICAgICB0aGlzLl9NQVRSSVhfQSA9IDB4OTkwOGIwZGY7XG4gICAgICAgIHRoaXMuX1VQUEVSX01BU0sgPSAweDgwMDAwMDAwO1xuICAgICAgICB0aGlzLl9MT1dFUl9NQVNLID0gMHg3ZmZmZmZmZjtcbiAgICAgICAgdGhpcy5fbXQgPSBuZXcgQXJyYXkodGhpcy5fTik7XG4gICAgICAgIHRoaXMuX210aSA9IHRoaXMuTiArIDE7XG5cbiAgICAgICAgdGhpcy5zZWVkID0gX3NlZWQgfHwgbmV3IERhdGUoKS5nZXRUaW1lKCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIHNldCBzZWVkKF9zZWVkKSB7XG4gICAgICAgIHRoaXMuX3NlZWQgPSBfc2VlZDtcbiAgICAgICAgbGV0IG10ID0gdGhpcy5fbXQ7XG5cbiAgICAgICAgbXRbMF0gPSBfc2VlZCA+Pj4gMDtcbiAgICAgICAgZm9yICh0aGlzLl9tdGkgPSAxOyB0aGlzLl9tdGkgPCB0aGlzLl9OOyB0aGlzLl9tdGkgKz0gMSkge1xuICAgICAgICAgICAgbGV0IG10aSA9IHRoaXMuX210aTtcbiAgICAgICAgICAgIGxldCBzID0gbXRbbXRpIC0gMV0gXiAobXRbbXRpIC0gMV0gPj4+IDMwKTtcbiAgICAgICAgICAgIG10W210aV0gPSAoKCgocyAmIDB4ZmZmZjAwMDApID4+PiAxNikgKiAxODEyNDMzMjUzKSA8PCAxNikgKyAocyAmIDB4MDAwMGZmZmYpICogMTgxMjQzMzI1MyArIG10aTtcbiAgICAgICAgICAgIG10W210aV0gPj4+PSAwO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgc2VlZCBvZiB0aGUgcmFuZG9tIG51bWJlciBnZW5lcmF0b3IuXG4gICAgICogQHJldHVybnMge051bWJlcn0gLSBUaGUgc2VlZC5cbiAgICAgKi9cbiAgICBnZXQgc2VlZCgpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX3NlZWQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyBhIGZsb2F0IGJldHdlZW4gMCBhbmQgMS5cbiAgICAgKiBAcmV0dXJucyB7TnVtYmVyfSAtIEEgcmFuZG9tIG51bWJlciBiZXR3ZWVuIFswLCAxXVxuICAgICAqL1xuICAgIGdldCByYW5kb20oKSB7XG4gICAgICAgIHJldHVybiB0aGlzLnJhbmRvbV9pbnQgKiAoMS4wIC8gNDI5NDk2NzI5Ni4wKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIGFuIGludGVnZXIgYmV0d2VlbiAwIGFuZCBNQVhfSU5URUdFUi5cbiAgICAgKiBAcmV0dXJucyB7SW50ZWdlcn0gLSBBIHJhbmRvbSBpbnRlZ2VyLlxuICAgICAqL1xuICAgIGdldCByYW5kb21faW50KCkge1xuICAgICAgICBsZXQgeSxcbiAgICAgICAgICAgIG1hZzAxID0gbmV3IEFycmF5KDB4MCwgdGhpcy5fTUFUUklYX0EpO1xuICAgICAgICBpZiAodGhpcy5fbXRpID49IHRoaXMuX04pIHtcbiAgICAgICAgICAgIGxldCBraztcblxuICAgICAgICAgICAgLyogaWYgKHRoaXMuX210aSA9PSB0aGlzLl9OICsgMSkge1xuICAgICAgICAgICAgICAgIHRoaXMuc2VlZCA9IDU0ODk7XG4gICAgICAgICAgICB9ICovXG5cbiAgICAgICAgICAgIGxldCBOX00gPSB0aGlzLl9OIC0gdGhpcy5fTTtcbiAgICAgICAgICAgIGxldCBNX04gPSB0aGlzLl9NIC0gdGhpcy5fTjtcblxuICAgICAgICAgICAgZm9yIChrayA9IDA7IGtrIDwgTl9NOyArK2trKSB7XG4gICAgICAgICAgICAgICAgeSA9ICh0aGlzLl9tdFtra10gJiB0aGlzLl9VUFBFUl9NQVNLKSB8ICh0aGlzLl9tdFtrayArIDFdICYgdGhpcy5fTE9XRVJfTUFTSyk7XG4gICAgICAgICAgICAgICAgdGhpcy5fbXRba2tdID0gdGhpcy5fbXRba2sgKyB0aGlzLl9NXSBeICh5ID4+PiAxKSBeIG1hZzAxW3kgJiAweDFdO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yICg7IGtrIDwgdGhpcy5fTiAtIDE7ICsra2spIHtcbiAgICAgICAgICAgICAgICB5ID0gKHRoaXMuX210W2trXSAmIHRoaXMuX1VQUEVSX01BU0spIHwgKHRoaXMuX210W2trICsgMV0gJiB0aGlzLl9MT1dFUl9NQVNLKTtcbiAgICAgICAgICAgICAgICB0aGlzLl9tdFtra10gPSB0aGlzLl9tdFtrayArIE1fTl0gXiAoeSA+Pj4gMSkgXiBtYWcwMVt5ICYgMHgxXTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgeSA9ICh0aGlzLl9tdFt0aGlzLl9OIC0gMV0gJiB0aGlzLl9VUFBFUl9NQVNLKSB8ICh0aGlzLl9tdFswXSAmIHRoaXMuX0xPV0VSX01BU0spO1xuICAgICAgICAgICAgdGhpcy5fbXRbdGhpcy5fTiAtIDFdID0gdGhpcy5fbXRbdGhpcy5fTSAtIDFdIF4gKHkgPj4+IDEpIF4gbWFnMDFbeSAmIDB4MV07XG5cbiAgICAgICAgICAgIHRoaXMuX210aSA9IDA7XG4gICAgICAgIH1cblxuICAgICAgICB5ID0gdGhpcy5fbXRbKHRoaXMuX210aSArPSAxKV07XG4gICAgICAgIHkgXj0geSA+Pj4gMTE7XG4gICAgICAgIHkgXj0gKHkgPDwgNykgJiAweDlkMmM1NjgwO1xuICAgICAgICB5IF49ICh5IDw8IDE1KSAmIDB4ZWZjNjAwMDA7XG4gICAgICAgIHkgXj0geSA+Pj4gMTg7XG5cbiAgICAgICAgcmV0dXJuIHkgPj4+IDA7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyBzYW1wbGVzIGZyb20gYW4gaW5wdXQgTWF0cml4IG9yIEFycmF5LlxuICAgICAqIEBwYXJhbSB7TWF0cml4fEFycmF5fEZsb2F0NjRBcnJheX0gQSAtIFRoZSBpbnB1dCBNYXRyaXggb3IgQXJyYXkuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG4gLSBUaGUgbnVtYmVyIG9mIHNhbXBsZXMuXG4gICAgICogQHJldHVybnMge0FycmF5fSAtIEEgcmFuZG9tIHNlbGVjdGlvbiBmb3JtIHtAbGluayBBfSBvZiB7QGxpbmsgbn0gc2FtcGxlcy5cbiAgICAgKi9cbiAgICBjaG9pY2UoQSwgbikge1xuICAgICAgICBpZiAoQSBpbnN0YW5jZW9mIE1hdHJpeCkge1xuICAgICAgICAgICAgbGV0IHJvd3MgPSBBLnNoYXBlWzBdO1xuICAgICAgICAgICAgaWYgKG4gPiByb3dzKSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKFwibiBiaWdnZXIgdGhhbiBBIVwiKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGxldCBzYW1wbGUgPSBuZXcgQXJyYXkobik7XG4gICAgICAgICAgICBsZXQgaW5kZXhfbGlzdCA9IGxpbnNwYWNlKDAsIHJvd3MgLSAxKTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSAwLCBsID0gaW5kZXhfbGlzdC5sZW5ndGg7IGkgPCBuOyArK2ksIC0tbCkge1xuICAgICAgICAgICAgICAgIGxldCByYW5kb21faW5kZXggPSB0aGlzLnJhbmRvbV9pbnQgJSBsO1xuICAgICAgICAgICAgICAgIHNhbXBsZVtpXSA9IGluZGV4X2xpc3Quc3BsaWNlKHJhbmRvbV9pbmRleCwgMSlbMF07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gc2FtcGxlLm1hcCgoZCkgPT4gQS5yb3coZCkpO1xuICAgICAgICB9IGVsc2UgaWYgKEFycmF5LmlzQXJyYXkoQSkgfHwgQSBpbnN0YW5jZW9mIEZsb2F0NjRBcnJheSkge1xuICAgICAgICAgICAgbGV0IHJvd3MgPSBBLmxlbmd0aDtcbiAgICAgICAgICAgIGlmIChuID4gcm93cykge1xuICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcIm4gYmlnZ2VyIHRoYW4gQSFcIik7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBsZXQgc2FtcGxlID0gbmV3IEFycmF5KG4pO1xuICAgICAgICAgICAgbGV0IGluZGV4X2xpc3QgPSBsaW5zcGFjZSgwLCByb3dzIC0gMSk7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMCwgbCA9IGluZGV4X2xpc3QubGVuZ3RoOyBpIDwgbjsgKytpLCAtLWwpIHtcbiAgICAgICAgICAgICAgICBsZXQgcmFuZG9tX2luZGV4ID0gdGhpcy5yYW5kb21faW50ICUgbDtcbiAgICAgICAgICAgICAgICBzYW1wbGVbaV0gPSBpbmRleF9saXN0LnNwbGljZShyYW5kb21faW5kZXgsIDEpWzBdO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIHNhbXBsZS5tYXAoKGQpID0+IEFbZF0pO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHN0YXRpY1xuICAgICAqIFJldHVybnMgc2FtcGxlcyBmcm9tIGFuIGlucHV0IE1hdHJpeCBvciBBcnJheS5cbiAgICAgKiBAcGFyYW0ge01hdHJpeHxBcnJheXxGbG9hdDY0QXJyYXl9IEEgLSBUaGUgaW5wdXQgTWF0cml4IG9yIEFycmF5LlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBuIC0gVGhlIG51bWJlciBvZiBzYW1wbGVzLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBzZWVkIC0gVGhlIHNlZWQgZm9yIHRoZSByYW5kb20gbnVtYmVyIGdlbmVyYXRvci5cbiAgICAgKiBAcmV0dXJucyB7QXJyYXl9IC0gQSByYW5kb20gc2VsZWN0aW9uIGZvcm0ge0BsaW5rIEF9IG9mIHtAbGluayBufSBzYW1wbGVzLlxuICAgICAqL1xuICAgIHN0YXRpYyBjaG9pY2UoQSwgbiwgc2VlZCA9IDEyMTIpIHtcbiAgICAgICAgY29uc3QgUiA9IG5ldyBSYW5kb21pemVyKHNlZWQpO1xuICAgICAgICByZXR1cm4gUi5jaG9pY2UoQSwgbik7XG4gICAgICAgIC8qIGxldCByb3dzID0gQS5zaGFwZVswXTtcbiAgICAgICAgaWYgKG4gPiByb3dzKSB7XG4gICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoXCJuIGJpZ2dlciB0aGFuIEEhXCIpO1xuICAgICAgICB9XG4gICAgICAgIGxldCByYW5kID0gbmV3IFJhbmRvbWl6ZXIoc2VlZCk7XG4gICAgICAgIGxldCBzYW1wbGUgPSBuZXcgQXJyYXkobik7XG4gICAgICAgIGxldCBpbmRleF9saXN0ID0gbGluc3BhY2UoMCwgcm93cyAtIDEpO1xuICAgICAgICBmb3IgKGxldCBpID0gMCwgbCA9IGluZGV4X2xpc3QubGVuZ3RoOyBpIDwgbjsgKytpLCAtLWwpIHtcbiAgICAgICAgICAgIGxldCByYW5kb21faW5kZXggPSByYW5kLnJhbmRvbV9pbnQgJSBsO1xuICAgICAgICAgICAgc2FtcGxlW2ldID0gaW5kZXhfbGlzdC5zcGxpY2UocmFuZG9tX2luZGV4LCAxKVswXTtcbiAgICAgICAgfVxuICAgICAgICAvL3JldHVybiByZXN1bHQ7XG4gICAgICAgIC8vcmV0dXJuIG5ldyBNYXRyaXgobiwgY29scywgKHJvdywgY29sKSA9PiBBLmVudHJ5KHNhbXBsZVtyb3ddLCBjb2wpKVxuICAgICAgICByZXR1cm4gc2FtcGxlLm1hcCgoZCkgPT4gQS5yb3coZCkpOyAqL1xuICAgIH1cbn1cbiIsIi8qKlxuICogUmV0dXJucyBtYXhpbXVtIGluIEFycmF5IHtAbGluayB2YWx1ZXN9LlxuICogQG1lbWJlcm9mIG1vZHVsZTp1dGlsc1xuICogQGFsaWFzIG1heFxuICogQHBhcmFtIHtBcnJheX0gdmFsdWVzIFxuICogQHJldHVybnMge051bWJlcn1cbiAqL1xuZXhwb3J0IGRlZmF1bHQgZnVuY3Rpb24gKHZhbHVlcykge1xuICAgIGxldCBtYXg7XG4gICAgZm9yIChjb25zdCB2YWx1ZSBvZiB2YWx1ZXMpIHtcbiAgICAgICAgaWYgKHZhbHVlICE9IG51bGwgJiYgKG1heCA8IHZhbHVlIHx8IChtYXggPT09IHVuZGVmaW5lZCAmJiB2YWx1ZSA+PSB2YWx1ZSkpKSB7XG4gICAgICAgICAgICBtYXggPSB2YWx1ZTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gbWF4O1xufSIsIi8qKlxuICogUmV0dXJucyBtYXhpbXVtIGluIEFycmF5IHtAbGluayB2YWx1ZXN9LlxuICogQG1lbWJlcm9mIG1vZHVsZTp1dGlsc1xuICogQGFsaWFzIG1pblxuICogQHBhcmFtIHtBcnJheX0gdmFsdWVzXG4gKiBAcmV0dXJucyB7TnVtYmVyfVxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAodmFsdWVzKSB7XG4gICAgbGV0IG1pbjtcbiAgICBmb3IgKGNvbnN0IHZhbHVlIG9mIHZhbHVlcykge1xuICAgICAgICBpZiAodmFsdWUgIT0gbnVsbCAmJiAobWluID4gdmFsdWUgfHwgKG1pbiA9PT0gdW5kZWZpbmVkICYmIHZhbHVlIDw9IHZhbHVlKSkpIHtcbiAgICAgICAgICAgIG1pbiA9IHZhbHVlO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBtaW47XG59IiwiLyoqXG4gKiBAY2xhc3NcbiAqIEBhbGlhcyBIZWFwXG4gKi9cbmV4cG9ydCBjbGFzcyBIZWFwIHtcbiAgICAvKipcbiAgICAgKiBBIGhlYXAgaXMgYSBkYXRhc3RydWN0dXJlIGhvbGRpbmcgaXRzIGVsZW1lbnRzIGluIGEgc3BlY2lmaWMgd2F5LCBzbyB0aGF0IHRoZSB0b3AgZWxlbWVudCB3b3VsZCBiZSB0aGUgZmlyc3QgZW50cnkgb2YgYW4gb3JkZXJlZCBsaXN0LlxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGF0YXN0cnVjdHVyZVxuICAgICAqIEBhbGlhcyBIZWFwXG4gICAgICogQHBhcmFtIHtBcnJheT19IGVsZW1lbnRzIC0gQ29udGFpbnMgdGhlIGVsZW1lbnRzIGZvciB0aGUgSGVhcC4ge0BsaW5rIGVsZW1lbnRzfSBjYW4gYmUgbnVsbC5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbYWNjZXNzb3IgPSAoZCkgPT4gZF0gLSBGdW5jdGlvbiByZXR1cm5zIHRoZSB2YWx1ZSBvZiB0aGUgZWxlbWVudC5cbiAgICAgKiBAcGFyYW0geyhcIm1pblwifFwibWF4XCJ8RnVuY3Rpb24pfSBbY29tcGFyYXRvciA9IFwibWluXCJdIC0gRnVuY3Rpb24gcmV0dXJuaW5nIHRydWUgb3IgZmFsc2UgZGVmaW5pbmcgdGhlIHdpc2hlZCBvcmRlciBvZiB0aGUgSGVhcCwgb3IgU3RyaW5nIGZvciBwcmVkZWZpbmVkIGZ1bmN0aW9uLiAoXCJtaW5cIiBmb3IgYSBNaW4tSGVhcCwgXCJtYXhcIiBmb3IgYSBNYXhfaGVhcClcbiAgICAgKiBAcmV0dXJucyB7SGVhcH1cbiAgICAgKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9CaW5hcnlfaGVhcH1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihlbGVtZW50cyA9IG51bGwsIGFjY2Vzc29yID0gZCA9PiBkLCBjb21wYXJhdG9yID0gXCJtaW5cIikge1xuICAgICAgICBpZiAoZWxlbWVudHMpIHtcbiAgICAgICAgICAgIHJldHVybiBIZWFwLmhlYXBpZnkoZWxlbWVudHMsIGFjY2Vzc29yLCBjb21wYXJhdG9yKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRoaXMuX2FjY2Vzc29yID0gYWNjZXNzb3I7XG4gICAgICAgICAgICB0aGlzLl9jb250YWluZXIgPSBbXTtcbiAgICAgICAgICAgIGlmIChjb21wYXJhdG9yID09IFwibWluXCIpIHtcbiAgICAgICAgICAgICAgICB0aGlzLl9jb21wYXJhdG9yID0gKGEsIGIpID0+IGEgPCBiO1xuICAgICAgICAgICAgfSBlbHNlIGlmIChjb21wYXJhdG9yID09IFwibWF4XCIpIHtcbiAgICAgICAgICAgICAgICB0aGlzLl9jb21wYXJhdG9yID0gKGEsIGIpID0+IGEgPiBiO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICB0aGlzLl9jb21wYXJhdG9yID0gY29tcGFyYXRvcjtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiB0aGlzXG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDcmVhdGVzIGEgSGVhcCBmcm9tIGFuIEFycmF5XG4gICAgICogQHBhcmFtIHtBcnJheXxTZXR9IGVsZW1lbnRzIC0gQ29udGFpbnMgdGhlIGVsZW1lbnRzIGZvciB0aGUgSGVhcC5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9uPX0gW2FjY2Vzc29yID0gKGQpID0+IGRdIC0gRnVuY3Rpb24gcmV0dXJucyB0aGUgdmFsdWUgb2YgdGhlIGVsZW1lbnQuXG4gICAgICogQHBhcmFtIHsoU3RyaW5nPXxGdW5jdGlvbil9IFtjb21wYXJhdG9yID0gXCJtaW5cIl0gLSBGdW5jdGlvbiByZXR1cm5pbmcgdHJ1ZSBvciBmYWxzZSBkZWZpbmluZyB0aGUgd2lzaGVkIG9yZGVyIG9mIHRoZSBIZWFwLCBvciBTdHJpbmcgZm9yIHByZWRlZmluZWQgZnVuY3Rpb24uIChcIm1pblwiIGZvciBhIE1pbi1IZWFwLCBcIm1heFwiIGZvciBhIE1heF9oZWFwKVxuICAgICAqIEByZXR1cm5zIHtIZWFwfVxuICAgICAqL1xuICAgIHN0YXRpYyBoZWFwaWZ5KGVsZW1lbnRzLCBhY2Nlc3NvciA9IGQgPT4gZCwgY29tcGFyYXRvciA9IFwibWluXCIpIHtcbiAgICAgICAgY29uc3QgaGVhcCA9IG5ldyBIZWFwKG51bGwsIGFjY2Vzc29yLCBjb21wYXJhdG9yKTtcbiAgICAgICAgY29uc3QgY29udGFpbmVyID0gaGVhcC5fY29udGFpbmVyO1xuICAgICAgICBmb3IgKGNvbnN0IGUgb2YgZWxlbWVudHMpIHtcbiAgICAgICAgICAgIGNvbnRhaW5lci5wdXNoKHtcbiAgICAgICAgICAgICAgICBcImVsZW1lbnRcIjogZSxcbiAgICAgICAgICAgICAgICBcInZhbHVlXCI6IGFjY2Vzc29yKGUpLFxuICAgICAgICAgICAgfSk7XG4gICAgICAgIH1cbiAgICAgICAgZm9yIChsZXQgaSA9IE1hdGguZmxvb3IoKGVsZW1lbnRzLmxlbmd0aCAvIDIpIC0gMSk7IGkgPj0gMDsgLS1pKSB7XG4gICAgICAgICAgICBoZWFwLl9oZWFwaWZ5X2Rvd24oaSk7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIGhlYXA7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogU3dhcHMgZWxlbWVudHMgb2YgY29udGFpbmVyIGFycmF5LlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IGluZGV4X2EgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IGluZGV4X2IgXG4gICAgICovXG4gICAgX3N3YXAoaW5kZXhfYSwgaW5kZXhfYikge1xuICAgICAgICBjb25zdCBjb250YWluZXIgPSB0aGlzLl9jb250YWluZXI7XG4gICAgICAgIFtjb250YWluZXJbaW5kZXhfYl0sIGNvbnRhaW5lcltpbmRleF9hXV0gPSBbY29udGFpbmVyW2luZGV4X2FdLCBjb250YWluZXJbaW5kZXhfYl1dO1xuICAgICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKi9cbiAgICBfaGVhcGlmeV91cCgpIHtcbiAgICAgICAgY29uc3QgY29udGFpbmVyID0gdGhpcy5fY29udGFpbmVyO1xuICAgICAgICBsZXQgaW5kZXggPSBjb250YWluZXIubGVuZ3RoIC0gMTtcbiAgICAgICAgd2hpbGUgKGluZGV4ID4gMCkge1xuICAgICAgICAgICAgbGV0IHBhcmVudEluZGV4ID0gTWF0aC5mbG9vcigoaW5kZXggLSAxKSAvIDIpO1xuICAgICAgICAgICAgaWYgKCF0aGlzLl9jb21wYXJhdG9yKGNvbnRhaW5lcltpbmRleF0udmFsdWUsIGNvbnRhaW5lcltwYXJlbnRJbmRleF0udmFsdWUpKSB7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdGhpcy5fc3dhcChwYXJlbnRJbmRleCwgaW5kZXgpXG4gICAgICAgICAgICBpbmRleCA9IHBhcmVudEluZGV4O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUHVzaGVzIHRoZSBlbGVtZW50IHRvIHRoZSBoZWFwLlxuICAgICAqIEBwYXJhbSB7fSBlbGVtZW50XG4gICAgICogQHJldHVybnMge0hlYXB9XG4gICAgICovXG4gICAgcHVzaChlbGVtZW50KSB7XG4gICAgICAgIGNvbnN0IHZhbHVlID0gdGhpcy5fYWNjZXNzb3IoZWxlbWVudCk7XG4gICAgICAgIC8vY29uc3Qgbm9kZSA9IG5ldyBOb2RlKGVsZW1lbnQsIHZhbHVlKTtcbiAgICAgICAgY29uc3Qgbm9kZSA9IHtcImVsZW1lbnRcIjogZWxlbWVudCwgXCJ2YWx1ZVwiOiB2YWx1ZX07XG4gICAgICAgIHRoaXMuX2NvbnRhaW5lci5wdXNoKG5vZGUpO1xuICAgICAgICB0aGlzLl9oZWFwaWZ5X3VwKCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzdGFydF9pbmRleCA9IDBdIFxuICAgICAqL1xuICAgIF9oZWFwaWZ5X2Rvd24oc3RhcnRfaW5kZXg9MCkge1xuICAgICAgICBjb25zdCBjb250YWluZXIgPSB0aGlzLl9jb250YWluZXI7XG4gICAgICAgIGNvbnN0IGNvbXBhcmF0b3IgPSB0aGlzLl9jb21wYXJhdG9yO1xuICAgICAgICBjb25zdCBsZW5ndGggPSBjb250YWluZXIubGVuZ3RoO1xuICAgICAgICBsZXQgbGVmdCA9IDIgKiBzdGFydF9pbmRleCArIDE7XG4gICAgICAgIGxldCByaWdodCA9IDIgKiBzdGFydF9pbmRleCArIDI7XG4gICAgICAgIGxldCBpbmRleCA9IHN0YXJ0X2luZGV4O1xuICAgICAgICBpZiAoaW5kZXggPiBsZW5ndGgpIHRocm93IFwiaW5kZXggaGlnaGVyIHRoYW4gbGVuZ3RoXCJcbiAgICAgICAgaWYgKGxlZnQgPCBsZW5ndGggJiYgY29tcGFyYXRvcihjb250YWluZXJbbGVmdF0udmFsdWUsIGNvbnRhaW5lcltpbmRleF0udmFsdWUpKSB7XG4gICAgICAgICAgICBpbmRleCA9IGxlZnQ7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHJpZ2h0IDwgbGVuZ3RoICYmIGNvbXBhcmF0b3IoY29udGFpbmVyW3JpZ2h0XS52YWx1ZSwgY29udGFpbmVyW2luZGV4XS52YWx1ZSkpIHtcbiAgICAgICAgICAgIGluZGV4ID0gcmlnaHQ7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKGluZGV4ICE9PSBzdGFydF9pbmRleCkge1xuICAgICAgICAgICAgdGhpcy5fc3dhcChzdGFydF9pbmRleCwgaW5kZXgpO1xuICAgICAgICAgICAgdGhpcy5faGVhcGlmeV9kb3duKGluZGV4KTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJlbW92ZXMgYW5kIHJldHVybnMgdGhlIHRvcCBlbnRyeSBvZiB0aGUgaGVhcC5cbiAgICAgKiBAcmV0dXJucyB7T2JqZWN0fSBPYmplY3QgY29uc2lzdHMgb2YgdGhlIGVsZW1lbnQgYW5kIGl0cyB2YWx1ZSAoY29tcHV0ZWQgYnkge0BsaW5rIGFjY2Vzc29yfSkuXG4gICAgICovXG4gICAgcG9wKCkge1xuICAgICAgICBjb25zdCBjb250YWluZXIgPSB0aGlzLl9jb250YWluZXI7XG4gICAgICAgIGlmIChjb250YWluZXIubGVuZ3RoID09PSAwKSB7XG4gICAgICAgICAgICByZXR1cm4gbnVsbDtcbiAgICAgICAgfSBlbHNlIGlmIChjb250YWluZXIubGVuZ3RoID09PSAxKSB7XG4gICAgICAgICAgICByZXR1cm4gY29udGFpbmVyLnBvcCgpO1xuICAgICAgICB9XG4gICAgICAgIHRoaXMuX3N3YXAoMCwgY29udGFpbmVyLmxlbmd0aCAtIDEpO1xuICAgICAgICBjb25zdCBpdGVtID0gY29udGFpbmVyLnBvcCgpO1xuICAgICAgICB0aGlzLl9oZWFwaWZ5X2Rvd24oKTtcbiAgICAgICAgcmV0dXJuIGl0ZW07XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyB0aGUgdG9wIGVudHJ5IG9mIHRoZSBoZWFwIHdpdGhvdXQgcmVtb3ZpbmcgaXQuXG4gICAgICogQHJldHVybnMge09iamVjdH0gT2JqZWN0IGNvbnNpc3RzIG9mIHRoZSBlbGVtZW50IGFuZCBpdHMgdmFsdWUgKGNvbXB1dGVkIGJ5IHtAbGluayBhY2Nlc3Nvcn0pLlxuICAgICAqL1xuICAgIGdldCBmaXJzdCgpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX2NvbnRhaW5lci5sZW5ndGggPiAwID8gdGhpcy5fY29udGFpbmVyWzBdIDogbnVsbDtcbiAgICB9XG5cblxuICAgIC8qKlxuICAgICAqIFlpZWxkcyB0aGUgcmF3IGRhdGFcbiAgICAgKiBAeWllbGRzIHtPYmplY3R9IE9iamVjdCBjb25zaXN0cyBvZiB0aGUgZWxlbWVudCBhbmQgaXRzIHZhbHVlIChjb21wdXRlZCBieSB7QGxpbmsgYWNjZXNzb3J9KS5cbiAgICAgKi9cbiAgICAqIGl0ZXJhdGUoKSB7XG4gICAgICAgIGZvciAobGV0IGkgPSAwLCBuID0gdGhpcy5fY29udGFpbmVyLmxlbmd0aDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgeWllbGQgdGhpcy5fY29udGFpbmVyW2ldLmVsZW1lbnQ7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRoZSBoZWFwIGFzIG9yZGVyZWQgYXJyYXkuXG4gICAgICogQHJldHVybnMge0FycmF5fSBBcnJheSBjb25zaXN0aW5nIHRoZSBlbGVtZW50cyBvcmRlcmVkIGJ5IHtAbGluayBjb21wYXJhdG9yfS5cbiAgICAgKi9cbiAgICB0b0FycmF5KCkge1xuICAgICAgICByZXR1cm4gdGhpcy5kYXRhKClcbiAgICAgICAgICAgIC5zb3J0KChhLGIpID0+IHRoaXMuX2NvbXBhcmF0b3IoYSwgYikgPyAtMSA6IDApXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogUmV0dXJucyBlbGVtZW50cyBvZiBjb250YWluZXIgYXJyYXkuXG4gICAgICogQHJldHVybnMge0FycmF5fSBBcnJheSBjb25zaXN0aW5nIHRoZSBlbGVtZW50cy5cbiAgICAgKi9cbiAgICBkYXRhKCkge1xuICAgICAgICByZXR1cm4gdGhpcy5fY29udGFpbmVyXG4gICAgICAgICAgICAubWFwKGQgPT4gZC5lbGVtZW50KVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgdGhlIGNvbnRhaW5lciBhcnJheS5cbiAgICAgKiBAcmV0dXJucyB7QXJyYXl9IFRoZSBjb250YWluZXIgYXJyYXkuXG4gICAgICovXG4gICAgcmF3X2RhdGEoKSB7XG4gICAgICAgIHJldHVybiB0aGlzLl9jb250YWluZXI7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogVGhlIHNpemUgb2YgdGhlIGhlYXAuXG4gICAgICogQHJldHVybnMge051bWJlcn1cbiAgICAgKi9cbiAgICBnZXQgbGVuZ3RoKCkge1xuICAgICAgICByZXR1cm4gdGhpcy5fY29udGFpbmVyLmxlbmd0aDtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIGZhbHNlIGlmIHRoZSB0aGUgaGVhcCBoYXMgZW50cmllcywgdHJ1ZSBpZiB0aGUgaGVhcCBoYXMgbm8gZW50cmllcy5cbiAgICAgKiBAcmV0dXJucyB7Qm9vbGVhbn1cbiAgICAgKi9cbiAgICBnZXQgZW1wdHkoKSB7XG4gICAgICAgIHJldHVybiB0aGlzLmxlbmd0aCA9PT0gMDtcbiAgICB9XG59IiwiLyoqXG4gKiBAY2xhc3NcbiAqIEBhbGlhcyBEaXNqb2ludFNldFxuICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvRGlzam9pbnQtc2V0X2RhdGFfc3RydWN0dXJlfVxuICovXG5leHBvcnQgY2xhc3MgRGlzam9pbnRTZXQge1xuICAgIC8qKlxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBhbGlhcyBEaXNqb2ludFNldFxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGF0YXN0cnVjdHVyZVxuICAgICAqIEBwYXJhbSB7QXJyYXk9fSBlbGVtZW50cyBcbiAgICAgKiBAcmV0dXJucyB7RGlzam9pbnRTZXR9XG4gICAgICovXG4gICAgY29uc3RydWN0b3IoZWxlbWVudHMgPSBudWxsKSB7XG4gICAgICAgIHRoaXMuX2xpc3QgPSBuZXcgU2V0KCk7XG4gICAgICAgIGlmIChlbGVtZW50cykge1xuICAgICAgICAgICAgZm9yIChjb25zdCBlIG9mIGVsZW1lbnRzKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5tYWtlX3NldChlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICBtYWtlX3NldCh4KSB7XG4gICAgICAgIGNvbnN0IGxpc3QgPSB0aGlzLl9saXN0O1xuICAgICAgICBpZiAoIWxpc3QuaGFzKHgpKSB7XG4gICAgICAgICAgICBsaXN0LmFkZCh4KTtcbiAgICAgICAgICAgIHguX19kaXNqb2ludF9zZXQgPSB7fTtcbiAgICAgICAgICAgIHguX19kaXNqb2ludF9zZXQucGFyZW50ID0geDtcbiAgICAgICAgICAgIHguX19kaXNqb2ludF9zZXQuY2hpbGRyZW4gPSBuZXcgU2V0KFt4XSk7XG4gICAgICAgICAgICB4Ll9fZGlzam9pbnRfc2V0LnNpemUgPSAxO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIGZpbmQoeCkge1xuICAgICAgICBjb25zdCBsaXN0ID0gdGhpcy5fbGlzdDtcbiAgICAgICAgaWYgKGxpc3QuaGFzKHgpKSB7XG4gICAgICAgICAgICBpZiAoeC5fX2Rpc2pvaW50X3NldC5wYXJlbnQgIT09IHgpIHtcbiAgICAgICAgICAgICAgICB4Ll9fZGlzam9pbnRfc2V0LmNoaWxkcmVuLmFkZCguLi54KTtcbiAgICAgICAgICAgICAgICB4Ll9fZGlzam9pbnRfc2V0LnBhcmVudCA9IHRoaXMuZmluZCh4Ll9fZGlzam9pbnRfc2V0LnBhcmVudCk7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHguX19kaXNqb2ludF9zZXQucGFyZW50O1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICByZXR1cm4geDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHJldHVybiBudWxsO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgdW5pb24oeCwgeSkge1xuICAgICAgICBsZXQgbm9kZV94ID0gdGhpcy5maW5kKHgpO1xuICAgICAgICBsZXQgbm9kZV95ID0gdGhpcy5maW5kKHkpO1xuXG4gICAgICAgIGlmIChub2RlX3ggPT09IG5vZGVfeSkgcmV0dXJuIHRoaXM7XG4gICAgICAgIGlmIChub2RlX3guX19kaXNqb2ludF9zZXQuc2l6ZSA8IG5vZGVfeS5fX2Rpc2pvaW50X3NldC5zaXplKSBbbm9kZV94LCBub2RlX3ldID0gW25vZGVfeSwgbm9kZV94XTtcblxuICAgICAgICBub2RlX3kuX19kaXNqb2ludF9zZXQucGFyZW50ID0gbm9kZV94O1xuICAgICAgICAvLyBrZWVwIHRyYWNrIG9mIGNoaWxkcmVuP1xuICAgICAgICBub2RlX3kuX19kaXNqb2ludF9zZXQuY2hpbGRyZW4uZm9yRWFjaChub2RlX3guX19kaXNqb2ludF9zZXQuY2hpbGRyZW4uYWRkLCBub2RlX3guX19kaXNqb2ludF9zZXQuY2hpbGRyZW4pO1xuICAgICAgICBub2RlX3guX19kaXNqb2ludF9zZXQuc2l6ZSArPSBub2RlX3kuX19kaXNqb2ludF9zZXQuc2l6ZTtcblxuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG59IiwiaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IEhlYXAgfSBmcm9tIFwiLi4vZGF0YXN0cnVjdHVyZS9pbmRleC5qc1wiO1xuLyoqXG4gKiBAY2xhc3NcbiAqIEBhbGlhcyBCYWxsVHJlZVxuICovXG5leHBvcnQgY2xhc3MgQmFsbFRyZWUge1xuICAgIC8qKlxuICAgICAqIEdlbmVyYXRlcyBhIEJhbGxUcmVlIHdpdGggZ2l2ZW4ge0BsaW5rIGVsZW1lbnRzfS5cbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmtublxuICAgICAqIEBhbGlhcyBCYWxsVHJlZVxuICAgICAqIEBwYXJhbSB7QXJyYXk9fSBlbGVtZW50cyAtIEVsZW1lbnRzIHdoaWNoIHNob3VsZCBiZSBhZGRlZCB0byB0aGUgQmFsbFRyZWVcbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSBtZXRyaWMgdG8gdXNlOiAoYSwgYikgPT4gZGlzdGFuY2VcbiAgICAgKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9CYWxsX3RyZWV9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9naXRodWIuY29tL2ludmlzYWwvbm9vYmpzL2Jsb2IvbWFzdGVyL3NyYy90cmVlL0JhbGxUcmVlLmpzfVxuICAgICAqIEByZXR1cm5zIHtCYWxsVHJlZX1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihlbGVtZW50cyA9IG51bGwsIG1ldHJpYyA9IGV1Y2xpZGVhbikge1xuICAgICAgICB0aGlzLl9Ob2RlID0gY2xhc3Mge1xuICAgICAgICAgICAgY29uc3RydWN0b3IocGl2b3QsIGNoaWxkMT1udWxsLCBjaGlsZDI9bnVsbCwgcmFkaXVzPW51bGwpIHtcbiAgICAgICAgICAgICAgICB0aGlzLnBpdm90ID0gcGl2b3Q7XG4gICAgICAgICAgICAgICAgdGhpcy5jaGlsZDEgPSBjaGlsZDE7XG4gICAgICAgICAgICAgICAgdGhpcy5jaGlsZDIgPSBjaGlsZDI7XG4gICAgICAgICAgICAgICAgdGhpcy5yYWRpdXMgPSByYWRpdXM7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdGhpcy5fTGVhZiA9IGNsYXNzIHtcbiAgICAgICAgICAgIGNvbnN0cnVjdG9yKHBvaW50cykge1xuICAgICAgICAgICAgICAgIHRoaXMucG9pbnRzID0gcG9pbnRzO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHRoaXMuX21ldHJpYyA9IG1ldHJpYztcbiAgICAgICAgaWYgKGVsZW1lbnRzKSB7XG4gICAgICAgICAgICB0aGlzLmFkZChlbGVtZW50cyk7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHtBcnJheTwqPn0gZWxlbWVudHMgLSBuZXcgZWxlbWVudHMuXG4gICAgICogQHJldHVybnMge0JhbGxUcmVlfVxuICAgICAqL1xuICAgIGFkZChlbGVtZW50cykge1xuICAgICAgICBlbGVtZW50cyA9IGVsZW1lbnRzLm1hcCgoZWxlbWVudCwgaW5kZXgpID0+IHtcbiAgICAgICAgICAgIHJldHVybiB7aW5kZXg6IGluZGV4LCBlbGVtZW50OiBlbGVtZW50fVxuICAgICAgICB9KVxuICAgICAgICB0aGlzLl9yb290ID0gdGhpcy5fY29uc3RydWN0KGVsZW1lbnRzKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge0FycmF5PCo+fSBlbGVtZW50cyBcbiAgICAgKiBAcmV0dXJucyB7Tm9kZX0gcm9vdCBvZiBiYWxsdHJlZS5cbiAgICAgKi9cbiAgICBfY29uc3RydWN0KGVsZW1lbnRzKSB7XG4gICAgICAgIGlmIChlbGVtZW50cy5sZW5ndGggPT09IDEpIHtcbiAgICAgICAgICAgIHJldHVybiBuZXcgdGhpcy5fTGVhZihlbGVtZW50cyk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBsZXQgYyA9IHRoaXMuX2dyZWF0ZXN0X3NwcmVhZChlbGVtZW50cyk7XG4gICAgICAgICAgICBsZXQgc29ydGVkX2VsZW1lbnRzID0gZWxlbWVudHMuc29ydCgoYSwgYikgPT4gYS5lbGVtZW50W2NdIC0gYi5lbGVtZW50W2NdKTtcbiAgICAgICAgICAgIGxldCBuID0gc29ydGVkX2VsZW1lbnRzLmxlbmd0aDtcbiAgICAgICAgICAgIGxldCBwX2luZGV4ID0gTWF0aC5mbG9vcihuIC8gMik7XG4gICAgICAgICAgICBsZXQgcCA9IGVsZW1lbnRzW3BfaW5kZXhdO1xuICAgICAgICAgICAgbGV0IEwgPSBzb3J0ZWRfZWxlbWVudHMuc2xpY2UoMCwgcF9pbmRleCk7XG4gICAgICAgICAgICBsZXQgUiA9IHNvcnRlZF9lbGVtZW50cy5zbGljZShwX2luZGV4LCBuKTtcbiAgICAgICAgICAgIGxldCByYWRpdXMgPSBNYXRoLm1heCguLi5lbGVtZW50cy5tYXAoZCA9PiB0aGlzLl9tZXRyaWMocC5lbGVtZW50LCBkLmVsZW1lbnQpKSk7XG4gICAgICAgICAgICBsZXQgQlxuICAgICAgICAgICAgaWYgKEwubGVuZ3RoID4gMCAmJiBSLmxlbmd0aCA+IDApIHsgICAgICAgICBcbiAgICAgICAgICAgICAgICBCID0gbmV3IHRoaXMuX05vZGUocCwgdGhpcy5fY29uc3RydWN0KEwpLCB0aGlzLl9jb25zdHJ1Y3QoUiksIHJhZGl1cyk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIEIgPSBuZXcgdGhpcy5fTGVhZihlbGVtZW50cyk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gQjtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtOb2RlfSBCIFxuICAgICAqIEByZXR1cm5zIHtOdW1iZXJ9XG4gICAgICovXG4gICAgX2dyZWF0ZXN0X3NwcmVhZChCKSB7XG4gICAgICAgIGxldCBkID0gQlswXS5lbGVtZW50Lmxlbmd0aDtcbiAgICAgICAgbGV0IHN0YXJ0ID0gbmV3IEFycmF5KGQpO1xuXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZDsgKytpKSB7XG4gICAgICAgICAgICBzdGFydFtpXSA9IFtJbmZpbml0eSwgLUluZmluaXR5XTtcbiAgICAgICAgfVxuXG4gICAgICAgIGxldCBzcHJlYWQgPSBCLnJlZHVjZSgoYWNjLCBjdXJyZW50KSA9PiB7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGQ7ICsraSkge1xuICAgICAgICAgICAgICAgIGFjY1tpXVswXSA9IE1hdGgubWluKGFjY1tpXVswXSwgY3VycmVudC5lbGVtZW50W2ldKTtcbiAgICAgICAgICAgICAgICBhY2NbaV1bMV0gPSBNYXRoLm1heChhY2NbaV1bMV0sIGN1cnJlbnQuZWxlbWVudFtpXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4gYWNjO1xuICAgICAgICB9LCBzdGFydCk7XG4gICAgICAgIHNwcmVhZCA9IHNwcmVhZC5tYXAoZCA9PiBkWzFdIC0gZFswXSk7XG4gICAgICAgIFxuICAgICAgICBsZXQgYyA9IDA7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZDsgKytpKSB7XG4gICAgICAgICAgICBjID0gc3ByZWFkW2ldID4gc3ByZWFkW2NdID8gaSA6IGM7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIGM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHsqfSB0IC0gcXVlcnkgZWxlbWVudC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2sgPSA1XSAtIG51bWJlciBvZiBuZWFyZXN0IG5laWdoYm9ycyB0byByZXR1cm4uXG4gICAgICogQHJldHVybnMge0hlYXB9IC0gSGVhcCBjb25zaXN0cyBvZiB0aGUge0BsaW5rIGt9IG5lYXJlc3QgbmVpZ2hib3JzLlxuICAgICAqL1xuICAgIHNlYXJjaCh0LCBrID0gNSkge1xuICAgICAgICByZXR1cm4gdGhpcy5fc2VhcmNoKHQsIGssIG5ldyBIZWFwKG51bGwsIGQgPT4gdGhpcy5fbWV0cmljKGQuZWxlbWVudCwgdCksIFwibWF4XCIpLCB0aGlzLl9yb290KTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7Kn0gdCAtIHF1ZXJ5IGVsZW1lbnQuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtrID0gNV0gLSBudW1iZXIgb2YgbmVhcmVzdCBuZWlnaGJvcnMgdG8gcmV0dXJuLlxuICAgICAqIEBwYXJhbSB7SGVhcH0gUSAtIEhlYXAgY29uc2lzdHMgb2YgdGhlIGN1cnJlbnRseSBmb3VuZCB7QGxpbmsga30gbmVhcmVzdCBuZWlnaGJvcnMuXG4gICAgICogQHBhcmFtIHtOb2RlfExlYWZ9IEIgXG4gICAgICovXG4gICAgX3NlYXJjaCh0LCBrLCBRLCBCKSB7XG4gICAgICAgIC8vIEIgaXMgTm9kZVxuICAgICAgICBpZiAoUS5sZW5ndGggPj0gayAmJiBCLnBpdm90ICYmIEIucmFkaXVzICYmIHRoaXMuX21ldHJpYyh0LCBCLnBpdm90LmVsZW1lbnQpIC0gQi5yYWRpdXMgPj0gUS5maXJzdC52YWx1ZSkge1xuICAgICAgICAgICAgcmV0dXJuIFE7XG4gICAgICAgIH0gXG4gICAgICAgIGlmIChCLmNoaWxkMSkgdGhpcy5fc2VhcmNoKHQsIGssIFEsIEIuY2hpbGQxKTtcbiAgICAgICAgaWYgKEIuY2hpbGQyKSB0aGlzLl9zZWFyY2godCwgaywgUSwgQi5jaGlsZDIpO1xuICAgICAgICBcbiAgICAgICAgLy8gQiBpcyBsZWFmXG4gICAgICAgIGlmIChCLnBvaW50cykge1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDAsIG4gPSBCLnBvaW50cy5sZW5ndGg7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgICAgICAgICBsZXQgcCA9IEIucG9pbnRzW2ldO1xuICAgICAgICAgICAgICAgIGlmIChrID4gUS5sZW5ndGgpIHtcbiAgICAgICAgICAgICAgICAgICAgUS5wdXNoKHApO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIFEucHVzaChwKTtcbiAgICAgICAgICAgICAgICAgICAgUS5wb3AoKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIFE7XG4gICAgfVxufSIsImltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBIZWFwIH0gZnJvbSBcIi4uL2RhdGFzdHJ1Y3R1cmUvaW5kZXguanNcIjtcbmltcG9ydCB7IGRpc3RhbmNlX21hdHJpeCwgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIEtOTlxuICovXG5leHBvcnQgY2xhc3MgS05OIHtcbiAgICAvKipcbiAgICAgKiBHZW5lcmF0ZXMgYSBLTk4gbGlzdCB3aXRoIGdpdmVuIHtAbGluayBlbGVtZW50c30uXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTprbm5cbiAgICAgKiBAYWxpYXMgS05OXG4gICAgICogQHBhcmFtIHtBcnJheT19IGVsZW1lbnRzIC0gRWxlbWVudHMgd2hpY2ggc2hvdWxkIGJlIGFkZGVkIHRvIHRoZSBLTk4gbGlzdFxuICAgICAqIEBwYXJhbSB7RnVuY3Rpb258XCJwcmVjb21wdXRlZFwifSBbbWV0cmljID0gZXVjbGlkZWFuXSBtZXRyaWMgaXMgZWl0aGVyIHByZWNvbXB1dGVkIG9yIGEgZnVuY3Rpb24gdG8gdXNlOiAoYSwgYikgPT4gZGlzdGFuY2VcbiAgICAgKiBAcmV0dXJucyB7S05OfVxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKGVsZW1lbnRzPW51bGwsIG1ldHJpYz1ldWNsaWRlYW4pIHtcbiAgICAgICAgdGhpcy5fbWV0cmljID0gbWV0cmljO1xuICAgICAgICB0aGlzLl9lbGVtZW50cyA9IGVsZW1lbnRzIGluc3RhbmNlb2YgTWF0cml4ID8gZWxlbWVudHMgOiBNYXRyaXguZnJvbShlbGVtZW50cyk7XG4gICAgICAgIGNvbnN0IE4gPSB0aGlzLl9lbGVtZW50cy5zaGFwZVswXTtcbiAgICAgICAgaWYgKG1ldHJpYyA9PT0gXCJwcmVjb21wdXRlZFwiKSB7XG4gICAgICAgICAgICB0aGlzLl9EID0gdGhpcy5fZWxlbWVudHMuY2xvbmUoKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRoaXMuX0QgPSBkaXN0YW5jZV9tYXRyaXgodGhpcy5fZWxlbWVudHMsIG1ldHJpYyk7XG4gICAgICAgIH1cbiAgICAgICAgdGhpcy5LTk4gPSBbXTtcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgTjsgKytyb3cpIHtcbiAgICAgICAgICAgIGNvbnN0IGRpc3RhbmNlcyA9IHRoaXMuX0Qucm93KHJvdyk7XG4gICAgICAgICAgICBjb25zdCBIID0gbmV3IEhlYXAobnVsbCwgZCA9PiBkLnZhbHVlLCBcIm1pblwiKTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgSC5wdXNoKHtcbiAgICAgICAgICAgICAgICAgICAgdmFsdWU6IGRpc3RhbmNlc1tqXSxcbiAgICAgICAgICAgICAgICAgICAgaW5kZXg6IGosXG4gICAgICAgICAgICAgICAgfSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB0aGlzLktOTi5wdXNoKEgpO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHtBcnJheXxOdW1iZXJ9IHQgLSBxdWVyeSBlbGVtZW50IG9yIGluZGV4LlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbayA9IDVdIC0gbnVtYmVyIG9mIG5lYXJlc3QgbmVpZ2hib3JzIHRvIHJldHVybi5cbiAgICAgKiBAcmV0dXJucyB7SGVhcH0gLSBIZWFwIGNvbnNpc3RzIG9mIHRoZSB7QGxpbmsga30gbmVhcmVzdCBuZWlnaGJvcnMuXG4gICAgICovXG4gICAgc2VhcmNoKHQsIGsgPSA1KSB7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgS05OID0gdGhpcy5LTk47XG4gICAgICAgIGxldCBIO1xuICAgICAgICBpZiAoQXJyYXkuaXNBcnJheSh0KSkge1xuICAgICAgICAgICAgaWYgKHRoaXMuX21ldHJpYyA9PSBcInByZWNvbXB1dGVkXCIpIHtcbiAgICAgICAgICAgICAgICB0aHJvdyBcIlNlYXJjaCBieSBxdWVyeSBlbGVtZW50IGlzIG9ubHkgcG9zc2libGUgd2hlbiBub3QgdXNpbmcgYSBwcmVjb21wdXRlZCBkaXN0YW5jZSBtYXRyaXghXCJcbiAgICAgICAgICAgIH0gXG4gICAgICAgICAgICBjb25zdCBlbGVtZW50cyA9IHRoaXMuX2VsZW1lbnRzO1xuICAgICAgICAgICAgY29uc3QgTiA9IEtOTi5sZW5ndGg7XG4gICAgICAgICAgICBsZXQgbmVhcmVzdF9lbGVtZW50X2luZGV4ID0gbnVsbDtcbiAgICAgICAgICAgIGxldCBuZWFyZXN0X2Rpc3QgPSBJbmZpbml0eTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgZWxlbWVudCA9IGVsZW1lbnRzLnJvdyhpKTtcbiAgICAgICAgICAgICAgICBjb25zdCBkaXN0ID0gbWV0cmljKHQsIGVsZW1lbnQpO1xuICAgICAgICAgICAgICAgIGlmIChkaXN0IDwgbmVhcmVzdF9kaXN0KSB7XG4gICAgICAgICAgICAgICAgICAgIG5lYXJlc3RfZWxlbWVudF9pbmRleCA9IGk7XG4gICAgICAgICAgICAgICAgICAgIG5lYXJlc3RfZGlzdCA9IGRpc3Q7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgSCA9IEtOTltuZWFyZXN0X2VsZW1lbnRfaW5kZXhdO1xuICAgICAgICB9IGVsc2UgaWYgKE51bWJlci5pc0ludGVnZXIodCkpIHtcbiAgICAgICAgICAgIEggPSBLTk5bdF1cbiAgICAgICAgfVxuXG4gICAgICAgIGxldCByZXN1bHQgPSBbXVxuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGs7ICsraSkge1xuICAgICAgICAgICAgcmVzdWx0LnB1c2goSC5wb3AoKSlcbiAgICAgICAgfVxuICAgICAgICByZXN1bHQuZm9yRWFjaChyZXMgPT4gSC5wdXNoKHJlcy5lbGVtZW50KSlcbiAgICAgICAgcmV0dXJuIHJlc3VsdFxuICAgIH0gICAgXG59XG4iLCJpbXBvcnQgeyBNYXRyaXgsIG5vcm0gfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBldWNsaWRlYW4gfSBmcm9tIFwiLi4vbWV0cmljcy9pbmRleC5qc1wiO1xuaW1wb3J0IHsgbmV1bWFpcl9zdW0gfSBmcm9tIFwiLi4vbnVtZXJpY2FsL2luZGV4LmpzXCI7XG5cbi8qKlxuICogQ29tcHV0ZXMgdGhlIFFSIERlY29tcG9zaXRpb24gb2YgdGhlIE1hdHJpeCB7QGxpbmsgQX0gdXNpbmcgR3JhbS1TY2htaWR0IHByb2Nlc3MuXG4gKiBAbWVtYmVyb2YgbW9kdWxlOmxpbmVhcl9hbGdlYnJhXG4gKiBAYWxpYXMgcXJcbiAqIEBwYXJhbSB7TWF0cml4fSBBXG4gKiBAcmV0dXJucyB7e1I6IE1hdHJpeCwgUTogTWF0cml4fX1cbiAqIEBzZWUge0BsaW5rIGh0dHBzOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL1FSX2RlY29tcG9zaXRpb24jVXNpbmdfdGhlX0dyYW0lRTIlODAlOTNTY2htaWR0X3Byb2Nlc3N9XG4gKi9cbmV4cG9ydCBkZWZhdWx0IGZ1bmN0aW9uKEEpIHtcbiAgICBjb25zdCBbcm93cywgY29sc10gPSBBLnNoYXBlO1xuICAgIGNvbnN0IFEgPSBuZXcgTWF0cml4KHJvd3MsIGNvbHMsIFwiaWRlbnRpdHlcIik7XG4gICAgY29uc3QgUiA9IG5ldyBNYXRyaXgoY29scywgY29scywgMCk7XG5cbiAgICBmb3IgKGxldCBqID0gMDsgaiA8IGNvbHM7ICsraikge1xuICAgICAgICBsZXQgdiA9IEEuY29sKGopO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGo7ICsraSkge1xuICAgICAgICAgICAgY29uc3QgcSA9IFEuY29sKGkpO1xuICAgICAgICAgICAgY29uc3QgcV9kb3RfdiA9IG5ldW1haXJfc3VtKHEubWFwKChxXywgaykgPT4gcV8gKiB2W2tdKSk7XG4gICAgICAgICAgICBSLnNldF9lbnRyeShpLGosIHFfZG90X3YpO1xuICAgICAgICAgICAgdiA9IHYubWFwKCh2XywgaykgPT4gdl8gLSBxX2RvdF92ICogcVtrXSk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3Qgdl9ub3JtID0gbm9ybSh2LCBldWNsaWRlYW4pO1xuICAgICAgICBmb3IgKGxldCBrID0gMDsgayA8IHJvd3M7ICsraykge1xuICAgICAgICAgICAgUS5zZXRfZW50cnkoaywgaiwgdltrXSAvIHZfbm9ybSk7XG4gICAgICAgIH1cbiAgICAgICAgUi5zZXRfZW50cnkoaixqLCB2X25vcm0pXG4gICAgfVxuICAgIHJldHVybiB7XCJSXCI6IFIsIFwiUVwiOiBRfTtcbn1cblxuIiwiaW1wb3J0IHsgcXIgfSBmcm9tIFwiLi9pbmRleC5qc1wiO1xuaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgUmFuZG9taXplciB9IGZyb20gXCIuLi91dGlsL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBuZXVtYWlyX3N1bSB9IGZyb20gXCIuLi9udW1lcmljYWwvaW5kZXguanNcIjtcblxuLyoqXG4gKiBDb21wdXRlcyB0aGUge0BsaW5rIGt9IGJpZ2dlc3QgRWlnZW52ZWN0b3JzIGFuZCBFaWdlbnZhbHVlcyBmcm9tIE1hdHJpeCB7QGxpbmsgQX0gd2l0aCB0aGUgUVItQWxnb3JpdGhtLlxuICogQHBhcmFtIHtNYXRyaXh9IEEgLSBUaGUgTWF0cml4XG4gKiBAcGFyYW0ge051bWJlcn0gayAtIFRoZSBudW1iZXIgb2YgZWlnZW52ZWN0b3JzIGFuZCBlaWdlbnZhbHVlcyB0byBjb21wdXRlLlxuICogQHBhcmFtIHtOdW1iZXJ9IFttYXhfaXRlcmF0aW9ucz0xMDBdIC0gVGhlIG51bWJlciBvZiBtYXhpdW11bSBpdGVyYXRpb25zIHRoZSBhbGdvcml0aG0gc2hvdWxkIHJ1bi5cbiAqIEBwYXJhbSB7TnVtYmVyfFJhbmRvbWl6ZXJ9IFtzZWVkPTEyMTJdIC0gVGhlIHNlZWQgdmFsdWUgb3IgYSByYW5kb21pemVyIHVzZWQgaW4gdGhlIGFsZ29yaXRobS5cbiAqIEByZXR1cm5zIHt7ZWlnZW52YWx1ZXM6IEFycmF5LCBlaWdlbnZlY3RvcnM6IEFycmF5fX0gLSBUaGUge0BsaW5rIGt9IGJpZ2dlc3QgZWlnZW52ZWN0b3JzIGFuZCBlaWdlbnZhbHVlcyBvZiBNYXRyaXgge0BsaW5rIEF9LlxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbihBLCBrID0gMiwgbWF4X2l0ZXJhdGlvbnM9MTAwLCBzZWVkPTEyMTIpIHtcbiAgICBjb25zdCByYW5kb21pemVyID0gc2VlZCBpbnN0YW5jZW9mIFJhbmRvbWl6ZXIgPyBzZWVkIDogbmV3IFJhbmRvbWl6ZXIoc2VlZCk7XG4gICAgaWYgKCEoQSBpbnN0YW5jZW9mIE1hdHJpeCkpIEEgPSBNYXRyaXguZnJvbShBKTtcbiAgICBjb25zdCBuID0gQS5zaGFwZVswXVxuICAgIGxldCB7IFE6IFEsIFI6IFIgfSA9IHFyKG5ldyBNYXRyaXgobiwgaywgKCkgPT4gcmFuZG9taXplci5yYW5kb20pKTtcbiAgICB3aGlsZSAobWF4X2l0ZXJhdGlvbnMtLSkge1xuICAgICAgICBjb25zdCBvbGRSID0gUi5jbG9uZSgpO1xuICAgICAgICBjb25zdCBaID0gQS5kb3QoUSk7XG4gICAgICAgIGNvbnN0IFFSID0gcXIoWik7IFxuICAgICAgICBRID0gUVIuUTtcbiAgICAgICAgUiA9IFFSLlI7XG4gICAgICAgIGlmIChuZXVtYWlyX3N1bShSLnN1YihvbGRSKS5kaWFnKSAvIG4gPCAxZS0xMikge1xuICAgICAgICAgICAgbWF4X2l0ZXJhdGlvbnMgPSAwO1xuICAgICAgICB9ICAgICAgICBcbiAgICB9XG5cbiAgICBjb25zdCBlaWdlbnZhbHVlcyA9IFIuZGlhZztcbiAgICBjb25zdCBlaWdlbnZlY3RvcnMgPSBRLnRyYW5zcG9zZSgpLnRvMmRBcnJheTtcbiAgICByZXR1cm4ge1xuICAgICAgICBcImVpZ2VudmFsdWVzXCI6IGVpZ2VudmFsdWVzLFxuICAgICAgICBcImVpZ2VudmVjdG9yc1wiOiBlaWdlbnZlY3RvcnMsXG4gICAgfTtcbn1cblxuIiwiaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IFJhbmRvbWl6ZXIgfSBmcm9tIFwiLi4vdXRpbC9pbmRleC5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIERSXG4gKiBAYm9ycm93cyBEUiNwYXJhbWV0ZXIgYXMgRFIjcGFyYVxuICogQGJvcnJvd3MgRFIjcGFyYW1ldGVyIGFzIERSI3BcbiAqL1xuZXhwb3J0IGNsYXNzIERSIHtcbiAgICAvL3N0YXRpYyBwYXJhbWV0ZXJfbGlzdCA9IFtdO1xuICAgIGdldCBwYXJhbWV0ZXJfbGlzdCgpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX3BhcmFtZXRlcl9saXN0O1xuICAgIH1cblxuICAgIHNldCBwYXJhbWV0ZXJfbGlzdChsaXN0KSB7XG4gICAgICAgIHRoaXMuX3BhcmFtZXRlcl9saXN0ID0gbGlzdDtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuICAgIC8qKlxuICAgICAqXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTpkaW1lbnNpb25hbGl0eV9yZWR1Y3Rpb25cbiAgICAgKiBAYWxpYXMgRFJcbiAgICAgKiBAcGFyYW0ge01hdHJpeHxBcnJheTxBcnJheTxOdW1iZXI+Pn0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuXG4gICAgICogQHBhcmFtIHtudW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHBhcmFtIHtmdW5jdGlvbn0gW21ldHJpYyA9IGV1Y2xpZGVhbl0gLSB0aGUgbWV0cmljIHdoaWNoIGRlZmluZXMgdGhlIGRpc3RhbmNlIGJldHdlZW4gdHdvIHBvaW50cy5cbiAgICAgKiBAcGFyYW0ge3NlZWR9IFtzZWVkID0gMTIxMl0gLSB0aGUgc2VlZCB2YWx1ZSBmb3IgdGhlIHJhbmRvbSBudW1iZXIgZ2VuZXJhdG9yLlxuICAgICAqIEByZXR1cm5zIHtEUn1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBkID0gMiwgbWV0cmljID0gZXVjbGlkZWFuLCBzZWVkID0gMTIxMikge1xuICAgICAgICBpZiAoQXJyYXkuaXNBcnJheShYKSkge1xuICAgICAgICAgICAgdGhpcy5fdHlwZSA9IFwiYXJyYXlcIjtcbiAgICAgICAgICAgIHRoaXMuWCA9IE1hdHJpeC5mcm9tKFgpO1xuICAgICAgICB9IGVsc2UgaWYgKFggaW5zdGFuY2VvZiBNYXRyaXgpIHtcbiAgICAgICAgICAgIHRoaXMuX3R5cGUgPSBcIm1hdHJpeFwiO1xuICAgICAgICAgICAgdGhpcy5YID0gWDtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcIm5vIHZhbGlkIHR5cGUgZm9yIFhcIik7XG4gICAgICAgIH1cbiAgICAgICAgW3RoaXMuX04sIHRoaXMuX0RdID0gdGhpcy5YLnNoYXBlO1xuICAgICAgICB0aGlzLl9kID0gZDtcbiAgICAgICAgdGhpcy5fbWV0cmljID0gbWV0cmljO1xuICAgICAgICB0aGlzLl9zZWVkID0gc2VlZDtcbiAgICAgICAgdGhpcy5fcmFuZG9taXplciA9IG5ldyBSYW5kb21pemVyKHNlZWQpO1xuICAgICAgICB0aGlzLl9pc19pbml0aWFsaXplZCA9IGZhbHNlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBTZXQgYW5kIGdldCBwYXJhbWV0ZXJzXG4gICAgICogQHBhcmFtIHtTdHJpbmd9IG5hbWUgLSBuYW1lIG9mIHRoZSBwYXJhbWV0ZXIuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFt2YWx1ZSA9IG51bGxdIC0gdmFsdWUgb2YgdGhlIHBhcmFtZXRlciB0byBzZXQsIGlmIDxjb2RlPnZhbHVlID09IG51bGw8L2NvZGU+IHRoZW4gcmV0dXJuIGFjdHVhbCBwYXJhbWV0ZXIgdmFsdWUuXG4gICAgICogQG1lbWJlcm9mIERSXG4gICAgICovXG4gICAgcGFyYW1ldGVyKG5hbWUsIHZhbHVlID0gbnVsbCkge1xuICAgICAgICBpZiAoIXRoaXMucGFyYW1ldGVyX2xpc3QuaW5jbHVkZXMobmFtZSkpIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgJHtuYW1lfSBpcyBub3QgYSB2YWxpZCBwYXJhbWV0ZXIhYCk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHZhbHVlKSB7XG4gICAgICAgICAgICB0aGlzW2BfJHtuYW1lfWBdID0gdmFsdWU7XG4gICAgICAgICAgICB0aGlzLl9pc19pbml0aWFsaXplZCA9IGZhbHNlO1xuICAgICAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICByZXR1cm4gdGhpc1tgXyR7bmFtZX1gXTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHBhcmEobmFtZSwgdmFsdWUgPSBudWxsKSB7XG4gICAgICAgIHJldHVybiB0aGlzLnBhcmFtZXRlcihuYW1lLCB2YWx1ZSk7XG4gICAgfVxuXG4gICAgcChuYW1lLCB2YWx1ZSA9IG51bGwpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMucGFyYW1ldGVyKG5hbWUsIHZhbHVlKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fSBSZXR1cm5zIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqL1xuICAgIHRyYW5zZm9ybSgpIHtcbiAgICAgICAgdGhpcy5jaGVja19pbml0KCk7XG4gICAgICAgIHJldHVybiB0aGlzLnByb2plY3Rpb247XG4gICAgfVxuXG4gICAgKmdlbmVyYXRvcigpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMudHJhbnNmb3JtKCk7XG4gICAgfVxuXG4gICAgY2hlY2tfaW5pdCgpIHtcbiAgICAgICAgaWYgKCF0aGlzLl9pc19pbml0aWFsaXplZCAmJiB0eXBlb2YgdGhpcy5pbml0ID09PSBcImZ1bmN0aW9uXCIpIHtcbiAgICAgICAgICAgIHRoaXMuaW5pdCgpO1xuICAgICAgICAgICAgdGhpcy5faXNfaW5pdGlhbGl6ZWQgPSB0cnVlO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHJldHVybnMge01hdHJpeH0gUmV0dXJucyB0aGUgcHJvamVjdGlvbi5cbiAgICAgKi9cbiAgICBnZXQgcHJvamVjdGlvbigpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX3R5cGUgPT09IFwibWF0cml4XCIgPyB0aGlzLlkgOiB0aGlzLlkudG8yZEFycmF5O1xuICAgIH1cblxuICAgIGFzeW5jIHRyYW5zZm9ybV9hc3luYyguLi5hcmdzKSB7XG4gICAgICAgIGNvbnN0IGRyID0gbmV3IHRoaXMoLi4uYXJncyk7XG4gICAgICAgIHJldHVybiBkci50cmFuc2Zvcm0oKTtcbiAgICB9XG5cbiAgICBzdGF0aWMgdHJhbnNmb3JtKC4uLmFyZ3MpIHtcbiAgICAgICAgbGV0IGRyID0gbmV3IHRoaXMoLi4uYXJncyk7XG4gICAgICAgIHJldHVybiBkci50cmFuc2Zvcm0oKTtcbiAgICB9XG5cbiAgICBzdGF0aWMgYXN5bmMgdHJhbnNmb3JtX2FzeW5jKC4uLmFyZ3MpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMudHJhbnNmb3JtKC4uLmFyZ3MpO1xuICAgIH1cblxuICAgIHN0YXRpYyAqZ2VuZXJhdG9yKC4uLmFyZ3MpIHtcbiAgICAgICAgY29uc3QgZHIgPSBuZXcgdGhpcyguLi5hcmdzKTtcbiAgICAgICAgY29uc3QgZ2VuID0gZHIuZ2VuZXJhdG9yKCk7XG4gICAgICAgIGZvciAoY29uc3QgcmVzIG9mIGdlbikge1xuICAgICAgICAgICAgeWllbGQgcmVzO1xuICAgICAgICB9XG4gICAgfVxufVxuIiwiaW1wb3J0IHsgc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uIH0gZnJvbSBcIi4uL2xpbmVhcl9hbGdlYnJhL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBNYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgUENBXG4gKiBAYXVnbWVudHMgRFJcbiAqL1xuZXhwb3J0IGNsYXNzIFBDQSBleHRlbmRzIERSIHtcbiAgICAvKipcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBQQ0FcbiAgICAgKiBAcGFyYW0ge01hdHJpeHxBcnJheTxBcnJheTxOdW1iZXI+Pn0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHJldHVybnMge1BDQX1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBkID0gMikge1xuICAgICAgICBzdXBlcihYLCBkKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogVHJhbnNmb3JtcyB0aGUgaW5wdXRkYXRhIHtAbGluayBYfSB0byBkaW1lbnNpb25hbGl0eSB7QGxpbmsgZH0uIElmIHBhcmFtZXRlciB7QGxpbmsgQX0gaXMgZ2l2ZW4sIHRoZW4gcHJvamVjdCB7QGxpbmsgQX0gd2l0aCB0aGUgcHJpbmNpcGFsIGNvbXBvbmVudHMgb2Yge0BsaW5rIFh9LlxuICAgICAqIEBwYXJhbSB7bnVsbHxNYXRyaXh8QXJyYXl9IFtBID0gbnVsbF0gLSBJZiBnaXZlbiwgdGhlIGRhdGEgdG8gcHJvamVjdC5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fEFycmF5fSAtIFRoZSBwcm9qZWN0ZWQgZGF0YS5cbiAgICAgKi9cbiAgICB0cmFuc2Zvcm0oQSA9IG51bGwpIHtcbiAgICAgICAgY29uc3QgViA9IHRoaXMucHJpbmNpcGFsX2NvbXBvbmVudHMoKTtcbiAgICAgICAgaWYgKEEgPT0gbnVsbCkge1xuICAgICAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgICAgIHRoaXMuWSA9IFguZG90KFYpO1xuICAgICAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICAgICAgfSBlbHNlIGlmIChBcnJheS5pc0FycmF5KEEpKSB7XG4gICAgICAgICAgICByZXR1cm4gZHJ1aWQuTWF0cml4LmZyb20oQSkuZG90KFYpLmFzQXJyYXk7XG4gICAgICAgIH0gZWxzZSBpZiAoQSBpbnN0YW5jZW9mIE1hdHJpeCkge1xuICAgICAgICAgICAgcmV0dXJuIEEuZG90KFYpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKFwiTm8gdmFsaWQgdHlwZSBmb3IgQSFcIik7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyB0aGUge0BsaW5rIGR9IHByaW5jaXBhbCBjb21wb25lbnRzIG9mIE1hdHJpeCB7QGxpbmsgWH0uXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBwcmluY2lwYWxfY29tcG9uZW50cygpIHtcbiAgICAgICAgaWYgKHRoaXMuVikge1xuICAgICAgICAgICAgcmV0dXJuIHRoaXMuVjtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCBYID0gdGhpcy5YO1xuICAgICAgICBjb25zdCBtZWFucyA9IE1hdHJpeC5mcm9tKFgubWVhbkNvbHMpO1xuICAgICAgICBjb25zdCBYX2NlbnQgPSBYLnN1YihtZWFucyk7XG4gICAgICAgIGNvbnN0IEMgPSBYX2NlbnQudHJhbnNwb3NlKCkuZG90KFhfY2VudCk7XG4gICAgICAgIGNvbnN0IHsgZWlnZW52ZWN0b3JzOiBWIH0gPSBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb24oQywgdGhpcy5fZCk7XG4gICAgICAgIHRoaXMuViA9IE1hdHJpeC5mcm9tKFYpLnRyYW5zcG9zZSgpO1xuICAgICAgICByZXR1cm4gdGhpcy5WO1xuICAgIH1cbn1cbiIsImltcG9ydCB7IHNpbXVsdGFuZW91c19wb3dlcml0ZXJhdGlvbn0gZnJvbSBcIi4uL2xpbmVhcl9hbGdlYnJhL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBkaXN0YW5jZV9tYXRyaXgsIE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgTURTXG4gKiBAZXh0ZW5kcyBEUlxuICovXG5leHBvcnQgY2xhc3MgTURTIGV4dGVuZHMgRFJ7XG4gICAgLyoqXG4gICAgICogQ2xhc3NpY2FsIE1EUy5cbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBNRFNcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHBhcmFtIHtGdW5jdGlvbnxcInByZWNvbXB1dGVkXCJ9IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gdGhlIG1ldHJpYyB3aGljaCBkZWZpbmVzIHRoZSBkaXN0YW5jZSBiZXR3ZWVuIHR3byBwb2ludHMuICBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW3NlZWQgPSAxMjEyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBkPTIsIG1ldHJpYz1ldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBUcmFuc2Zvcm1zIHRoZSBpbnB1dGRhdGEge0BsaW5rIFh9IHRvIGRpbWVuc2lvbmFsaXR5IHtAbGluayBkfS5cbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fEFycmF5fVxuICAgICAqL1xuICAgIHRyYW5zZm9ybSgpIHtcbiAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgY29uc3Qgcm93cyA9IFguc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgQSA9IG1ldHJpYyA9PT0gXCJwcmVjb21wdXRlZFwiID8gWCA6IGRpc3RhbmNlX21hdHJpeChYLCBtZXRyaWMpOyBcbiAgICAgICAgY29uc3QgYWlfID0gQS5tZWFuQ29scztcbiAgICAgICAgY29uc3QgYV9qID0gQS5tZWFuUm93cztcbiAgICAgICAgY29uc3QgYV9fID0gQS5tZWFuO1xuXG4gICAgICAgIHRoaXMuX2RfWCA9IEE7XG4gICAgICAgIGNvbnN0IEIgPSBuZXcgTWF0cml4KHJvd3MsIHJvd3MsIChpLCBqKSA9PiAoQS5lbnRyeShpLCBqKSAtIGFpX1tpXSAtIGFfaltqXSArIGFfXykpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICBjb25zdCB7IGVpZ2VudmVjdG9yczogViB9ID0gc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uKEIsIHRoaXMuX2QpO1xuICAgICAgICB0aGlzLlkgPSBNYXRyaXguZnJvbShWKS50cmFuc3Bvc2UoKVxuICAgICAgICBcbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcmV0dXJucyB7TnVtYmVyfSAtIHRoZSBzdHJlc3Mgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICovXG4gICAgc3RyZXNzKCkge1xuICAgICAgICBjb25zdCBOID0gdGhpcy5YLnNoYXBlWzBdO1xuICAgICAgICBjb25zdCBZID0gdGhpcy5ZO1xuICAgICAgICBjb25zdCBkX1ggPSB0aGlzLl9kX1g7IC8qbmV3IE1hdHJpeCgpO1xuICAgICAgICBkX1guc2hhcGUgPSBbTiwgTiwgKGksIGopID0+IHtcbiAgICAgICAgICAgIHJldHVybiBpIDwgaiA/IG1ldHJpYyhYLnJvdyhpKSwgWC5yb3coaikpIDogZF9YLmVudHJ5KGosIGkpO1xuICAgICAgICB9XSovXG4gICAgICAgIGNvbnN0IGRfWSA9IG5ldyBNYXRyaXgoKTtcbiAgICAgICAgZF9ZLnNoYXBlID0gW04sIE4sIChpLCBqKSA9PiB7XG4gICAgICAgICAgICByZXR1cm4gaSA8IGogPyBldWNsaWRlYW4oWS5yb3coaSksIFkucm93KGopKSA6IGRfWS5lbnRyeShqLCBpKTtcbiAgICAgICAgfV1cbiAgICAgICAgbGV0IHRvcF9zdW0gPSAwO1xuICAgICAgICBsZXQgYm90dG9tX3N1bSA9IDA7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gaSArIDE7IGogPCBOOyArK2opIHtcbiAgICAgICAgICAgICAgICB0b3Bfc3VtICs9IE1hdGgucG93KGRfWC5lbnRyeShpLCBqKSAtIGRfWS5lbnRyeShpLCBqKSwgMik7XG4gICAgICAgICAgICAgICAgYm90dG9tX3N1bSArPSBNYXRoLnBvdyhkX1guZW50cnkoaSwgaiksIDIpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiBNYXRoLnNxcnQodG9wX3N1bSAvIGJvdHRvbV9zdW0pO1xuICAgIH1cbn0iLCJpbXBvcnQgeyBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb259IGZyb20gXCIuLi9saW5lYXJfYWxnZWJyYS9pbmRleC5qc1wiO1xuaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IEhlYXAgfSBmcm9tIFwiLi4vZGF0YXN0cnVjdHVyZS9pbmRleC5qc1wiO1xuaW1wb3J0IHsgRFIgfSBmcm9tIFwiLi9EUi5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIElTT01BUFxuICogQGV4dGVuZHMgRFJcbiAqL1xuZXhwb3J0IGNsYXNzIElTT01BUCBleHRlbmRzIERSIHtcbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBJU09NQVBcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuIFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBuZWlnaGJvcnMgLSB0aGUgbnVtYmVyIG9mIG5laWdoYm9ycyB7QGxpbmsgSVNPTUFQfSBzaG91bGQgdXNlIHRvIHByb2plY3QgdGhlIGRhdGEuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uIFxuICAgICAqIEBwYXJhbSB7RnVuY3Rpb259IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gdGhlIG1ldHJpYyB3aGljaCBkZWZpbmVzIHRoZSBkaXN0YW5jZSBiZXR3ZWVuIHR3byBwb2ludHMuIFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbc2VlZCA9IDEyMTJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKFgsIG5laWdoYm9ycywgZCA9IDIsIG1ldHJpYyA9IGV1Y2xpZGVhbiwgc2VlZD0xMjEyKSB7XG4gICAgICAgIHN1cGVyKFgsIGQsIG1ldHJpYywgc2VlZCk7XG4gICAgICAgIHN1cGVyLnBhcmFtZXRlcl9saXN0ID0gW1wia1wiXTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJrXCIsIE1hdGgubWluKG5laWdoYm9ycyA/PyBNYXRoLm1heChNYXRoLmZsb29yKHRoaXMuWC5zaGFwZVswXSAvIDEwKSwgMiksIHRoaXMuX04gLTEpKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQ29tcHV0ZXMgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHJldHVybnMge01hdHJpeH0gUmV0dXJucyB0aGUgcHJvamVjdGlvbi5cbiAgICAgKi9cbiAgICB0cmFuc2Zvcm0oKSB7XG4gICAgICAgIHRoaXMuY2hlY2tfaW5pdCgpO1xuICAgICAgICBjb25zdCBYID0gdGhpcy5YO1xuICAgICAgICBjb25zdCByb3dzID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgbWV0cmljID0gdGhpcy5fbWV0cmljO1xuICAgICAgICAvLyBUT0RPOiBtYWtlIGtubiBleHRlcm4gYW5kIHBhcmFtZXRlciBmb3IgY29uc3RydWN0b3Igb3IgdHJhbnNmb3JtP1xuICAgICAgICBjb25zdCBEID0gbmV3IE1hdHJpeCgpO1xuICAgICAgICBELnNoYXBlID0gW3Jvd3MsIHJvd3MsIChpLGopID0+IGkgPD0gaiA/IG1ldHJpYyhYLnJvdyhpKSwgWC5yb3coaikpIDogRC5lbnRyeShqLGkpXVxuICAgICAgICBjb25zdCBrTmVhcmVzdE5laWdoYm9ycyA9IFtdO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHJvd3M7ICsraSkge1xuICAgICAgICAgICAgY29uc3Qgcm93ID0gW107XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IHJvd3M7ICsraikge1xuICAgICAgICAgICAgICAgIHJvdy5wdXNoKHtcbiAgICAgICAgICAgICAgICAgICAgXCJpbmRleFwiOiBqLFxuICAgICAgICAgICAgICAgICAgICBcImRpc3RhbmNlXCI6IEQuZW50cnkoaSwgaiksXG4gICAgICAgICAgICAgICAgfSlcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGNvbnN0IEggPSBuZXcgSGVhcChyb3csIGQgPT4gZC5kaXN0YW5jZSwgXCJtaW5cIik7XG4gICAgICAgICAgICBrTmVhcmVzdE5laWdoYm9ycy5wdXNoKEgudG9BcnJheSgpLnNsaWNlKDEsIHRoaXMuX2sgKyAxKSlcbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgLypEID0gZGlqa3N0cmEoa05lYXJlc3ROZWlnaGJvcnMpOyovXG4gICAgICAgIC8vIGNvbXB1dGUgc2hvcnRlc3QgcGF0aHNcbiAgICAgICAgLy8gVE9ETzogbWFrZSBleHRlcm5cbiAgICAgICAgLyoqIEBzZWUge0BsaW5rIGh0dHBzOi8vZW4ud2lraXBlZGlhLm9yZy93aWtpL0Zsb3lkJUUyJTgwJTkzV2Fyc2hhbGxfYWxnb3JpdGhtfSAqL1xuICAgICAgICBjb25zdCBHID0gbmV3IE1hdHJpeChyb3dzLCByb3dzLCAoaSxqKSA9PiB7XG4gICAgICAgICAgICBjb25zdCBvdGhlciA9IGtOZWFyZXN0TmVpZ2hib3JzW2ldLmZpbmQobiA9PiBuLmluZGV4ID09PSBqKTtcbiAgICAgICAgICAgIHJldHVybiBvdGhlciA/IG90aGVyLmRpc3RhbmNlIDogSW5maW5pdHlcbiAgICAgICAgfSk7XG5cbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCByb3dzOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgcm93czsgKytqKSB7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgayA9IDA7IGsgPCByb3dzOyArK2spIHtcbiAgICAgICAgICAgICAgICAgICAgRy5zZXRfZW50cnkoaSwgaiwgTWF0aC5taW4oRy5lbnRyeShpLCBqKSwgRy5lbnRyeShpLCBrKSArIEcuZW50cnkoaywgaikpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgXG4gICAgICAgIGxldCBhaV8gPSBuZXcgRmxvYXQ2NEFycmF5KHJvd3MpO1xuICAgICAgICBsZXQgYV9qID0gbmV3IEZsb2F0NjRBcnJheShyb3dzKTtcbiAgICAgICAgbGV0IGFfXyA9IDA7XG4gICAgICAgIGxldCBBID0gbmV3IE1hdHJpeChyb3dzLCByb3dzLCAoaSxqKSA9PiB7XG4gICAgICAgICAgICBsZXQgdmFsID0gRy5lbnRyeShpLCBqKTtcbiAgICAgICAgICAgIHZhbCA9IHZhbCA9PT0gSW5maW5pdHkgPyAwIDogdmFsO1xuICAgICAgICAgICAgYWlfW2ldICs9IHZhbDtcbiAgICAgICAgICAgIGFfaltqXSArPSB2YWw7XG4gICAgICAgICAgICBhX18gKz0gdmFsO1xuICAgICAgICAgICAgcmV0dXJuIHZhbDtcbiAgICAgICAgfSk7XG4gICAgICAgIFxuICAgICAgICBhaV8gPSBhaV8ubWFwKHYgPT4gdiAvIHJvd3MpO1xuICAgICAgICBhX2ogPSBhX2oubWFwKHYgPT4gdiAvIHJvd3MpO1xuICAgICAgICBhX18gLz0gKHJvd3MgKiogMik7XG4gICAgICAgIGNvbnN0IEIgPSBuZXcgTWF0cml4KHJvd3MsIHJvd3MsIChpLGopID0+IChBLmVudHJ5KGksaikgLSBhaV9baV0gLSBhX2pbal0gKyBhX18pKTtcbiAgICAgICAgICAgICBcbiAgICAgICAgLy8gY29tcHV0ZSBkIGVpZ2VudmVjdG9yc1xuICAgICAgICBjb25zdCB7IGVpZ2VudmVjdG9yczogViB9ID0gc2ltdWx0YW5lb3VzX3Bvd2VyaXRlcmF0aW9uKEIsIHRoaXMuX2QpO1xuICAgICAgICB0aGlzLlkgPSBNYXRyaXguZnJvbShWKS50cmFuc3Bvc2UoKTtcbiAgICAgICAgLy8gcmV0dXJuIGVtYmVkZGluZ1xuICAgICAgICByZXR1cm4gdGhpcy5wcm9qZWN0aW9uO1xuICAgIH1cblxuXG59IiwiaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IERSIH0gZnJvbSBcIi4vRFIuanNcIjtcbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgRkFTVE1BUFxuICogQGV4dGVuZHMgRFJcbiAqL1xuZXhwb3J0IGNsYXNzIEZBU1RNQVAgZXh0ZW5kcyBEUntcbiAgICAvKipcbiAgICAgKiBGYXN0TWFwOiBhIGZhc3QgYWxnb3JpdGhtIGZvciBpbmRleGluZywgZGF0YS1taW5pbmcgYW5kIHZpc3VhbGl6YXRpb24gb2YgdHJhZGl0aW9uYWwgYW5kIG11bHRpbWVkaWEgZGF0YXNldHNcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBGQVNUTUFQXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IFggLSB0aGUgaGlnaC1kaW1lbnNpb25hbCBkYXRhLiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2QgPSAyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzLiAgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzZWVkID0gMTIxMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHJldHVybnMge0ZBU1RNQVB9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9kb2kub3JnLzEwLjExNDUvMjIzNzg0LjIyMzgxMn1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBkPTIsIG1ldHJpYz1ldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDaG9vc2VzIHR3byBwb2ludHMgd2hpY2ggYXJlIHRoZSBtb3N0IGRpc3RhbnQgaW4gdGhlIGFjdHVhbCBwcm9qZWN0aW9uLlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtmdW5jdGlvbn0gZGlzdCBcbiAgICAgKiBAcmV0dXJucyB7QXJyYXl9IEFuIGFycmF5IGNvbnNpc3Rpbmcgb2YgZmlyc3QgaW5kZXgsIHNlY29uZCBpbmRleCwgYW5kIGRpc3RhbmNlIGJldHdlZW4gdGhlIHR3byBwb2ludHMuXG4gICAgICovXG4gICAgX2Nob29zZV9kaXN0YW50X29iamVjdHMoZGlzdCkge1xuICAgICAgICBjb25zdCBYID0gdGhpcy5YO1xuICAgICAgICBjb25zdCBOID0gWC5zaGFwZVswXTtcbiAgICAgICAgbGV0IGFfaW5kZXggPSB0aGlzLl9yYW5kb21pemVyLnJhbmRvbV9pbnQgJSBOIC0gMTtcbiAgICAgICAgbGV0IGJfaW5kZXggPSBudWxsO1xuICAgICAgICBsZXQgbWF4X2Rpc3QgPSAtSW5maW5pdHk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBkX2FpID0gZGlzdChhX2luZGV4LCBpKVxuICAgICAgICAgICAgaWYgKGRfYWkgPiBtYXhfZGlzdCkge1xuICAgICAgICAgICAgICAgIG1heF9kaXN0ID0gZF9haTtcbiAgICAgICAgICAgICAgICBiX2luZGV4ID0gaTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBtYXhfZGlzdCA9IC1JbmZpbml0eTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGNvbnN0IGRfYmkgPSBkaXN0KGJfaW5kZXgsIGkpXG4gICAgICAgICAgICBpZiAoZF9iaSA+IG1heF9kaXN0KSB7XG4gICAgICAgICAgICAgICAgbWF4X2Rpc3QgPSBkX2JpO1xuICAgICAgICAgICAgICAgIGFfaW5kZXggPSBpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiBbYV9pbmRleCwgYl9pbmRleCwgbWF4X2Rpc3RdO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENvbXB1dGVzIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9IFRoZSB7QGxpbmsgZH0tZGltZW5zaW9uYWwgcHJvamVjdGlvbiBvZiB0aGUgZGF0YSBtYXRyaXgge0BsaW5rIFh9LlxuICAgICAqL1xuICAgIHRyYW5zZm9ybSgpIHtcbiAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgY29uc3QgTiA9IFguc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IGQgPSB0aGlzLl9kO1xuICAgICAgICBjb25zdCBtZXRyaWMgPSB0aGlzLl9tZXRyaWM7XG4gICAgICAgIGNvbnN0IFkgPSBuZXcgTWF0cml4KE4sIGQsIDApO1xuICAgICAgICBsZXQgZGlzdCA9IChhLCBiKSA9PiBtZXRyaWMoWC5yb3coYSksIFgucm93KGIpKTtcblxuICAgICAgICBmb3IgKGxldCBfY29sID0gMDsgX2NvbCA8IGQ7ICsrX2NvbCkge1xuICAgICAgICAgICAgbGV0IG9sZF9kaXN0ID0gZGlzdDtcbiAgICAgICAgICAgIC8vIGNob29zZSBwaXZvdCBvYmplY3RzXG4gICAgICAgICAgICBjb25zdCBbYV9pbmRleCwgYl9pbmRleCwgZF9hYl0gPSB0aGlzLl9jaG9vc2VfZGlzdGFudF9vYmplY3RzKGRpc3QpO1xuICAgICAgICAgICAgLy8gcmVjb3JkIGlkIG9mIHBpdm90IG9iamVjdHNcbiAgICAgICAgICAgIC8vUEFbMF0ucHVzaChhX2luZGV4KTtcbiAgICAgICAgICAgIC8vUEFbMV0ucHVzaChiX2luZGV4KTtcbiAgICAgICAgICAgIC8qIGlmIChkX2FiID09PSAwKSB7XG4gICAgICAgICAgICAgICAgLy8gYmVjYXVzZSBhbGwgaW50ZXItb2JqZWN0IGRpc3RhbmNlcyBhcmUgemVyb3NcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgICAgICAgICBZLnNldF9lbnRyeShpLCBfY29sLCAwKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2UgeyAqL1xuICAgICAgICAgICAgaWYgKGRfYWIgIT09IDApIHtcbiAgICAgICAgICAgICAgICAvLyBwcm9qZWN0IHRoZSBvYmplY3RzIG9uIHRoZSBsaW5lIChPX2EsIE9fYilcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBkX2FpID0gZGlzdChhX2luZGV4LCBpKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgZF9iaSA9IGRpc3QoYl9pbmRleCwgaSk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHlfaSA9IChkX2FpICoqIDIgKyBkX2FiICoqIDIgLSBkX2JpICoqIDIpIC8gKDIgKiBkX2FiKTtcbiAgICAgICAgICAgICAgICAgICAgWS5zZXRfZW50cnkoaSwgX2NvbCwgeV9pKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgLy8gY29uc2lkZXIgdGhlIHByb2plY3Rpb25zIG9mIHRoZSBvYmplY3RzIG9uIGFcbiAgICAgICAgICAgICAgICAvLyBoeXBlcnBsYW5lIHBlcnBlbmRpY2x1YXIgdG8gdGhlIGxpbmUgKGEsIGIpO1xuICAgICAgICAgICAgICAgIC8vIHRoZSBkaXN0YW5jZSBmdW5jdGlvbiBEJygpIGJldHdlZW4gdHdvIFxuICAgICAgICAgICAgICAgIC8vIHByb2plY3Rpb25zIGlzIGdpdmVuIGJ5IEVxLjRcbiAgICAgICAgICAgICAgICBkaXN0ID0gKGEsIGIpID0+IE1hdGguc3FydChvbGRfZGlzdChhLCBiKSAqKiAyIC0gKFkuZW50cnkoYSwgX2NvbCkgLSBZLmVudHJ5KGIsIF9jb2wpKSAqKiAyKVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIC8vIHJldHVybiBlbWJlZGRpbmdcbiAgICAgICAgdGhpcy5ZID0gWTtcbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG59IiwiaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IHNpbXVsdGFuZW91c19wb3dlcml0ZXJhdGlvbn0gZnJvbSBcIi4uL2xpbmVhcl9hbGdlYnJhL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgTERBXG4gKiBAZXh0ZW5kcyBEUlxuICovXG5leHBvcnQgY2xhc3MgTERBIGV4dGVuZHMgRFIge1xuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGltZW5zaW9uYWxpdHlfcmVkdWN0aW9uXG4gICAgICogQGFsaWFzIExEQVxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gdGhlIGhpZ2gtZGltZW5zaW9uYWwgZGF0YS5cbiAgICAgKiBAcGFyYW0ge0FycmF5fSBsYWJlbHMgLSB0aGUgbGFiZWwgLyBjbGFzcyBvZiBlYWNoIGRhdGEgcG9pbnQuXG4gICAgICogQHBhcmFtIHtudW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHBhcmFtIHtmdW5jdGlvbn0gW21ldHJpYyA9IGV1Y2xpZGVhbl0gLSB0aGUgbWV0cmljIHdoaWNoIGRlZmluZXMgdGhlIGRpc3RhbmNlIGJldHdlZW4gdHdvIHBvaW50cy4gIFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbc2VlZCA9IDEyMTJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKFgsIGxhYmVscywgZCA9IDIsIG1ldHJpYyA9IGV1Y2xpZGVhbiwgc2VlZD0xMjEyKSB7XG4gICAgICAgIHN1cGVyKFgsIGQsIG1ldHJpYywgc2VlZCk7XG4gICAgICAgIHN1cGVyLnBhcmFtZXRlcl9saXN0ID0gW1wibGFiZWxzXCJdO1xuICAgICAgICB0aGlzLnBhcmFtZXRlcihcImxhYmVsc1wiLCBsYWJlbHMpO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBUcmFuc2Zvcm1zIHRoZSBpbnB1dGRhdGEge0BsaW5rIFh9IHRvIGRpbWVuaW9uYWxpdHkge0BsaW5rIGR9LlxuICAgICAqL1xuICAgIHRyYW5zZm9ybSgpIHtcbiAgICAgICAgbGV0IFggPSB0aGlzLlg7XG4gICAgICAgIGxldCBbIHJvd3MsIGNvbHMgXSA9IFguc2hhcGU7XG4gICAgICAgIGxldCBsYWJlbHMgPSB0aGlzLl9sYWJlbHM7XG4gICAgICAgIGxldCB1bmlxdWVfbGFiZWxzID0ge307XG4gICAgICAgIGxldCBsYWJlbF9pZCA9IDA7XG4gICAgICAgIGxhYmVscy5mb3JFYWNoKChsLCBpKSA9PiB7XG4gICAgICAgICAgICBpZiAobCBpbiB1bmlxdWVfbGFiZWxzKSB7XG4gICAgICAgICAgICAgICAgdW5pcXVlX2xhYmVsc1tsXS5jb3VudCsrO1xuICAgICAgICAgICAgICAgIHVuaXF1ZV9sYWJlbHNbbF0ucm93cy5wdXNoKFgucm93KGkpKTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgdW5pcXVlX2xhYmVsc1tsXSA9IHtcbiAgICAgICAgICAgICAgICAgICAgXCJpZFwiOiBsYWJlbF9pZCsrLFxuICAgICAgICAgICAgICAgICAgICBcImNvdW50XCI6IDEsXG4gICAgICAgICAgICAgICAgICAgIFwicm93c1wiOiBbWC5yb3coaSldXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSlcbiAgICAgICAgXG4gICAgICAgIC8vIGNyZWF0ZSBYX21lYW4gYW5kIHZlY3RvciBtZWFucztcbiAgICAgICAgbGV0IFhfbWVhbiA9IFgubWVhbjtcbiAgICAgICAgbGV0IFZfbWVhbiA9IG5ldyBNYXRyaXgobGFiZWxfaWQsIGNvbHMpXG4gICAgICAgIGZvciAobGV0IGxhYmVsIGluIHVuaXF1ZV9sYWJlbHMpIHtcbiAgICAgICAgICAgIGxldCBWID0gTWF0cml4LmZyb20odW5pcXVlX2xhYmVsc1tsYWJlbF0ucm93cyk7XG4gICAgICAgICAgICBsZXQgdl9tZWFuID0gVi5tZWFuQ29scztcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgY29sczsgKytqKSB7XG4gICAgICAgICAgICAgICAgVl9tZWFuLnNldF9lbnRyeSh1bmlxdWVfbGFiZWxzW2xhYmVsXS5pZCwgaiwgdl9tZWFuW2pdKTtcbiAgICAgICAgICAgIH0gICAgICAgICAgIFxuICAgICAgICB9XG4gICAgICAgIC8vIHNjYXR0ZXJfYmV0d2VlblxuICAgICAgICBsZXQgU19iID0gbmV3IE1hdHJpeChjb2xzLCBjb2xzKTtcbiAgICAgICAgZm9yIChsZXQgbGFiZWwgaW4gdW5pcXVlX2xhYmVscykge1xuICAgICAgICAgICAgbGV0IHYgPSBWX21lYW4ucm93KHVuaXF1ZV9sYWJlbHNbbGFiZWxdLmlkKTtcbiAgICAgICAgICAgIGxldCBtID0gbmV3IE1hdHJpeChjb2xzLCAxLCAoaikgPT4gdltqXSAtIFhfbWVhbik7XG4gICAgICAgICAgICBsZXQgTiA9IHVuaXF1ZV9sYWJlbHNbbGFiZWxdLmNvdW50O1xuICAgICAgICAgICAgU19iID0gU19iLmFkZChtLmRvdChtLnRyYW5zcG9zZSgpKS5tdWx0KE4pKTtcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIHNjYXR0ZXJfd2l0aGluXG4gICAgICAgIGxldCBTX3cgPSBuZXcgTWF0cml4KGNvbHMsIGNvbHMpO1xuICAgICAgICBmb3IgKGxldCBsYWJlbCBpbiB1bmlxdWVfbGFiZWxzKSB7XG4gICAgICAgICAgICBsZXQgdiA9IFZfbWVhbi5yb3codW5pcXVlX2xhYmVsc1tsYWJlbF0uaWQpO1xuICAgICAgICAgICAgbGV0IG0gPSBuZXcgTWF0cml4KGNvbHMsIDEsIChqKSA9PiB2W2pdKVxuICAgICAgICAgICAgbGV0IFIgPSB1bmlxdWVfbGFiZWxzW2xhYmVsXS5yb3dzO1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDAsIG4gPSB1bmlxdWVfbGFiZWxzW2xhYmVsXS5jb3VudDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgICAgIGxldCByb3dfdiA9IG5ldyBNYXRyaXgoY29scywgMSwgKGosXykgPT4gUltpXVtqXSAtIG0uZW50cnkoaiwgMCkpO1xuICAgICAgICAgICAgICAgIFNfdyA9IFNfdy5hZGQocm93X3YuZG90KHJvd192LnRyYW5zcG9zZSgpKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBsZXQgeyBlaWdlbnZlY3RvcnM6IFYgfSA9IHNpbXVsdGFuZW91c19wb3dlcml0ZXJhdGlvbihTX3cuaW52ZXJzZSgpLmRvdChTX2IpLCB0aGlzLl9kKVxuICAgICAgICBWID0gTWF0cml4LmZyb20oVikudHJhbnNwb3NlKClcbiAgICAgICAgdGhpcy5ZID0gWC5kb3QoVilcblxuICAgICAgICAvLyByZXR1cm4gZW1iZWRkaW5nXG4gICAgICAgIHJldHVybiB0aGlzLnByb2plY3Rpb247XG4gICAgfVxufSIsImltcG9ydCB7IE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb259IGZyb20gXCIuLi9saW5lYXJfYWxnZWJyYS9pbmRleC5qc1wiO1xuaW1wb3J0IHsga19uZWFyZXN0X25laWdoYm9ycyB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IG5ldW1haXJfc3VtIH0gZnJvbSBcIi4uL251bWVyaWNhbC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgRFIgfSBmcm9tIFwiLi9EUi5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIExMRVxuICogQGV4dGVuZHMgRFJcbiAqL1xuZXhwb3J0IGNsYXNzIExMRSBleHRlbmRzIERSIHtcbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBMTEVcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG5laWdoYm9ycyAtIHRoZSBsYWJlbCAvIGNsYXNzIG9mIGVhY2ggZGF0YSBwb2ludC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2QgPSAyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzLiAgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzZWVkID0gMTIxMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICovXG4gICAgY29uc3RydWN0b3IoWCwgbmVpZ2hib3JzLCBkPTIsIG1ldHJpYz1ldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICBzdXBlci5wYXJhbWV0ZXJfbGlzdCA9IFtcImtcIl07XG4gICAgICAgIHRoaXMucGFyYW1ldGVyKFwia1wiLCBNYXRoLm1pbihuZWlnaGJvcnMgPz8gTWF0aC5tYXgoTWF0aC5mbG9vcih0aGlzLl9OIC8gMTApLCAyKSwgdGhpcy5fTiAtIDEpKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogVHJhbnNmb3JtcyB0aGUgaW5wdXRkYXRhIHtAbGluayBYfSB0byBkaW1lbmlvbmFsaXR5IHtAbGluayBkfS5cbiAgICAgKi9cbiAgICB0cmFuc2Zvcm0oKSB7XG4gICAgICAgIGNvbnN0IFggPSB0aGlzLlg7XG4gICAgICAgIGNvbnN0IGQgPSB0aGlzLl9kO1xuICAgICAgICBjb25zdCByb3dzID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgY29scyA9IHRoaXMuX0Q7XG4gICAgICAgIGNvbnN0IGsgPSB0aGlzLnBhcmFtZXRlcihcImtcIik7XG4gICAgICAgIGNvbnN0IG5OID0ga19uZWFyZXN0X25laWdoYm9ycyhYLCBrLCBudWxsLCB0aGlzLl9tZXRyaWMpO1xuICAgICAgICBjb25zdCBPID0gbmV3IE1hdHJpeChrLCAxLCAxKTtcbiAgICAgICAgY29uc3QgVyA9IG5ldyBNYXRyaXgocm93cywgcm93cyk7XG5cbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIGNvbnN0IG5OX3JvdyA9IG5OW3Jvd107XG4gICAgICAgICAgICBjb25zdCBaID0gbmV3IE1hdHJpeChrLCBjb2xzLCAoaSwgaikgPT4gWC5lbnRyeShuTl9yb3dbaV0uaiwgaikgLSBYLmVudHJ5KHJvdywgaikpO1xuICAgICAgICAgICAgY29uc3QgQyA9IFouZG90KFouVCk7XG4gICAgICAgICAgICBpZiAoIGsgPiBjb2xzICkge1xuICAgICAgICAgICAgICAgIGNvbnN0IENfdHJhY2UgPSBuZXVtYWlyX3N1bShDLmRpYWcpIC8gMTAwMDtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IGs7ICsraikge1xuICAgICAgICAgICAgICAgICAgICBDLnNldF9lbnRyeShqLCBqLCBDLmVudHJ5KGosIGopICsgQ190cmFjZSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgLy8gcmVjb25zdHJ1Y3Q7XG4gICAgICAgICAgICBsZXQgdyA9IE1hdHJpeC5zb2x2ZV9DRyhDLCBPLCB0aGlzLl9yYW5kb21pemVyKTtcbiAgICAgICAgICAgIHcgPSB3LmRpdmlkZSh3LnN1bSk7XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IGs7ICsraikge1xuICAgICAgICAgICAgICAgIFcuc2V0X2VudHJ5KHJvdywgbk5fcm93W2pdLmosIHcuZW50cnkoaiwgMCkpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIC8vIGNvbXAgZW1iZWRkaW5nXG4gICAgICAgIGNvbnN0IEkgPSBuZXcgTWF0cml4KHJvd3MsIHJvd3MsIFwiaWRlbnRpdHlcIik7XG4gICAgICAgIGNvbnN0IElXID0gSS5zdWIoVyk7XG4gICAgICAgIGNvbnN0IE0gPSBJVy5ULmRvdChJVyk7XG4gICAgICAgIGNvbnN0IHsgZWlnZW52ZWN0b3JzOiBWIH0gPSBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb24oTS5ULmludmVyc2UoKSwgZCArIDEpO1xuICAgICAgICB0aGlzLlkgPSBNYXRyaXguZnJvbShWLnNsaWNlKDEsIDEgKyBkKSkuVDtcblxuICAgICAgICAvLyByZXR1cm4gZW1iZWRkaW5nXG4gICAgICAgIHJldHVybiB0aGlzLnByb2plY3Rpb247XG4gICAgfVxufSIsImltcG9ydCB7IE1hdHJpeCwga19uZWFyZXN0X25laWdoYm9ycyB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb259IGZyb20gXCIuLi9saW5lYXJfYWxnZWJyYS9pbmRleC5qc1wiO1xuaW1wb3J0IHsgRFIgfSBmcm9tIFwiLi9EUi5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIExUU0FcbiAqIEBleHRlbmRzIERSXG4gKi9cbmV4cG9ydCBjbGFzcyBMVFNBIGV4dGVuZHMgRFIge1xuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGltZW5zaW9uYWxpdHlfcmVkdWN0aW9uXG4gICAgICogQGFsaWFzIExUU0FcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG5laWdoYm9ycyAtIHRoZSBsYWJlbCAvIGNsYXNzIG9mIGVhY2ggZGF0YSBwb2ludC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2QgPSAyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzLiAgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzZWVkID0gMTIxMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9lcHVicy5zaWFtLm9yZy9kb2kvYWJzLzEwLjExMzcvUzEwNjQ4Mjc1MDI0MTkxNTR9XG4gICAgICovXG4gICAgY29uc3RydWN0b3IoWCwgbmVpZ2hib3JzLCBkPTIsIG1ldHJpYz1ldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICBzdXBlci5wYXJhbWV0ZXJfbGlzdCA9IFtcImtcIl07XG4gICAgICAgIHRoaXMucGFyYW1ldGVyKFwia1wiLCBNYXRoLm1pbihuZWlnaGJvcnMgPz8gTWF0aC5tYXgoTWF0aC5mbG9vcih0aGlzLl9OIC8gMTApLCAyKSwgdGhpcy5fTiAtIDEpKTtcbiAgICAgICAgaWYgKHRoaXMuX0QgPD0gZCkgdGhyb3cgYERpbWVuc2lvbmFsaXR5IG9mIFggKEQgPSAke3RoaXMuX0R9KSBtdXN0IGJlIGdyZWF0ZXIgdGhhbiB0aGUgcmVxdWlyZWQgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHJlc3VsdCAoZCA9ICR7ZH0pIWA7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFRyYW5zZm9ybXMgdGhlIGlucHV0ZGF0YSB7QGxpbmsgWH0gdG8gZGltZW5pb25hbGl0eSB7QGxpbmsgZH0uXG4gICAgICovXG4gICAgdHJhbnNmb3JtKCkge1xuICAgICAgICBjb25zdCBYID0gdGhpcy5YO1xuICAgICAgICBjb25zdCBkID0gdGhpcy5fZDtcbiAgICAgICAgY29uc3QgWyByb3dzLCBEIF0gPSBYLnNoYXBlO1xuICAgICAgICBjb25zdCBrID0gdGhpcy5wYXJhbWV0ZXIoXCJrXCIpO1xuICAgICAgICAvLyAxLjEgZGV0ZXJtaW5lIGsgbmVhcmVzdCBuZWlnaGJvcnNcbiAgICAgICAgY29uc3Qgbk4gPSBrX25lYXJlc3RfbmVpZ2hib3JzKFgsIGssIG51bGwsIHRoaXMuX21ldHJpYyk7XG4gICAgICAgIC8vIGNlbnRlciBtYXRyaXhcbiAgICAgICAgY29uc3QgTyA9IG5ldyBNYXRyaXgoRCwgRCwgXCJjZW50ZXJcIik7XG4gICAgICAgIGNvbnN0IEIgPSBuZXcgTWF0cml4KHJvd3MsIHJvd3MsIDApO1xuICAgICAgICBcbiAgICAgICAgZm9yIChsZXQgcm93ID0gMDsgcm93IDwgcm93czsgKytyb3cpIHtcbiAgICAgICAgICAgIC8vIDEuMiBjb21wdXRlIHRoZSBkIGxhcmdlc3QgZWlnZW52ZWN0b3JzIG9mIHRoZSBjb3JyZWxhdGlvbiBtYXRyaXhcbiAgICAgICAgICAgIGNvbnN0IElfaSA9IFtyb3csIC4uLm5OW3Jvd10ubWFwKG4gPT4gbi5qKV1cbiAgICAgICAgICAgIGxldCBYX2kgPSBNYXRyaXguZnJvbShJX2kubWFwKG4gPT4gWC5yb3cobikpKTtcbiAgICAgICAgICAgIC8vIGNlbnRlciBYX2lcbiAgICAgICAgICAgIFhfaSA9IFhfaS5kb3QoTylcbiAgICAgICAgICAgIC8vIGNvcnJlbGF0aW9uIG1hdHJpeFxuICAgICAgICAgICAgY29uc3QgQyA9IFhfaS5kb3QoWF9pLnRyYW5zcG9zZSgpKTtcbiAgICAgICAgICAgIGNvbnN0IHsgZWlnZW52ZWN0b3JzOiBnIH0gPSBzaW11bHRhbmVvdXNfcG93ZXJpdGVyYXRpb24oQywgZCk7XG4gICAgICAgICAgICAvL2cucHVzaChsaW5zcGFjZSgwLCBrKS5tYXAoXyA9PiAxIC8gTWF0aC5zcXJ0KGsgKyAxKSkpO1xuICAgICAgICAgICAgY29uc3QgR19pX3QgPSBNYXRyaXguZnJvbShnKTtcbiAgICAgICAgICAgIC8vIDIuIENvbnN0cnVjdGluZyBhbGlnbm1lbnQgbWF0cml4XG4gICAgICAgICAgICBjb25zdCBXX2kgPSBHX2lfdC50cmFuc3Bvc2UoKS5kb3QoR19pX3QpLmFkZCgxIC8gTWF0aC5zcXJ0KGsgKyAxKSk7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGsgKyAxOyArK2kpIHtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IGsgKyAxOyArK2opIHtcbiAgICAgICAgICAgICAgICAgICAgQi5zZXRfZW50cnkoSV9pW2ldLCBJX2lbal0sIEIuZW50cnkoSV9pW2ldLCBJX2lbal0pIC0gKGkgPT09IGogPyAxIDogMCApICsgV19pLmVudHJ5KGksIGopKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICAvLyAzLiBBbGlnbmluZyBnbG9iYWwgY29vcmRpbmF0ZXNcbiAgICAgICAgY29uc3QgeyBlaWdlbnZlY3RvcnM6IFkgfSA9IHNpbXVsdGFuZW91c19wb3dlcml0ZXJhdGlvbihCLCBkICsgMSk7XG4gICAgICAgIHRoaXMuWSA9IE1hdHJpeC5mcm9tKFkuc2xpY2UoMSkpLnRyYW5zcG9zZSgpO1xuXG4gICAgICAgIC8vIHJldHVybiBlbWJlZGRpbmdcbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG59IiwiaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgZXVjbGlkZWFuIH0gZnJvbSBcIi4uL21ldHJpY3MvaW5kZXguanNcIjtcbmltcG9ydCB7IERSIH0gZnJvbSBcIi4vRFIuanNcIjtcblxuLyoqXG4gKiBAY2xhc3NcbiAqIEBhbGlhcyBUU05FXG4gKiBAZXh0ZW5kcyBEUlxuICovXG5leHBvcnQgY2xhc3MgVFNORSBleHRlbmRzIERSIHtcbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBUU05FXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IFggLSB0aGUgaGlnaC1kaW1lbnNpb25hbCBkYXRhLiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW3BlcnBsZXhpdHkgPSA1MF0gLSBwZXJwbGV4aXR5LlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbZXBzaWxvbiA9IDEwXSAtIGxlYXJuaW5nIHBhcmFtZXRlci5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2QgPSAyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzLiAgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzZWVkID0gMTIxMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHJldHVybnMge1RTTkV9XG4gICAgICovXG4gICAgY29uc3RydWN0b3IoWCwgcGVycGxleGl0eT01MCwgZXBzaWxvbj0xMCwgZD0yLCBtZXRyaWM9ZXVjbGlkZWFuLCBzZWVkPTEyMTIpIHtcbiAgICAgICAgc3VwZXIoWCwgZCwgbWV0cmljLCBzZWVkKTtcbiAgICAgICAgc3VwZXIucGFyYW1ldGVyX2xpc3QgPSBbXCJwZXJwbGV4aXR5XCIsIFwiZXBzaWxvblwiXTtcbiAgICAgICAgWyB0aGlzLl9OLCB0aGlzLl9EIF0gPSB0aGlzLlguc2hhcGU7XG4gICAgICAgIHRoaXMucGFyYW1ldGVyKFwicGVycGxleGl0eVwiLCBNYXRoLm1pbihwZXJwbGV4aXR5LCB0aGlzLl9OIC0gMSkpO1xuICAgICAgICB0aGlzLnBhcmFtZXRlcihcImVwc2lsb25cIiwgZXBzaWxvbik7XG4gICAgICAgIHRoaXMuX2l0ZXIgPSAwO1xuICAgICAgICB0aGlzLlkgPSBuZXcgTWF0cml4KHRoaXMuX04sIHRoaXMuX2QsICgpID0+IHRoaXMuX3JhbmRvbWl6ZXIucmFuZG9tKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IGRpc3RhbmNlX21hdHJpeCAtIGFjY2VwdHMgYSBwcmVjb21wdXRlZCBkaXN0YW5jZSBtYXRyaXhcbiAgICAgKiBAcmV0dXJucyB7VFNORX1cbiAgICAgKi9cbiAgICBpbml0KGRpc3RhbmNlX21hdHJpeD1udWxsKSB7XG4gICAgICAgIC8vIGluaXRcbiAgICAgICAgY29uc3QgSHRhcmdldCA9IE1hdGgubG9nKHRoaXMuX3BlcnBsZXhpdHkpO1xuICAgICAgICBjb25zdCBOID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgRCA9IHRoaXMuX0Q7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgbGV0IERlbHRhO1xuICAgICAgICBpZiAoZGlzdGFuY2VfbWF0cml4KSB7XG4gICAgICAgICAgICBEZWx0YSA9IGRpc3RhbmNlX21hdHJpeDtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIERlbHRhID0gbmV3IE1hdHJpeChOLCBOKTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgWF9pID0gWC5yb3coaSk7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IGkgKyAxOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IGRpc3RhbmNlID0gbWV0cmljKFhfaSwgWC5yb3coaikpXG4gICAgICAgICAgICAgICAgICAgIERlbHRhLnNldF9lbnRyeShpLCBqLCBkaXN0YW5jZSk7XG4gICAgICAgICAgICAgICAgICAgIERlbHRhLnNldF9lbnRyeShqLCBpLCBkaXN0YW5jZSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgIH0gXG4gICAgICAgICAgICBcbiAgICAgICAgY29uc3QgUCA9IG5ldyBNYXRyaXgoTiwgTiwgXCJ6ZXJvc1wiKTtcblxuICAgICAgICB0aGlzLl95c3RlcCA9IG5ldyBNYXRyaXgoTiwgRCwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgdGhpcy5fZ2FpbnMgPSBuZXcgTWF0cml4KE4sIEQsIDEpO1xuXG4gICAgICAgIC8vIHNlYXJjaCBmb3IgZml0dGluZyBzaWdtYVxuICAgICAgICBsZXQgcHJvdyA9IG5ldyBBcnJheShOKS5maWxsKDApO1xuICAgICAgICBjb25zdCB0b2wgPSAxZS00O1xuICAgICAgICBjb25zdCBtYXh0cmllcyA9IDUwO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgbGV0IGJldGFtaW4gPSAtSW5maW5pdHk7XG4gICAgICAgICAgICBsZXQgYmV0YW1heCA9IEluZmluaXR5O1xuICAgICAgICAgICAgbGV0IGJldGEgPSAxO1xuICAgICAgICAgICAgbGV0IGRvbmUgPSBmYWxzZTtcblxuICAgICAgICAgICAgbGV0IG51bSA9IDA7XG4gICAgICAgICAgICB3aGlsZSghZG9uZSkge1xuICAgICAgICAgICAgICAgIGxldCBwc3VtID0gMDtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IE47ICsraikge1xuICAgICAgICAgICAgICAgICAgICBsZXQgcGogPSBNYXRoLmV4cCgtRGVsdGEuZW50cnkoaSwgaikgKiBiZXRhKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGkgPT09IGopIHBqID0gMDtcbiAgICAgICAgICAgICAgICAgICAgcHJvd1tqXSA9IHBqO1xuICAgICAgICAgICAgICAgICAgICBwc3VtICs9IHBqO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBsZXQgSGhlcmUgPSAwO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgICAgIGxldCBwaiA9IChwc3VtID09PSAwKSA/IDAgOiBwcm93W2pdIC8gcHN1bTtcbiAgICAgICAgICAgICAgICAgICAgcHJvd1tqXSA9IHBqO1xuICAgICAgICAgICAgICAgICAgICBpZiAocGogPiAxZS03KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBIaGVyZSAtPSBwaiAqIE1hdGgubG9nKHBqKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBpZiAoSGhlcmUgPiBIdGFyZ2V0KSB7XG4gICAgICAgICAgICAgICAgICAgIGJldGFtaW4gPSBiZXRhO1xuICAgICAgICAgICAgICAgICAgICBiZXRhID0gKGJldGFtYXggPT09IEluZmluaXR5KSA/IChiZXRhICogMikgOiAoKGJldGEgKyBiZXRhbWF4KSAvIDIpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGJldGFtYXggPSBiZXRhO1xuICAgICAgICAgICAgICAgICAgICBiZXRhID0gKGJldGFtaW4gPT09IC1JbmZpbml0eSkgPyAoYmV0YSAvIDIpIDogKChiZXRhICsgYmV0YW1pbikgLyAyKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgKytudW07XG4gICAgICAgICAgICAgICAgaWYgKE1hdGguYWJzKEhoZXJlIC0gSHRhcmdldCkgPCB0b2wpIGRvbmUgPSB0cnVlO1xuICAgICAgICAgICAgICAgIGlmIChudW0gPj0gbWF4dHJpZXMpIGRvbmUgPSB0cnVlO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IE47ICsraikge1xuICAgICAgICAgICAgICAgIFAuc2V0X2VudHJ5KGksIGosIHByb3dbal0pO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgLy9jb21wdXRlIHByb2JhYmlsaXRpZXNcbiAgICAgICAgY29uc3QgUG91dCA9IG5ldyBNYXRyaXgoTiwgTiwgXCJ6ZXJvc1wiKVxuICAgICAgICBjb25zdCBOMiA9IE4gKiAyO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgZm9yIChsZXQgaiA9IGk7IGogPCBOOyArK2opIHtcbiAgICAgICAgICAgICAgICBjb25zdCBwID0gTWF0aC5tYXgoKFAuZW50cnkoaSwgaikgKyBQLmVudHJ5KGosIGkpKSAvIE4yLCAxZS0xMDApO1xuICAgICAgICAgICAgICAgIFBvdXQuc2V0X2VudHJ5KGksIGosIHApO1xuICAgICAgICAgICAgICAgIFBvdXQuc2V0X2VudHJ5KGosIGksIHApO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHRoaXMuX1AgPSBQb3V0O1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2l0ZXJhdGlvbnM9NTAwXSAtIG51bWJlciBvZiBpdGVyYXRpb25zLlxuICAgICAqIEB5aWVsZHMge01hdHJpeHxBcnJheTxBcnJheT59IC0gdGhlIHByb2plY3Rpb24uXG4gICAgICovXG4gICAgdHJhbnNmb3JtKGl0ZXJhdGlvbnM9NTAwKSB7XG4gICAgICAgIHRoaXMuY2hlY2tfaW5pdCgpO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGl0ZXJhdGlvbnM7ICsraSkge1xuICAgICAgICAgICAgdGhpcy5uZXh0KCk7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2l0ZXJhdGlvbnM9NTAwXSAtIG51bWJlciBvZiBpdGVyYXRpb25zLlxuICAgICAqIEB5aWVsZHMge01hdHJpeHxBcnJheTxBcnJheT59IC0gdGhlIHByb2plY3Rpb24uXG4gICAgICovXG4gICAgKiBnZW5lcmF0b3IoaXRlcmF0aW9ucz01MDApIHtcbiAgICAgICAgdGhpcy5jaGVja19pbml0KCk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgaXRlcmF0aW9uczsgKytpKSB7XG4gICAgICAgICAgICB0aGlzLm5leHQoKTtcbiAgICAgICAgICAgIHlpZWxkIHRoaXMucHJvamVjdGlvbjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcy5wcm9qZWN0aW9uO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIHBlcmZvcm1zIGEgb3B0aW1pemF0aW9uIHN0ZXBcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9XG4gICAgICovXG4gICAgbmV4dCgpIHtcbiAgICAgICAgY29uc3QgaXRlciA9ICsrdGhpcy5faXRlcjtcbiAgICAgICAgY29uc3QgUCA9IHRoaXMuX1A7XG4gICAgICAgIGNvbnN0IHlzdGVwID0gdGhpcy5feXN0ZXA7XG4gICAgICAgIGNvbnN0IGdhaW5zID0gdGhpcy5fZ2FpbnM7XG4gICAgICAgIGNvbnN0IE4gPSB0aGlzLl9OO1xuICAgICAgICBjb25zdCBlcHNpbG9uID0gdGhpcy5fZXBzaWxvbjtcbiAgICAgICAgY29uc3QgZGltID0gdGhpcy5fZDtcbiAgICAgICAgbGV0IFkgPSB0aGlzLlk7XG5cbiAgICAgICAgLy9jYWxjIGNvc3QgZ3JhZGllbnQ7XG4gICAgICAgIGNvbnN0IHBtdWwgPSBpdGVyIDwgMTAwID8gNCA6IDE7XG4gICAgICAgIFxuICAgICAgICAvLyBjb21wdXRlIFEgZGlzdCAodW5ub3JtYWxpemVkKVxuICAgICAgICBjb25zdCBRdSA9IG5ldyBNYXRyaXgoTiwgTiwgXCJ6ZXJvc1wiKVxuICAgICAgICBsZXQgcXN1bSA9IDA7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gaSArIDE7IGogPCBOOyArK2opIHtcbiAgICAgICAgICAgICAgICBsZXQgZHN1bSA9IDA7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgZCA9IDA7IGQgPCBkaW07ICsrZCkge1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBkaGVyZSA9IFkuZW50cnkoaSwgZCkgLSBZLmVudHJ5KGosIGQpO1xuICAgICAgICAgICAgICAgICAgICBkc3VtICs9IGRoZXJlICogZGhlcmU7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGNvbnN0IHF1ID0gMSAvICgxICsgZHN1bSk7XG4gICAgICAgICAgICAgICAgUXUuc2V0X2VudHJ5KGksIGosIHF1KTtcbiAgICAgICAgICAgICAgICBRdS5zZXRfZW50cnkoaiwgaSwgcXUpO1xuICAgICAgICAgICAgICAgIHFzdW0gKz0gMiAqIHF1O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgLy8gbm9ybWFsaXplIFEgZGlzdFxuICAgICAgICBjb25zdCBRID0gbmV3IE1hdHJpeChOLCBOLCAwKVxuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgZm9yIChsZXQgaiA9IGkgKyAxOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgdmFsID0gTWF0aC5tYXgoUXUuZW50cnkoaSwgaikgLyBxc3VtLCAxZS0xMDApO1xuICAgICAgICAgICAgICAgIFEuc2V0X2VudHJ5KGksIGosIHZhbCk7XG4gICAgICAgICAgICAgICAgUS5zZXRfZW50cnkoaiwgaSwgdmFsKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGNvbnN0IGdyYWQgPSBuZXcgTWF0cml4KE4sIGRpbSwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgcHJlbXVsdCA9IDQgKiAocG11bCAqIFAuZW50cnkoaSwgaikgLSBRLmVudHJ5KGksIGopKSAqIFF1LmVudHJ5KGksIGopO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGQgPSAwOyBkIDwgZGltOyArK2QpIHtcbiAgICAgICAgICAgICAgICAgICAgZ3JhZC5zZXRfZW50cnkoaSwgZCwgZ3JhZC5lbnRyeShpLCBkKSArIHByZW11bHQgKiAoWS5lbnRyeShpLCBkKSAtIFkuZW50cnkoaiwgZCkpKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICAvLyBwZXJmb3JtIGdyYWRpZW50IHN0ZXBcbiAgICAgICAgbGV0IHltZWFuID0gbmV3IEZsb2F0NjRBcnJheShkaW0pO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgZm9yIChsZXQgZCA9IDA7IGQgPCBkaW07ICsrZCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IGdpZCA9IGdyYWQuZW50cnkoaSwgZCk7XG4gICAgICAgICAgICAgICAgY29uc3Qgc2lkID0geXN0ZXAuZW50cnkoaSwgZCk7XG4gICAgICAgICAgICAgICAgY29uc3QgZ2FpbmlkID0gZ2FpbnMuZW50cnkoaSwgZCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgbGV0IG5ld2dhaW4gPSBNYXRoLnNpZ24oZ2lkKSA9PT0gTWF0aC5zaWduKHNpZCkgPyBnYWluaWQgKiAuOCA6IGdhaW5pZCArIC4yO1xuICAgICAgICAgICAgICAgIGlmIChuZXdnYWluIDwgLjAxKSBuZXdnYWluID0gLjAxO1xuICAgICAgICAgICAgICAgIGdhaW5zLnNldF9lbnRyeShpLCBkLCBuZXdnYWluKTtcblxuICAgICAgICAgICAgICAgIGNvbnN0IG1vbXZhbCA9IGl0ZXIgPCAyNTAgPyAuNSA6IC44O1xuICAgICAgICAgICAgICAgIGNvbnN0IG5ld3NpZCA9IG1vbXZhbCAqIHNpZCAtIGVwc2lsb24gKiBuZXdnYWluICogZ2lkO1xuICAgICAgICAgICAgICAgIHlzdGVwLnNldF9lbnRyeShpLCBkLCBuZXdzaWQpO1xuXG4gICAgICAgICAgICAgICAgWS5zZXRfZW50cnkoaSwgZCwgWS5lbnRyeShpLCBkKSArIG5ld3NpZCk7XG4gICAgICAgICAgICAgICAgeW1lYW5bZF0gKz0gWS5lbnRyeShpLCBkKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBkID0gMDsgZCA8IDI7ICsrZCkge1xuICAgICAgICAgICAgICAgIFkuc2V0X2VudHJ5KGksIGQsIFkuZW50cnkoaSwgZCkgLSB5bWVhbltkXSAvIE4pXG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4gdGhpcy5ZO1xuICAgIH1cbn0gIiwiLyoqXG4gKlxuICogQG1lbWJlcm9mIG1vZHVsZTpvcHRpbWl6YXRpb25cbiAqIEBhbGlhcyBwb3dlbGxcbiAqIEBwYXJhbSB7RnVuY3Rpb259IGZcbiAqIEBwYXJhbSB7QXJyYXl9IHgwXG4gKiBAcGFyYW0ge051bWJlcn0gW21heF9pdGVyID0gMzAwXVxuICogQHJldHVybnMge0FycmF5fVxuICogQHNlZSBodHRwOi8vb3B0aW1pemF0aW9uLWpzLmdpdGh1Yi5pby9vcHRpbWl6YXRpb24tanMvb3B0aW1pemF0aW9uLmpzLmh0bWwjbGluZTQzOFxuICovXG5leHBvcnQgZGVmYXVsdCBmdW5jdGlvbiAoZiwgeDAsIG1heF9pdGVyID0gMzAwKSB7XG4gICAgY29uc3QgZXBzaWxvbiA9IDFlLTI7XG4gICAgY29uc3QgbiA9IHgwLmxlbmd0aDtcbiAgICBsZXQgYWxwaGEgPSAxZS0zO1xuICAgIGxldCBwZnggPSAxMDAwMDtcbiAgICBsZXQgeCA9IHgwLnNsaWNlKCk7XG4gICAgbGV0IGZ4ID0gZih4KTtcbiAgICBsZXQgY29udmVyZ2VuY2UgPSBmYWxzZTtcblxuICAgIHdoaWxlIChtYXhfaXRlci0tID49IDAgJiYgIWNvbnZlcmdlbmNlKSB7XG4gICAgICAgIGNvbnZlcmdlbmNlID0gdHJ1ZTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgICAgIHhbaV0gKz0gMWUtNjtcbiAgICAgICAgICAgIGxldCBmeGkgPSBmKHgpO1xuICAgICAgICAgICAgeFtpXSAtPSAxZS02O1xuICAgICAgICAgICAgbGV0IGR4ID0gKGZ4aSAtIGZ4KSAvIDFlLTY7XG4gICAgICAgICAgICBpZiAoTWF0aC5hYnMoZHgpID4gZXBzaWxvbikge1xuICAgICAgICAgICAgICAgIGNvbnZlcmdlbmNlID0gZmFsc2U7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB4W2ldIC09IGFscGhhICogZHg7XG4gICAgICAgICAgICBmeCA9IGYoeCk7XG4gICAgICAgIH1cbiAgICAgICAgYWxwaGEgKj0gcGZ4ID49IGZ4ID8gMS4wNSA6IDAuNDtcbiAgICAgICAgcGZ4ID0gZng7XG4gICAgfVxuICAgIHJldHVybiB4O1xufVxuIiwiaW1wb3J0IHsgTWF0cml4IH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgZXVjbGlkZWFuLCBldWNsaWRlYW5fc3F1YXJlZCB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBCYWxsVHJlZSB9IGZyb20gXCIuLi9rbm4vaW5kZXguanNcIjtcbmltcG9ydCB7IG5ldW1haXJfc3VtIH0gZnJvbSBcIi4uL251bWVyaWNhbC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgbGluc3BhY2UgfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBwb3dlbGwgfSBmcm9tIFwiLi4vb3B0aW1pemF0aW9uL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5pbXBvcnQgeyBtYXggfSBmcm9tIFwiLi4vdXRpbC9pbmRleC5qc1wiO1xuaW1wb3J0IHsgS05OIH0gZnJvbSBcIi4uL2tubi9pbmRleC5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIFVNQVBcbiAqIEBleHRlbmRzIERSXG4gKi9cbmV4cG9ydCBjbGFzcyBVTUFQIGV4dGVuZHMgRFIge1xuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTpkaW1lbnNpb25hbGl0eV9yZWR1Y3Rpb25cbiAgICAgKiBAYWxpYXMgVU1BUFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gdGhlIGhpZ2gtZGltZW5zaW9uYWwgZGF0YS4gXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtuX25laWdoYm9ycyA9IDE1XSAtIHNpemUgb2YgdGhlIGxvY2FsIG5laWdoYm9yaG9vZC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2xvY2FsX2Nvbm5lY3Rpdml0eSA9IDFdIC0gbnVtYmVyIG9mIG5lYXJlc3QgbmVpZ2hib3JzIGNvbm5lY3RlZCBpbiB0aGUgbG9jYWwgbmVpZ2hib3Job29kLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbbWluX2Rpc3QgPSAxXSAtIGNvbnRyb2xzIGhvdyB0aWdodGx5IHBvaW50cyBnZXQgcGFja2VkIHRvZ2V0aGVyLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbZCA9IDJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEBwYXJhbSB7RnVuY3Rpb259IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gdGhlIG1ldHJpYyB3aGljaCBkZWZpbmVzIHRoZSBkaXN0YW5jZSBiZXR3ZWVuIHR3byBwb2ludHMgaW4gdGhlIGhpZ2gtZGltZW5zaW9uYWwgc3BhY2UuICBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW3NlZWQgPSAxMjEyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcmV0dXJucyB7VU1BUH1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBuX25laWdoYm9ycz0xNSwgbG9jYWxfY29ubmVjdGl2aXR5PTEsIG1pbl9kaXN0PTEsIGQ9MiwgbWV0cmljPWV1Y2xpZGVhbiwgc2VlZD0xMjEyKSB7XG4gICAgICAgIHN1cGVyKFgsIGQsIG1ldHJpYywgc2VlZClcbiAgICAgICAgc3VwZXIucGFyYW1ldGVyX2xpc3QgPSBbXCJuX25laWdoYm9yc1wiLCBcImxvY2FsX2Nvbm5lY3Rpdml0eVwiLCBcIm1pbl9kaXN0XCJdO1xuICAgICAgICBbIHRoaXMuX04sIHRoaXMuX0QgXSA9IHRoaXMuWC5zaGFwZTtcbiAgICAgICAgbl9uZWlnaGJvcnMgPSBNYXRoLm1pbih0aGlzLl9OIC0gMSwgbl9uZWlnaGJvcnMpO1xuICAgICAgICB0aGlzLnBhcmFtZXRlcihcIm5fbmVpZ2hib3JzXCIsIG5fbmVpZ2hib3JzKTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJsb2NhbF9jb25uZWN0aXZpdHlcIiwgTWF0aC5taW4obG9jYWxfY29ubmVjdGl2aXR5LCBuX25laWdoYm9ycyAtIDEpKTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJtaW5fZGlzdFwiLCBtaW5fZGlzdCk7XG4gICAgICAgIHRoaXMuX2l0ZXIgPSAwO1xuICAgICAgICB0aGlzLl9zcHJlYWQgPSAxO1xuICAgICAgICB0aGlzLl9zZXRfb3BfbWl4X3JhdGlvID0gMTtcbiAgICAgICAgdGhpcy5fcmVwdWxzaW9uX3N0cmVuZ3RoID0gMTtcbiAgICAgICAgdGhpcy5fbmVnYXRpdmVfc2FtcGxlX3JhdGUgPSA1O1xuICAgICAgICB0aGlzLl9uX2Vwb2NocyA9IDM1MDtcbiAgICAgICAgdGhpcy5faW5pdGlhbF9hbHBoYSA9IDE7XG4gICAgICAgIHRoaXMuWSA9IG5ldyBNYXRyaXgodGhpcy5fTiwgdGhpcy5fZCwgKCkgPT4gdGhpcy5fcmFuZG9taXplci5yYW5kb20pO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBzcHJlYWQgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG1pbl9kaXN0IFxuICAgICAqIEByZXR1cm5zIHtBcnJheX1cbiAgICAgKi9cbiAgICBfZmluZF9hYl9wYXJhbXMoc3ByZWFkLCBtaW5fZGlzdCkge1xuICAgICAgICBjb25zdCBjdXJ2ZSA9ICh4LCBhLCBiKSA9PiAxIC8gKDEgKyBhICogTWF0aC5wb3coeCwgMiAqIGIpKTtcbiAgICAgICAgY29uc3QgeHYgPSBsaW5zcGFjZSgwLCBzcHJlYWQgKiAzLCAzMDApO1xuICAgICAgICBjb25zdCB5diA9IGxpbnNwYWNlKDAsIHNwcmVhZCAqIDMsIDMwMCk7XG4gICAgICAgIFxuICAgICAgICBmb3IgKGxldCBpID0gMCwgbiA9IHh2Lmxlbmd0aDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgY29uc3QgeHZfaSA9IHh2W2ldO1xuICAgICAgICAgICAgeXZbaV0gPSAoeHZfaSA8IG1pbl9kaXN0ID8gMSA6IE1hdGguZXhwKC0oeHZfaSAtIG1pbl9kaXN0KSAvIHNwcmVhZCkpO1xuICAgICAgICB9XG4gICAgICBcbiAgICAgICAgY29uc3QgZXJyID0gKHApID0+IHtcbiAgICAgICAgICAgIGNvbnN0IGVycm9yID0gbGluc3BhY2UoMSwgMzAwKS5tYXAoKF8sIGkpID0+IHl2W2ldIC0gY3VydmUoeHZbaV0sIHBbMF0sIHBbMV0pKTtcbiAgICAgICAgICAgIHJldHVybiBNYXRoLnNxcnQobmV1bWFpcl9zdW0oZXJyb3IubWFwKGUgPT4gZSAqIGUpKSk7XG4gICAgICAgIH1cbiAgICAgIFxuICAgICAgICByZXR1cm4gcG93ZWxsKGVyciwgWzEsIDFdKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7QXJyYXk8QXJyYXk+fSBkaXN0YW5jZXMgXG4gICAgICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSBzaWdtYXMgXG4gICAgICogQHBhcmFtIHtBcnJheTxOdW1iZXI+fSByaG9zIFxuICAgICAqIEByZXR1cm5zIHtBcnJheX1cbiAgICAgKi9cbiAgICBfY29tcHV0ZV9tZW1iZXJzaGlwX3N0cmVuZ3RocyhkaXN0YW5jZXMsIHNpZ21hcywgcmhvcykge1xuICAgICAgICBmb3IgKGxldCBpID0gMCwgbiA9IGRpc3RhbmNlcy5sZW5ndGg7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwLCBtID0gZGlzdGFuY2VzW2ldLmxlbmd0aDsgaiA8IG07ICsraikge1xuICAgICAgICAgICAgICAgIGNvbnN0IHYgPSBkaXN0YW5jZXNbaV1bal0udmFsdWUgLSByaG9zW2ldO1xuICAgICAgICAgICAgICAgIGRpc3RhbmNlc1tpXVtqXS52YWx1ZSA9IHYgPiAwID8gTWF0aC5leHAoLXYgLyBzaWdtYXNbaV0pIDogMTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gZGlzdGFuY2VzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtLTk58QmFsbFRyZWV9IGtubiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gayBcbiAgICAgKiBAcmV0dXJucyB7T2JqZWN0fVxuICAgICAqL1xuICAgIF9zbW9vdGhfa25uX2Rpc3Qoa25uLCBrKSB7XG4gICAgICAgIGNvbnN0IFNNT09USF9LX1RPTEVSQU5DRSA9IDFlLTU7XG4gICAgICAgIGNvbnN0IE1JTl9LX0RJU1RfU0NBTEUgPSAxZS0zO1xuICAgICAgICBjb25zdCBuX2l0ZXIgPSA2NDtcbiAgICAgICAgY29uc3QgbG9jYWxfY29ubmVjdGl2aXR5ID0gdGhpcy5fbG9jYWxfY29ubmVjdGl2aXR5O1xuICAgICAgICBjb25zdCB0YXJnZXQgPSBNYXRoLmxvZzIoayk7XG4gICAgICAgIGNvbnN0IHJob3MgPSBbXTtcbiAgICAgICAgY29uc3Qgc2lnbWFzID0gW107XG4gICAgICAgIGNvbnN0IFggPSB0aGlzLlg7XG4gICAgICAgIGNvbnN0IE4gPSBYLnNoYXBlWzBdO1xuICAgICAgICAvL2NvbnN0IGRpc3RhbmNlcyA9IFsuLi5YXS5tYXAoeF9pID0+IGtubi5zZWFyY2goeF9pLCBrKS5yYXdfZGF0YSgpLnJldmVyc2UoKSk7XG5cbiAgICAgICAgY29uc3QgZGlzdGFuY2VzID0gW107XG4gICAgICAgIGlmICh0aGlzLl9tZXRyaWMgPT09IFwicHJlY29tcHV0ZWRcIikge1xuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgICAgICBkaXN0YW5jZXMucHVzaChrbm4uc2VhcmNoKGksIGspLnJldmVyc2UoKSlcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgZm9yIChjb25zdCB4X2kgb2YgWCkge1xuICAgICAgICAgICAgICAgIGRpc3RhbmNlcy5wdXNoKGtubi5zZWFyY2goeF9pLCBrKS5yYXdfZGF0YSgpLnJldmVyc2UoKSlcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBsZXQgbG8gPSAwO1xuICAgICAgICAgICAgbGV0IGhpID0gSW5maW5pdHk7XG4gICAgICAgICAgICBsZXQgbWlkID0gMTtcblxuICAgICAgICAgICAgY29uc3Qgc2VhcmNoX3Jlc3VsdCA9IGRpc3RhbmNlc1tpXVxuICAgICAgICAgICAgY29uc3Qgbm9uX3plcm9fZGlzdCA9IHNlYXJjaF9yZXN1bHQuZmlsdGVyKGQgPT4gZC52YWx1ZSA+IDApO1xuICAgICAgICAgICAgY29uc3Qgbm9uX3plcm9fZGlzdF9sZW5ndGggPSBub25femVyb19kaXN0Lmxlbmd0aDtcbiAgICAgICAgICAgIGlmIChub25femVyb19kaXN0X2xlbmd0aCA+PSBsb2NhbF9jb25uZWN0aXZpdHkpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBpbmRleCA9IE1hdGguZmxvb3IobG9jYWxfY29ubmVjdGl2aXR5KTtcbiAgICAgICAgICAgICAgICBjb25zdCBpbnRlcnBvbGF0aW9uID0gbG9jYWxfY29ubmVjdGl2aXR5IC0gaW5kZXg7XG4gICAgICAgICAgICAgICAgaWYgKGluZGV4ID4gMCkge1xuICAgICAgICAgICAgICAgICAgICByaG9zLnB1c2gobm9uX3plcm9fZGlzdFtpbmRleCAtIDFdKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGludGVycG9sYXRpb24gPiBTTU9PVEhfS19UT0xFUkFOQ0UpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJob3NbaV0udmFsdWUgKz0gaW50ZXJwb2xhdGlvbiAqIChub25femVyb19kaXN0W2luZGV4XS52YWx1ZSAtIG5vbl96ZXJvX2Rpc3RbaW5kZXggLSAxXSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICByaG9zW2ldLnZhbHVlID0gaW50ZXJwb2xhdGlvbiAqIG5vbl96ZXJvX2Rpc3RbMF0udmFsdWU7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfSBlbHNlIGlmIChub25femVyb19kaXN0X2xlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICByaG9zW2ldID0gbm9uX3plcm9fZGlzdFtub25femVyb19kaXN0X2xlbmd0aCAtIDFdLnZhbHVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgeCA9IDA7IHggPCBuX2l0ZXI7ICsreCkge1xuICAgICAgICAgICAgICAgIGxldCBwc3VtID0gMDtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IGs7ICsraikge1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBkID0gc2VhcmNoX3Jlc3VsdFtqXS52YWx1ZSAtIHJob3NbaV07XG4gICAgICAgICAgICAgICAgICAgIHBzdW0gKz0gKGQgPiAwID8gTWF0aC5leHAoLShkIC8gbWlkKSkgOiAxKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgaWYgKE1hdGguYWJzKHBzdW0gLSB0YXJnZXQpIDwgU01PT1RIX0tfVE9MRVJBTkNFKSB7XG4gICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBpZiAocHN1bSA+IHRhcmdldCkge1xuICAgICAgICAgICAgICAgICAgICBbaGksIG1pZF0gPSBbbWlkLCAobG8gKyBoaSkgLyAyXTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBpZiAoaGkgPT09IEluZmluaXR5KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBbbG8sIG1pZF0gPSBbbWlkLCBtaWQgKiAyXTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIFtsbywgbWlkXSA9IFttaWQsIChsbyArIGhpKSAvIDJdO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgc2lnbWFzW2ldID0gbWlkO1xuXG4gICAgICAgICAgICBjb25zdCBtZWFuX2l0aGQgPSBzZWFyY2hfcmVzdWx0LnJlZHVjZSgoYSwgYikgPT4gYSArIGIudmFsdWUsIDApIC8gc2VhcmNoX3Jlc3VsdC5sZW5ndGg7XG4gICAgICAgICAgICAvL2xldCBtZWFuX2QgPSBudWxsO1xuICAgICAgICAgICAgaWYgKHJob3NbaV0gPiAwKSB7XG4gICAgICAgICAgICAgICAgaWYgKHNpZ21hc1tpXSA8IE1JTl9LX0RJU1RfU0NBTEUgKiBtZWFuX2l0aGQpIHtcbiAgICAgICAgICAgICAgICAgICAgc2lnbWFzW2ldID0gTUlOX0tfRElTVF9TQ0FMRSAqIG1lYW5faXRoZDtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIGNvbnN0IG1lYW5fZCA9IGRpc3RhbmNlcy5yZWR1Y2UoKGFjYywgcmVzKSA9PiBhY2MgKyByZXMucmVkdWNlKChhLCBiKSA9PiBhICsgYi52YWx1ZSwgMCkgLyByZXMubGVuZ3RoKTtcbiAgICAgICAgICAgICAgICBpZiAoc2lnbWFzW2ldID4gTUlOX0tfRElTVF9TQ0FMRSAqIG1lYW5fZCkge1xuICAgICAgICAgICAgICAgICAgICBzaWdtYXNbaV0gPSBNSU5fS19ESVNUX1NDQUxFICogbWVhbl9kO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgXCJkaXN0YW5jZXNcIjogZGlzdGFuY2VzLCBcbiAgICAgICAgICAgIFwic2lnbWFzXCI6IHNpZ21hcywgXG4gICAgICAgICAgICBcInJob3NcIjogcmhvc1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gbl9uZWlnaGJvcnMgXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBfZnV6enlfc2ltcGxpY2lhbF9zZXQoWCwgbl9uZWlnaGJvcnMpIHtcbiAgICAgICAgY29uc3QgTiA9IFguc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3Qga25uID0gbWV0cmljID09PSBcInByZWNvbXB1dGVkXCIgPyBuZXcgS05OKFgsIFwicHJlY29tcHV0ZWRcIikgOiBuZXcgQmFsbFRyZWUoWC50bzJkQXJyYXksIG1ldHJpYyk7XG4gICAgICAgIGxldCB7IGRpc3RhbmNlcywgc2lnbWFzLCByaG9zIH0gPSB0aGlzLl9zbW9vdGhfa25uX2Rpc3Qoa25uLCBuX25laWdoYm9ycyk7XG4gICAgICAgIGRpc3RhbmNlcyA9IHRoaXMuX2NvbXB1dGVfbWVtYmVyc2hpcF9zdHJlbmd0aHMoZGlzdGFuY2VzLCBzaWdtYXMsIHJob3MpO1xuICAgICAgICBjb25zdCByZXN1bHQgPSBuZXcgTWF0cml4KE4sIE4sIFwiemVyb3NcIik7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBkaXN0YW5jZXNfaSA9IGRpc3RhbmNlc1tpXTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgZGlzdGFuY2VzX2kubGVuZ3RoOyArK2opIHtcbiAgICAgICAgICAgICAgICByZXN1bHQuc2V0X2VudHJ5KGksIGRpc3RhbmNlc19pW2pdLmVsZW1lbnQuaW5kZXgsIGRpc3RhbmNlc19pW2pdLnZhbHVlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjb25zdCB0cmFuc3Bvc2VkX3Jlc3VsdCA9IHJlc3VsdC5UO1xuICAgICAgICBjb25zdCBwcm9kX21hdHJpeCA9IHJlc3VsdC5tdWx0KHRyYW5zcG9zZWRfcmVzdWx0KTtcbiAgICAgICAgcmV0dXJuIHJlc3VsdFxuICAgICAgICAgICAgLmFkZCh0cmFuc3Bvc2VkX3Jlc3VsdClcbiAgICAgICAgICAgIC5zdWIocHJvZF9tYXRyaXgpXG4gICAgICAgICAgICAubXVsdCh0aGlzLl9zZXRfb3BfbWl4X3JhdGlvKVxuICAgICAgICAgICAgLmFkZChwcm9kX21hdHJpeC5tdWx0KDEgLSB0aGlzLl9zZXRfb3BfbWl4X3JhdGlvKSk7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gbl9lcG9jaHMgXG4gICAgICogQHJldHVybnMge0FycmF5fVxuICAgICAqL1xuICAgIF9tYWtlX2Vwb2Noc19wZXJfc2FtcGxlKG5fZXBvY2hzKSB7XG4gICAgICAgIGNvbnN0IHdlaWdodHMgPSB0aGlzLl93ZWlnaHRzO1xuICAgICAgICBjb25zdCByZXN1bHQgPSBuZXcgRmxvYXQzMkFycmF5KHdlaWdodHMubGVuZ3RoKS5maWxsKC0xKTtcbiAgICAgICAgY29uc3Qgd2VpZ2h0c19tYXggPSBtYXgod2VpZ2h0cyk7XG4gICAgICAgIGNvbnN0IG5fc2FtcGxlcyA9IHdlaWdodHMubWFwKHcgPT4gbl9lcG9jaHMgKiAodyAvIHdlaWdodHNfbWF4KSk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgcmVzdWx0Lmxlbmd0aDsgKytpKSBcbiAgICAgICAgICBpZiAobl9zYW1wbGVzW2ldID4gMCkgcmVzdWx0W2ldID0gTWF0aC5yb3VuZChuX2Vwb2NocyAvIG5fc2FtcGxlc1tpXSk7XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gZ3JhcGggXG4gICAgICogQHJldHVybnMge09iamVjdH1cbiAgICAgKi9cbiAgICBfdG9jb28oZ3JhcGgpIHtcbiAgICAgICAgY29uc3Qgcm93cyA9IFtdO1xuICAgICAgICBjb25zdCBjb2xzID0gW107XG4gICAgICAgIGNvbnN0IGRhdGEgPSBbXTtcbiAgICAgICAgY29uc3QgWyByb3dzX24sIGNvbHNfbiBdID0gZ3JhcGguc2hhcGU7XG4gICAgICAgIGZvciAobGV0IHJvdyA9IDA7IHJvdyA8IHJvd3NfbjsgKytyb3cpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGNvbCA9IDA7IGNvbCA8IGNvbHNfbjsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBlbnRyeSA9IGdyYXBoLmVudHJ5KHJvdywgY29sKTtcbiAgICAgICAgICAgICAgICBpZiAoZW50cnkgIT09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgcm93cy5wdXNoKHJvdyk7XG4gICAgICAgICAgICAgICAgICAgIGNvbHMucHVzaChjb2wpO1xuICAgICAgICAgICAgICAgICAgICBkYXRhLnB1c2goZW50cnkpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgXCJyb3dzXCI6IHJvd3MsIFxuICAgICAgICAgICAgXCJjb2xzXCI6IGNvbHMsIFxuICAgICAgICAgICAgXCJkYXRhXCI6IGRhdGFcbiAgICAgICAgfTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyBhbGwgbmVjZXNzYXJ5IFxuICAgICAqIEByZXR1cm5zIHtVTUFQfVxuICAgICAqL1xuICAgIGluaXQoKSB7XG4gICAgICAgIGNvbnN0IFsgYSwgYiBdID0gdGhpcy5fZmluZF9hYl9wYXJhbXModGhpcy5fc3ByZWFkLCB0aGlzLl9taW5fZGlzdCk7XG4gICAgICAgIHRoaXMuX2EgPSBhO1xuICAgICAgICB0aGlzLl9iID0gYjtcbiAgICAgICAgdGhpcy5fZ3JhcGggPSB0aGlzLl9mdXp6eV9zaW1wbGljaWFsX3NldCh0aGlzLlgsIHRoaXMuX25fbmVpZ2hib3JzKTtcbiAgICAgICAgY29uc3QgeyByb3dzLCBjb2xzLCBkYXRhOiB3ZWlnaHRzIH0gPSB0aGlzLl90b2Nvbyh0aGlzLl9ncmFwaCk7XG4gICAgICAgIHRoaXMuX2hlYWQgPSByb3dzO1xuICAgICAgICB0aGlzLl90YWlsID0gY29scztcbiAgICAgICAgdGhpcy5fd2VpZ2h0cyA9IHdlaWdodHM7XG4gICAgICAgIHRoaXMuX2Vwb2Noc19wZXJfc2FtcGxlID0gdGhpcy5fbWFrZV9lcG9jaHNfcGVyX3NhbXBsZSh0aGlzLl9uX2Vwb2Nocyk7XG4gICAgICAgIHRoaXMuX2Vwb2Noc19wZXJfbmVnYXRpdmVfc2FtcGxlID0gdGhpcy5fZXBvY2hzX3Blcl9zYW1wbGUubWFwKGQgPT4gZCAqIHRoaXMuX25lZ2F0aXZlX3NhbXBsZV9yYXRlKTtcbiAgICAgICAgdGhpcy5fZXBvY2hfb2ZfbmV4dF9zYW1wbGUgPSB0aGlzLl9lcG9jaHNfcGVyX3NhbXBsZS5zbGljZSgpO1xuICAgICAgICB0aGlzLl9lcG9jaF9vZl9uZXh0X25lZ2F0aXZlX3NhbXBsZSA9IHRoaXMuX2Vwb2Noc19wZXJfbmVnYXRpdmVfc2FtcGxlLnNsaWNlKCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIHNldCBsb2NhbF9jb25uZWN0aXZpdHkodmFsdWUpIHtcbiAgICAgICAgdGhpcy5fbG9jYWxfY29ubmVjdGl2aXR5ID0gdmFsdWU7XG4gICAgfVxuXG4gICAgZ2V0IGxvY2FsX2Nvbm5lY3Rpdml0eSgpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX2xvY2FsX2Nvbm5lY3Rpdml0eTtcbiAgICB9XG5cbiAgICBzZXQgbWluX2Rpc3QodmFsdWUpIHtcbiAgICAgICAgdGhpcy5fbWluX2Rpc3QgPSB2YWx1ZTtcbiAgICB9XG5cbiAgICBnZXQgbWluX2Rpc3QoKSB7XG4gICAgICAgIHJldHVybiB0aGlzLl9taW5fZGlzdDtcbiAgICB9XG5cbiAgICBncmFwaCgpIHtcbiAgICAgICAgdGhpcy5jaGVja19pbml0KCk7XG4gICAgICAgIHJldHVybiB7IGNvbHM6IHRoaXMuX2hlYWQsIHJvd3M6IHRoaXMuX3RhaWwsIHdlaWdodHM6IHRoaXMuX3dlaWdodHMgfTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2l0ZXJhdGlvbnM9MzUwXSAtIG51bWJlciBvZiBpdGVyYXRpb25zLlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh8QXJyYXl9XG4gICAgICovXG4gICAgdHJhbnNmb3JtKGl0ZXJhdGlvbnM9MzUwKSB7XG4gICAgICAgIGlmICh0aGlzLl9uX2Vwb2NocyAhPSBpdGVyYXRpb25zKSB7XG4gICAgICAgICAgICB0aGlzLl9uX2Vwb2NocyA9IGl0ZXJhdGlvbnM7XG4gICAgICAgICAgICB0aGlzLmluaXQoKTtcbiAgICAgICAgfVxuICAgICAgICB0aGlzLmNoZWNrX2luaXQoKTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBpdGVyYXRpb25zOyArK2kpIHtcbiAgICAgICAgICAgIHRoaXMubmV4dCgpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzLnByb2plY3Rpb247XG4gICAgfVxuXG5cbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2l0ZXJhdGlvbnM9MzUwXSAtIG51bWJlciBvZiBpdGVyYXRpb25zLlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh8QXJyYXl9XG4gICAgICovXG4gICAgKiBnZW5lcmF0b3IoaXRlcmF0aW9ucz0zNTApIHtcbiAgICAgICAgaWYgKHRoaXMuX25fZXBvY2hzICE9IGl0ZXJhdGlvbnMpIHtcbiAgICAgICAgICAgIHRoaXMuX25fZXBvY2hzID0gaXRlcmF0aW9ucztcbiAgICAgICAgICAgIHRoaXMuaW5pdCgpO1xuICAgICAgICB9XG4gICAgICAgIHRoaXMuY2hlY2tfaW5pdCgpO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGl0ZXJhdGlvbnM7ICsraSkge1xuICAgICAgICAgICAgdGhpcy5uZXh0KCk7XG4gICAgICAgICAgICB5aWVsZCB0aGlzLnByb2plY3Rpb247XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSB4IFxuICAgICAqIEByZXR1cm5zIHtOdW1iZXJ9XG4gICAgICovXG4gICAgX2NsaXAoeCkge1xuICAgICAgICBpZiAoeCA+IDQpIHJldHVybiA0O1xuICAgICAgICBpZiAoeCA8IC00KSByZXR1cm4gLTQ7XG4gICAgICAgIHJldHVybiB4O1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIHBlcmZvcm1zIHRoZSBvcHRpbWl6YXRpb24gc3RlcC5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBoZWFkX2VtYmVkZGluZyBcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gdGFpbF9lbWJlZGRpbmcgXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IGhlYWQgXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IHRhaWwgXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBfb3B0aW1pemVfbGF5b3V0KGhlYWRfZW1iZWRkaW5nLCB0YWlsX2VtYmVkZGluZywgaGVhZCwgdGFpbCkge1xuICAgICAgICBjb25zdCB7IFxuICAgICAgICAgICAgX2Q6IGRpbSwgXG4gICAgICAgICAgICBfYWxwaGE6IGFscGhhLCBcbiAgICAgICAgICAgIF9yZXB1bHNpb25fc3RyZW5ndGg6IHJlcHVsc2lvbl9zdHJlbmd0aCwgXG4gICAgICAgICAgICBfYTogYSwgXG4gICAgICAgICAgICBfYjogYixcbiAgICAgICAgICAgIF9lcG9jaHNfcGVyX3NhbXBsZTogZXBvY2hzX3Blcl9zYW1wbGUsXG4gICAgICAgICAgICBfZXBvY2hzX3Blcl9uZWdhdGl2ZV9zYW1wbGU6IGVwb2Noc19wZXJfbmVnYXRpdmVfc2FtcGxlLFxuICAgICAgICAgICAgX2Vwb2NoX29mX25leHRfbmVnYXRpdmVfc2FtcGxlOiBlcG9jaF9vZl9uZXh0X25lZ2F0aXZlX3NhbXBsZSxcbiAgICAgICAgICAgIF9lcG9jaF9vZl9uZXh0X3NhbXBsZTogZXBvY2hfb2ZfbmV4dF9zYW1wbGUsXG4gICAgICAgICAgICBfY2xpcDogY2xpcFxuICAgICAgICB9ID0gdGhpcztcbiAgICAgICAgY29uc3QgdGFpbF9sZW5ndGggPSB0YWlsLmxlbmd0aDtcblxuICAgICAgICBmb3IgKGxldCBpID0gMCwgbiA9IGVwb2Noc19wZXJfc2FtcGxlLmxlbmd0aDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgaWYgKGVwb2NoX29mX25leHRfc2FtcGxlW2ldIDw9IHRoaXMuX2l0ZXIpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBqID0gaGVhZFtpXTtcbiAgICAgICAgICAgICAgICBjb25zdCBrID0gdGFpbFtpXTtcbiAgICAgICAgICAgICAgICBjb25zdCBjdXJyZW50ID0gaGVhZF9lbWJlZGRpbmcucm93KGopO1xuICAgICAgICAgICAgICAgIGNvbnN0IG90aGVyID0gdGFpbF9lbWJlZGRpbmcucm93KGspO1xuICAgICAgICAgICAgICAgIGNvbnN0IGRpc3QgPSBldWNsaWRlYW5fc3F1YXJlZChjdXJyZW50LCBvdGhlcik7XG4gICAgICAgICAgICAgICAgbGV0IGdyYWRfY29lZmYgPSAwO1xuICAgICAgICAgICAgICAgIGlmIChkaXN0ID4gMCkge1xuICAgICAgICAgICAgICAgICAgICBncmFkX2NvZWZmID0gKC0yICogYSAqIGIgKiBNYXRoLnBvdyhkaXN0LCBiIC0gMSkpIC8gKGEgKiBNYXRoLnBvdyhkaXN0LCBiKSArIDEpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBmb3IgKGxldCBkID0gMDsgZCA8IGRpbTsgKytkKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IGdyYWRfZCA9IGNsaXAoZ3JhZF9jb2VmZiAqIChjdXJyZW50W2RdIC0gb3RoZXJbZF0pKSAqIGFscGhhO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBjID0gY3VycmVudFtkXSArIGdyYWRfZDtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgbyA9IG90aGVyW2RdIC0gZ3JhZF9kO1xuICAgICAgICAgICAgICAgICAgICBjdXJyZW50W2RdID0gYztcbiAgICAgICAgICAgICAgICAgICAgb3RoZXJbZF0gPSBvO1xuICAgICAgICAgICAgICAgICAgICBoZWFkX2VtYmVkZGluZy5zZXRfZW50cnkoaiwgZCwgYyk7XG4gICAgICAgICAgICAgICAgICAgIHRhaWxfZW1iZWRkaW5nLnNldF9lbnRyeShrLCBkLCBvKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZXBvY2hfb2ZfbmV4dF9zYW1wbGVbaV0gKz0gZXBvY2hzX3Blcl9zYW1wbGVbaV07XG4gICAgICAgICAgICAgICAgY29uc3Qgbl9uZWdfc2FtcGxlcyA9ICh0aGlzLl9pdGVyIC0gZXBvY2hfb2ZfbmV4dF9uZWdhdGl2ZV9zYW1wbGVbaV0pIC8gZXBvY2hzX3Blcl9uZWdhdGl2ZV9zYW1wbGVbaV07XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgcCA9IDA7IHAgPCBuX25lZ19zYW1wbGVzOyArK3ApIHtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgayA9IE1hdGguZmxvb3IodGhpcy5fcmFuZG9taXplci5yYW5kb20gKiB0YWlsX2xlbmd0aCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IG90aGVyID0gdGFpbF9lbWJlZGRpbmcucm93KHRhaWxba10pO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBkaXN0ID0gZXVjbGlkZWFuX3NxdWFyZWQoY3VycmVudCwgb3RoZXIpO1xuICAgICAgICAgICAgICAgICAgICBsZXQgZ3JhZF9jb2VmZiA9IDA7XG4gICAgICAgICAgICAgICAgICAgIGlmIChkaXN0ID4gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgZ3JhZF9jb2VmZiA9ICgyICogcmVwdWxzaW9uX3N0cmVuZ3RoICogYikgLyAoKC4wMSArIGRpc3QpICogKGEgKiBNYXRoLnBvdyhkaXN0LCBiKSArIDEpKTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChqID09PSBrKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBmb3IgKGxldCBkID0gMDsgZCA8IGRpbTsgKytkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjb25zdCBncmFkX2QgPSBjbGlwKGdyYWRfY29lZmYgKiAoY3VycmVudFtkXSAtIG90aGVyW2RdKSkgKiBhbHBoYTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGNvbnN0IGMgPSBjdXJyZW50W2RdICsgZ3JhZF9kO1xuICAgICAgICAgICAgICAgICAgICAgICAgY29uc3QgbyA9IG90aGVyW2RdIC0gZ3JhZF9kO1xuICAgICAgICAgICAgICAgICAgICAgICAgY3VycmVudFtkXSA9IGM7XG4gICAgICAgICAgICAgICAgICAgICAgICBvdGhlcltkXSA9IG87XG4gICAgICAgICAgICAgICAgICAgICAgICBoZWFkX2VtYmVkZGluZy5zZXRfZW50cnkoaiwgZCwgYyk7XG4gICAgICAgICAgICAgICAgICAgICAgICB0YWlsX2VtYmVkZGluZy5zZXRfZW50cnkodGFpbFtrXSwgZCwgbyk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZXBvY2hfb2ZfbmV4dF9uZWdhdGl2ZV9zYW1wbGVbaV0gKz0gKG5fbmVnX3NhbXBsZXMgKiBlcG9jaHNfcGVyX25lZ2F0aXZlX3NhbXBsZVtpXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIGhlYWRfZW1iZWRkaW5nO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICBuZXh0KCkge1xuICAgICAgICBsZXQgaXRlciA9ICsrdGhpcy5faXRlcjtcbiAgICAgICAgbGV0IFkgPSB0aGlzLlk7XG5cbiAgICAgICAgdGhpcy5fYWxwaGEgPSAodGhpcy5faW5pdGlhbF9hbHBoYSAqICgxIC0gaXRlciAvIHRoaXMuX25fZXBvY2hzKSk7XG4gICAgICAgIHRoaXMuWSA9IHRoaXMuX29wdGltaXplX2xheW91dChZLCBZLCB0aGlzLl9oZWFkLCB0aGlzLl90YWlsKTtcblxuICAgICAgICByZXR1cm4gdGhpcy5ZO1xuICAgIH1cbn0gIiwiaW1wb3J0IHsgTWF0cml4LCBsaW5zcGFjZSB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBQQ0EgfSBmcm9tIFwiLi9QQ0EuanNcIjtcbmltcG9ydCB7IEJhbGxUcmVlIH0gZnJvbSBcIi4uL2tubi9pbmRleC5qc1wiO1xuaW1wb3J0IHsgRFIgfSBmcm9tIFwiLi9EUi5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIFRyaU1hcFxuICogQGV4dGVuZHMgRFJcbiAqL1xuZXhwb3J0IGNsYXNzIFRyaU1hcCBleHRlbmRzIERSe1xuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGltZW5zaW9uYWxpdHlfcmVkdWN0aW9uXG4gICAgICogQGFsaWFzIFRyaU1hcFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gdGhlIGhpZ2gtZGltZW5zaW9uYWwgZGF0YS4gXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFt3ZWlnaHRfYWRqID0gNTAwXSAtIHNjYWxpbmcgZmFjdG9yLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbYyA9IDVdIC0gbnVtYmVyIG9mIHRyaXBsZXRzIG11bHRpcGxpZXIuXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtkID0gMl0gLSB0aGUgZGltZW5zaW9uYWxpdHkgb2YgdGhlIHByb2plY3Rpb24uXG4gICAgICogQHBhcmFtIHtGdW5jdGlvbn0gW21ldHJpYyA9IGV1Y2xpZGVhbl0gLSB0aGUgbWV0cmljIHdoaWNoIGRlZmluZXMgdGhlIGRpc3RhbmNlIGJldHdlZW4gdHdvIHBvaW50cy4gIFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbc2VlZCA9IDEyMTJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEByZXR1cm5zIHtUcmlNYXB9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9hcnhpdi5vcmcvcGRmLzE5MTAuMDAyMDR2MS5wZGZ9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9naXRodWIuY29tL2VhbWlkL3RyaW1hcH1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCB3ZWlnaHRfYWRqID0gNTAwLCBjID0gNSwgZCA9IDIsIG1ldHJpYyA9IGV1Y2xpZGVhbiwgc2VlZD0xMjEyKSB7XG4gICAgICAgIHN1cGVyKFgsIGQsIG1ldHJpYywgc2VlZCk7XG4gICAgICAgIHN1cGVyLnBhcmFtZXRlcl9saXN0ID0gW1wid2VpZ2h0X2FkalwiLCBcImNcIl07XG4gICAgICAgIHRoaXMucGFyYW1ldGVyKFwid2VpZ2h0X2FkalwiLCB3ZWlnaHRfYWRqKTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJjXCIsIGMpXG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBbcGNhID0gbnVsbF0gLSBJbml0aWFsIEVtYmVkZGluZyAoaWYgbnVsbCB0aGVuIFBDQSBnZXRzIHVzZWQpLiBcbiAgICAgKiBAcGFyYW0ge0tOTn0gW2tubiA9IG51bGxdIC0gS05OIE9iamVjdCAoaWYgbnVsbCB0aGVuIEJhbGxUcmVlIGdldHMgdXNlZCkuIFxuICAgICAqL1xuICAgIGluaXQocGNhID0gbnVsbCwga25uID0gbnVsbCkge1xuICAgICAgICBjb25zdCBYID0gdGhpcy5YO1xuICAgICAgICBjb25zdCBOID0gWC5zaGFwZVswXTtcbiAgICAgICAgY29uc3QgZCA9IHRoaXMuX2Q7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgYyA9IHRoaXMuX2M7XG4gICAgICAgIHRoaXMubl9pbmxpZXJzID0gMiAqIGM7XG4gICAgICAgIHRoaXMubl9vdXRsaWVycyA9IDEgKiBjO1xuICAgICAgICB0aGlzLm5fcmFuZG9tID0gMSAqIGM7XG4gICAgICAgIHRoaXMuWSA9IHBjYSB8fCBuZXcgUENBKFgsIGQpLnRyYW5zZm9ybSgpLy8ubXVsdCguMDEpO1xuICAgICAgICB0aGlzLmtubiA9IGtubiB8fCBuZXcgQmFsbFRyZWUoWC50bzJkQXJyYXksIG1ldHJpYyk7XG4gICAgICAgIGNvbnN0IHt0cmlwbGV0cywgd2VpZ2h0c30gPSB0aGlzLl9nZW5lcmF0ZV90cmlwbGV0cyh0aGlzLm5faW5saWVycywgdGhpcy5uX291dGxpZXJzLCB0aGlzLm5fcmFuZG9tKTtcbiAgICAgICAgdGhpcy50cmlwbGV0cyA9IHRyaXBsZXRzO1xuICAgICAgICB0aGlzLndlaWdodHMgPSB3ZWlnaHRzO1xuICAgICAgICB0aGlzLmxyID0gMTAwMCAqIE4gLyB0cmlwbGV0cy5zaGFwZVswXTtcbiAgICAgICAgdGhpcy5DID0gSW5maW5pdHk7XG4gICAgICAgIHRoaXMudG9sID0gMWUtNztcbiAgICAgICAgdGhpcy52ZWwgPSBuZXcgTWF0cml4KE4sIGQsIDApO1xuICAgICAgICB0aGlzLmdhaW4gPSBuZXcgTWF0cml4KE4sIGQsIDEpO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBHZW5lcmF0ZXMge0BsaW5rIG5faW5saWVyc30geCB7QGxpbmsgbl9vdXRsaWVyc30geCB7QGxpbmsgbl9yYW5kb219IHRyaXBsZXRzLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBuX2lubGllcnMgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG5fb3V0bGllcnMgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG5fcmFuZG9tIFxuICAgICAqL1xuICAgIF9nZW5lcmF0ZV90cmlwbGV0cyhuX2lubGllcnMsIG5fb3V0bGllcnMsIG5fcmFuZG9tKSB7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3Qgd2VpZ2h0X2FkaiA9IHRoaXMuX3dlaWdodF9hZGo7XG4gICAgICAgIGNvbnN0IFggPSB0aGlzLlg7XG4gICAgICAgIGNvbnN0IE4gPSBYLnNoYXBlWzBdO1xuICAgICAgICBjb25zdCBrbm4gPSB0aGlzLmtubjtcbiAgICAgICAgY29uc3Qgbl9leHRyYSA9IE1hdGgubWluKG5faW5saWVycyArIDIwLCBOKTtcbiAgICAgICAgY29uc3QgbmJycyA9IG5ldyBNYXRyaXgoTiwgbl9leHRyYSk7XG4gICAgICAgIGNvbnN0IGtubl9kaXN0YW5jZXMgPSBuZXcgTWF0cml4KE4sIG5fZXh0cmEpO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAga25uLnNlYXJjaChYLnJvdyhpKSwgbl9leHRyYSArIDEpXG4gICAgICAgICAgICAgICAgLnJhd19kYXRhKClcbiAgICAgICAgICAgICAgICAuZmlsdGVyKGQgPT4gZC52YWx1ZSAhPSAwKVxuICAgICAgICAgICAgICAgIC5zb3J0KChhLCBiKSA9PiBhLnZhbHVlIC0gYi52YWx1ZSlcbiAgICAgICAgICAgICAgICAuZm9yRWFjaCgoZCwgaikgPT4ge1xuICAgICAgICAgICAgICAgICAgICBuYnJzLnNldF9lbnRyeShpLCBqLCBkLmVsZW1lbnQuaW5kZXgpXG4gICAgICAgICAgICAgICAgICAgIGtubl9kaXN0YW5jZXMuc2V0X2VudHJ5KGksIGosIGQudmFsdWUpXG4gICAgICAgICAgICAgICAgfSk7XG4gICAgICAgIH1cbiAgICAgICAgLy8gc2NhbGUgcGFyYW1ldGVyXG4gICAgICAgIGNvbnN0IHNpZyA9IG5ldyBGbG9hdDY0QXJyYXkoTik7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBzaWdbaV0gPSBNYXRoLm1heChcbiAgICAgICAgICAgICAgICAgICAoa25uX2Rpc3RhbmNlcy5lbnRyeShpLCAzKSArXG4gICAgICAgICAgICAgICAgICAgIGtubl9kaXN0YW5jZXMuZW50cnkoaSwgNCkgK1xuICAgICAgICAgICAgICAgICAgICBrbm5fZGlzdGFuY2VzLmVudHJ5KGksIDUpICtcbiAgICAgICAgICAgICAgICAgICAga25uX2Rpc3RhbmNlcy5lbnRyeShpLCA2KSkgLyA0LFxuICAgICAgICAgICAgICAgICAgICAxZS0xMCk7XG4gICAgICAgIH1cbiAgICAgICAgXG4gICAgICAgIGNvbnN0IFAgPSB0aGlzLl9maW5kX3Aoa25uX2Rpc3RhbmNlcywgc2lnLCBuYnJzKTtcbiAgICAgICAgXG4gICAgICAgIGxldCB0cmlwbGV0cyA9IHRoaXMuX3NhbXBsZV9rbm5fdHJpcGxldHMoUCwgbmJycywgbl9pbmxpZXJzLCBuX291dGxpZXJzKTtcbiAgICAgICAgbGV0IG5fdHJpcGxldHMgPSB0cmlwbGV0cy5zaGFwZVswXTtcbiAgICAgICAgY29uc3Qgb3V0bGllcl9kaXN0YW5jZXMgPSBuZXcgRmxvYXQ2NEFycmF5KG5fdHJpcGxldHMpO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IG5fdHJpcGxldHM7ICsraSkge1xuICAgICAgICAgICAgY29uc3QgaiA9IHRyaXBsZXRzLmVudHJ5KGksIDApO1xuICAgICAgICAgICAgY29uc3QgayA9IHRyaXBsZXRzLmVudHJ5KGksIDIpO1xuICAgICAgICAgICAgb3V0bGllcl9kaXN0YW5jZXNbaV0gPSBtZXRyaWMoWC5yb3coaiksIFgucm93KGspKTtcbiAgICAgICAgfVxuICAgICAgICBsZXQgd2VpZ2h0cyA9IHRoaXMuX2ZpbmRfd2VpZ2h0cyh0cmlwbGV0cywgUCwgbmJycywgb3V0bGllcl9kaXN0YW5jZXMsIHNpZyk7XG4gICAgICAgIFxuICAgICAgICBpZiAobl9yYW5kb20gPiAwKSB7XG4gICAgICAgICAgICBjb25zdCB7cmFuZG9tX3RyaXBsZXRzLCByYW5kb21fd2VpZ2h0c30gPSB0aGlzLl9zYW1wbGVfcmFuZG9tX3RyaXBsZXRzKFgsIG5fcmFuZG9tLCBzaWcpO1xuICAgICAgICAgICAgdHJpcGxldHMgPSB0cmlwbGV0cy5jb25jYXQocmFuZG9tX3RyaXBsZXRzLCBcInZlcnRpY2FsXCIpO1xuICAgICAgICAgICAgd2VpZ2h0cyA9IEZsb2F0NjRBcnJheS5mcm9tKFsuLi53ZWlnaHRzLCAuLi5yYW5kb21fd2VpZ2h0c10pXG4gICAgICAgIH1cbiAgICAgICAgbl90cmlwbGV0cyA9IHRyaXBsZXRzLnNoYXBlWzBdO1xuICAgICAgICBsZXQgbWF4X3dlaWdodCA9IC1JbmZpbml0eTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuX3RyaXBsZXRzOyArK2kpIHtcbiAgICAgICAgICAgIGlmIChpc05hTih3ZWlnaHRzW2ldKSkge3dlaWdodHNbaV0gPSAwO31cbiAgICAgICAgICAgIGlmIChtYXhfd2VpZ2h0IDwgd2VpZ2h0c1tpXSkgbWF4X3dlaWdodCA9IHdlaWdodHNbaV07XG4gICAgICAgIH1cbiAgICAgICAgbGV0IG1heF93ZWlnaHRfMiA9IC1JbmZpbml0eTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuX3RyaXBsZXRzOyArK2kpIHtcbiAgICAgICAgICAgIHdlaWdodHNbaV0gLz0gbWF4X3dlaWdodDtcbiAgICAgICAgICAgIHdlaWdodHNbaV0gKz0gLjAwMDE7XG4gICAgICAgICAgICB3ZWlnaHRzW2ldID0gTWF0aC5sb2coMSArIHdlaWdodF9hZGogKiB3ZWlnaHRzW2ldKTtcbiAgICAgICAgICAgIGlmIChtYXhfd2VpZ2h0XzIgPCB3ZWlnaHRzW2ldKSBtYXhfd2VpZ2h0XzIgPSB3ZWlnaHRzW2ldO1xuICAgICAgICB9XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbl90cmlwbGV0czsgKytpKSB7XG4gICAgICAgICAgICB3ZWlnaHRzW2ldIC89IG1heF93ZWlnaHRfMjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgXCJ0cmlwbGV0c1wiOiB0cmlwbGV0cyxcbiAgICAgICAgICAgIFwid2VpZ2h0c1wiOiB3ZWlnaHRzLFxuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQ2FsY3VsYXRlcyB0aGUgc2ltaWxhcml0eSBtYXRyaXggUFxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IGtubl9kaXN0YW5jZXMgLSBtYXRyaXggb2YgcGFpcndpc2Uga25uIGRpc3RhbmNlc1xuICAgICAqIEBwYXJhbSB7RmxvYXQ2NEFycmF5fSBzaWcgLSBzY2FsaW5nIGZhY3RvciBmb3IgdGhlIGRpc3RhbmNlc1xuICAgICAqIEBwYXJhbSB7TWF0cml4fSBuYnJzIC0gbmVhcmVzdCBuZWlnaGJvcnNcbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fSBwYWlyd2lzZSBzaW1pbGFyaXR5IG1hdHJpeFxuICAgICAqL1xuICAgIF9maW5kX3Aoa25uX2Rpc3RhbmNlcywgc2lnLCBuYnJzKSB7XG4gICAgICAgIGNvbnN0IFtOLCBuX25laWdoYm9yc10gPSBrbm5fZGlzdGFuY2VzLnNoYXBlO1xuICAgICAgICByZXR1cm4gbmV3IE1hdHJpeChOLCBuX25laWdoYm9ycywgKGksIGopID0+IHtcbiAgICAgICAgICAgIHJldHVybiBNYXRoLmV4cCgtKChrbm5fZGlzdGFuY2VzLmVudHJ5KGksIGopICoqIDIpIC8gc2lnW2ldIC8gc2lnW25icnMuZW50cnkoaSwgaildKSk7XG4gICAgICAgIH0pO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFNhbXBsZSBuZWFyZXN0IG5laWdoYm9ycyB0cmlwbGV0cyBiYXNlZCBvbiB0aGUgc2ltaWxhcml0eSB2YWx1ZXMgZ2l2ZW4gaW4gUC5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBQIC0gTWF0cml4IG9mIHBhaXJ3aXNlIHNpbWlsYXJpdGllcyBiZXR3ZWVuIGVhY2ggcG9pbnQgYW5kIGl0cyBuZWlnaGJvcnMgZ2l2ZW4gaW4gbWF0cml4IG5icnMuXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IG5icnMgLSBOZWFyZXN0IG5laWdoYm9ycyBpbmRpY2VzIGZvciBlYWNoIHBvaW50LiBUaGUgc2ltaWxhcml0eSB2YWx1ZXMgYXJlIGdpdmVuIGluIG1hdHJpeCB7QGxpbmsgUH0uIFJvdyBpIGNvcnJlc3BvbmRzIHRvIHRoZSBpLXRoIHBvaW50LlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBuX2lubGllcnMgLSBOdW1iZXIgb2YgaW5saWVyIHBvaW50cy5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gbl9vdXRsaWVycyAtIE51bWJlciBvZiBvdXRsaWVyIHBvaW50cy5cbiAgICAgKiBcbiAgICAgKi9cbiAgICBfc2FtcGxlX2tubl90cmlwbGV0cyhQLCBuYnJzLCBuX2lubGllcnMsIG5fb3V0bGllcnMpIHtcbiAgICAgICAgY29uc3QgTiA9IG5icnMuc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IHRyaXBsZXRzID0gbmV3IE1hdHJpeChOICogbl9pbmxpZXJzICogbl9vdXRsaWVycywgMyk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBsZXQgbl9pID0gaSAqIG5faW5saWVycyAqIG5fb3V0bGllcnM7XG4gICAgICAgICAgICBjb25zdCBzb3J0X2luZGljZXMgPSB0aGlzLl9fYXJnc29ydChQLnJvdyhpKS5tYXAoZCA9PiAtZCkpO1xuICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBuX2lubGllcnM7ICsraikge1xuICAgICAgICAgICAgICAgIGxldCBuX2ogPSBqICogbl9vdXRsaWVycztcbiAgICAgICAgICAgICAgICBjb25zdCBzaW0gPSBuYnJzLmVudHJ5KGksIHNvcnRfaW5kaWNlc1tqXSk7XG4gICAgICAgICAgICAgICAgY29uc3Qgc2FtcGxlcyA9IHRoaXMuX3JlamVjdGlvbl9zYW1wbGUobl9vdXRsaWVycywgTiwgc29ydF9pbmRpY2VzLnNsaWNlKDAsIGogKyAxKSk7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgayA9IDA7IGsgPCBuX291dGxpZXJzOyArK2spIHtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgaW5kZXggPSBuX2kgKyBuX2ogKyBrO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBvdXQgPSBzYW1wbGVzW2tdO1xuICAgICAgICAgICAgICAgICAgICB0cmlwbGV0cy5zZXRfZW50cnkoaW5kZXgsIDAsIGkpO1xuICAgICAgICAgICAgICAgICAgICB0cmlwbGV0cy5zZXRfZW50cnkoaW5kZXgsIDEsIHNpbSk7XG4gICAgICAgICAgICAgICAgICAgIHRyaXBsZXRzLnNldF9lbnRyeShpbmRleCwgMiwgb3V0KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRyaXBsZXRzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFNob3VsZCBkbyB0aGUgc2FtZSBhcyBucC5hcmdzb3J0KClcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7QXJyYXl9IEEgXG4gICAgICovXG4gICAgX19hcmdzb3J0KEEpIHtcbiAgICAgICAgcmV0dXJuIEFcbiAgICAgICAgICAgIC5tYXAoKGQsIGkpID0+IHtyZXR1cm4ge2Q6IGQsIGk6IGl9O30pXG4gICAgICAgICAgICAuc29ydCgoYSwgYikgPT4gYS5kIC0gYi5kKVxuICAgICAgICAgICAgLm1hcCgoZCkgPT4gZC5pKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBTYW1wbGVzIHtAbGluayBuX3NhbXBsZXN9IGludGVnZXJzIGZyb20gYSBnaXZlbiBpbnRlcnZhbCBbMCwge0BsaW5rIG1heF9pbnR9XSB3aGlsZSByZWplY3Rpb24gdGhlIHZhbHVlcyB0aGF0IGFyZSBpbiB0aGUge0BsaW5rIHJlamVjdHN9LlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHsqfSBuX3NhbXBsZXMgXG4gICAgICogQHBhcmFtIHsqfSBtYXhfaW50IFxuICAgICAqIEBwYXJhbSB7Kn0gcmVqZWN0cyBcbiAgICAgKi9cbiAgICBfcmVqZWN0aW9uX3NhbXBsZShuX3NhbXBsZXMsIG1heF9pbnQsIHJlamVjdHMpIHtcbiAgICAgICAgY29uc3QgcmFuZG9taXplciA9IHRoaXMuX3JhbmRvbWl6ZXI7XG4gICAgICAgIGNvbnN0IGludGVydmFsID0gbGluc3BhY2UoMCwgbWF4X2ludCAtIDEpLmZpbHRlcihkID0+IHJlamVjdHMuaW5kZXhPZihkKSA8IDApO1xuICAgICAgICByZXR1cm4gcmFuZG9taXplci5jaG9pY2UoaW50ZXJ2YWwsIE1hdGgubWluKG5fc2FtcGxlcywgaW50ZXJ2YWwubGVuZ3RoIC0gMikpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENhbGN1bGF0ZXMgdGhlIHdlaWdodHMgZm9yIHRoZSBzYW1wbGVkIG5lYXJlc3QgbmVpZ2hib3JzIHRyaXBsZXRzXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gdHJpcGxldHMgLSBTYW1wbGVkIFRyaXBsZXRzLlxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBQIC0gUGFpcndpc2Ugc2ltaWxhcml0eSBtYXRyaXguXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IG5icnMgLSBuZWFyZXN0IE5laWdoYm9yc1xuICAgICAqIEBwYXJhbSB7RmxvYXQ2NEFycmF5fSBvdXRsaWVyX2Rpc3RhbmNlcyAtIE1hdHJpeCBvZiBwYWlyd2lzZSBvdXRsaWVyIGRpc3RhbmNlc1xuICAgICAqIEBwYXJhbSB7RmxvYXQ2NEFycmF5fSBzaWcgLSBzY2FsaW5nIGZhY3RvciBmb3IgdGhlIGRpc3RhbmNlcy5cbiAgICAgKi9cbiAgICBfZmluZF93ZWlnaHRzKHRyaXBsZXRzLCBQLCBuYnJzLCBvdXRsaWVyX2Rpc3RhbmNlcywgc2lnKSB7XG4gICAgICAgIGNvbnN0IG5fdHJpcGxldHMgPSB0cmlwbGV0cy5zaGFwZVswXTtcbiAgICAgICAgY29uc3Qgd2VpZ2h0cyA9IG5ldyBGbG9hdDY0QXJyYXkobl90cmlwbGV0cyk7XG4gICAgICAgIGZvciAobGV0IHQgPSAwOyB0IDwgbl90cmlwbGV0czsgKyt0KSB7XG4gICAgICAgICAgICBjb25zdCBpID0gdHJpcGxldHMuZW50cnkodCwgMCk7XG4gICAgICAgICAgICBjb25zdCBzaW0gPSBuYnJzLnJvdyhpKS5pbmRleE9mKHRyaXBsZXRzLmVudHJ5KHQsIDEpKTtcbiAgICAgICAgICAgIGNvbnN0IHBfc2ltID0gUC5lbnRyeShpLCBzaW0pO1xuICAgICAgICAgICAgbGV0IHBfb3V0ID0gTWF0aC5leHAoLShvdXRsaWVyX2Rpc3RhbmNlc1t0XSAqKiAyIC8gKHNpZ1tpXSAqIHNpZ1t0cmlwbGV0cy5lbnRyeSh0LCAyKV0pKSk7XG4gICAgICAgICAgICBpZiAocF9vdXQgPCAxZS0yMCkgcF9vdXQgPSAxZS0yMDtcbiAgICAgICAgICAgIHdlaWdodHNbdF0gPSBwX3NpbSAvIHBfb3V0O1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB3ZWlnaHRzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFNhbXBsZSB1bmlmb3JtbHkgcmFub20gdHJpcGxldHNcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gRGF0YSBtYXRyaXguXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG5fcmFuZG9tIC0gTnVtYmVyIG9mIHJhbmRvbSB0cmlwbGV0cyBwZXIgcG9pbnRcbiAgICAgKiBAcGFyYW0ge0Zsb2F0NjRBcnJheX0gc2lnIC0gU2NhbGluZyBmYWN0b3IgZm9yIHRoZSBkaXN0YW5jZXNcbiAgICAgKi9cbiAgICBfc2FtcGxlX3JhbmRvbV90cmlwbGV0cyhYLCBuX3JhbmRvbSwgc2lnKSB7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgcmFuZG9taXplciA9IHRoaXMuX3JhbmRvbWl6ZXI7XG4gICAgICAgIGNvbnN0IE4gPSBYLnNoYXBlWzBdO1xuICAgICAgICBjb25zdCByYW5kb21fdHJpcGxldHMgPSBuZXcgTWF0cml4KE4gKiBuX3JhbmRvbSwgMyk7XG4gICAgICAgIGNvbnN0IHJhbmRvbV93ZWlnaHRzID0gbmV3IEZsb2F0NjRBcnJheShOICogbl9yYW5kb20pO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgY29uc3Qgbl9pID0gaSAqIG5fcmFuZG9tO1xuICAgICAgICAgICAgY29uc3QgaW5kaWNlcyA9IFsuLi5saW5zcGFjZSgwLCBpIC0gMSksIC4uLmxpbnNwYWNlKGkgKyAxLCBOIC0gMSldXG4gICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IG5fcmFuZG9tOyArK2opIHtcbiAgICAgICAgICAgICAgICBsZXQgW3NpbSwgb3V0XSA9IHJhbmRvbWl6ZXIuY2hvaWNlKGluZGljZXMsIDIpO1xuICAgICAgICAgICAgICAgIGxldCBwX3NpbSA9IE1hdGguZXhwKC0oKG1ldHJpYyhYLnJvdyhpKSwgWC5yb3coc2ltKSkgKiogMikgLyAoc2lnW2ldICogc2lnW3NpbV0pKSk7XG4gICAgICAgICAgICAgICAgaWYgKHBfc2ltIDwgMWUtMjApIHBfc2ltID0gMWUtMjA7XG4gICAgICAgICAgICAgICAgbGV0IHBfb3V0ID0gTWF0aC5leHAoLSgobWV0cmljKFgucm93KGkpLCBYLnJvdyhvdXQpKSAqKiAyKSAvIChzaWdbaV0gKiBzaWdbb3V0XSkpKTsgXG4gICAgICAgICAgICAgICAgaWYgKHBfb3V0IDwgMWUtMjApIHBfb3V0ID0gMWUtMjA7XG5cbiAgICAgICAgICAgICAgICBpZiAocF9zaW0gPCBwX291dCkge1xuICAgICAgICAgICAgICAgICAgICBbc2ltLCBvdXRdID0gW291dCwgc2ltXTtcbiAgICAgICAgICAgICAgICAgICAgW3Bfc2ltLCBwX291dF0gPSBbcF9vdXQsIHBfc2ltXTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgY29uc3QgaW5kZXggPSBuX2kgKyBqO1xuICAgICAgICAgICAgICAgIHJhbmRvbV90cmlwbGV0cy5zZXRfZW50cnkoaW5kZXgsIDAsIGkpO1xuICAgICAgICAgICAgICAgIHJhbmRvbV90cmlwbGV0cy5zZXRfZW50cnkoaW5kZXgsIDEsIHNpbSk7XG4gICAgICAgICAgICAgICAgcmFuZG9tX3RyaXBsZXRzLnNldF9lbnRyeShpbmRleCwgMiwgb3V0KTtcbiAgICAgICAgICAgICAgICByYW5kb21fd2VpZ2h0c1tpbmRleF0gPSBwX3NpbSAvIHBfb3V0O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICBcInJhbmRvbV90cmlwbGV0c1wiOiByYW5kb21fdHJpcGxldHMsXG4gICAgICAgICAgICBcInJhbmRvbV93ZWlnaHRzXCI6IHJhbmRvbV93ZWlnaHRzLFxuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQ29tcHV0ZXMgdGhlIGdyYWRpZW50IGZvciB1cGRhdGluZyB0aGUgZW1iZWRkaW5nLlxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBZIC0gVGhlIGVtYmVkZGluZ1xuICAgICAqL1xuICAgIF9ncmFkKFkpIHtcbiAgICAgICAgY29uc3Qgbl9pbmxpZXJzID0gdGhpcy5uX2lubGllcnM7XG4gICAgICAgIGNvbnN0IG5fb3V0bGllcnMgPSB0aGlzLm5fb3V0bGllcnM7XG4gICAgICAgIGNvbnN0IHRyaXBsZXRzID0gdGhpcy50cmlwbGV0cztcbiAgICAgICAgY29uc3Qgd2VpZ2h0cyA9IHRoaXMud2VpZ2h0cztcbiAgICAgICAgY29uc3QgW04sIGRpbV0gPSBZLnNoYXBlO1xuICAgICAgICBjb25zdCBuX3RyaXBsZXRzID0gdHJpcGxldHMuc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IGdyYWQgPSBuZXcgTWF0cml4KE4sIGRpbSwgMCk7XG4gICAgICAgIGxldCB5X2lqID0gbmV3IEFycmF5KGRpbSkuZmlsbCgwKTtcbiAgICAgICAgbGV0IHlfaWsgPSBuZXcgQXJyYXkoZGltKS5maWxsKDApO1xuICAgICAgICBsZXQgZF9paiA9IDE7XG4gICAgICAgIGxldCBkX2lrID0gMTtcbiAgICAgICAgbGV0IG5fdmlvbCA9IDA7XG4gICAgICAgIGxldCBsb3NzID0gMDtcbiAgICAgICAgY29uc3Qgbl9rbm5fdHJpcGxldHMgPSBOICogbl9pbmxpZXJzICogbl9vdXRsaWVycztcblxuICAgICAgICBmb3IgKGxldCB0ID0gMDsgdCA8IG5fdHJpcGxldHM7ICsrdCkge1xuICAgICAgICAgICAgY29uc3QgW2ksIGosIGtdID0gdHJpcGxldHMucm93KHQpO1xuICAgICAgICAgICAgLy8gdXBkYXRlIHlfaWosIHlfaWssIGRfaWosIGRfaWtcbiAgICAgICAgICAgIGlmICh0ICUgbl9vdXRsaWVycyA9PSAwIHx8IHQgPj0gbl9rbm5fdHJpcGxldHMpIHtcbiAgICAgICAgICAgICAgICBkX2lqID0gMVxuICAgICAgICAgICAgICAgIGRfaWsgPSAxXG4gICAgICAgICAgICAgICAgZm9yIChsZXQgZCA9IDA7IGQgPCBkaW07ICsrZCkge1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBZX2lkID0gWS5lbnRyeShpLCBkKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgWV9qZCA9IFkuZW50cnkoaiwgZCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IFlfa2QgPSBZLmVudHJ5KGssIGQpO1xuICAgICAgICAgICAgICAgICAgICB5X2lqW2RdID0gWV9pZCAtIFlfamQ7XG4gICAgICAgICAgICAgICAgICAgIHlfaWtbZF0gPSBZX2lkIC0gWV9rZDtcbiAgICAgICAgICAgICAgICAgICAgZF9paiArPSAoeV9paltkXSAqKiAyKTtcbiAgICAgICAgICAgICAgICAgICAgZF9payArPSAoeV9pa1tkXSAqKiAyKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAvLyB1cGRhdGUgeV9payBhbmQgZF9payBvbmx5XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIGRfaWsgPSAxO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGQgPSAwOyBkIDwgZGltOyArK2QpIHtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgWV9pZCA9IFkuZW50cnkoaSwgZCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IFlfa2QgPSBZLmVudHJ5KGssIGQpO1xuICAgICAgICAgICAgICAgICAgICB5X2lrW2RdID0gWV9pZCAtIFlfa2Q7XG4gICAgICAgICAgICAgICAgICAgIGRfaWsgKz0gKHlfaWtbZF0gKiogMik7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBpZiAoZF9paiA+IGRfaWspICsrbl92aW9sO1xuICAgICAgICAgICAgbG9zcyArPSB3ZWlnaHRzW3RdIC8gKDEgKyBkX2lrIC8gZF9paik7XG4gICAgICAgICAgICBjb25zdCB3ID0gKHdlaWdodHNbdF0gLyAoZF9paiArIGRfaWspKSAqKiAyO1xuICAgICAgICAgICAgZm9yIChsZXQgZCA9IDA7IGQgPCBkaW07ICsrZCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IGdzID0geV9paltkXSAqIGRfaWsgKiB3O1xuICAgICAgICAgICAgICAgIGNvbnN0IGdvID0geV9pa1tkXSAqIGRfaWogKiB3O1xuICAgICAgICAgICAgICAgIGdyYWQuc2V0X2VudHJ5KGksIGQsIGdyYWQuZW50cnkoaSwgZCkgKyBncyAtIGdvKTtcbiAgICAgICAgICAgICAgICBncmFkLnNldF9lbnRyeShqLCBkLCBncmFkLmVudHJ5KGosIGQpIC0gZ3MpO1xuICAgICAgICAgICAgICAgIGdyYWQuc2V0X2VudHJ5KGssIGQsIGdyYWQuZW50cnkoaywgZCkgKyBnbyk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHtcbiAgICAgICAgICAgIFwiZ3JhZFwiOiBncmFkLFxuICAgICAgICAgICAgXCJsb3NzXCI6IGxvc3MsXG4gICAgICAgICAgICBcIm5fdmlvbFwiOiBuX3Zpb2wsXG4gICAgICAgIH07XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IG1heF9pdGVyYXRpb24gXG4gICAgICovXG4gICAgdHJhbnNmb3JtKG1heF9pdGVyYXRpb24gPSA0MDApIHtcbiAgICAgICAgdGhpcy5jaGVja19pbml0KCk7XG4gICAgICAgIGZvciAobGV0IGl0ZXIgPSAwOyBpdGVyIDwgbWF4X2l0ZXJhdGlvbjsgKytpdGVyKSB7XG4gICAgICAgICAgICB0aGlzLl9uZXh0KGl0ZXIpXG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAeWllbGRzIHtNYXRyaXh9XG4gICAgICogQHJldHVybnMge01hdHJpeH1cbiAgICAgKi9cbiAgICAqIGdlbmVyYXRvcigpIHtcbiAgICAgICAgdGhpcy5jaGVja19pbml0KCk7XG4gICAgICAgIGZvciAobGV0IGl0ZXIgPSAwOyBpdGVyIDwgODAwOyArK2l0ZXIpIHtcbiAgICAgICAgICAgIHRoaXMuX25leHQoaXRlcik7XG4gICAgICAgICAgICB5aWVsZCB0aGlzLnByb2plY3Rpb247XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBEb2VzIHRoZSBpdGVyYXRpb24gc3RlcC5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBpdGVyIFxuICAgICAqL1xuICAgIF9uZXh0KGl0ZXIpIHtcbiAgICAgICAgY29uc3QgZ2FtbWEgPSBpdGVyID4gMTUwID8gLjUgOiAuMztcbiAgICAgICAgY29uc3Qgb2xkX0MgPSB0aGlzLkM7XG4gICAgICAgIGNvbnN0IHZlbCA9IHRoaXMudmVsO1xuICAgICAgICBjb25zdCBZID0gdGhpcy5ZLmFkZCh2ZWwubXVsdChnYW1tYSkpO1xuICAgICAgICBjb25zdCB7Z3JhZCwgbG9zcywgbl92aW9sfSA9IHRoaXMuX2dyYWQoWSk7XG4gICAgICAgIHRoaXMuQyA9IGxvc3M7XG4gICAgICAgIHRoaXMuWSA9IHRoaXMuX3VwZGF0ZV9lbWJlZGRpbmcoWSwgaXRlciwgZ3JhZCk7XG4gICAgICAgIHRoaXMubHIgKj0gKG9sZF9DID4gbG9zcyArIHRoaXMudG9sKSAgPyAxLjAxIDogLjk7XG4gICAgICAgIHJldHVybiB0aGlzLlk7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogVXBkYXRlcyB0aGUgZW1iZWRkaW5nLlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IFkgXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IGl0ZXIgXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IGdyYWQgXG4gICAgICovXG4gICAgX3VwZGF0ZV9lbWJlZGRpbmcoWSwgaXRlciwgZ3JhZCkge1xuICAgICAgICBjb25zdCBbTiwgZGltXSA9IFkuc2hhcGU7XG4gICAgICAgIGNvbnN0IGdhbW1hID0gaXRlciA+IDE1MCA/IC45IDogLjU7IC8vIG1vbWVudCBwYXJhbWV0ZXJcbiAgICAgICAgY29uc3QgbWluX2dhaW4gPSAuMDE7XG4gICAgICAgIGNvbnN0IGdhaW4gPSB0aGlzLmdhaW47XG4gICAgICAgIGNvbnN0IHZlbCA9IHRoaXMudmVsO1xuICAgICAgICBjb25zdCBsciA9IHRoaXMubHI7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBmb3IgKGxldCBkID0gMDsgZCA8IGRpbTsgKytkKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgbmV3X2dhaW4gPSAoTWF0aC5zaWduKHZlbC5lbnRyeShpLCBkKSkgIT0gTWF0aC5zaWduKGdyYWQuZW50cnkoaSwgZCkpKSA/IGdhaW4uZW50cnkoaSwgZCkgKyAuMiA6IE1hdGgubWF4KGdhaW4uZW50cnkoaSwgZCkgKiAuOCwgbWluX2dhaW4pO1xuICAgICAgICAgICAgICAgIGdhaW4uc2V0X2VudHJ5KGksIGQsIG5ld19nYWluKTtcbiAgICAgICAgICAgICAgICB2ZWwuc2V0X2VudHJ5KGksIGQsIGdhbW1hICogdmVsLmVudHJ5KGksIGQpIC0gbHIgKiBnYWluLmVudHJ5KGksIGQpICogZ3JhZC5lbnRyeShpLCBkKSk7XG4gICAgICAgICAgICAgICAgWS5zZXRfZW50cnkoaSwgZCwgWS5lbnRyeShpLCBkKSArIHZlbC5lbnRyeShpLCBkKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIFk7XG4gICAgfVxufSIsImltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBNYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIEhpZXJhcmNoaWNhbF9DbHVzdGVyaW5nXG4gKi9cbmV4cG9ydCBjbGFzcyBIaWVyYXJjaGljYWxfQ2x1c3RlcmluZyB7XG4gICAgLyoqXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTpjbHVzdGVyaW5nXG4gICAgICogQGFsaWFzIEhpZXJhcmNoaWNhbF9DbHVzdGVyaW5nXG4gICAgICogQHRvZG8gbmVlZHMgcmVzdHJ1Y3R1cmluZy5cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gLSBEYXRhIG9yIGRpc3RhbmNlIG1hdHJpeCBpZiBtZXRyaWMgaXMgJ3ByZWNvbXB1dGVkJ1xuICAgICAqIEBwYXJhbSB7KFwic2luZ2xlXCJ8XCJjb21wbGV0ZVwifFwiYXZlcmFnZVwiKX0gW2xpbmthZ2UgPSBcImNvbXBsZXRlXCJdXG4gICAgICogQHBhcmFtIHtGdW5jdGlvbnxcInByZWNvbXB1dGVkXCJ9IFttZXRyaWMgPSBldWNsaWRlYW5dXG4gICAgICogQHJldHVybnMge0hpZXJhcmNoaWNhbF9DbHVzdGVyaW5nfVxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKG1hdHJpeCwgbGlua2FnZSA9IFwiY29tcGxldGVcIiwgbWV0cmljID0gZXVjbGlkZWFuKSB7XG4gICAgICAgIHRoaXMuX2lkID0gMDtcbiAgICAgICAgdGhpcy5fbWF0cml4ID0gbWF0cml4IGluc3RhbmNlb2YgTWF0cml4ID8gbWF0cml4IDogTWF0cml4LmZyb20obWF0cml4KTtcbiAgICAgICAgdGhpcy5fbWV0cmljID0gbWV0cmljO1xuICAgICAgICB0aGlzLl9saW5rYWdlID0gbGlua2FnZTtcbiAgICAgICAgaWYgKG1ldHJpYyA9PT0gXCJwcmVjb21wdXRlZFwiICYmIHRoaXMuX21hdHJpeC5zaGFwZVswXSAhPT0gdGhpcy5fbWF0cml4LnNoYXBlWzFdKSB7XG4gICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoXCJJZiBtZXRyaWMgaXMgJ3ByZWNvbXB1dGVkJywgdGhlbiBtYXRyaXggaGFzIHRvIGJlIHNxdWFyZSFcIik7XG4gICAgICAgIH1cbiAgICAgICAgdGhpcy5pbml0KCk7XG4gICAgICAgIHRoaXMucm9vdCA9IHRoaXMuZG8oKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICpcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gdmFsdWUgLSB2YWx1ZSB3aGVyZSB0byBjdXQgdGhlIHRyZWUuXG4gICAgICogQHBhcmFtIHsoXCJkaXN0YW5jZVwifFwiZGVwdGhcIil9IFt0eXBlID0gXCJkaXN0YW5jZVwiXSAtIHR5cGUgb2YgdmFsdWUuXG4gICAgICogQHJldHVybnMge0FycmF5PEFycmF5Pn0gLSBBcnJheSBvZiBjbHVzdGVycyB3aXRoIHRoZSBpbmRpY2VzIG9mIHRoZSByb3dzIGluIGdpdmVuIHtAbGluayBtYXRyaXh9LlxuICAgICAqL1xuICAgIGdldF9jbHVzdGVycyh2YWx1ZSwgdHlwZSA9IFwiZGlzdGFuY2VcIikge1xuICAgICAgICBsZXQgY2x1c3RlcnMgPSBbXTtcbiAgICAgICAgbGV0IGFjY2Vzc29yO1xuICAgICAgICBzd2l0Y2ggKHR5cGUpIHtcbiAgICAgICAgICAgIGNhc2UgXCJkaXN0YW5jZVwiOlxuICAgICAgICAgICAgICAgIGFjY2Vzc29yID0gKGQpID0+IGQuZGlzdDtcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJkZXB0aFwiOlxuICAgICAgICAgICAgICAgIGFjY2Vzc29yID0gKGQpID0+IGQuZGVwdGg7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBkZWZhdWx0OlxuICAgICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcImludmFsaWQgdHlwZVwiKTtcbiAgICAgICAgfVxuICAgICAgICB0aGlzLl90cmF2ZXJzZSh0aGlzLnJvb3QsIGFjY2Vzc29yLCB2YWx1ZSwgY2x1c3RlcnMpO1xuICAgICAgICByZXR1cm4gY2x1c3RlcnM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge30gbm9kZVxuICAgICAqIEBwYXJhbSB7Kn0gZlxuICAgICAqIEBwYXJhbSB7Kn0gdmFsdWVcbiAgICAgKiBAcGFyYW0geyp9IHJlc3VsdFxuICAgICAqL1xuICAgIF90cmF2ZXJzZShub2RlLCBmLCB2YWx1ZSwgcmVzdWx0KSB7XG4gICAgICAgIGlmIChmKG5vZGUpIDw9IHZhbHVlKSB7XG4gICAgICAgICAgICByZXN1bHQucHVzaChub2RlLmxlYXZlcygpKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRoaXMuX3RyYXZlcnNlKG5vZGUubGVmdCwgZiwgdmFsdWUsIHJlc3VsdCk7XG4gICAgICAgICAgICB0aGlzLl90cmF2ZXJzZShub2RlLnJpZ2h0LCBmLCB2YWx1ZSwgcmVzdWx0KTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIGNvbXB1dGVzIHRoZSB0cmVlLlxuICAgICAqL1xuICAgIGluaXQoKSB7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgQSA9IHRoaXMuX21hdHJpeDtcbiAgICAgICAgY29uc3QgbiA9ICh0aGlzLl9uID0gQS5zaGFwZVswXSk7XG4gICAgICAgIGNvbnN0IGRfbWluID0gKHRoaXMuX2RfbWluID0gbmV3IEZsb2F0NjRBcnJheShuKSk7XG4gICAgICAgIGxldCBkaXN0YW5jZV9tYXRyaXg7XG4gICAgICAgIGlmIChtZXRyaWMgIT09IFwicHJlY29tcHV0ZWRcIikge1xuICAgICAgICAgICAgZGlzdGFuY2VfbWF0cml4ID0gbmV3IE1hdHJpeChuLCBuLCAwKTsgLy9uZXcgQXJyYXkobik7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgICAgIGRfbWluW2ldID0gMDtcbiAgICAgICAgICAgICAgICAvL2Rpc3RhbmNlX21hdHJpeFtpXSA9IG5ldyBGbG9hdDY0QXJyYXkobik7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBuOyArK2opIHtcbiAgICAgICAgICAgICAgICAgICAgZGlzdGFuY2VfbWF0cml4LnNldF9lbnRyeShpLCBqLCBpID09PSBqID8gSW5maW5pdHkgOiBtZXRyaWMoQS5yb3coaSksIEEucm93KGopKSk7XG4gICAgICAgICAgICAgICAgICAgIGlmIChkaXN0YW5jZV9tYXRyaXguZW50cnkoaSwgZF9taW5baV0pID4gZGlzdGFuY2VfbWF0cml4LmVudHJ5KGksIGopKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBkX21pbltpXSA9IGo7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBkaXN0YW5jZV9tYXRyaXggPSB0aGlzLl9tYXRyaXguY2xvbmUoKTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbjsgKytpKSB7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBuOyArK2opIHtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGkgPT09IGopIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRpc3RhbmNlX21hdHJpeC5zZXRfZW50cnkoaSwgaiwgSW5maW5pdHkpO1xuICAgICAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGRpc3RhbmNlX21hdHJpeC5lbnRyeShpLCBkX21pbltpXSkgPiBkaXN0YW5jZV9tYXRyaXguZW50cnkoaSwgaikpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRfbWluW2ldID0gajtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICB0aGlzLl9kaXN0YW5jZV9tYXRyaXggPSBkaXN0YW5jZV9tYXRyaXg7XG4gICAgICAgIGNvbnN0IGNsdXN0ZXJzID0gKHRoaXMuX2NsdXN0ZXJzID0gbmV3IEFycmF5KG4pKTtcbiAgICAgICAgY29uc3QgY19zaXplID0gKHRoaXMuX2Nfc2l6ZSA9IG5ldyBVaW50MTZBcnJheShuKSk7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbjsgKytpKSB7XG4gICAgICAgICAgICBjbHVzdGVyc1tpXSA9IFtdO1xuICAgICAgICAgICAgY2x1c3RlcnNbaV1bMF0gPSBuZXcgQ2x1c3Rlcih0aGlzLl9pZCsrLCBudWxsLCBudWxsLCAwLCBBLnJvdyhpKSwgaSwgMSwgMCk7XG4gICAgICAgICAgICBjX3NpemVbaV0gPSAxO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIGNvbXB1dGVzIHRoZSB0cmVlLlxuICAgICAqL1xuICAgIGRvKCkge1xuICAgICAgICBjb25zdCBuID0gdGhpcy5fbjtcbiAgICAgICAgY29uc3QgZF9taW4gPSB0aGlzLl9kX21pbjtcbiAgICAgICAgY29uc3QgRCA9IHRoaXMuX2Rpc3RhbmNlX21hdHJpeDtcbiAgICAgICAgY29uc3QgY2x1c3RlcnMgPSB0aGlzLl9jbHVzdGVycztcbiAgICAgICAgY29uc3QgY19zaXplID0gdGhpcy5fY19zaXplO1xuICAgICAgICBjb25zdCBsaW5rYWdlID0gdGhpcy5fbGlua2FnZTtcbiAgICAgICAgbGV0IHJvb3QgPSBudWxsO1xuICAgICAgICBmb3IgKGxldCBwID0gMCwgcF9tYXggPSBuIC0gMTsgcCA8IHBfbWF4OyArK3ApIHtcbiAgICAgICAgICAgIGxldCBjMSA9IDA7XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgICAgIGxldCBEX2lfbWluID0gRC5lbnRyeShpLCBkX21pbltpXSk7XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgaiA9IGkgKyAxOyBqIDwgbjsgKytqKSB7XG4gICAgICAgICAgICAgICAgICAgIGlmIChEX2lfbWluID4gRC5lbnRyeShpLCBqKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgZF9taW5baV0gPSBqO1xuICAgICAgICAgICAgICAgICAgICAgICAgRF9pX21pbiA9IEQuZW50cnkoaSwgZF9taW5baV0pO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBuOyArK2kpIHtcbiAgICAgICAgICAgICAgICBpZiAoRC5lbnRyeShpLCBkX21pbltpXSkgPCBELmVudHJ5KGMxLCBkX21pbltjMV0pKSB7XG4gICAgICAgICAgICAgICAgICAgIGMxID0gaTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBsZXQgYzIgPSBkX21pbltjMV07XG4gICAgICAgICAgICBsZXQgYzFfY2x1c3RlciA9IGNsdXN0ZXJzW2MxXVswXTtcbiAgICAgICAgICAgIGxldCBjMl9jbHVzdGVyID0gY2x1c3RlcnNbYzJdWzBdO1xuICAgICAgICAgICAgbGV0IGMxX2NsdXN0ZXJfaW5kaWNlcyA9IGMxX2NsdXN0ZXIuaXNMZWFmID8gW2MxX2NsdXN0ZXIuaW5kZXhdIDogYzFfY2x1c3Rlci5pbmRleDtcbiAgICAgICAgICAgIGxldCBjMl9jbHVzdGVyX2luZGljZXMgPSBjMl9jbHVzdGVyLmlzTGVhZiA/IFtjMl9jbHVzdGVyLmluZGV4XSA6IGMyX2NsdXN0ZXIuaW5kZXg7XG4gICAgICAgICAgICBsZXQgaW5kaWNlcyA9IGMxX2NsdXN0ZXJfaW5kaWNlcy5jb25jYXQoYzJfY2x1c3Rlcl9pbmRpY2VzKTtcbiAgICAgICAgICAgIGxldCBuZXdfY2x1c3RlciA9IG5ldyBDbHVzdGVyKHRoaXMuX2lkKyssIGMxX2NsdXN0ZXIsIGMyX2NsdXN0ZXIsIEQuZW50cnkoYzEsIGMyKSwgbnVsbCwgaW5kaWNlcyk7XG4gICAgICAgICAgICBjMV9jbHVzdGVyLnBhcmVudCA9IG5ld19jbHVzdGVyO1xuICAgICAgICAgICAgYzJfY2x1c3Rlci5wYXJlbnQgPSBuZXdfY2x1c3RlcjtcbiAgICAgICAgICAgIGNsdXN0ZXJzW2MxXS51bnNoaWZ0KG5ld19jbHVzdGVyKTtcbiAgICAgICAgICAgIGNfc2l6ZVtjMV0gKz0gY19zaXplW2MyXTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgbjsgKytqKSB7XG4gICAgICAgICAgICAgICAgY29uc3QgRF9jMV9qID0gRC5lbnRyeShjMSwgaik7XG4gICAgICAgICAgICAgICAgY29uc3QgRF9jMl9qID0gRC5lbnRyeShjMiwgaik7XG4gICAgICAgICAgICAgICAgbGV0IHZhbHVlO1xuICAgICAgICAgICAgICAgIHN3aXRjaCAobGlua2FnZSkge1xuICAgICAgICAgICAgICAgICAgICBjYXNlIFwic2luZ2xlXCI6XG4gICAgICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IE1hdGgubWluKERfYzFfaiwgRF9jMl9qKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICBjYXNlIFwiY29tcGxldGVcIjpcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhbHVlID0gTWF0aC5tYXgoRF9jMV9qLCBEX2MyX2opO1xuICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgIGNhc2UgXCJhdmVyYWdlXCI6XG4gICAgICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IChjX3NpemVbYzFdICogRF9jMV9qICsgY19zaXplW2MyXSAqIERfYzJfaikgLyAoY19zaXplW2MxXSArIGNfc2l6ZVtqXSk7XG4gICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgRC5zZXRfZW50cnkoaiwgYzEsIHZhbHVlKTtcbiAgICAgICAgICAgICAgICBELnNldF9lbnRyeShjMSwgaiwgdmFsdWUpO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBELnNldF9lbnRyeShjMSwgYzEsIEluZmluaXR5KTtcbiAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbjsgKytpKSB7XG4gICAgICAgICAgICAgICAgRC5zZXRfZW50cnkoaSwgYzIsIEluZmluaXR5KTtcbiAgICAgICAgICAgICAgICBELnNldF9lbnRyeShjMiwgaSwgSW5maW5pdHkpO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvKiBmb3IgKGxldCBqID0gMDsgaiA8IG47ICsraikge1xuICAgICAgICAgICAgICAgIGlmIChkX21pbltqXSA9PT0gYzIpIHtcbiAgICAgICAgICAgICAgICAgICAgZF9taW5bal0gPSBjMTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgaWYgKEQuZW50cnkoYzEsIGopIDwgRC5lbnRyeShjMSwgZF9taW5bYzFdKSkge1xuICAgICAgICAgICAgICAgICAgICBkX21pbltjMV0gPSBqO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0gKi9cbiAgICAgICAgICAgIHJvb3QgPSBuZXdfY2x1c3RlcjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gcm9vdDtcbiAgICB9XG59XG5cbmNsYXNzIENsdXN0ZXIge1xuICAgIGNvbnN0cnVjdG9yKGlkLCBsZWZ0LCByaWdodCwgZGlzdCwgY2VudHJvaWQsIGluZGV4LCBzaXplLCBkZXB0aCkge1xuICAgICAgICB0aGlzLmlkID0gaWQ7XG4gICAgICAgIHRoaXMubGVmdCA9IGxlZnQ7XG4gICAgICAgIHRoaXMucmlnaHQgPSByaWdodDtcbiAgICAgICAgdGhpcy5kaXN0ID0gZGlzdDtcbiAgICAgICAgdGhpcy5pbmRleCA9IGluZGV4O1xuICAgICAgICB0aGlzLnNpemUgPSBzaXplID8/IGxlZnQuc2l6ZSArIHJpZ2h0LnNpemU7XG4gICAgICAgIHRoaXMuZGVwdGggPSBkZXB0aCA/PyAxICsgTWF0aC5tYXgobGVmdC5kZXB0aCwgcmlnaHQuZGVwdGgpO1xuICAgICAgICB0aGlzLmNlbnRyb2lkID0gY2VudHJvaWQgPz8gdGhpcy5fY2FsY3VsYXRlX2NlbnRyb2lkKGxlZnQsIHJpZ2h0KTtcbiAgICAgICAgdGhpcy5wYXJlbnQgPSBudWxsO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICBfY2FsY3VsYXRlX2NlbnRyb2lkKGxlZnQsIHJpZ2h0KSB7XG4gICAgICAgIGNvbnN0IGxfc2l6ZSA9IGxlZnQuc2l6ZTtcbiAgICAgICAgY29uc3Qgcl9zaXplID0gcmlnaHQuc2l6ZTtcbiAgICAgICAgY29uc3QgbF9jZW50cm9pZCA9IGxlZnQuY2VudHJvaWQ7XG4gICAgICAgIGNvbnN0IHJfY2VudHJvaWQgPSByaWdodC5jZW50cm9pZDtcbiAgICAgICAgY29uc3Qgc2l6ZSA9IHRoaXMuc2l6ZTtcbiAgICAgICAgY29uc3QgbiA9IGxlZnQuY2VudHJvaWQubGVuZ3RoO1xuICAgICAgICBjb25zdCBuZXdfY2VudHJvaWQgPSBuZXcgRmxvYXQ2NEFycmF5KG4pO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgbmV3X2NlbnRyb2lkW2ldID0gKGxfc2l6ZSAqIGxfY2VudHJvaWRbaV0gKyByX3NpemUgKiByX2NlbnRyb2lkW2ldKSAvIHNpemU7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIG5ld19jZW50cm9pZDtcbiAgICB9XG5cbiAgICBnZXQgaXNMZWFmKCkge1xuICAgICAgICByZXR1cm4gdGhpcy5kZXB0aCA9PT0gMDtcbiAgICB9XG5cbiAgICBsZWF2ZXMoKSB7XG4gICAgICAgIGlmICh0aGlzLmlzTGVhZikgcmV0dXJuIFt0aGlzXTtcbiAgICAgICAgY29uc3QgbGVmdCA9IHRoaXMubGVmdDtcbiAgICAgICAgY29uc3QgcmlnaHQgPSB0aGlzLnJpZ2h0O1xuICAgICAgICByZXR1cm4gKGxlZnQuaXNMZWFmID8gW2xlZnRdIDogbGVmdC5sZWF2ZXMoKSkuY29uY2F0KHJpZ2h0LmlzTGVhZiA/IFtyaWdodF0gOiByaWdodC5sZWF2ZXMoKSk7XG4gICAgfVxuXG4gICAgZGVzY2VuZGFudHMoKSB7XG4gICAgICAgIGlmICh0aGlzLmlzTGVhZikgcmV0dXJuIFt0aGlzXTtcbiAgICAgICAgY29uc3QgbGVmdF9kZXNjZW5kYW50cyA9IHRoaXMubGVmdC5kZXNjZW5kYW50cygpO1xuICAgICAgICBjb25zdCByaWdodF9kZXNjZW5kYW50cyA9IHRoaXMucmlnaHQuZGVzY2VuZGFudHMoKTtcbiAgICAgICAgcmV0dXJuIGxlZnRfZGVzY2VuZGFudHMuY29uY2F0KHJpZ2h0X2Rlc2NlbmRhbnRzKS5jb25jYXQoW3RoaXNdKTtcbiAgICB9XG59XG4iLCJpbXBvcnQgeyBldWNsaWRlYW4gfSBmcm9tIFwiLi4vbWV0cmljcy9pbmRleC5qc1wiO1xuaW1wb3J0IHsgUmFuZG9taXplciB9IGZyb20gXCIuLi91dGlsL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBIZWFwIH0gZnJvbSBcIi4uL2RhdGFzdHJ1Y3R1cmUvaW5kZXguanNcIjtcbmltcG9ydCB7IGxpbnNwYWNlIH0gZnJvbSBcIi4uL21hdHJpeC9pbmRleC5qc1wiO1xuXG4vKipcbiAqIEBjbGFzc1xuICogQGFsaWFzIEtNZWFuc1xuICovXG5leHBvcnQgY2xhc3MgS01lYW5zIHtcbiAgICAvKipcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmNsdXN0ZXJpbmdcbiAgICAgKiBAYWxpYXMgS01lYW5zXG4gICAgICogQHRvZG8gbmVlZHMgcmVzdHJ1Y3R1cmluZy4gXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IG1hdHJpeCBcbiAgICAgKiBAcGFyYW0ge051bWJlcnN9IEsgXG4gICAgICogQHBhcmFtIHtGdW5jdGlvbn0gW21ldHJpYyA9IGV1Y2xpZGVhbl0gXG4gICAgICogQHBhcmFtIHtOdW1iZXJ9IFtzZWVkID0gMTk4N11cbiAgICAgKiBAcGFyYW0ge0Jvb2xlYW59IFtpbml0ID0gdHJ1ZV1cbiAgICAgKiBAcmV0dXJucyB7S01lYW5zfVxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKG1hdHJpeCwgSywgbWV0cmljID0gZXVjbGlkZWFuLCBzZWVkPTE5ODcsIGluaXQgPSB0cnVlKSB7XG4gICAgICAgIHRoaXMuX21ldHJpYyA9IG1ldHJpYztcbiAgICAgICAgdGhpcy5fbWF0cml4ID0gbWF0cml4O1xuICAgICAgICB0aGlzLl9LID0gSztcbiAgICAgICAgY29uc3QgW04sIERdID0gbWF0cml4LnNoYXBlO1xuICAgICAgICB0aGlzLl9OID0gTjtcbiAgICAgICAgdGhpcy5fRCA9IEQ7XG4gICAgICAgIGlmIChLID4gTikgSyA9IE47XG4gICAgICAgIHRoaXMuX3JhbmRvbWl6ZXIgPSBuZXcgUmFuZG9taXplcihzZWVkKTtcbiAgICAgICAgdGhpcy5fY2x1c3RlcnMgPSBuZXcgQXJyYXkoTikuZmlsbCh1bmRlZmluZWQpO1xuICAgICAgICB0aGlzLl9jbHVzdGVyX2NlbnRyb2lkcyA9IHRoaXMuX2dldF9yYW5kb21fY2VudHJvaWRzKEspO1xuICAgICAgICBpZiAoaW5pdCkgdGhpcy5pbml0KEssIHRoaXMuX2NsdXN0ZXJfY2VudHJvaWRzKTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHJldHVybnMge0FycmF5PEFycmF5Pn0gLSBBcnJheSBvZiBjbHVzdGVycyB3aXRoIHRoZSBpbmRpY2VzIG9mIHRoZSByb3dzIGluIGdpdmVuIHtAbGluayBtYXRyaXh9LiBcbiAgICAgKi9cbiAgICBnZXRfY2x1c3RlcnMoKSB7XG4gICAgICAgIGNvbnN0IEsgPSB0aGlzLl9LO1xuICAgICAgICBjb25zdCBjbHVzdGVycyA9IHRoaXMuX2NsdXN0ZXJzO1xuICAgICAgICBjb25zdCByZXN1bHQgPSBuZXcgQXJyYXkoSykuZmlsbCgpLm1hcCgoKSA9PiBuZXcgQXJyYXkoKSk7XG4gICAgICAgIGNsdXN0ZXJzLmZvckVhY2goKGMsIGkpID0+IHJlc3VsdFtjXS5wdXNoKGkpKTtcbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7QXJyYXl9IHBvaW50cyBcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBjYW5kaWRhdGVzIFxuICAgICAqL1xuICAgIF9mdXJ0aGVzdF9wb2ludChwb2ludHMsIGNhbmRpZGF0ZXMpIHtcbiAgICAgICAgY29uc3QgQSA9IHRoaXMuX21hdHJpeDtcbiAgICAgICAgY29uc3QgbWV0cmljID0gdGhpcy5fbWV0cmljO1xuICAgICAgICBsZXQgaSA9IHBvaW50cy5sZW5ndGg7XG4gICAgICAgIGxldCBIID0gSGVhcC5oZWFwaWZ5KFxuICAgICAgICAgICAgY2FuZGlkYXRlcywgXG4gICAgICAgICAgICAoZCkgPT4ge1xuICAgICAgICAgICAgICAgIGNvbnN0IEFkID0gQS5yb3coZClcbiAgICAgICAgICAgICAgICBsZXQgc3VtID0gMDtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IGk7ICsraikge1xuICAgICAgICAgICAgICAgICAgICBzdW0gKz0gbWV0cmljKEFkLCBwb2ludHNbal0pXG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHJldHVybiBzdW07XG4gICAgICAgICAgICB9LCBcbiAgICAgICAgICAgIFwibWF4XCJcbiAgICAgICAgKVxuICAgICAgICByZXR1cm4gSC5wb3AoKS5lbGVtZW50O1xuICAgIH1cblxuICAgIF9nZXRfcmFuZG9tX2NlbnRyb2lkcyhLKSB7XG4gICAgICAgIGNvbnN0IE4gPSB0aGlzLl9OO1xuICAgICAgICBjb25zdCByYW5kb21pemVyID0gdGhpcy5fcmFuZG9taXplcjtcbiAgICAgICAgY29uc3QgQSA9IHRoaXMuX21hdHJpeDtcbiAgICAgICAgY29uc3QgY2x1c3Rlcl9jZW50cm9pZHMgPSBuZXcgQXJyYXkoSykuZmlsbCgpXG4gICAgICAgIGNvbnN0IGluZGljZXMgPSBsaW5zcGFjZSgwLCBOIC0gMSk7XG4gICAgICAgIGNvbnN0IHJhbmRvbV9wb2ludCA9IHJhbmRvbWl6ZXIucmFuZG9tX2ludCAlIChOIC0gMSk7XG4gICAgICAgIGNsdXN0ZXJfY2VudHJvaWRzWzBdID0gQS5yb3cocmFuZG9tX3BvaW50KTtcbiAgICAgICAgY29uc3QgaW5pdF9wb2ludHMgPSBbcmFuZG9tX3BvaW50XTtcbiAgICAgICAgY29uc3Qgc2FtcGxlX3NpemUgPSBNYXRoLmZsb29yKChOIC0gSykgLyBLKTsvLyAvIEtcbiAgICAgICAgZm9yIChsZXQgaSA9IDE7IGkgPCBLOyArK2kpIHtcbiAgICAgICAgICAgIC8vIHNhbXBsaW5nICsga21lYW5zKysgaW1wcm92ZW1lbnQ/XG4gICAgICAgICAgICBjb25zdCBzYW1wbGUgPSByYW5kb21pemVyLmNob2ljZShpbmRpY2VzLmZpbHRlcihkID0+IGluaXRfcG9pbnRzLmluZGV4T2YoZCkgPT0gLTEpLCBzYW1wbGVfc2l6ZSk7XG4gICAgICAgICAgICBjb25zdCBmdXJ0aGVzdF9wb2ludCA9IHRoaXMuX2Z1cnRoZXN0X3BvaW50KGNsdXN0ZXJfY2VudHJvaWRzLnNsaWNlKDAsIGkpLCBzYW1wbGUpO1xuICAgICAgICAgICAgaW5pdF9wb2ludHMucHVzaChmdXJ0aGVzdF9wb2ludCk7XG4gICAgICAgICAgICBjbHVzdGVyX2NlbnRyb2lkc1tpXSA9IEEucm93KGZ1cnRoZXN0X3BvaW50KTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gY2x1c3Rlcl9jZW50cm9pZHM7XG4gICAgfVxuXG4gICAgX2l0ZXJhdGlvbihjbHVzdGVyX2NlbnRyb2lkcykge1xuICAgICAgICBjb25zdCBLID0gY2x1c3Rlcl9jZW50cm9pZHMubGVuZ3RoO1xuICAgICAgICBjb25zdCBOID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgRCA9IHRoaXMuX0Q7XG4gICAgICAgIGNvbnN0IEEgPSB0aGlzLl9tYXRyaXg7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgY2x1c3RlcnMgPSB0aGlzLl9jbHVzdGVycztcbiAgICAgICAgbGV0IGNsdXN0ZXJzX2NoYW5nZWQgPSBmYWxzZTtcbiAgICAgICAgLy8gZmluZCBuZWFyZXN0IGNsdXN0ZXIgY2VudHJvaWQuXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBBaSA9IEEucm93KGkpXG4gICAgICAgICAgICBsZXQgbWluX2Rpc3QgPSBJbmZpbml0eTtcbiAgICAgICAgICAgIGxldCBtaW5fY2x1c3RlciA9IG51bGw7XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IEs7ICsraikge1xuICAgICAgICAgICAgICAgIGxldCBkID0gbWV0cmljKGNsdXN0ZXJfY2VudHJvaWRzW2pdLCBBaSk7XG4gICAgICAgICAgICAgICAgaWYgKGQgPCBtaW5fZGlzdCkge1xuICAgICAgICAgICAgICAgICAgICBtaW5fZGlzdCA9IGQ7XG4gICAgICAgICAgICAgICAgICAgIG1pbl9jbHVzdGVyID0gajsgXG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKGNsdXN0ZXJzW2ldICE9PSBtaW5fY2x1c3Rlcikge1xuICAgICAgICAgICAgICAgIGNsdXN0ZXJzX2NoYW5nZWQgPSB0cnVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgY2x1c3RlcnNbaV0gPSBtaW5fY2x1c3RlcjtcbiAgICAgICAgfVxuICAgICAgICAvLyB1cGRhdGUgY2x1c3RlciBjZW50cm9pZFxuICAgICAgICAvLyByZXNldCBjbHVzdGVyIGNlbnRyb2lkcyB0byAwXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgSzsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBjZW50cm9pZCA9IGNsdXN0ZXJfY2VudHJvaWRzW2ldO1xuICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBEOyArK2opIHtcbiAgICAgICAgICAgICAgICBjZW50cm9pZFtqXSA9IDA7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgLy8gY29tcHV0ZSBjZW50cm9pZFxuICAgICAgICB0aGlzLl9jb21wdXRlX2NlbnRyb2lkKGNsdXN0ZXJfY2VudHJvaWRzKTtcblxuICAgICAgICByZXR1cm4geyAgIFxuICAgICAgICAgICAgXCJjbHVzdGVyc19jaGFuZ2VkXCI6IGNsdXN0ZXJzX2NoYW5nZWQsXG4gICAgICAgICAgICBcImNsdXN0ZXJfY2VudHJvaWRzXCI6IGNsdXN0ZXJfY2VudHJvaWRzXG4gICAgICAgIH07XG4gICAgfVxuXG4gICAgX2NvbXB1dGVfY2VudHJvaWQoY2x1c3Rlcl9jZW50cm9pZHMpIHtcbiAgICAgICAgY29uc3QgSyA9IGNsdXN0ZXJfY2VudHJvaWRzLmxlbmd0aDtcbiAgICAgICAgY29uc3QgTiA9IHRoaXMuX047XG4gICAgICAgIGNvbnN0IEQgPSB0aGlzLl9EO1xuICAgICAgICBjb25zdCBBID0gdGhpcy5fbWF0cml4O1xuICAgICAgICBjb25zdCBjbHVzdGVycyA9IHRoaXMuX2NsdXN0ZXJzO1xuICAgICAgICBjb25zdCBjbHVzdGVyX2NvdW50ZXIgPSBuZXcgQXJyYXkoSykuZmlsbCgwKTtcblxuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgY29uc3QgQWkgPSBBLnJvdyhpKTtcbiAgICAgICAgICAgIGNvbnN0IGNpID0gY2x1c3RlcnNbaV07XG4gICAgICAgICAgICBjbHVzdGVyX2NvdW50ZXJbY2ldKys7XG4gICAgICAgICAgICBjb25zdCBjZW50cm9pZCA9IGNsdXN0ZXJfY2VudHJvaWRzW2NpXTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgRDsgKytqKSB7XG4gICAgICAgICAgICAgICAgY2VudHJvaWRbal0gKz0gQWlbal07XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBLOyArK2kpIHtcbiAgICAgICAgICAgIGNvbnN0IG4gPSBjbHVzdGVyX2NvdW50ZXJbaV07XG4gICAgICAgICAgICBjbHVzdGVyX2NlbnRyb2lkc1tpXSA9IGNsdXN0ZXJfY2VudHJvaWRzW2ldLm1hcChjID0+IGMgLyBuKTtcbiAgICAgICAgfVxuICAgICAgICBcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyB7QGxpbmsgS30gY2x1c3RlcnMgb3V0IG9mIHRoZSB7QGxpbmsgbWF0cml4fS5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gSyAtIG51bWJlciBvZiBjbHVzdGVycy5cbiAgICAgKi9cbiAgICBpbml0KEssIGNsdXN0ZXJfY2VudHJvaWRzKSB7XG4gICAgICAgIGlmICghSykgSyA9IHRoaXMuX0s7XG4gICAgICAgIGlmICghY2x1c3Rlcl9jZW50cm9pZHMpIGNsdXN0ZXJfY2VudHJvaWRzID0gdGhpcy5fZ2V0X3JhbmRvbV9jZW50cm9pZHMoSyk7XG4gICAgICAgIGxldCBjbHVzdGVyc19jaGFuZ2VkID0gZmFsc2U7XG4gICAgICAgIGRvIHtcbiAgICAgICAgICAgIGNvbnN0IGl0ZXJhdGlvbl9yZXN1bHQgPSB0aGlzLl9pdGVyYXRpb24oY2x1c3Rlcl9jZW50cm9pZHMpXG4gICAgICAgICAgICBjbHVzdGVyX2NlbnRyb2lkcyA9IGl0ZXJhdGlvbl9yZXN1bHQuY2x1c3Rlcl9jZW50cm9pZHM7XG4gICAgICAgICAgICBjbHVzdGVyc19jaGFuZ2VkID0gaXRlcmF0aW9uX3Jlc3VsdC5jbHVzdGVyc19jaGFuZ2VkO1xuICAgICAgICB9IHdoaWxlIChjbHVzdGVyc19jaGFuZ2VkKVxuICAgIH1cbiAgICBcbn1cbiIsImltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBSYW5kb21pemVyIH0gZnJvbSBcIi4uL3V0aWwvaW5kZXguanNcIjtcbmltcG9ydCB7IGxpbnNwYWNlLCBNYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBtaW4gfSBmcm9tIFwiLi4vdXRpbC9pbmRleC5qc1wiO1xuLyoqXG4gKiBAY2xhc3NcbiAqIEBhbGlhcyBLTWVkb2lkc1xuICovXG5leHBvcnQgY2xhc3MgS01lZG9pZHMge1xuICAgIC8qKlxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6Y2x1c3RlcmluZ1xuICAgICAqIEBhbGlhcyBLTWVkb2lkc1xuICAgICAqIEB0b2RvIG5lZWRzIHJlc3RydWN0dXJpbmcuIFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBtYXRyaXggLSBkYXRhIG1hdHJpeFxuICAgICAqIEBwYXJhbSB7TnVtYmVyc30gSyAtIG51bWJlciBvZiBjbHVzdGVyc1xuICAgICAqIEBwYXJhbSB7bnVtYmVyfSBbbWF4X2l0ZXI9bnVsbF0gLSBtYXhpbXVtIG51bWJlciBvZiBpdGVyYXRpb25zLiBEZWZhdWx0IGlzIDEwICogTWF0aC5sb2cxMChOKVxuICAgICAqIEBwYXJhbSB7RnVuY3Rpb259IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gbWV0cmljIGRlZmluaW5nIHRoZSBkaXNzaW1pbGFyaXR5IFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbc2VlZCA9IDEyMTJdIC0gc2VlZCB2YWx1ZSBmb3IgcmFuZG9tIG51bWJlciBnZW5lcmF0b3JcbiAgICAgKiBAcmV0dXJucyB7S01lZG9pZHN9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9saW5rLnNwcmluZ2VyLmNvbS9jaGFwdGVyLzEwLjEwMDcvOTc4LTMtMDMwLTMyMDQ3LThfMTZ9IEZhc3RlciBrLU1lZG9pZHMgQ2x1c3RlcmluZzogSW1wcm92aW5nIHRoZSBQQU0sIENMQVJBLCBhbmQgQ0xBUkFOUyBBbGdvcml0aG1zXG4gICAgICovXG4gICAgY29uc3RydWN0b3IobWF0cml4LCBLLCBtYXhfaXRlcj1udWxsLCBtZXRyaWMgPSBldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICB0aGlzLl9tZXRyaWMgPSBtZXRyaWM7XG4gICAgICAgIHRoaXMuX21hdHJpeCA9IG1hdHJpeDtcbiAgICAgICAgdGhpcy5fQSA9IHRoaXMuX21hdHJpeC50bzJkQXJyYXk7XG4gICAgICAgIHRoaXMuX0sgPSBLO1xuICAgICAgICBjb25zdCBbTiwgRF0gPSBtYXRyaXguc2hhcGU7XG4gICAgICAgIHRoaXMuX04gPSBOO1xuICAgICAgICB0aGlzLl9EID0gRDtcbiAgICAgICAgdGhpcy5fbWF4X2l0ZXIgPSBtYXhfaXRlciB8fCAxMCAqIE1hdGgubG9nMTAoTikgXG4gICAgICAgIHRoaXMuX2Rpc3RhbmNlX21hdHJpeCA9IG5ldyBNYXRyaXgoTiwgTiwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgLyogZm9yIChsZXQgaSA9IDE7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSBpICsgMTsgaiA8IE47ICsraikge1xuICAgICAgICAgICAgICAgIGxldCBkaXN0ID0gbWV0cmljKHRoaXMuX0FbaV0sIHRoaXMuX0Fbal0pO1xuICAgICAgICAgICAgICAgIHRoaXMuX2Rpc3RhbmNlX21hdHJpeC5zZXRfZW50cnkoaSwgaiwgZGlzdCk7XG4gICAgICAgICAgICAgICAgdGhpcy5fZGlzdGFuY2VfbWF0cml4LnNldF9lbnRyeShqLCBpLCBkaXN0KVxuICAgICAgICAgICAgfVxuICAgICAgICB9ICovXG4gICAgICAgIGlmIChLID4gTikgSyA9IE47XG4gICAgICAgIHRoaXMuX3JhbmRvbWl6ZXIgPSBuZXcgUmFuZG9taXplcihzZWVkKTtcbiAgICAgICAgdGhpcy5fY2x1c3RlcnMgPSBuZXcgQXJyYXkoTikuZmlsbCh1bmRlZmluZWQpO1xuICAgICAgICB0aGlzLl9jbHVzdGVyX21lZG9pZHMgPSB0aGlzLl9nZXRfcmFuZG9tX21lZG9pZHMoSyk7XG4gICAgICAgIC8vaWYgKGluaXQpIHRoaXMuaW5pdChLLCB0aGlzLl9jbHVzdGVyX21lZG9pZHMpO1xuICAgICAgICB0aGlzLl9pc19pbml0aWFsaXplZCA9IGZhbHNlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcmV0dXJucyB7QXJyYXk8QXJyYXk+fSAtIEFycmF5IG9mIGNsdXN0ZXJzIHdpdGggdGhlIGluZGljZXMgb2YgdGhlIHJvd3MgaW4gZ2l2ZW4ge0BsaW5rIG1hdHJpeH0uIFxuICAgICAqL1xuICAgIGdldF9jbHVzdGVycygpIHtcbiAgICAgICAgY29uc3QgSyA9IHRoaXMuX0s7XG4gICAgICAgIGNvbnN0IEEgPSB0aGlzLl9BO1xuICAgICAgICBpZiAoIXRoaXMuX2lzX2luaXRpYWxpemVkKSB7XG4gICAgICAgICAgICB0aGlzLmluaXQoSywgdGhpcy5fY2x1c3Rlcl9tZWRvaWRzKTtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCByZXN1bHQgPSBuZXcgQXJyYXkoSykuZmlsbCgpLm1hcCgoKSA9PiBuZXcgQXJyYXkoKSk7XG4gICAgICAgIEEuZm9yRWFjaCgoeF9qLCBqKSA9PiB7XG4gICAgICAgICAgICByZXN1bHRbdGhpcy5fbmVhcmVzdF9tZWRvaWQoeF9qLCBqKS5pbmRleF9uZWFyZXN0XS5wdXNoKGopO1xuICAgICAgICB9KVxuICAgICAgICByZXN1bHQubWVkb2lkcyA9IHRoaXMuX2NsdXN0ZXJfbWVkb2lkcztcbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG5cbiAgICBhc3luYyogZ2VuZXJhdG9yKCkge1xuICAgICAgICBjb25zdCBtYXhfaXRlciA9IHRoaXMuX21heF9pdGVyO1xuICAgICAgICB5aWVsZCB0aGlzLmdldF9jbHVzdGVycygpXG4gICAgICAgIGxldCBmaW5pc2ggPSBmYWxzZTtcbiAgICAgICAgbGV0IGkgPSAwXG4gICAgICAgIGRvIHtcbiAgICAgICAgICAgIGZpbmlzaCA9IHRoaXMuX2l0ZXJhdGlvbigpO1xuICAgICAgICAgICAgeWllbGQgdGhpcy5nZXRfY2x1c3RlcnMoKTtcbiAgICAgICAgfSB3aGlsZSAoIWZpbmlzaCAmJiArK2kgPCBtYXhfaXRlcilcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBBbGdvcml0aG0gMS4gRmFzdFBBTTE6IEltcHJvdmVkIFNXQVAgYWxnb3JpdGhtXG4gICAgICovXG4gICAgLyogX2l0ZXJhdGlvbl8xKCkge1xuICAgICAgICBjb25zdCBBID0gdGhpcy5fQTtcbiAgICAgICAgY29uc3QgTiA9IHRoaXMuX047XG4gICAgICAgIGNvbnN0IEsgPSB0aGlzLl9LO1xuICAgICAgICBjb25zdCBtZWRvaWRzID0gdGhpcy5fY2x1c3Rlcl9tZWRvaWRzO1xuICAgICAgICBsZXQgRGVsdGFURCA9IDA7XG4gICAgICAgIGxldCBtMCA9IG51bGw7XG4gICAgICAgIGxldCB4MCA9IG51bGw7XG4gICAgICAgIEEuZm9yRWFjaCgoeF9qLCBqKSA9PiB7XG4gICAgICAgICAgICBpZiAobWVkb2lkcy5maW5kSW5kZXgobSA9PiBtID09PSBqKSA8IDApIHtcbiAgICAgICAgICAgICAgICBjb25zdCBuZWFyZXN0X21lZG9pZCA9IHRoaXMuX25lYXJlc3RfbWVkb2lkKHhfaiwgaik7XG4gICAgICAgICAgICAgICAgY29uc3QgZF9qID0gbmVhcmVzdF9tZWRvaWQuZGlzdGFuY2VfbmVhcmVzdDsgLy8gZGlzdGFuY2UgdG8gY3VycmVudCBtZWRvaWRcbiAgICAgICAgICAgICAgICBjb25zdCBkZWx0YVREID0gbmV3IEFycmF5KEspLmZpbGwoLWRfaik7IC8vIGNoYW5nZSBpZiBtYWtpbmcgaiBhIG1lZG9pZFxuICAgICAgICAgICAgICAgIEEuZm9yRWFjaCgoeF9vLCBvKSA9PiB7XG4gICAgICAgICAgICAgICAgICAgIC8vIGRpc2FuY2UgdG8gbmV3IG1lZG9pZFxuICAgICAgICAgICAgICAgICAgICBjb25zdCBkX29qID0gdGhpcy5fZ2V0X2Rpc3RhbmNlKG8sIGosIHhfbywgeF9qKTtcbiAgICAgICAgICAgICAgICAgICAgY29uc3Qge1xuICAgICAgICAgICAgICAgICAgICAgICAgXCJpbmRleF9uZWFyZXN0XCI6IG4sXG4gICAgICAgICAgICAgICAgICAgICAgICBcImRpc3RhbmNlX25lYXJlc3RcIjogZF9uLFxuICAgICAgICAgICAgICAgICAgICAgICAgXCJkaXN0YW5jZV9zZWNvbmRcIjogZF9zLFxuICAgICAgICAgICAgICAgICAgICB9ID0gdGhpcy5fbmVhcmVzdF9tZWRvaWQoeF9vLCBvKTsgXG4gICAgICAgICAgICAgICAgICAgIHRoaXMuX2NsdXN0ZXJzW29dID0gbjsgLy8gY2FjaGVkIHZhbHVlc1xuICAgICAgICAgICAgICAgICAgICBkZWx0YVREW25dICs9IE1hdGgubWluKGRfb2osIGRfcykgLSBkX247IC8vIGxvc3MgY2hhbmdlXG4gICAgICAgICAgICAgICAgICAgIGlmIChkX29qIDwgZF9uKSB7IC8vIHJlYXNzaWdubWVudCBjaGVja1xuICAgICAgICAgICAgICAgICAgICAgICAgZGVsdGFURC5mb3JFYWNoKChkX2ksIGkpID0+IHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAobiAhPT0gaSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBkZWx0YVREW2ldID0gZF9pICsgZF9vaiAtIGRfbjsgLy8gdXBkYXRlIGxvc3MgY2hhbmdlXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgfSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9KTtcbiAgICAgICAgICAgICAgICAvLyBjaG9vc2UgYmVzdCBtZWRvaWQgaTtcbiAgICAgICAgICAgICAgICBjb25zdCBpID0gZGVsdGFURFxuICAgICAgICAgICAgICAgICAgICAubWFwKChkLCBpKSA9PiBbZCwgaV0pXG4gICAgICAgICAgICAgICAgICAgIC5zb3J0KChkMSwgZDIpID0+IGQxWzBdIC0gZDJbMF0pWzBdWzFdO1xuICAgICAgICAgICAgICAgIGNvbnN0IGRlbHRhVERfaSA9IGRlbHRhVERbaV07XG4gICAgICAgICAgICAgICAgLy8gc3RvcmVcbiAgICAgICAgICAgICAgICBpZiAoZGVsdGFURF9pIDwgRGVsdGFURCkge1xuICAgICAgICAgICAgICAgICAgICBEZWx0YVREID0gZGVsdGFURF9pO1xuICAgICAgICAgICAgICAgICAgICBtMCA9IGk7XG4gICAgICAgICAgICAgICAgICAgIHgwID0gajtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pO1xuXG4gICAgICAgIGlmIChEZWx0YVREID49IDApIHtcbiAgICAgICAgICAgIHJldHVybiB0cnVlIC8vIGJyZWFrIGxvb3AgaWYgRGVsdGFURCA+PSAwXG4gICAgICAgIH1cbiAgICAgICAgLy8gc3dhcCByb2xlcyBvZiBtZWRvaWQgbSBhbmQgbm9uLW1lZG9pZCB4O1xuICAgICAgICBtZWRvaWRzW20wXSA9IHgwO1xuICAgICAgICB0aGlzLl9jbHVzdGVyX21lZG9pZHMgPSBtZWRvaWRzO1xuICAgICAgICByZXR1cm4gZmFsc2VcbiAgICB9ICovXG5cbiAgICAvKiogQWxnb3JpdGhtIDIuIEZhc3RQQU0yOiBTV0FQIHdpdGggbXVsdGlwbGUgY2FuZGlkYXRlc1xuICAgICAqIFxuICAgICAqL1xuICAgIF9pdGVyYXRpb24oKSB7XG4gICAgICAgIGNvbnN0IEEgPSB0aGlzLl9BO1xuICAgICAgICBjb25zdCBLID0gdGhpcy5fSztcbiAgICAgICAgY29uc3QgbWVkb2lkcyA9IHRoaXMuX2NsdXN0ZXJfbWVkb2lkcztcbiAgICAgICAgY29uc3QgY2FjaGUgPSBBLm1hcCgoeF9vLCBvKSA9PiB0aGlzLl9uZWFyZXN0X21lZG9pZCh4X28sIG8pKTtcbiAgICAgICAgLy8gZW1wdHkgYmVzdCBjYW5kaWRhdGVzIGFycmF5XG4gICAgICAgIGNvbnN0IERlbHRhVEQgPSBuZXcgQXJyYXkoSykuZmlsbCgwKTtcbiAgICAgICAgY29uc3QgeHMgPSBuZXcgQXJyYXkoSykuZmlsbChudWxsKTtcbiAgICAgICAgQS5mb3JFYWNoKCh4X2osIGopID0+IHtcbiAgICAgICAgICAgIGlmIChtZWRvaWRzLmZpbmRJbmRleChtID0+IG0gPT09IGopIDwgMCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IGRfaiA9IGNhY2hlW2pdLmRpc3RhbmNlX25lYXJlc3Q7IC8vIGRpc3RhbmNlIHRvIGN1cnJlbnQgbWVkb2lkXG4gICAgICAgICAgICAgICAgY29uc3QgZGVsdGFURCA9IG5ldyBBcnJheShLKS5maWxsKC1kX2opOyAvLyBjaGFuZ2UgaWYgbWFraW5nIGogYSBtZWRvaWRcbiAgICAgICAgICAgICAgICBBLmZvckVhY2goKHhfbywgbykgPT4ge1xuICAgICAgICAgICAgICAgICAgICBpZiAoaiA9PT0gbykgcmV0dXJuO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBkX29qID0gdGhpcy5fZ2V0X2Rpc3RhbmNlKG8sIGosIHhfbywgeF9qKTsgLy8gZGlzdGFuY2UgdG8gbmV3IG1lZG9pZFxuICAgICAgICAgICAgICAgICAgICBjb25zdCB7XCJpbmRleF9uZWFyZXN0XCI6IG4sIFwiZGlzdGFuY2VfbmVhcmVzdFwiOiBkX24sIFwiZGlzdGFuY2Vfc2Vjb25kXCI6IGRfc30gPSBjYWNoZVtvXTsgLy8gY2FjaGVkXG4gICAgICAgICAgICAgICAgICAgIGRlbHRhVERbbl0gKz0gTWF0aC5taW4oZF9vaiwgZF9zKSAtIGRfbjsgLy8gbG9zcyBjaGFuZ2UgZm9yIHhfb1xuICAgICAgICAgICAgICAgICAgICAvLyBSZWFzc2lnbm1lbnQgY2hlY2tcbiAgICAgICAgICAgICAgICAgICAgaWYgKGRfb2ogPCBkX24pIHsgXG4gICAgICAgICAgICAgICAgICAgICAgICAvLyB1cGRhdGUgbG9zcyBjaGFuZ2VcbiAgICAgICAgICAgICAgICAgICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgSzsgKytpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGkgIT09IG4pIGRlbHRhVERbaV0gKz0gZF9vaiAtIGRfbjtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH0pO1xuICAgICAgICAgICAgICAgIC8vIHJlbWVtYmVyIGJlc3Qgc3dhcCBmb3IgaTtcbiAgICAgICAgICAgICAgICBkZWx0YVREXG4gICAgICAgICAgICAgICAgICAgIC5tYXAoKGQsIGkpID0+IFtkLCBpXSlcbiAgICAgICAgICAgICAgICAgICAgLmZpbHRlcigoW2QsIGldKSA9PiBkIDwgRGVsdGFURFtpXSlcbiAgICAgICAgICAgICAgICAgICAgLmZvckVhY2goKFtkLCBpXSkgPT4ge1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGQgPCBEZWx0YVREW2ldKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgRGVsdGFURFtpXSA9IGQ7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgeHNbaV0gPSBqO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9KVxuICAgICAgICAgICAgfVxuICAgICAgICB9KVxuICAgICAgICAvLyBzdG9wIGlmIG5vIGltcHJvdmVtZW50cyB3ZXJlIGZvdW5kXG4gICAgICAgIGlmIChtaW4oRGVsdGFURCkgPj0gMCkgcmV0dXJuIHRydWU7IFxuXG4gICAgICAgIC8vIGV4ZWN1dGUgYWxsIGltcHJvdmVtZW50c1xuICAgICAgICB3aGlsZSAobWluKERlbHRhVEQpIDwgMCkge1xuICAgICAgICAgICAgLy8gc3dhcCByb2xlcyBvZiBtZWRvaWQgbV9pIGFuZCBub25fbWVkb2lkIHhzX2lcbiAgICAgICAgICAgIGNvbnN0IGkgPSBEZWx0YVREXG4gICAgICAgICAgICAgICAgLm1hcCgoZCwgaSkgPT4gW2QsIGldKVxuICAgICAgICAgICAgICAgIC5zb3J0KChbYV0sIFtiXSkgPT4gYSAtIGIpWzBdWzFdO1xuICAgICAgICAgICAgaWYgKG1lZG9pZHMuZmlsdGVyKG0gPT4gbSA9PSB4c1tpXSkubGVuZ3RoID09IDApIHtcbiAgICAgICAgICAgICAgICBtZWRvaWRzW2ldID0geHNbaV07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICAvLyBkaXNhYmxlIHRoZSBzd2FwIGp1c3QgcGVyZm9ybWVkXG4gICAgICAgICAgICBEZWx0YVREW2ldID0gMDsgXG4gICAgICAgICAgICAvLyByZWNvbXB1dGUgVEQgZm9yIHJlbWFpbmluZyBzd2FwIGNhbmRpZGF0ZXNcbiAgICAgICAgICAgIERlbHRhVERcbiAgICAgICAgICAgICAgICAubWFwKChkX2osIGopID0+IFtkX2osIGpdKVxuICAgICAgICAgICAgICAgIC5maWx0ZXIoKFtkX2pdKSA9PiBkX2ogPCAwKVxuICAgICAgICAgICAgICAgIC5mb3JFYWNoKChbXywgal0pID0+IHtcbiAgICAgICAgICAgICAgICAgICAgY29uc3QgeF9qID0gQVtqXTtcbiAgICAgICAgICAgICAgICAgICAgbGV0IHN1bSA9IDA7XG4gICAgICAgICAgICAgICAgICAgIEEuZm9yRWFjaCgoeF9vLCBvKSA9PiB7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAobWVkb2lkcy5maW5kSW5kZXgobSA9PiBtICE9IGogJiYgbSA9PSBvKSA+PSAwKSByZXR1cm47XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoaSA9PSBqKSByZXR1cm47XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoY2FjaGVbb10uaW5kZXhfbmVhcmVzdCA9PT0gbWVkb2lkc1tqXSlcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdW0gKz0gKE1hdGgubWluKHRoaXMuX2dldF9kaXN0YW5jZShvLCBqLCB4X28sIHhfaiksIGNhY2hlW29dLmRpc3RhbmNlX3NlY29uZCkgLSBjYWNoZVtvXS5kaXN0YW5jZV9uZWFyZXN0KTsgXG4gICAgICAgICAgICAgICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdW0gKz0gKE1hdGgubWluKHRoaXMuX2dldF9kaXN0YW5jZShvLCBqLCB4X28sIHhfaikgLSBjYWNoZVtvXS5kaXN0YW5jZV9uZWFyZXN0LCAwKSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH0pO1xuICAgICAgICAgICAgICAgICAgICBEZWx0YVREW2pdID0gc3VtO1xuICAgICAgICAgICAgICAgIH0pXG4gICAgICAgIH1cbiAgICAgICAgdGhpcy5fY2x1c3Rlcl9tZWRvaWRzID0gbWVkb2lkcztcbiAgICAgICAgcmV0dXJuIGZhbHNlO1xuICAgIH1cblxuICAgIF9nZXRfZGlzdGFuY2UoaSwgaiwgeF9pPW51bGwsIHhfaj1udWxsKSB7XG4gICAgICAgIGlmIChpID09PSBqKSByZXR1cm4gMDtcbiAgICAgICAgY29uc3QgRCA9IHRoaXMuX2Rpc3RhbmNlX21hdHJpeDtcbiAgICAgICAgY29uc3QgQSA9IHRoaXMuX0E7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgbGV0IGRfaWogPSBELmVudHJ5KGksIGopO1xuICAgICAgICBpZiAoZF9paiA9PT0gMCkge1xuICAgICAgICAgICAgZF9paiA9IG1ldHJpYyh4X2kgfHwgQVtpXSwgeF9qIHx8IEFbal0pO1xuICAgICAgICAgICAgRC5zZXRfZW50cnkoaSwgaiwgZF9paik7XG4gICAgICAgICAgICBELnNldF9lbnRyeShqLCBpLCBkX2lqKTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gZF9pajtcbiAgICB9XG5cbiAgICBfbmVhcmVzdF9tZWRvaWQoeF9qLCBqKSB7XG4gICAgICAgIGNvbnN0IG1lZG9pZHMgPSB0aGlzLl9jbHVzdGVyX21lZG9pZHM7XG4gICAgICAgIGNvbnN0IEEgPSB0aGlzLl9BO1xuICAgICAgICBjb25zdCBbbmVhcmVzdCwgc2Vjb25kXSA9IG1lZG9pZHNcbiAgICAgICAgICAgIC5tYXAoKG0sIGkpID0+IHtcbiAgICAgICAgICAgICAgICBjb25zdCB4X20gPSBBW21dOyBcbiAgICAgICAgICAgICAgICByZXR1cm4gW3RoaXMuX2dldF9kaXN0YW5jZShqLCBtLCB4X2osIHhfbSksIGldO1xuICAgICAgICAgICAgfSlcbiAgICAgICAgICAgIC5zb3J0KChtMSwgbTIpID0+IG0xWzBdIC0gbTJbMF0pO1xuICAgICAgICBcbiAgICAgICAgcmV0dXJuIHsgXG4gICAgICAgICAgICBcImRpc3RhbmNlX25lYXJlc3RcIjogbmVhcmVzdFswXSwgXG4gICAgICAgICAgICBcImluZGV4X25lYXJlc3RcIjogbmVhcmVzdFsxXSxcbiAgICAgICAgICAgIFwiZGlzdGFuY2Vfc2Vjb25kXCI6IHNlY29uZFswXSxcbiAgICAgICAgICAgIFwiaW5kZXhfc2Vjb25kXCI6IHNlY29uZFsxXSxcbiAgICAgICAgfTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBDb21wdXRlcyB7QGxpbmsgS30gY2x1c3RlcnMgb3V0IG9mIHRoZSB7QGxpbmsgbWF0cml4fS5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gSyAtIG51bWJlciBvZiBjbHVzdGVycy5cbiAgICAgKi9cbiAgICBpbml0KEssIGNsdXN0ZXJfbWVkb2lkcykge1xuICAgICAgICBpZiAoIUspIEsgPSB0aGlzLl9LO1xuICAgICAgICBpZiAoIWNsdXN0ZXJfbWVkb2lkcykgY2x1c3Rlcl9tZWRvaWRzID0gdGhpcy5fZ2V0X3JhbmRvbV9tZWRvaWRzKEspO1xuICAgICAgICBjb25zdCBtYXhfaXRlciA9IHRoaXMuX21heF9pdGVyO1xuICAgICAgICBsZXQgZmluaXNoID0gZmFsc2U7XG4gICAgICAgIGxldCBpID0gMFxuICAgICAgICBkbyB7XG4gICAgICAgICAgICBmaW5pc2ggPSB0aGlzLl9pdGVyYXRpb24oKTtcbiAgICAgICAgfSB3aGlsZSAoIWZpbmlzaCAmJiArK2kgPCBtYXhfaXRlcilcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQWxnb3JpdGhtIDMuIEZhc3RQQU0gTEFCOiBMaW5lYXIgQXBwcm94aW1hdGUgQlVJTEQgaW5pdGlhbGl6YXRpb24uXG4gICAgICogQHBhcmFtIHtudW1iZXJ9IEsgLSBudW1iZXIgb2YgY2x1c3RlcnNcbiAgICAgKiBcbiAgICAgKi9cbiAgICBfZ2V0X3JhbmRvbV9tZWRvaWRzKEspIHtcbiAgICAgICAgY29uc3QgTiA9IHRoaXMuX047XG4gICAgICAgIGNvbnN0IEEgPSB0aGlzLl9BO1xuICAgICAgICBjb25zdCBpbmRpY2VzID0gbGluc3BhY2UoMCwgTiAtIDEpO1xuICAgICAgICBjb25zdCByYW5kb21pemVyID0gdGhpcy5fcmFuZG9taXplcjtcbiAgICAgICAgY29uc3QgbiA9IE1hdGgubWluKE4sIDEwICsgTWF0aC5jZWlsKE1hdGguc3FydChOKSkpO1xuICAgICAgICBjb25zdCBURCA9IG5ldyBBcnJheShuKS5maWxsKEluZmluaXR5KTtcbiAgICAgICAgY29uc3QgbWVkb2lkcyA9IFtdO1xuICAgICAgICAvLyBmaXJzdCBtZWRvaWRcbiAgICAgICAgbGV0IFREMCA9IEluZmluaXR5O1xuICAgICAgICBsZXQgUyA9IHJhbmRvbWl6ZXIuY2hvaWNlKGluZGljZXMsIG4pO1xuICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IG47ICsraikge1xuICAgICAgICAgICAgY29uc3QgU19qID0gU1tqXTtcbiAgICAgICAgICAgIGNvbnN0IHhfaiA9IEFbU19qXTtcbiAgICAgICAgICAgIGZvciAobGV0IG8gPSAwOyBvIDwgbjsgKytvKSB7XG4gICAgICAgICAgICAgICAgaWYgKG8gPT09IGopIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIGNvbnN0IHhfbyA9IEFbU1tvXV07XG4gICAgICAgICAgICAgICAgVERbal0gKz0gdGhpcy5fZ2V0X2Rpc3RhbmNlKGosIG8sIHhfaiwgeF9vKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChURFtqXSA8IFREMCkge1xuICAgICAgICAgICAgICAgIFREMCA9IFREW2pdOyAvLyBzbWFsbGVzdCBkaXN0YW5jZSBzdW1cbiAgICAgICAgICAgICAgICBtZWRvaWRzLnB1c2goU19qKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICAvLyBvdGhlciBtZWRvaWRzXG4gICAgICAgIGZvciAobGV0IGkgPSAxOyBpIDwgSzsgKytpKSB7XG4gICAgICAgICAgICBsZXQgRGVsdGFURCA9IEluZmluaXR5O1xuICAgICAgICAgICAgUyA9IHJhbmRvbWl6ZXIuY2hvaWNlKGluZGljZXMuZmlsdGVyKGluZGV4ID0+IG1lZG9pZHMuZmluZEluZGV4KGQgPT4gZCA9PT0gaW5kZXgpIDwgMCksIG4pO1xuICAgICAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBuOyArK2opIHtcbiAgICAgICAgICAgICAgICBsZXQgZGVsdGFURCA9IDA7XG4gICAgICAgICAgICAgICAgY29uc3QgU19qID0gU1tqXTtcbiAgICAgICAgICAgICAgICBjb25zdCB4X2ogPSBBW1Nfal07XG4gICAgICAgICAgICAgICAgZm9yIChsZXQgbyA9IDA7IG8gPCBuOyArK28pIHtcbiAgICAgICAgICAgICAgICAgICAgaWYgKG8gPT09IGopIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCBTX28gPSBTW29dO1xuICAgICAgICAgICAgICAgICAgICBjb25zdCB4X28gPSBBW1Nfb107XG4gICAgICAgICAgICAgICAgICAgIGxldCBkZWx0YSA9IHRoaXMuX2dldF9kaXN0YW5jZShTX2osIFNfbywgeF9qLCB4X28pIC0gbWluKG1lZG9pZHMubWFwKG0gPT4gdGhpcy5fZ2V0X2Rpc3RhbmNlKFNfbywgbSwgeF9vKSkpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoZGVsdGEgPCAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBkZWx0YVREID0gZGVsdGFURCArIGRlbHRhO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIC8vIGJlc3QgcmVkdWN0aW9uXG4gICAgICAgICAgICAgICAgaWYgKGRlbHRhVEQgPCBEZWx0YVREKSB7XG4gICAgICAgICAgICAgICAgICAgIERlbHRhVEQgPSBkZWx0YVREO1xuICAgICAgICAgICAgICAgICAgICBtZWRvaWRzLnB1c2goU19qKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBURDAgKz0gRGVsdGFURDtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gbWVkb2lkcy5zbGljZSgwLCBLKTtcbiAgICB9XG4gICAgXG59XG4iLCJpbXBvcnQgeyBldWNsaWRlYW4gfSBmcm9tIFwiLi4vbWV0cmljcy9pbmRleC5qc1wiO1xuaW1wb3J0IHsgSGVhcCB9IGZyb20gXCIuLi9kYXRhc3RydWN0dXJlL2luZGV4LmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgT1BUSUNTXG4gKi9cbmV4cG9ydCBjbGFzcyBPUFRJQ1Mge1xuICAgIC8qKlxuICAgICAqICoqTyoqcmRlcmluZyAqKlAqKm9pbnRzICoqVCoqbyAqKkkqKmRlbnRpZnkgdGhlICoqQyoqbHVzdGVyaW5nICoqUyoqdHJ1Y3R1cmUuXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTpjbHVzdGVyaW5nXG4gICAgICogQGFsaWFzIE9QVElDU1xuICAgICAqIEB0b2RvIG5lZWRzIHJlc3RydWN0dXJpbmcuIFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBtYXRyaXggLSB0aGUgZGF0YS5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gZXBzaWxvbiAtIHRoZSBtaW5pbXVtIGRpc3RhbmNlIHdoaWNoIGRlZmluZXMgd2hldGhlciBhIHBvaW50IGlzIGEgbmVpZ2hib3Igb3Igbm90LlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBtaW5fcG9pbnRzIC0gdGhlIG1pbmltdW0gbnVtYmVyIG9mIHBvaW50cyB3aGljaCBhIHBvaW50IG5lZWRzIHRvIGNyZWF0ZSBhIGNsdXN0ZXIuIChTaG91bGQgYmUgaGlnaGVyIHRoYW4gMSwgZWxzZSBlYWNoIHBvaW50IGNyZWF0ZXMgYSBjbHVzdGVyLilcbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBkaXN0YW5jZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzIG9mIHRoZSB7QGxpbmsgbWF0cml4fS5cbiAgICAgKiBAcmV0dXJucyB7T1BUSUNTfVxuICAgICAqIEBzZWUge0BsaW5rIGh0dHBzOi8vd3d3LmRicy5pZmkubG11LmRlL1B1Ymxpa2F0aW9uZW4vUGFwZXJzL09QVElDUy5wZGZ9XG4gICAgICogQHNlZSB7QGxpbmsgaHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvT1BUSUNTX2FsZ29yaXRobX1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihtYXRyaXgsIGVwc2lsb24sIG1pbl9wb2ludHMsIG1ldHJpYyA9IGV1Y2xpZGVhbikge1xuICAgICAgICB0aGlzLl9tYXRyaXggPSBtYXRyaXg7XG4gICAgICAgIHRoaXMuX2Vwc2lsb24gPSBlcHNpbG9uO1xuICAgICAgICB0aGlzLl9taW5fcG9pbnRzID0gbWluX3BvaW50cztcbiAgICAgICAgdGhpcy5fbWV0cmljID0gbWV0cmljO1xuXG4gICAgICAgIHRoaXMuX29yZGVyZWRfbGlzdCA9IFtdO1xuICAgICAgICB0aGlzLl9jbHVzdGVycyA9IFtdO1xuICAgICAgICB0aGlzLl9EQiA9IG5ldyBBcnJheShtYXRyaXguc2hhcGVbMF0pLmZpbGwoKTtcbiAgICAgICAgdGhpcy5pbml0KCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENvbXB1dGVzIHRoZSBjbHVzdGVyaW5nLlxuICAgICAqL1xuICAgIGluaXQoKSB7XG4gICAgICAgIGNvbnN0IG9yZGVyZWRfbGlzdCA9IHRoaXMuX29yZGVyZWRfbGlzdDtcbiAgICAgICAgY29uc3QgbWF0cml4ID0gdGhpcy5fbWF0cml4O1xuICAgICAgICBjb25zdCBOID0gbWF0cml4LnNoYXBlWzBdO1xuICAgICAgICBjb25zdCBEQiA9IHRoaXMuX0RCO1xuICAgICAgICBjb25zdCBjbHVzdGVycyA9IHRoaXMuX2NsdXN0ZXJzO1xuICAgICAgICBsZXQgY2x1c3Rlcl9pbmRleCA9IHRoaXMuX2NsdXN0ZXJfaW5kZXggPSAwO1xuXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBEQltpXSA9IHtcbiAgICAgICAgICAgICAgICBcImVsZW1lbnRcIjogbWF0cml4LnJvdyhpKSxcbiAgICAgICAgICAgICAgICBcImluZGV4XCI6IGksXG4gICAgICAgICAgICAgICAgXCJyZWFjaGFiaWxpdHlfZGlzdGFuY2VcIjogdW5kZWZpbmVkLFxuICAgICAgICAgICAgICAgIFwicHJvY2Vzc2VkXCI6IGZhbHNlLFxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGZvciAoY29uc3QgcCBvZiBEQikge1xuICAgICAgICAgICAgaWYgKHAucHJvY2Vzc2VkKSBjb250aW51ZTtcbiAgICAgICAgICAgIHAubmVpZ2hib3JzID0gdGhpcy5fZ2V0X25laWdoYm9ycyhwKTtcbiAgICAgICAgICAgIHAucHJvY2Vzc2VkID0gdHJ1ZTtcbiAgICAgICAgICAgIGNsdXN0ZXJzLnB1c2goW3AuaW5kZXhdKVxuICAgICAgICAgICAgY2x1c3Rlcl9pbmRleCA9IGNsdXN0ZXJzLmxlbmd0aCAtIDE7XG4gICAgICAgICAgICBvcmRlcmVkX2xpc3QucHVzaChwKTtcbiAgICAgICAgICAgIGlmICh0aGlzLl9jb3JlX2Rpc3RhbmNlKHApICE9IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IHNlZWRzID0gbmV3IEhlYXAobnVsbCwgZCA9PiBkLnJlYWNoYWJpbGl0eV9kaXN0YW5jZSwgXCJtaW5cIilcbiAgICAgICAgICAgICAgICB0aGlzLl91cGRhdGUocCwgc2VlZHMpO1xuICAgICAgICAgICAgICAgIHRoaXMuX2V4cGFuZF9jbHVzdGVyKHNlZWRzLCBjbHVzdGVyc1tjbHVzdGVyX2luZGV4XSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge09iamVjdH0gcCAtIGEgcG9pbnQgb2Yge0BsaW5rIG1hdHJpeH0uXG4gICAgICogQHJldHVybnMge0FycmF5fSBBbiBhcnJheSBjb25zaXN0aW5nIG9mIHRoZSB7QGxpbmsgZXBzaWxvbn0tbmVpZ2hib3Job29kIG9mIHtAbGluayBwfS5cbiAgICAgKi9cbiAgICBfZ2V0X25laWdoYm9ycyhwKSB7XG4gICAgICAgIGlmIChcIm5laWdoYm9yc1wiIGluIHApIHJldHVybiBwLm5laWdoYm9ycztcbiAgICAgICAgY29uc3QgREIgPSB0aGlzLl9EQjtcbiAgICAgICAgY29uc3QgbWV0cmljID0gdGhpcy5fbWV0cmljO1xuICAgICAgICBjb25zdCBlcHNpbG9uID0gdGhpcy5fZXBzaWxvbjtcbiAgICAgICAgY29uc3QgbmVpZ2hib3JzID0gW107XG4gICAgICAgIGZvciAoY29uc3QgcSBvZiBEQikge1xuICAgICAgICAgICAgaWYgKHEuaW5kZXggPT0gcC5pbmRleCkgY29udGludWU7XG4gICAgICAgICAgICBpZiAobWV0cmljKHAuZWxlbWVudCwgcS5lbGVtZW50KSA8IGVwc2lsb24pIHtcbiAgICAgICAgICAgICAgICBuZWlnaGJvcnMucHVzaChxKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gbmVpZ2hib3JzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtPYmplY3R9IHAgLSBhIHBvaW50IG9mIHtAbGluayBtYXRyaXh9LlxuICAgICAqIEByZXR1cm5zIHtOdW1iZXJ9IFRoZSBkaXN0YW5jZSB0byB0aGUge0BsaW5rIG1pbl9wb2ludHN9LXRoIG5lYXJlc3QgcG9pbnQgb2Yge0BsaW5rIHB9LCBvciB1bmRlZmluZWQgaWYgdGhlIHtAbGluayBlcHNpbG9ufS1uZWlnaGJvcmhvb2QgaGFzIGZld2VyIGVsZW1lbnRzIHRoYW4ge0BsaW5rIG1pbl9wb2ludHN9LlxuICAgICAqL1xuICAgIF9jb3JlX2Rpc3RhbmNlKHApIHtcbiAgICAgICAgY29uc3QgbWluX3BvaW50cyA9IHRoaXMuX21pbl9wb2ludHM7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgaWYgKHAubmVpZ2hib3JzICYmIHAubmVpZ2hib3JzLmxlbmd0aCA8PSBtaW5fcG9pbnRzKSB7XG4gICAgICAgICAgICByZXR1cm4gdW5kZWZpbmVkO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBtZXRyaWMocC5lbGVtZW50LCBwLm5laWdoYm9yc1ttaW5fcG9pbnRzXS5lbGVtZW50KTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBVcGRhdGVzIHRoZSByZWFjaGFiaWxpdHkgZGlzdGFuY2Ugb2YgdGhlIHBvaW50cy5cbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7T2JqZWN0fSBwIFxuICAgICAqIEBwYXJhbSB7SGVhcH0gc2VlZHMgXG4gICAgICovXG4gICAgX3VwZGF0ZShwLCBzZWVkcykge1xuICAgICAgICBjb25zdCBtZXRyaWMgPSB0aGlzLl9tZXRyaWM7XG4gICAgICAgIGNvbnN0IGNvcmVfZGlzdGFuY2UgPSB0aGlzLl9jb3JlX2Rpc3RhbmNlKHApO1xuICAgICAgICBjb25zdCBuZWlnaGJvcnMgPSB0aGlzLl9nZXRfbmVpZ2hib3JzKHApOy8vcC5uZWlnaGJvcnM7XG4gICAgICAgIGZvciAoY29uc3QgcSBvZiBuZWlnaGJvcnMpIHtcbiAgICAgICAgICAgIGlmIChxLnByb2Nlc3NlZCkgY29udGludWU7XG4gICAgICAgICAgICBjb25zdCBuZXdfcmVhY2hhYmlsaXR5X2Rpc3RhbmNlID0gTWF0aC5tYXgoY29yZV9kaXN0YW5jZSwgbWV0cmljKHAuZWxlbWVudCwgcS5lbGVtZW50KSk7XG4gICAgICAgICAgICAvL2lmIChxLnJlYWNoYWJpbGl0eV9kaXN0YW5jZSA9PSB1bmRlZmluZWQpIHsgLy8gcSBpcyBub3QgaW4gc2VlZHNcbiAgICAgICAgICAgIGlmIChzZWVkcy5yYXdfZGF0YSgpLmZpbmRJbmRleChkID0+IGQuZWxlbWVudCA9PSBxKSA8IDApIHtcbiAgICAgICAgICAgICAgICBxLnJlYWNoYWJpbGl0eV9kaXN0YW5jZSA9IG5ld19yZWFjaGFiaWxpdHlfZGlzdGFuY2U7XG4gICAgICAgICAgICAgICAgc2VlZHMucHVzaChxKTtcbiAgICAgICAgICAgIH0gZWxzZSB7IC8vIHEgaXMgaW4gc2VlZHNcbiAgICAgICAgICAgICAgICBpZiAobmV3X3JlYWNoYWJpbGl0eV9kaXN0YW5jZSA8IHEucmVhY2hhYmlsaXR5X2Rpc3RhbmNlKSB7XG4gICAgICAgICAgICAgICAgICAgIHEucmVhY2hhYmlsaXR5X2Rpc3RhbmNlID0gbmV3X3JlYWNoYWJpbGl0eV9kaXN0YW5jZTtcbiAgICAgICAgICAgICAgICAgICAgc2VlZHMgPSBIZWFwLmhlYXBpZnkoc2VlZHMuZGF0YSgpLCBkID0+IGQucmVhY2hhYmlsaXR5X2Rpc3RhbmNlLCBcIm1pblwiKTsgLy8gc2VlZHMgY2hhbmdlIGtleSA9L1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEV4cGFuZHMgdGhlIHtAbGluayBjbHVzdGVyfSB3aXRoIHBvaW50cyBpbiB7QGxpbmsgc2VlZHN9LlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtIZWFwfSBzZWVkcyBcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBjbHVzdGVyIFxuICAgICAqL1xuICAgIF9leHBhbmRfY2x1c3RlcihzZWVkcywgY2x1c3Rlcikge1xuICAgICAgICBjb25zdCBvcmRlcmVkX2xpc3QgPSB0aGlzLl9vcmRlcmVkX2xpc3Q7XG4gICAgICAgIHdoaWxlICghc2VlZHMuZW1wdHkpIHtcbiAgICAgICAgICAgIGNvbnN0IHEgPSBzZWVkcy5wb3AoKS5lbGVtZW50O1xuICAgICAgICAgICAgcS5uZWlnaGJvcnMgPSB0aGlzLl9nZXRfbmVpZ2hib3JzKHEpO1xuICAgICAgICAgICAgcS5wcm9jZXNzZWQgPSB0cnVlO1xuICAgICAgICAgICAgY2x1c3Rlci5wdXNoKHEuaW5kZXgpO1xuICAgICAgICAgICAgb3JkZXJlZF9saXN0LnB1c2gocSk7XG4gICAgICAgICAgICBpZiAodGhpcy5fY29yZV9kaXN0YW5jZShxKSAhPSB1bmRlZmluZWQpIHtcbiAgICAgICAgICAgICAgICB0aGlzLl91cGRhdGUocSwgc2VlZHMpO1xuICAgICAgICAgICAgICAgIHRoaXMuX2V4cGFuZF9jbHVzdGVyKHNlZWRzLCBjbHVzdGVyKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFJldHVybnMgYW4gYXJyYXkgb2YgY2x1c3RlcnMuXG4gICAgICogQHJldHVybnMge0FycmF5PEFycmF5Pn0gQXJyYXkgb2YgY2x1c3RlcnMgd2l0aCB0aGUgaW5kaWNlcyBvZiB0aGUgcm93cyBpbiBnaXZlbiB7QGxpbmsgbWF0cml4fS5cbiAgICAgKi9cbiAgICBnZXRfY2x1c3RlcnMoKSB7XG4gICAgICAgIGNvbnN0IGNsdXN0ZXJzID0gW107XG4gICAgICAgIGNvbnN0IG91dGxpZXJzID0gW107XG4gICAgICAgIGNvbnN0IG1pbl9wb2ludHMgPSB0aGlzLl9taW5fcG9pbnRzO1xuICAgICAgICBmb3IgKGNvbnN0IGNsdXN0ZXIgb2YgdGhpcy5fY2x1c3RlcnMpIHtcbiAgICAgICAgICAgIGlmIChjbHVzdGVyLmxlbmd0aCA8IG1pbl9wb2ludHMpIHtcbiAgICAgICAgICAgICAgICBvdXRsaWVycy5wdXNoKC4uLmNsdXN0ZXIpO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBjbHVzdGVycy5wdXNoKGNsdXN0ZXIpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGNsdXN0ZXJzLnB1c2gob3V0bGllcnMpO1xuICAgICAgICByZXR1cm4gY2x1c3RlcnM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHJldHVybnMge0FycmF5fSBSZXR1cm5zIGFuIGFycmF5LCB3aGVyZSB0aGUgaXRoIGVudHJ5IGRlZmluZXMgdGhlIGNsdXN0ZXIgYWZmaXJtYXRpb24gb2YgdGhlIGl0aCBwb2ludCBvZiB7QGxpbmsgbWF0cml4fS4gKC0xIHN0YW5kcyBmb3Igb3V0bGllcilcbiAgICAgKi9cbiAgICBnZXRfY2x1c3Rlcl9hZmZpcm1hdGlvbigpIHtcbiAgICAgICAgY29uc3QgTiA9IHRoaXMuX21hdHJpeC5zaGFwZVswXTtcbiAgICAgICAgY29uc3QgcmVzdWx0ID0gbmV3IEFycmF5KE4pLmZpbGwoKTtcbiAgICAgICAgY29uc3QgY2x1c3RlcnMgPSB0aGlzLmdldF9jbHVzdGVycygpO1xuICAgICAgICBmb3IgKGxldCBpID0gMCwgbiA9IGNsdXN0ZXJzLmxlbmd0aDsgaSA8IG47ICsraSkge1xuICAgICAgICAgICAgY29uc3QgY2x1c3RlciA9IGNsdXN0ZXJzW2ldXG4gICAgICAgICAgICBmb3IgKGNvbnN0IGluZGV4IG9mIGNsdXN0ZXIpIHtcbiAgICAgICAgICAgICAgICByZXN1bHRbaW5kZXhdID0gKGkgPCBuIC0gMSkgPyBpIDogLTE7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG59XG4iLCJpbXBvcnQgeyBNYXRyaXggfSBmcm9tIFwiLi4vbWF0cml4L2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5pbXBvcnQgeyBNRFMgfSBmcm9tIFwiLi9NRFMuanNcIjtcbmltcG9ydCB7IEtNZWRvaWRzIH0gZnJvbSBcIi4uL2NsdXN0ZXJpbmcvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBCYWxsVHJlZSB9IGZyb20gXCIuLi9rbm4vaW5kZXguanNcIjtcbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgTFNQXG4gKiBAZXh0ZW5kcyBEUlxuICovXG5leHBvcnQgY2xhc3MgTFNQIGV4dGVuZHMgRFIge1xuICAgIC8qKlxuICAgICAqIFxuICAgICAqIEBjb25zdHJ1Y3RvclxuICAgICAqIEBtZW1iZXJvZiBtb2R1bGU6ZGltZW5zaW9uYWxpdHlfcmVkdWN0aW9uXG4gICAgICogQGFsaWFzIExTUFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gdGhlIGhpZ2gtZGltZW5zaW9uYWwgZGF0YS4gXG4gICAgICogQHBhcmFtIHtudW1iZXJ9IFtrID0gTWF0aC5tYXgoTWF0aC5mbG9vcihOIC8gMTApLCAyKV0gLSBudW1iZXIgb2YgbmVpZ2hib3JzIHRvIGNvbnNpZGVyLlxuICAgICAqIEBwYXJhbSB7bnVtYmVyfSBbY29udHJvbF9wb2ludHMgPSBNYXRoLmNlaWwoTWF0aC5zcXJ0KE4pKV0gLSBudW1iZXIgb2YgY29udHJvbHBvaW50c1xuICAgICAqIEBwYXJhbSB7bnVtYmVyfSBbZCA9IDJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEBwYXJhbSB7ZnVuY3Rpb259IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gdGhlIG1ldHJpYyB3aGljaCBkZWZpbmVzIHRoZSBkaXN0YW5jZSBiZXR3ZWVuIHR3byBwb2ludHMuICBcbiAgICAgKiBAcmV0dXJucyB7TFNQfVxuICAgICAqIEBzZWUge0BsaW5rIGh0dHBzOi8vaWVlZXhwbG9yZS5pZWVlLm9yZy9kb2N1bWVudC80Mzc4MzcwfVxuICAgICAqL1xuICAgIGNvbnN0cnVjdG9yKFgsIGssIGNvbnRyb2xfcG9pbnRzLCBkPTIsIG1ldHJpYz1ldWNsaWRlYW4sIHNlZWQ9MTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICBzdXBlci5wYXJhbWV0ZXJfbGlzdCA9IFtcImtcIiwgXCJjb250cm9sX3BvaW50c1wiXTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJrXCIsIE1hdGgubWluKGsgPz8gTWF0aC5tYXgoTWF0aC5mbG9vcih0aGlzLl9OIC8gMTApLCAyKSwgdGhpcy5fTiAtIDEpKTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJjb250cm9sX3BvaW50c1wiLCBNYXRoLm1pbihjb250cm9sX3BvaW50cyA/PyBNYXRoLmNlaWwoTWF0aC5zcXJ0KHRoaXMuX04pKSwgdGhpcy5fTiAtIDEpKTtcbiAgICAgICAgdGhpcy5faXNfaW5pdGlhbGl6ZWQgPSBmYWxzZTtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogXG4gICAgICogQHBhcmFtIHtEUn0gRFIgLSBtZXRob2QgdXNlZCBmb3IgcG9zaXRpb24gY29udHJvbCBwb2ludHMuXG4gICAgICogQHBhcmFtIHtEUl9wYXJhbWV0ZXJzfSBEUl9wYXJhbWV0ZXJzIC0gYXJyYXkgY29udGFpbmluZyBwYXJhbWV0ZXJzIGZvciB0aGUgRFIgbWV0aG9kIHdoaWNoIHByb2plY3RzIHRoZSBjb250cm9sIHBvaW50c1xuICAgICAqIEByZXR1cm5zIHtMU1B9IFxuICAgICAqL1xuICAgIGluaXQoRFI9TURTLCBEUl9wYXJhbWV0ZXJzPVtdLCBLTk49QmFsbFRyZWUpIHtcbiAgICAgICAgaWYgKHRoaXMuX2lzX2luaXRpYWxpemVkKSByZXR1cm4gdGhpcztcbiAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgY29uc3QgTiA9IHRoaXMuX047XG4gICAgICAgIGNvbnN0IEsgPSB0aGlzLnBhcmFtZXRlcihcImtcIik7XG4gICAgICAgIGNvbnN0IGQgPSB0aGlzLl9kO1xuICAgICAgICBjb25zdCBtZXRyaWMgPSB0aGlzLl9tZXRyaWM7XG4gICAgICAgIGNvbnN0IG5jID0gdGhpcy5wYXJhbWV0ZXIoXCJjb250cm9sX3BvaW50c1wiKTtcbiAgICAgICAgY29uc3QgY29udHJvbF9wb2ludHMgPSBuZXcgS01lZG9pZHMoWCwgbmMsIG51bGwsIG1ldHJpYykuZ2V0X2NsdXN0ZXJzKCkubWVkb2lkcztcbiAgICAgICAgY29uc3QgQyA9IG5ldyBNYXRyaXgobmMsIE4sIFwiemVyb3NcIilcbiAgICAgICAgY29udHJvbF9wb2ludHMuZm9yRWFjaCgoY19pLCBpKSA9PiB7XG4gICAgICAgICAgICBDLnNldF9lbnRyeShpLCBjX2ksIDEpO1xuICAgICAgICB9KVxuICAgICAgICBjb25zdCBZX0MgPSBuZXcgRFIoTWF0cml4LmZyb20oY29udHJvbF9wb2ludHMubWFwKGNfaSA9PiBYLnJvdyhjX2kpKSksIC4uLkRSX3BhcmFtZXRlcnMsIGQpLnRyYW5zZm9ybSgpO1xuICAgICAgICBcbiAgICAgICAgY29uc3QgWEEgPSBYLnRvMmRBcnJheTtcbiAgICAgICAgY29uc3Qga25uID0gbmV3IEtOTihYQSwgbWV0cmljKTtcbiAgICAgICAgY29uc3QgTCA9IG5ldyBNYXRyaXgoTiwgTiwgXCJJXCIpO1xuICAgICAgICBjb25zdCBhbHBoYSA9IC0xL0s7XG4gICAgICAgIFhBLmZvckVhY2goKHhfaSwgaSkgPT4ge1xuICAgICAgICAgICAgZm9yIChjb25zdCB7XCJpbmRleFwiOiBqfSBvZiBrbm4uc2VhcmNoKHhfaSwgSykuaXRlcmF0ZSgpKSB7XG4gICAgICAgICAgICAgICAgaWYgKGkgPT09IGopIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIEwuc2V0X2VudHJ5KGksIGosIGFscGhhKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSlcbiAgICAgICAgY29uc3QgQSA9IEwuY29uY2F0KEMsIFwidmVydGljYWxcIik7XG5cbiAgICAgICAgY29uc3QgeiA9IG5ldyBNYXRyaXgoTiwgZCwgXCJ6ZXJvc1wiKTtcbiAgICAgICAgY29uc3QgYiA9IHouY29uY2F0KFlfQywgXCJ2ZXJ0aWNhbFwiKTtcbiAgICAgICAgXG4gICAgICAgIHRoaXMuX0EgPSBBO1xuICAgICAgICB0aGlzLl9iID0gYjtcbiAgICAgICAgdGhpcy5faXNfaW5pdGlhbGl6ZWQgPSB0cnVlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cblxuICAgIC8qKlxuICAgICAqIENvbXB1dGVzIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEByZXR1cm5zIHtNYXRyaXh9IFJldHVybnMgdGhlIHByb2plY3Rpb24uXG4gICAgICovXG4gICAgdHJhbnNmb3JtKCkge1xuICAgICAgICB0aGlzLmNoZWNrX2luaXQoKTtcbiAgICAgICAgY29uc3QgQSA9IHRoaXMuX0E7XG4gICAgICAgIGNvbnN0IEFUID0gQS5UXG4gICAgICAgIGNvbnN0IGIgPSB0aGlzLl9iO1xuICAgICAgICBjb25zdCBBVEEgPSBBVC5kb3QoQSk7XG4gICAgICAgIGNvbnN0IEFUYiA9IEFULmRvdChiKTtcbiAgICAgICAgdGhpcy5ZID0gTWF0cml4LnNvbHZlX0NHKEFUQSwgQVRiLCB0aGlzLl9yYW5kb21pemVyKTtcbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG59ICIsImltcG9ydCB7IE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiB9IGZyb20gXCIuL0RSLmpzXCI7XG5pbXBvcnQgeyBEaXNqb2ludFNldCB9IGZyb20gXCIuLi9kYXRhc3RydWN0dXJlL2luZGV4LmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgVG9wb01hcFxuICogQG1lbWJlcm9mIG1vZHVsZTpkaW1lbnNpb25hbGl0eV9yZWR1Y3Rpb25cbiAqIEBleHRlbmRzIERSXG4gKi9cbmV4cG9ydCBjbGFzcyBUb3BvTWFwIGV4dGVuZHMgRFIge1xuICAgIC8qKlxuICAgICAqXG4gICAgICogQGNvbnN0cnVjdG9yXG4gICAgICogQG1lbWJlcm9mIG1vZHVsZTpkaW1lbnNpb25hbGl0eV9yZWR1Y3Rpb25cbiAgICAgKiBAYWxpYXMgVG9wb01hcFxuICAgICAqIEBwYXJhbSB7TWF0cml4fSBYIC0gdGhlIGhpZ2gtZGltZW5zaW9uYWwgZGF0YS5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW2QgPSAyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBbbWV0cmljID0gZXVjbGlkZWFuXSAtIHRoZSBtZXRyaWMgd2hpY2ggZGVmaW5lcyB0aGUgZGlzdGFuY2UgYmV0d2VlbiB0d28gcG9pbnRzLlxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbc2VlZCA9IDEyMTJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEByZXR1cm5zIHtUb3BvTWFwfVxuICAgICAqIEBzZWUge0BsaW5rIGh0dHBzOi8vYXJ4aXYub3JnL3BkZi8yMDA5LjAxNTEyLnBkZn1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBkID0gMiwgbWV0cmljID0gZXVjbGlkZWFuLCBzZWVkID0gMTIxMikge1xuICAgICAgICBzdXBlcihYLCBkLCBtZXRyaWMsIHNlZWQpO1xuICAgICAgICBzdXBlci5wYXJhbWV0ZXJfbGlzdCA9IFtdO1xuICAgICAgICBbdGhpcy5fTiwgdGhpcy5fRF0gPSB0aGlzLlguc2hhcGU7XG4gICAgICAgIHRoaXMuX2Rpc3RhbmNlX21hdHJpeCA9IG5ldyBNYXRyaXgodGhpcy5fTiwgdGhpcy5fTiwgMCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICovXG4gICAgX19sYXp5X2Rpc3RhbmNlX21hdHJpeChpLCBqLCBtZXRyaWMpIHtcbiAgICAgICAgY29uc3QgRCA9IHRoaXMuX2Rpc3RhbmNlX21hdHJpeDtcbiAgICAgICAgY29uc3QgWCA9IHRoaXMuWDtcbiAgICAgICAgY29uc3QgRF9paiA9IEQuZW50cnkoaSwgaik7XG4gICAgICAgIGlmIChEX2lqID09PSAwKSB7XG4gICAgICAgICAgICBsZXQgZGlzdCA9IG1ldHJpYyhYLnJvdyhpKSwgWC5yb3coaikpO1xuICAgICAgICAgICAgRC5zZXRfZW50cnkoaSwgaiwgZGlzdCk7XG4gICAgICAgICAgICBELnNldF9lbnRyeShqLCBpLCBkaXN0KTtcbiAgICAgICAgICAgIHJldHVybiBkaXN0O1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBEX2lqO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENvbXB1dGVzIHRoZSBtaW5pbXVtIHNwYW5uaW5nIHRyZWUsIHVzaW5nIGEgZ2l2ZW4gbWV0cmljXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge0Z1bmN0aW9ufSBtZXRyaWNcbiAgICAgKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lwZWRpYS5vcmcvd2lraS9LcnVza2FsJTI3c19hbGdvcml0aG19XG4gICAgICovXG4gICAgX21ha2VfbWluaW11bV9zcGFubmluZ190cmVlKG1ldHJpYyA9IGV1Y2xpZGVhbikge1xuICAgICAgICBjb25zdCBOID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgWCA9IFsuLi50aGlzLlhdO1xuXG4gICAgICAgIGxldCBkaXNqb2ludF9zZXQgPSBuZXcgRGlzam9pbnRTZXQoWCk7XG4gICAgICAgIGNvbnN0IEYgPSBbXTtcbiAgICAgICAgbGV0IEUgPSBbXTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSBpICsgMTsgaiA8IE47ICsraikge1xuICAgICAgICAgICAgICAgIEUucHVzaChbaSwgaiwgdGhpcy5fX2xhenlfZGlzdGFuY2VfbWF0cml4KGksIGosIG1ldHJpYyldKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBFID0gRS5zb3J0KChhLCBiKSA9PiBhWzJdIC0gYlsyXSk7XG5cbiAgICAgICAgZm9yIChjb25zdCBbdSwgdiwgd10gb2YgRSkge1xuICAgICAgICAgICAgY29uc3Qgc2V0X3UgPSBkaXNqb2ludF9zZXQuZmluZChYW3VdKTtcbiAgICAgICAgICAgIGNvbnN0IHNldF92ID0gZGlzam9pbnRfc2V0LmZpbmQoWFt2XSk7XG4gICAgICAgICAgICBpZiAoc2V0X3UgIT09IHNldF92KSB7XG4gICAgICAgICAgICAgICAgRi5wdXNoKFt1LCB2LCB3XSk7XG4gICAgICAgICAgICAgICAgZGlzam9pbnRfc2V0LnVuaW9uKHNldF91LCBzZXRfdik7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4gRi5zb3J0KChhLCBiKSA9PiBhWzJdIC0gYlsyXSk7XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogaW5pdGlhbGl6ZXMgVG9wb01hcC4gU2V0cyBhbGwgcHJvamN0ZWQgcG9pbnRzIHRvIHplcm8sIGFuZCBjb21wdXRlcyBhIG1pbmltdW0gc3Bhbm5pbmcgdHJlZS5cbiAgICAgKi9cbiAgICBpbml0KCkge1xuICAgICAgICB0aGlzLlkgPSBuZXcgTWF0cml4KHRoaXMuX04sIHRoaXMuX2QsIDApO1xuICAgICAgICB0aGlzLl9FbXN0ID0gdGhpcy5fbWFrZV9taW5pbXVtX3NwYW5uaW5nX3RyZWUodGhpcy5fbWV0cmljKTtcbiAgICAgICAgdGhpcy5faXNfaW5pdGlhbGl6ZWQgPSB0cnVlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBSZXR1cm5zIHRydWUgaWYgUG9pbnQgQyBpcyBsZWZ0IG9mIGxpbmUgQUIuXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBQb2ludEEgLSBQb2ludCBBIG9mIGxpbmUgQUJcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBQb2ludEIgLSBQb2ludCBCIG9mIGxpbmUgQUJcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBQb2ludEMgLSBQb2ludCBDXG4gICAgICogQHJldHVybnMge0Jvb2xlYW59XG4gICAgICovXG4gICAgX19odWxsX2Nyb3NzKFtheCwgYXldLCBbYngsIGJ5XSwgW3N4LCBzeV0pIHtcbiAgICAgICAgcmV0dXJuIChieCAtIGF4KSAqIChzeSAtIGF5KSAtIChieSAtIGF5KSAqIChzeCAtIGF4KSA8PSAwO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENvbXB1dGVzIHRoZSBjb252ZXggaHVsbCBvZiB0aGUgc2V0IG9mIFBvaW50cyBTXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBTIC0gU2V0IG9mIFBvaW50cy5cbiAgICAgKiBAc2VlIHtAbGluayBodHRwczovL2VuLndpa2lib29rcy5vcmcvd2lraS9BbGdvcml0aG1fSW1wbGVtZW50YXRpb24vR2VvbWV0cnkvQ29udmV4X2h1bGwvTW9ub3RvbmVfY2hhaW4jSmF2YVNjcmlwdH1cbiAgICAgKiBAcmV0dXJucyB7QXJyYXl9IGNvbnZleCBodWxsIG9mIFMuIFN0YXJ0cyBhdCB0aGUgYm90dG9tLW1vc3QgcG9pbnQgYW5kIGNvbnRpbnVlcyBjb3VudGVyLWNsb2Nrd2lzZS5cbiAgICAgKi9cbiAgICBfX2h1bGwoUykge1xuICAgICAgICBjb25zdCBwb2ludHMgPSBTLnNvcnQoKFt4MSwgeTFdLCBbeDIsIHkyXSkgPT4geTEgLSB5MiB8fCB4MSAtIHgyKTtcbiAgICAgICAgY29uc3QgTiA9IHBvaW50cy5sZW5ndGg7XG4gICAgICAgIGlmIChOIDw9IDIpIHJldHVybiBwb2ludHM7XG5cbiAgICAgICAgY29uc3QgbG93ZXIgPSBbXTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIHdoaWxlIChsb3dlci5sZW5ndGggPj0gMiAmJiB0aGlzLl9faHVsbF9jcm9zcyhsb3dlcltsb3dlci5sZW5ndGggLSAyXSwgbG93ZXJbbG93ZXIubGVuZ3RoIC0gMV0sIHBvaW50c1tpXSkpIHtcbiAgICAgICAgICAgICAgICBsb3dlci5wb3AoKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGxvd2VyLnB1c2gocG9pbnRzW2ldKTtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCB1cHBlciA9IFtdO1xuICAgICAgICBmb3IgKGxldCBpID0gTiAtIDE7IGkgPj0gMDsgLS1pKSB7XG4gICAgICAgICAgICB3aGlsZSAodXBwZXIubGVuZ3RoID49IDIgJiYgdGhpcy5fX2h1bGxfY3Jvc3ModXBwZXJbdXBwZXIubGVuZ3RoIC0gMl0sIHVwcGVyW3VwcGVyLmxlbmd0aCAtIDFdLCBwb2ludHNbaV0pKSB7XG4gICAgICAgICAgICAgICAgdXBwZXIucG9wKCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB1cHBlci5wdXNoKHBvaW50c1tpXSk7XG4gICAgICAgIH1cbiAgICAgICAgdXBwZXIucG9wKCk7XG4gICAgICAgIGxvd2VyLnBvcCgpO1xuICAgICAgICByZXR1cm4gbG93ZXIuY29uY2F0KHVwcGVyKTtcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBGaW5kcyB0aGUgYW5nbGUgdG8gcm90YXRlIFBvaW50IEEgYW5kIEIgdG8gbGllIG9uIGEgbGluZSBwYXJhbGxlbCB0byB0aGUgeC1heGlzLlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtBcnJheX0gUG9pbnRBXG4gICAgICogQHBhcmFtIHtBcnJheX0gUG9pbnRCXG4gICAgICogQHJldHVybiB7T2JqZWN0fSBPYmplY3QgY29udGFpbmluZyB0aGUgc2ludXMtIGFuZCBjb3NpbnVzLXZhbHVlcyBmb3IgYSByb3RhdGlvbi5cbiAgICAgKi9cbiAgICBfX2ZpbmRBbmdsZShbcDF4LCBwMXldLCBbcDJ4LCBwMnldKSB7XG4gICAgICAgIGNvbnN0IG4gPSBldWNsaWRlYW4oW3AxeCwgcDF5XSwgW3AyeCwgcDJ5XSk7XG4gICAgICAgIGlmIChuID09PSAwKVxuICAgICAgICAgICAgcmV0dXJuIHtcbiAgICAgICAgICAgICAgICBzaW46IDAsXG4gICAgICAgICAgICAgICAgY29zOiAxLFxuICAgICAgICAgICAgfTtcbiAgICAgICAgY29uc3QgdmVjID0gWyhwMnggLSBwMXgpIC8gbiwgKHAyeSAtIHAxeSkgLyBuXTtcbiAgICAgICAgY29uc3QgY29zID0gdmVjWzBdO1xuICAgICAgICBsZXQgc2luID0gTWF0aC5zcXJ0KDEgLSBjb3MgKiBjb3MpO1xuICAgICAgICBzaW4gPSB2ZWNbMV0gPj0gMCA/IC1zaW4gOiBzaW47XG4gICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICBzaW46IHNpbixcbiAgICAgICAgICAgIGNvczogY29zLFxuICAgICAgICB9O1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtBcnJheX0gaHVsbFxuICAgICAqIEBwYXJhbSB7QXJyYXl9IHBcbiAgICAgKiBAcGFyYW0ge0Jvb2x9IHRvcEVkZ2VcbiAgICAgKi9cbiAgICBfX2FsaWduX2h1bGwoaHVsbCwgcCwgdG9wRWRnZSkge1xuICAgICAgICBsZXQgdiA9IC0xO1xuICAgICAgICBsZXQgZDI7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgaHVsbC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgY29uc3QgZCA9IGV1Y2xpZGVhbihodWxsW2ldLCBwKTtcbiAgICAgICAgICAgIGlmICh2ID09PSAtMSkge1xuICAgICAgICAgICAgICAgIGQyID0gZDtcbiAgICAgICAgICAgICAgICB2ID0gaTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgaWYgKGQyID4gZCkge1xuICAgICAgICAgICAgICAgICAgICBkMiA9IGQ7XG4gICAgICAgICAgICAgICAgICAgIHYgPSBpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGxldCB2MTtcbiAgICAgICAgbGV0IHYyO1xuICAgICAgICBpZiAodG9wRWRnZSkge1xuICAgICAgICAgICAgdjEgPSBodWxsW3ZdO1xuICAgICAgICAgICAgdjIgPSBodWxsWyh2ICsgMSkgJSBodWxsLmxlbmd0aF07XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBpZiAodiA9PSAwKSB2ID0gaHVsbC5sZW5ndGggLSAxO1xuICAgICAgICAgICAgdjEgPSBodWxsW3ZdO1xuICAgICAgICAgICAgdjIgPSBodWxsWyh2IC0gMSkgJSBodWxsLmxlbmd0aF07XG4gICAgICAgIH1cblxuICAgICAgICBjb25zdCB0cmFuc2Zvcm1hdGlvbiA9IHtcbiAgICAgICAgICAgIHR4OiAtaHVsbFt2XVswXSxcbiAgICAgICAgICAgIHR5OiAtaHVsbFt2XVsxXSxcbiAgICAgICAgfTtcblxuICAgICAgICBpZiAoaHVsbC5sZW5ndGggPj0gMikge1xuICAgICAgICAgICAgY29uc3QgeyBzaW4sIGNvcyB9ID0gdGhpcy5fX2ZpbmRBbmdsZSh2MSwgdjIpO1xuICAgICAgICAgICAgdHJhbnNmb3JtYXRpb24uc2luID0gc2luO1xuICAgICAgICAgICAgdHJhbnNmb3JtYXRpb24uY29zID0gY29zO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdHJhbnNmb3JtYXRpb24uc2luID0gMDtcbiAgICAgICAgICAgIHRyYW5zZm9ybWF0aW9uLmNvcyA9IDE7XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4gdHJhbnNmb3JtYXRpb247XG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogQHByaXZhdGVcbiAgICAgKiBAcGFyYW0ge0FycmF5fSBQb2ludCAtIFRoZSBwb2ludCB3aGljaCBzaG91bGQgZ2V0IHRyYW5zZm9ybWVkLlxuICAgICAqIEBwYXJhbSB7T2JqZWN0fSBUcmFuc2Zvcm1hdGlvbiAtIGNvbnRhaW5zIHRoZSB2YWx1ZXMgZm9yIHRyYW5zbGF0aW9uIGFuZCByb3RhdGlvbi5cbiAgICAgKi9cbiAgICBfX3RyYW5zZm9ybShbcHgsIHB5XSwgeyB0eCwgdHksIHNpbiwgY29zIH0pIHtcbiAgICAgICAgbGV0IHggPSBweCArIHR4O1xuICAgICAgICBsZXQgeSA9IHB5ICsgdHk7XG4gICAgICAgIGxldCB4eCA9IHggKiBjb3MgLSB5ICogc2luO1xuICAgICAgICBsZXQgeXkgPSB4ICogc2luICsgeSAqIGNvcztcbiAgICAgICAgcmV0dXJuIFt4eCwgeXldO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIENhbGxzIHtAbGluayBfX3RyYW5zZm9ybX0gZm9yIGVhY2ggcG9pbnQgaW4gU2V0IENcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7QXJyYXl9IEMgLSBTZXQgb2YgcG9pbnRzLlxuICAgICAqIEBwYXJhbSB7T2JqZWN0fSB0IC0gVHJhbnNmb3JtIG9iamVjdC5cbiAgICAgKiBAcGFyYW0ge051bWJlcn0geU9mZnNldCAtIHZhbHVlIHRvIG9mZnNldCBzZXQgQy5cbiAgICAgKi9cbiAgICBfX3RyYW5zZm9ybV9jb21wb25lbnQoQywgdCwgeU9mZnNldCkge1xuICAgICAgICBjb25zdCBOID0gQy5sZW5ndGg7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgTjsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBjID0gQ1tpXTtcbiAgICAgICAgICAgIGNvbnN0IFtjeCwgY3ldID0gdGhpcy5fX3RyYW5zZm9ybShjLCB0KTtcbiAgICAgICAgICAgIGNbMF0gPSBjeDtcbiAgICAgICAgICAgIGNbMV0gPSBjeSArIHlPZmZzZXQ7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBAcHJpdmF0ZVxuICAgICAqIEBwYXJhbSB7QXJyYXl9IHUgLSBwb2ludCB1XG4gICAgICogQHBhcmFtIHtBcnJheX0gdiAtIHBvaW50IHZcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gdyAtIGVkZ2Ugd2VpZ2h0IHdcbiAgICAgKi9cbiAgICBfX2FsaWduX2NvbXBvbmVudHModSwgdiwgdykge1xuICAgICAgICBjb25zdCBwb2ludHNfdSA9IFsuLi51Ll9fZGlzam9pbnRfc2V0LmNoaWxkcmVuXTtcbiAgICAgICAgY29uc3QgcG9pbnRzX3YgPSBbLi4udi5fX2Rpc2pvaW50X3NldC5jaGlsZHJlbl07XG5cbiAgICAgICAgY29uc3QgaHVsbF91ID0gdGhpcy5fX2h1bGwocG9pbnRzX3UpO1xuICAgICAgICBjb25zdCBodWxsX3YgPSB0aGlzLl9faHVsbChwb2ludHNfdik7XG5cbiAgICAgICAgY29uc3QgdF91ID0gdGhpcy5fX2FsaWduX2h1bGwoaHVsbF91LCB1LCBmYWxzZSk7XG4gICAgICAgIGNvbnN0IHRfdiA9IHRoaXMuX19hbGlnbl9odWxsKGh1bGxfdiwgdiwgdHJ1ZSk7XG5cbiAgICAgICAgdGhpcy5fX3RyYW5zZm9ybV9jb21wb25lbnQocG9pbnRzX3UsIHRfdSwgMCk7XG4gICAgICAgIHRoaXMuX190cmFuc2Zvcm1fY29tcG9uZW50KHBvaW50c192LCB0X3YsIHcpO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIFRyYW5zZm9ybXMgdGhlIGlucHV0ZGF0YSB7QGxpbmsgWH0gdG8gZGltZW5zaW9uYWxpdHkgMi5cbiAgICAgKi9cbiAgICB0cmFuc2Zvcm0oKSB7XG4gICAgICAgIGlmICghdGhpcy5faXNfaW5pdGlhbGl6ZWQpIHRoaXMuaW5pdCgpO1xuICAgICAgICBjb25zdCBFbXN0ID0gdGhpcy5fRW1zdDtcbiAgICAgICAgY29uc3QgWSA9IFsuLi50aGlzLlldO1xuICAgICAgICBjb25zdCBjb21wb25lbnRzID0gbmV3IERpc2pvaW50U2V0KFxuICAgICAgICAgICAgWS5tYXAoKHksIGkpID0+IHtcbiAgICAgICAgICAgICAgICB5LmkgPSBpO1xuICAgICAgICAgICAgICAgIHJldHVybiB5O1xuICAgICAgICAgICAgfSlcbiAgICAgICAgKTtcblxuICAgICAgICBmb3IgKGNvbnN0IFt1LCB2LCB3XSBvZiBFbXN0KSB7XG4gICAgICAgICAgICBjb25zdCBjb21wb25lbnRfdSA9IGNvbXBvbmVudHMuZmluZChZW3VdKTtcbiAgICAgICAgICAgIGNvbnN0IGNvbXBvbmVudF92ID0gY29tcG9uZW50cy5maW5kKFlbdl0pO1xuICAgICAgICAgICAgaWYgKGNvbXBvbmVudF91ID09PSBjb21wb25lbnRfdikgY29udGludWU7XG4gICAgICAgICAgICB0aGlzLl9fYWxpZ25fY29tcG9uZW50cyhjb21wb25lbnRfdSwgY29tcG9uZW50X3YsIHcpO1xuICAgICAgICAgICAgY29tcG9uZW50cy51bmlvbihjb21wb25lbnRfdSwgY29tcG9uZW50X3YpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzLnByb2plY3Rpb247XG4gICAgfVxuXG4gICAgKmdlbmVyYXRvcigpIHtcbiAgICAgICAgaWYgKCF0aGlzLl9pc19pbml0aWFsaXplZCkgdGhpcy5pbml0KCk7XG4gICAgICAgIGNvbnN0IEVtc3QgPSB0aGlzLl9FbXN0O1xuICAgICAgICBjb25zdCBZID0gWy4uLnRoaXMuWV07XG4gICAgICAgIGNvbnN0IGNvbXBvbmVudHMgPSBuZXcgRGlzam9pbnRTZXQoXG4gICAgICAgICAgICBZLm1hcCgoeSwgaSkgPT4ge1xuICAgICAgICAgICAgICAgIHkuaSA9IGk7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHk7XG4gICAgICAgICAgICB9KVxuICAgICAgICApO1xuXG4gICAgICAgIGZvciAoY29uc3QgW3UsIHYsIHddIG9mIEVtc3QpIHtcbiAgICAgICAgICAgIGNvbnN0IGNvbXBvbmVudF91ID0gY29tcG9uZW50cy5maW5kKFlbdV0pO1xuICAgICAgICAgICAgY29uc3QgY29tcG9uZW50X3YgPSBjb21wb25lbnRzLmZpbmQoWVt2XSk7XG4gICAgICAgICAgICBpZiAoY29tcG9uZW50X3UgPT09IGNvbXBvbmVudF92KSBjb250aW51ZTtcbiAgICAgICAgICAgIHRoaXMuX19hbGlnbl9jb21wb25lbnRzKGNvbXBvbmVudF91LCBjb21wb25lbnRfdiwgdyk7XG4gICAgICAgICAgICBjb21wb25lbnRzLnVuaW9uKGNvbXBvbmVudF91LCBjb21wb25lbnRfdik7XG4gICAgICAgICAgICAvKiBsZXQgb2sgPSB0cnVlXG4gICAgICAgICAgICBZLmZvckVhY2goKFt4LCB5XSkgPT4gb2sgPSBvayAmJiAhaXNOYU4oeCkgJiYgIWlzTmFOKHkpKVxuICAgICAgICAgICAgaWYgKCFvaykge1xuICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKC4uLlkpIFxuICAgICAgICAgICAgICAgIHRocm93IFwiZXJyb3JcIiBcbiAgICAgICAgICAgIH0gKi9cbiAgICAgICAgICAgIHlpZWxkIHRoaXMucHJvamVjdGlvbjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcy5wcm9qZWN0aW9uO1xuICAgIH1cbn1cbiIsImltcG9ydCB7IE1hdHJpeCB9IGZyb20gXCIuLi9tYXRyaXgvaW5kZXguanNcIjtcbmltcG9ydCB7IGV1Y2xpZGVhbiB9IGZyb20gXCIuLi9tZXRyaWNzL2luZGV4LmpzXCI7XG5pbXBvcnQgeyBEUiBhcyBEaW1SZWR9IGZyb20gXCIuL0RSLmpzXCI7XG5cbi8qKlxuICogQGNsYXNzXG4gKiBAYWxpYXMgU0FNTU9OXG4gKiBAZXh0ZW5kcyBEUlxuICovXG5leHBvcnQgY2xhc3MgU0FNTU9OIGV4dGVuZHMgRGltUmVkIHtcbiAgICAvKipcbiAgICAgKiBcbiAgICAgKiBAY29uc3RydWN0b3JcbiAgICAgKiBAbWVtYmVyb2YgbW9kdWxlOmRpbWVuc2lvbmFsaXR5X3JlZHVjdGlvblxuICAgICAqIEBhbGlhcyBTQU1NT05cbiAgICAgKiBAcGFyYW0ge01hdHJpeH0gWCAtIHRoZSBoaWdoLWRpbWVuc2lvbmFsIGRhdGEuIFxuICAgICAqIEBwYXJhbSB7TnVtYmVyfSBbZCA9IDJdIC0gdGhlIGRpbWVuc2lvbmFsaXR5IG9mIHRoZSBwcm9qZWN0aW9uLlxuICAgICAqIEBwYXJhbSB7RnVuY3Rpb259IFttZXRyaWMgPSBldWNsaWRlYW5dIC0gdGhlIG1ldHJpYyB3aGljaCBkZWZpbmVzIHRoZSBkaXN0YW5jZSBiZXR3ZWVuIHR3byBwb2ludHMuICBcbiAgICAgKiBAcGFyYW0ge051bWJlcn0gW3NlZWQgPSAxMjEyXSAtIHRoZSBkaW1lbnNpb25hbGl0eSBvZiB0aGUgcHJvamVjdGlvbi5cbiAgICAgKiBAcmV0dXJucyB7U0FNTU9OfVxuICAgICAqIEBzZWUge0BsaW5rIGh0dHBzOi8vYXJ4aXYub3JnL3BkZi8yMDA5LjAxNTEyLnBkZn1cbiAgICAgKi9cbiAgICBjb25zdHJ1Y3RvcihYLCBtYWdpYz0wLjEsIGQ9MiwgbWV0cmljPWV1Y2xpZGVhbiwgc2VlZD0xMjEyKSB7XG4gICAgICAgIHN1cGVyKFgsIGQsIG1ldHJpYywgc2VlZClcbiAgICAgICAgc3VwZXIucGFyYW1ldGVyX2xpc3QgPSBbXCJtYWdpY1wiXTtcbiAgICAgICAgdGhpcy5wYXJhbWV0ZXIoXCJtYWdpY1wiLCBtYWdpYyk7XG4gICAgICAgIFsgdGhpcy5fTiwgdGhpcy5fRCBdID0gdGhpcy5YLnNoYXBlO1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICB9XG5cbiAgICAvKipcbiAgICAgKiBpbml0aWFsaXplcyBTQU1NT04uIFNldHMgYWxsIHByb2pjdGVkIHBvaW50cyB0byB6ZXJvLCBhbmQgY29tcHV0ZXMgYSBtaW5pbXVtIHNwYW5uaW5nIHRyZWUuXG4gICAgICovXG4gICAgaW5pdChEUj1cInJhbmRvbVwiLCBkaXN0YW5jZV9tYXRyaXg9bnVsbCkge1xuICAgICAgICBjb25zdCBOID0gdGhpcy5fTjtcbiAgICAgICAgY29uc3QgZCA9IHRoaXMuX2Q7XG5cbiAgICAgICAgaWYgKERSID09PSBcInJhbmRvbVwiKSB7XG4gICAgICAgICAgICBjb25zdCByYW5kb21pemVyID0gdGhpcy5fcmFuZG9taXplcjtcbiAgICAgICAgICAgIHRoaXMuWSA9IG5ldyBNYXRyaXgoTiwgZCwgKCkgPT4gcmFuZG9taXplci5yYW5kb20pO1xuICAgICAgICB9IGVsc2UgaWYgKERSIGluc3RhbmNlb2YgRGltUmVkKSB7XG4gICAgICAgICAgICB0aGlzLlkgPSBEUi50cmFuc2Zvcm0odGhpcy5YKTtcbiAgICAgICAgfVxuICAgICAgICB0aGlzLmRpc3RhbmNlX21hdHJpeCA9IGRpc3RhbmNlX21hdHJpeCB8fCB0aGlzLl9fZGlzdGFuY2VfbWF0cml4KHRoaXMuWCk7XG4gICAgICAgIHJldHVybiB0aGlzO1xuICAgIH1cblxuICAgIC8qKlxuICAgICAqIEBwcml2YXRlXG4gICAgICogQHBhcmFtIHtNYXRyaXh9IEFcbiAgICAgKiBAcmV0dXJucyB7TWF0cml4fSBcbiAgICAgKi9cbiAgICBfX2Rpc3RhbmNlX21hdHJpeChBKSB7XG4gICAgICAgIGNvbnN0IG1ldHJpYyA9IHRoaXMuX21ldHJpYztcbiAgICAgICAgY29uc3QgTiA9IEEuc2hhcGVbMF07XG4gICAgICAgIGNvbnN0IEQgPSBuZXcgTWF0cml4KE4sIE4pO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IE47ICsraSkge1xuICAgICAgICAgICAgY29uc3QgQV9pID0gQS5yb3coaSk7XG4gICAgICAgICAgICBmb3IgKGxldCBqID0gaTsgaiA8IE47ICsraikge1xuICAgICAgICAgICAgICAgIGxldCBkaXN0YW5jZSA9IChpID09PSBqID8gMCA6IG1ldHJpYyhBX2ksIEEucm93KGopKSk7XG4gICAgICAgICAgICAgICAgRC5zZXRfZW50cnkoaSwgaiwgZGlzdGFuY2UpO1xuICAgICAgICAgICAgICAgIEQuc2V0X2VudHJ5KGosIGksIGRpc3RhbmNlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gRDsgICAgICAgICAgICAgICAgXG4gICAgfVxuXG4gICAgLyoqXG4gICAgICogVHJhbnNmb3JtcyB0aGUgaW5wdXRkYXRhIHtAbGluayBYfSB0byBkaW1lbmlvbmFsaXR5IDIuXG4gICAgICovXG4gICAgdHJhbnNmb3JtKG1heF9pdGVyPTIwMCkge1xuICAgICAgICBpZiAoIXRoaXMuX2lzX2luaXRpYWxpemVkKSB0aGlzLmluaXQoKTtcbiAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBtYXhfaXRlcjsgKytqKSB7XG4gICAgICAgICAgICB0aGlzLl9zdGVwKClcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gdGhpcy5wcm9qZWN0aW9uO1xuICAgIH1cblxuICAgICogZ2VuZXJhdG9yKG1heF9pdGVyPTIwMCkge1xuICAgICAgICBpZiAoIXRoaXMuX2lzX2luaXRpYWxpemVkKSB0aGlzLmluaXQoKTtcblxuICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IG1heF9pdGVyOyArK2opIHtcbiAgICAgICAgICAgIHRoaXMuX3N0ZXAoKVxuICAgICAgICAgICAgeWllbGQgdGhpcy5wcm9qZWN0aW9uO1xuICAgICAgICB9XG5cbiAgICAgICAgcmV0dXJuIHRoaXMucHJvamVjdGlvbjtcbiAgICB9XG5cbiAgICBfc3RlcCgpIHtcbiAgICAgICAgY29uc3QgTUFHSUMgPSB0aGlzLnBhcmFtZXRlcihcIm1hZ2ljXCIpO1xuICAgICAgICBjb25zdCBEID0gdGhpcy5kaXN0YW5jZV9tYXRyaXg7XG4gICAgICAgIGNvbnN0IE4gPSB0aGlzLl9OO1xuICAgICAgICBjb25zdCBkID0gdGhpcy5fZDtcbiAgICAgICAgY29uc3QgbWV0cmljID0gdGhpcy5fbWV0cmljO1xuICAgICAgICBsZXQgWSA9IHRoaXMuWTtcbiAgICAgICAgXG4gICAgICAgIGxldCBHID0gbmV3IE1hdHJpeChOLCBkLCAwKTtcblxuICAgICAgICBsZXQgc3VtID0gbmV3IEZsb2F0NjRBcnJheShkKTtcbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGxldCBlMSA9IG5ldyBGbG9hdDY0QXJyYXkoZCk7XG4gICAgICAgICAgICBsZXQgZTIgPSBuZXcgRmxvYXQ2NEFycmF5KGQpO1xuICAgICAgICAgICAgY29uc3QgWWkgPSBZLnJvdyhpKTtcbiAgICAgICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgTjsgKytqKSB7XG4gICAgICAgICAgICAgICAgaWYgKGkgPT09IGopIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIGNvbnN0IFlqID0gWS5yb3coaik7XG4gICAgICAgICAgICAgICAgY29uc3QgZGVsdGEgPSBuZXcgRmxvYXQ2NEFycmF5KGQpO1xuICAgICAgICAgICAgICAgIGZvciAobGV0IGsgPSAwOyBrIDwgZDsgKytrKSB7XG4gICAgICAgICAgICAgICAgICAgIGRlbHRhW2tdID0gWWlba10gLSBZaltrXVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjb25zdCBkWSA9IG1ldHJpYyhZaSwgWWopO1xuICAgICAgICAgICAgICAgIGNvbnN0IGRYID0gRC5lbnRyeShpLCBqKTtcbiAgICAgICAgICAgICAgICBjb25zdCBkcSA9IGRYIC0gZFk7XG4gICAgICAgICAgICAgICAgY29uc3QgZHIgPSBNYXRoLm1heChkWCAqIGRZLCAxZS0yKTtcbiAgICAgICAgICAgICAgICBmb3IgKGxldCBrID0gMDsgayA8IGQ7ICsraykge1xuICAgICAgICAgICAgICAgICAgICBlMVtrXSArPSBkZWx0YVtrXSAqIGRxIC8gZHI7XG4gICAgICAgICAgICAgICAgICAgIGUyW2tdICs9IChkcSAtIE1hdGgucG93KGRlbHRhW2tdLCAyKSAqICgxICsgZHEgLyBkWSkgLyBkWSkgLyBkcjtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBmb3IgKGxldCBrID0gMDsgayA8IGQ7ICsraykge1xuICAgICAgICAgICAgICAgIGNvbnN0IHZhbCA9IFkuZW50cnkoaSwgaykgKyAoTUFHSUMgKiBlMVtrXSAvIE1hdGguYWJzKGUyW2tdKSB8fCAwKTtcbiAgICAgICAgICAgICAgICBHLnNldF9lbnRyeShpLCBrLCB2YWwpO1xuICAgICAgICAgICAgICAgIHN1bVtrXSArPSB2YWw7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgZm9yIChsZXQgayA9IDA7IGsgPCBkOyArK2spIHtcbiAgICAgICAgICAgIHN1bVtrXSAvPSBOO1xuICAgICAgICB9XG5cbiAgICAgICAgZm9yIChsZXQgaSA9IDA7IGkgPCBOOyArK2kpIHtcbiAgICAgICAgICAgIGZvciAobGV0IGsgPSAwOyBrIDwgZDsgKytrKSB7XG4gICAgICAgICAgICAgICAgWS5zZXRfZW50cnkoaSwgaywgRy5lbnRyeShpLCBrKSAtIHN1bVtrXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIFk7XG4gICAgfVxufSAiXSwibmFtZXMiOlsiZGlzdGFuY2VfbWF0cml4IiwiZG1hdHJpeCIsInNpbXVsdGFuZW91c19wb3dlcml0ZXJhdGlvbiIsIkRpbVJlZCIsIkRSIl0sIm1hcHBpbmdzIjoiOzs7Ozs7O0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNlLGtCQUFRLEVBQUUsQ0FBQyxFQUFFLENBQUMsRUFBRTtBQUMvQixJQUFJLE9BQU8sSUFBSSxDQUFDLElBQUksQ0FBQyxpQkFBaUIsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUM5Qzs7QUNYQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ2Usa0JBQVEsRUFBRSxRQUFRLEVBQUU7QUFDbkMsSUFBSSxJQUFJLENBQUMsR0FBRyxRQUFRLENBQUMsTUFBTSxDQUFDO0FBQzVCLElBQUksSUFBSSxHQUFHLEdBQUcsQ0FBQyxDQUFDO0FBQ2hCLElBQUksSUFBSSxZQUFZLEdBQUcsQ0FBQyxDQUFDO0FBQ3pCLElBQUksSUFBSSxDQUFDLEVBQUUsQ0FBQyxDQUFDO0FBQ2I7QUFDQSxJQUFJLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDaEMsUUFB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QUMsQ0FBQyxDQUFDO0FBQzNELEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksR0FBRyxDQUFDLEtBQUssRUFBRTtBQUNmLFFBQVEsT0FBTyxJQUFJLENBQUMsS0FBSyxFQUFFLENBQUMsTUFBTSxDQUFDLEtBQUssRUFBRSxDQUFDLENBQUMsRUFBRSxDQUFDLEtBQUssQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQzNELEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksR0FBRyxDQUFDLEtBQUssRUFBRTtBQUNmLFFBQVEsT0FBTyxJQUFJLENBQUMsS0FBSyxFQUFFLENBQUMsTUFBTSxDQUFDLEtBQUssRUFBRSxDQUFDLENBQUMsRUFBRSxDQUFDLEtBQUssQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQzNELEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxJQUFJLEtBQUssR0FBRztBQUNoQixRQUFRLE9BQU8sQ0FBQyxJQUFJLENBQUMsS0FBSyxFQUFFLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQztBQUN4QyxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxJQUFJLEtBQUssQ0FBQyxDQUFDLElBQUksRUFBRSxJQUFJLEVBQUUsS0FBSyxHQUFHLE1BQU0sQ0FBQyxDQUFDLEVBQUU7QUFDN0MsUUFBUSxJQUFJLENBQUMsS0FBSyxHQUFHLElBQUksQ0FBQztBQUMxQixRQUFRLElBQUksQ0FBQyxLQUFLLEdBQUcsSUFBSSxDQUFDO0FBQzFCLFFBQVEsSUFBSSxDQUFDLEtBQUssR0FBRyxJQUFJLFlBQVksQ0FBQyxJQUFJLEdBQUcsSUFBSSxDQUFDLENBQUM7QUFDbkQsUUFBUSxLQUFLLElBQUksR0FBRyxHQUFHLENBQUMsRUFBRSxHQUFHLEdBQUcsSUFB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QyxLQUFLLENBQUMsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxDQUFDLENBQUM7QUFDM0YsYUFBYTtBQUNiLFlBQVksQ0FBQyxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsR0FBRyxFQUFFLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsR0FBRyxFQUFFLEdBQUcsQ0FBQyxDQUFDLENBQUM7QUFDckUsU0FBUztBQUNUO0FBQ0E7QUFDQSxRQUFRLEtBQUssSUFBSSxHQUFHLEdBQUcsSUFBSSxHQUFHLENBQUMsRUFBRSxHQUFHLElBQUksQ0FBQyxFQUFFLEVBQUUsR0FBRyxFQUFFO0FBQ2xELFlBQVksS0FBSyxJQUFJLEdBQUcsR0FBRyxJQUFJLEdBQUcsQ0FBQyxFQUFFLEdBQUcsR0FBRyxHQUFHLEVBQUUsRUFBRSxHQUFHLEVBQUU7QUFDdkQsZ0JBQWdCLENBQUMsQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLEdBQUcsRUFBRSxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLEdBQUcsRUFBRSxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQzNGLGFBQWE7QUFDYixZQUFZLENBQUMsQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLEdBQUcsRUFBRSxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLEdBQUcsRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQ3JFLFNBQVM7QUFDVDtBQUNBLFFBQVEsT0FBTyxDQUFDLENBQUM7QUFDakIsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksT0FBTyxFQUFFLENBQUMsQ0FBQyxFQUFFO0FBQ2pCLFFBQVEsTUFBTSxJQUFJLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FB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LENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDcEMsWUFBWSxNQUFNLENBQUMsR0FBRyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQy9CLFlBQVksTUFBTSxPQUFPLEdBQUcsV0FBVyxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxLQUFLLEVBQUUsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3JFLFlBQVksQ0FBQyxDQUFDLFNBQVMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxFQUFFLE9BQU8sQ0FBQyxDQUFDO0FBQ3RDLFlBQVksQ0FBQyxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxLQUFLLEVBQUUsR0FBRyxPQUFPLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDdEQsU0FBUztBQUNULFFBQVEsTUFBTSxNQUFNLEdBQUcsSUFBSSxDQUFDLENBQUMsRUFBRSxTQUFTLENBQUMsQ0FBQztBQUMxQyxRQUFRLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxJQUFJLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDdkMsWUFBWSxDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUMsQ0FBQyxHQUFHLE1BQU0sQ0FBQyxDQUFDO0FBQzdDLFNBQVM7QUFDVCxRQUFRLENBQUMsQ0FBQyxTQUFTLENBQUMsQ0FBQyxDQUFDLENBQUMsRUFBRSxNQUFNLEVBQUM7QUFDaEMsS0FBSztBQUNMLElBQUksT0FBTyxDQUFDLEdBQUcsRUFBRSxDQUFDLEVBQUUsR0FBRyxFQUFFLENBQUMsQ0FBQyxDQUFDO0FBQzVCOztBQzNCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ2Usc0NBQVEsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxjQUFjLENBQUMsR0FBRyxF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QyxJQUFJLEVBQUU7QUFDcEUsUUFBUSxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxNQUFNLEVBQUUsSUFBSSxDQUFDLENBQUM7QUFDbEMsUUFBUSxLQUFLLENBQUMsY0FBYyxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUM7QUFDckMsUUFBUSxJQUFJLENBQUMsU0FBUyxDQUFDLEdBQUcsRUFBRSxJQUFJLENBQUMsR0FBRyxDQUFDLFNBQVMsSUFBSSxJQUFJLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxLQUFLLENBQUMsSUFBSSxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLEdBQUcsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLEVBQUUsSUFBSSxDQUFDLEVBQUUsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQzlHLFFBQVEsT0FBTyxJQUFJLENBQUM7QUFDcEIsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJLFNBQVMsR0FBRztBQUNoQixRQUFRLElBQUksQ0FBQyxVQUFVLEVBQUUsQ0FBQztBQUMxQixRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDLENBQUM7QUFDekIsUUFBUSxNQUFNLElBQUksR0FBRyxJQUFJLENBQUMsRUFBRSxDQUFDO0FBQzdCLFFBQVEsTUFBTSxNQUFNLEdBQUcsSUFBSSxDQUFDLE9BQU8sQ0FBQztBQUNwQztBQUNBLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxNQUFNLEVBQUUsQ0FBQztBQUMvQixRQUFRLENBQUMsQ0FBQyxLQUFLLEdBQUcsQ0FBQyxJQUFJLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUMsS0FBSyxDQUFDLElBQUksQ0FBQyxHQUFHLE1BQU0sQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsRUFBQztBQUMzRixR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QUFDLENBQUMsQ0FBQyxJQUFJLEdBQUcsQ0FBQztBQUMxQixZQUFZLEdBQUcsSUFBSSxHQUFHLENBQUM7QUFDdkIsWUFBWSxPQUFPLEdBQUcsQ0FBQztBQUN2QixTQUFTLENBQUMsQ0FBQztBQUNYO0FBQ0EsUUFBUSxHQUFHLEdBQUcsR0FBRyxDQUFDLEdBQUcsQ0FBQyxDQUFDLElBQUksQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDO0FBQ3JDLFFBQVEsR0FBRyxHQUFHLEdBQUcsQ0FBQyxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztBQUNyQyxRQUFRLEdBQUcsS0FBSyxJQUFJLElBQUksQ0FBQyxDQUFDLENBQUM7QUFDM0IsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxJQUFJLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUMsTUFBTSxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsR0FBRyxHQUFHLENBQUMsQ0FBQyxDQUFDLEdBQUcsR0FBRyxDQUFDLENBQUMsQ0FBQyxHQUFHLEdBQUcsQ0FBQyxDQUFDLENBQUM7QUFDMUY7QUFDQTtBQUNBLFFBQVEsTUFBTSxFQUFFLFlBQVksRUFBRSxDQUFDLEVBQUUsR0FBR0EsNkJBQTJCLENBQUMsQ0FBQyxFQUFFLElBQUksQ0FBQyxFQUFFLENBQUMsQ0FBQztBQUM1RSxRQUFRLElBQUksQ0FBQyxDQUFDLEdBQUcsTUFBTSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQUMsQ0FBQyxTQUFTLEVBQUUsQ0FBQztBQUM1QztBQUNBLFFBQVEsT0FBTyxJQUFJLENBQUMsVUFBVSxDQUFDO0FBQy9CLEtBQUs7QUFDTDtBQUNBO0FBQ0E7O0FDOUZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTyxNQUFNLE9BQU8sU0FBUyxFQUFFO0FBQy9CO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksV0FBVyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxFQUFFLE1BQU0sQ0FBQyxTQUFTLEVBQUUsSUFBSSxDQUFDLElBQUksRUFBRTtBQUNyRCxRQUFRLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLE1BQU0sRUFBRSxJQUFJLENBQUMsQ0FBQztBQUNsQyxRQUFRLE9BQU8sSUFBSSxDQUFDO0FBQ3BCLEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksdUJBQXVCLENBQUMsSUFBSSxFQUFFO0FBQ2xDLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQztBQUN6QixRQUFRLE1BQU0sQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDN0IsUUFBUSxJQUFJLE9BQU8sR0FBRyxJQUFJLENBQUMsV0FBVyxDQUFDLFVBQVUsR0FBRyxDQUFDLEdBQUcsQ0FBQyxDQUFDO0FBQzFELFFBQVEsSUFBSSxPQUFPLEdBQUcsSUFBSSxDQUFDO0FBQzNCLFFBQVEsSUFBSSxRQUFRLEdBQUcsQ0FBQyxRQUFRLENBQUM7QUFDakMsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3BDLFlBQVksTUFBTSxJQUFJLEdBQUcsSUFBSSxDQUFDLE9BQU8sRUFBRSxDQUFDLEVBQUM7QUFDekMsWUFBWSxJQUFJLElBQUksR0FBRyxRQUFRLEVBQUU7QUFDakMsZ0JBQWdCLFFBQVEsR0FBRyxJQUFJLENBQUM7QUFDaEMsZ0JBQWdCLE9BQU8sR0FBRyxDQUFDLENBQUM7QUFDNUIsYUFBYTtBQUNiLFNBQVM7QUFDVCxRQUFRLFFBQVEsR0FBRyxDQUFDLFFBQVEsQ0FBQztBQUM3QixRQUFRLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDcEMsWUFBWSxNQUFNLElBQUksR0FBRyxJQUFJLENBQUMsT0FBTyxFQUFFLENBQUMsRUFBQztBQUN6QyxZQUFZLElBQUksSUFBSSxHQUFHLFFBQVEsRUFBRTtBQUNqQyxnQkFBZ0IsUUFBUSxHQUFHLElBQUksQ0FBQztBQUNoQyxnQkFBZ0IsT0FBTyxHQUFHLENBQUMsQ0FBQztBQUM1QixhQUFhO0FBQ2IsU0FBUztBQUNULFFBQVEsT0FBTyxDQUFDLE9BQU8sRUFBRSxPQUFPLEVBQUUsUUFBUSxDQUFDLENBQUM7QUFDNUMsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJLFNBQVMsR0FBRztBQUNoQixRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDLENBQUM7QUFDekIsUUFBUSxNQUFNLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQzdCLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLEVBQUUsQ0FBQztBQUMxQixRQUFRLE1BQU0sTUFBTSxHQUFHLElBQUksQ0FBQyxPQUFPLENBQUM7QUFDcEMsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDO0FBQ3RDLFFBQVEsSUFBSSxJQUFJLEdBQUcsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxLQUFLLE1BQU0sQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUN4RDtBQUNBLFFBQVEsS0FBSyxJQUFJLElBQUksR0FBRyxDQUFDLEVBQUUsSUFBSSxHQUFHLENBQUMsRUFBRSxFQUFFLElBQUksRUFBRTtBQUM3QyxZQUFZLElBQUksUUFBUSxHQUFHLElBQUksQ0FBQztBQUNoQztBQUNBLFlBQVksTUFBTSxDQUFDLE9BQU8sRUFBRSxPQUFPLEVBQUUsSUFBSSxDQUFDLEdBQUcsSUFBSSxDQUFDLHVCQUF1QixDQUFDLElBQUksQ0FBQyxDQUFDO0FBQ2hGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFlBQVksSUFBSSxJQUFJLEtBQUssQ0FBQyxFQUFFO0FBQzVCO0FBQ0EsZ0JBQWdCLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDNUMsb0JBQW9CLE1BQU0sSUFBSSxHQUFHLElBQUksQ0FBQyxPQUFPLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDbEQsb0JBQW9CLE1BQU0sSUFBSSxHQUFHLElBQUksQ0FBQyxPQUFPLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDbEQsb0JBQW9CLE1BQU0sR0FBRyxHQUFHLENBQUMsSUFBSSxJQUFJLENBQUMsR0FBRyxJQUFJLElBQUksQ0FBQyxHQUFHLElBQUksSUFBSSxDQUFDLEtBQUssQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDO0FBQ2pGLG9CQUFvQixDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxJQUFJLEVBQUUsR0FBRyxDQUFDLENBQUM7QUFDOUMsaUJBQWlCO0FBQ2pCO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZ0JBQWdCLElBQUksR0FBRyxDQUFDLENBQUMsRUFBRSxDQUFDLEtBQUssSUFBSSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxJQUFJLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLElBQUksQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLElBQUksQ0FBQyxLQUFLLENBQUMsRUFBQztBQUM1RyxhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0EsUUFBUSxJQUFJLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQztBQUNuQixRQUFRLE9BQU8sSUFBSSxDQUFDLFVBQVUsQ0FBQztBQUMvQixLQUFLO0FBQ0w7O0FDL0ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTyxNQUFNLEdBQUcsU0FBUyxFQUFFLENBQUM7QUFDNUI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksV0FBVyxDQUFDLENBQUMsRUFBRSxNQUFNLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxNQUFNLEdBQUcsU0FBUyxFQUFFLElBQUksQ0FBQyxJQUFJLEVBQUU7QUFDakUsUUFBUSxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxNQUFNLEVBQUUsSUFBSSxDQUFDLENBQUM7QUFDbEMsUUFBUSxLQUFLLENBQUMsY0FBYyxHQUFHLENBQUMsUUFBUSxDQUFDLENBQUM7QUFDMUMsUUFBUSxJQUFJLENBQUMsU0FBUyxDQUFDLFFBQVEsRUFBRSxNQUFNLENBQUMsQ0FBQztBQUN6QyxRQUFRLE9BQU8sSUFBSSxDQUFDO0FBQ3BCLEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksU0FBUyxHQUFHO0FBQ2hCLFFBQVEsSUFBSSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQztBQUN2QixRQUFRLElBQUksRUFBRSxJQUFJLEVBQUUsSUFBSSxFQUFFLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQztBQUNyQyxRQUFRLElBQUksTUFBTSxHQUFHLElBQUksQ0FBQyxPQUFPLENBQUM7QUFDbEMsUUFBUSxJQUFJLGFBQWEsR0FBRyxFQUFFLENBQUM7QUFDL0IsUUFBUSxJQUFJLFFBQVEsR0FBRyxDQUFDLENBQUM7QUFDekIsUUFBUSxNQUFNLENBQUMsT0FBTyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsS0FBSztBQUNqQyxZQUFZLElBQUksQ0FBQyxJQUFJLGFBQWEsRUFBRTtBQUNwQyxnQkFBZ0IsYUFBYSxDQUFDLENBQUMsQ0FBQyxDQUFDLEtBQUssRUFBRSxDQUFDO0FBQ3pDLGdCQUFnQixhQUFhLENBQUMsQ0FBQyxDQUFDLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDckQsYUFBYSxNQUFNO0FBQ25CLGdCQUFnQixhQUFhLENBQUMsQ0FBQyxDQUFDLEdBQUc7QUFDbkMsb0JBQW9CLElBQUksRUFBRSxRQUFRLEVBQUU7QUFDcEMsb0JBQW9CLE9BQU8sRUFBRSxDQUFDO0FBQzlCLG9CQUFvQixNQUFNLEVBQUUsQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3RDLGlCQUFpQixDQUFDO0FBQ2xCLGFBQWE7QUFDYixTQUFTLEVBQUM7QUFDVjtBQUNBO0FBQ0EsUUFBUSxJQUFJLE1BQU0sR0FBRyxDQUFDLENBQUMsSUFBSSxDQUFDO0FBQzVCLFFBQVEsSUFBSSxNQUFNLEdBQUcsSUFBSSxNQUFNLENBQUMsUUFBUSxFQUFFLElBQUksRUFBQztBQUMvQyxRQUFRLEtBQUssSUFBSSxLQUFLLElBQUksYUFBYSxFQUFFO0FBQ3pDLFlBQVksSUFBSSxDQUFDLEdBQUcsTUFBTSxDQUFDLElBQUksQ0FBQyxhQUFhLENBQUMsS0FBSyxDQUFDLENBQUMsSUFBSSxDQUFDLENBQUM7QUFDM0QsWUFBWSxJQUFJLE1BQU0sR0FBRyxDQUFDLENBQUMsUUFBUSxDQUFDO0FBQ3BDLFlBQVksS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLElBQUksRUFBRSxFQUFFLENBQUMsRUFBRTtBQUMzQyxnQkFBZ0IsTUFBTSxDQUFDLFNBQVMsQ0FBQyxhQUFhLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRSxNQUFNLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUN4RSxhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0EsUUFBUSxJQUFJLEdBQUcsR0FBRyxJQUFJLE1BQU0sQ0FBQyxJQUFJLEVBQUUsSUFBSSxDQUFDLENBQUM7QUFDekMsUUFBUSxLQUFLLElBQUksS0FBSyxJQUFJLGFBQWEsRUFBRTtBQUN6QyxZQUFZLElBQUksQ0FBQyxHQUFHLE1BQU0sQ0FBQyxHQUFHLENBQUMsYUFBYSxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDO0FBQ3hELFlBQVksSUFBSSxDQUFDLEdBQUcsSUFBSSxNQUFNLENBQUMsSUFBSSxFQUFFLENBQUMsRUFBRSxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsQ0FBQyxDQUFDLEdBQUcsTUFBTSxDQUFDLENBQUM7QUFDOUQsWUFBWSxJQUFJLENBQUMsR0FBRyxhQUFhLENBQUMsS0FBSyxDQUFDLENBQUMsS0FBSyxDQUFDO0FBQy9DLFlBQVksR0FBRyxHQUFHLEdBQUcsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsU0FBUyxFQUFFLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUN4RCxTQUFTO0FBQ1Q7QUFDQTtBQUNBLFFBQVEsSUFBSSxHQUFHLEdBQUcsSUFBSSxNQUFNLENBQUMsSUFBSSxFQUFFLElBQUksQ0FBQyxDQUFDO0FBQ3pDLFFBQVEsS0FBSyxJQUFJLEtBQUssSUFBSSxhQUFhLEVBQUU7QUFDekMsWUFBWSxJQUFJLENBQUMsR0FBRyxNQUFNLENBQUMsR0FBRyxDQUFDLGFBQWEsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQztBQUN4RCxZQUFZLElBQUksQ0FBQyxHQUFHLElBQUksTUFBTSxDQUFDLElBQUksRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxFQUFDO0FBQ3BELFlBQVksSUFBSSxDQUFDLEdBQUcsYUFBYSxDQUFDLEtBQUssQ0FBQyxDQUFDLElBQUksQ0FBQztBQUM5QyxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxhQUFhLENBQUMsS0FBSyxDQUFDLENBQUMsS0FBSyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDeEUsZ0JBQWdCLElBQUksS0FBSyxHQUFHLElBQUksTUFBTSxDQUFDLElBQUksRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ2xGLGdCQUFnQixHQUFHLEdBQUcsR0FBRyxDQUFDLEdBQUcsQ0FBQyxLQUFLLENBQUMsR0FBRyxDQUFDLEtBQUssQ0FBQyxTQUFTLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDNUQsYUFBYTtBQUNiLFNBQVM7QUFDVDtBQUNBLFFBQVEsSUFBSSxFQUFFLFlBQVksRUFBRSxDQUFDLEVBQUUsR0FBR0EsNkJBQTJCLENBQUMsR0FBRyxDQUFDLE9BQU8sRUFBRSxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsRUFBRSxJQUFJLENBQUMsRUFBRSxFQUFDO0FBQzlGLFFBQVEsQ0FBQyxHQUFHLE1BQU0sQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUMsU0FBUyxHQUFFO0FBQ3RDLFFBQVEsSUFBSSxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsRUFBQztBQUN6QjtBQUNBO0FBQ0EsUUFBUSxPQUFPLElBQUksQ0FBQyxVQUFVLENBQUM7QUFDL0IsS0FBSztBQUNMOztBQ2xGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ08sTUFBTSxHQUFHLFNBQVMsRUFBRSxDQUFDO0FBQzVCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJLFdBQVcsQ0FBQyxDQUFDLEVBQUUsU0FBUyxFQUFFLENBQUMsQ0FBQyxDQUFDLEVBQUUsTUFBTSxDQUFDLFNBQVMsRUFBRSxJQUFJLENBQUMsSUFBSSxFQUFFO0FBQ2hFLFFBQVEsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsTUFBTSxFQUFFLElBQUksQ0FBQyxDQUFDO0FBQ2xDLFFBQVEsS0FBSyxDQUFDLGNBQWMsR0FBRyxDQUFDLEdBQUcsQ0FBQyxDQUFDO0FBQ3JDLFFBQVEsSUFBSSxDQUFDLFNBQVMsQ0FBQyxHQUFHLEVBQUUsSUFBSSxDQUFDLEdBQUcsQ0FBQyxTQUFTLElBQUksSUFBSSxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsS0FBSyxDQUFDLElBQUksQ0FBQyxFQUFFLEdBQUcsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLEVBQUUsSUFBSSxDQUFDLEVBQUUsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3ZHLFFBQVEsT0FBTyxJQUFJLENBQUM7QUFDcEIsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxTQUFTLEdBQUc7QUFDaEIsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQyxDQUFDO0FBQ3pCLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLEVBQUUsQ0FBQztBQUMxQixRQUFRLE1BQU0sSUFBSSxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDN0IsUUFBUSxNQUFNLElBQUksR0FBRyxJQUFJLENBQUMsRUFBRSxDQUFDO0FBQzdCLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLFNBQVMsQ0FBQyxHQUFHLENBQUMsQ0FBQztBQUN0QyxRQUFRLE1BQU0sRUFBRSxHQUFHLG1CQUFtQixDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsSUFBSSxFQUFFLElBQUksQ0FBQyxPQUFPLENBQUMsQ0FBQztBQUNqRSxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDdEMsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxJQUFJLEVBQUUsSUFBSSxDQUFDLENBQUM7QUFDekM7QUFDQSxRQUFRLEtBQUssSUFBSSxHQUFHLEdBQUcsQ0FBQyxFQUFFLEdBQUcsR0FBRyxJQUFJLEVBQUUsRUFBRSxHQUFHLEVBQUU7QUFDN0MsWUFBWSxNQUFNLE1BQU0sR0FBRyxFQUFFLENBQUMsR0FBRyxDQUFDLENBQUM7QUFDbkMsWUFBWSxNQUFNLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxDQUFDLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQyxFQUFFLENBQUMsS0FBSyxDQUFDLENBQUMsS0FBSyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxHQUFHLEVBQUUsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUMvRixZQUFZLE1BQU0sQ0FBQyxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ2pDLFlBQVksS0FBSyxDQUFDLEdBQUcsSUFBSSxHQUFHO0FBQzVCLGdCQUFnQixNQUFNLE9BQU8sR0FBRyxXQUFXLENBQUMsQ0FBQyxDQUFDLElBQUksQ0FBQyxHQUFHLElBQUksQ0FBQztBQUMzRCxnQkFBZ0IsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM1QyxvQkFBb0IsQ0FBQyxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLE9BQU8sQ0FBQyxDQUFDO0FBQy9ELGlCQUFpQjtBQUNqQixhQUFhO0FBQ2I7QUFDQSxZQUFZLElBQUksQ0FBQyxHQUFHLE1BQU0sQ0FBQyxRQUFRLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxJQUFJLENBQUMsV0FBVyxDQUFDLENBQUM7QUFDNUQsWUFBWSxDQUFDLEdBQUcsQ0FBQyxDQUFDLE1BQU0sQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUM7QUFDaEMsWUFBWSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3hDLGdCQUFnQixDQUFDLENBQUMsU0FBUyxDQUFDLEdBQUcsRUFBRSxNQUFNLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDN0QsYUFBYTtBQUNiLFNBQVM7QUFDVDtBQUNBLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxNQUFNLENBQUMsSUFBSSxFQUFFLElBQUksRUFBRSxVQUFVLENBQUMsQ0FBQztBQUNyRCxRQUFRLE1BQU0sRUFBRSxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDNUIsUUFBUSxNQUFNLENBQUMsR0FBRyxFQUFFLENBQUMsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsQ0FBQztBQUMvQixRQUFRLE1BQU0sRUFBRSxZQUFZLEVBQUUsQ0FBQyxFQUFFLEdBQUdBLDZCQUEyQixDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsT0FBTyxFQUFFLEVBQUUsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQ3RGLFFBQVEsSUFBSSxDQUFDLENBQUMsR0FBRyxNQUFNLENBQUMsSUFBSSxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUNsRDtBQUNBO0FBQ0EsUUFBUSxPQUFPLElBQUksQ0FBQyxVQUFVLENBQUM7QUFDL0IsS0FBSztBQUNMOztBQ2xFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ08sTUFBTSxJQUFJLFNBQVMsRUFBRSxDQUFDO0FBQzdCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksV0FBVyxDQUFDLENBQUMsRUFBRSxTQUFTLEVBQUUsQ0FBQyxDQUFDLENBQUMsRUFBRSxNQUFNLENBQUMsU0FBUyxFQUFFLElBQUksQ0FBQyxJQUFJLEVBQUU7QUFDaEUsUUFBUSxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxNQUFNLEVBQUUsSUFBSSxDQUFDLENBQUM7QUFDbEMsUUFBUSxLQUFLLENBQUMsY0FBYyxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUM7QUFDckMsUUFBUSxJQUFJLENBQUMsU0FBUyxDQUFDLEdBQUcsRUFBRSxJQUFJLENBQUMsR0FBRyxDQUFDLFNBQVMsSUFBSSxJQUFJLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxLQUFLLENBQUMsSUFBSSxDQUFDLEVBQUUsR0FBRyxFQUFFLENBQUMsRUFBRSxDQUFDLENBQUMsRUFBRSxJQUFJLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDdkcsUUFBUSxJQUFJLElBQUksQ0FBQyxFQUFFLElBQUksQ0FBQyxFQUFFLE1BQU0sQ0FBQyx5QkFBeUIsRUFBRSxJQUFJLENBQUMsRUFBRSxDQUFDLHNFQUFzRSxFQUFFLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQztBQUNsSixRQUFRLE9BQU8sSUFBSSxDQUFDO0FBQ3BCLEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksU0FBUyxHQUFHO0FBQ2hCLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQztBQUN6QixRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDMUIsUUFBUSxNQUFNLEVBQUUsSUFBSSxFQUFFLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUM7QUFDcEMsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLENBQUMsU0FBUyxDQUFDLEdBQUcsQ0FBQyxDQUFDO0FBQ3RDO0FBQ0EsUUFBUSxNQUFNLEVBQUUsR0FBRyxtQkFBbUIsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLElBQUksRUFBRSxJQUFJLENBQUMsT0FBTyxDQUFDLENBQUM7QUFDakU7QUFDQSxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsUUFBUSxDQUFDLENBQUM7QUFDN0MsUUFBUSxNQUFNLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxJQUFJLEVBQUUsSUFBSSxFQUFFLENBQUMsQ0FBQyxDQUFDO0FBQzVDO0FBQ0EsUUFBUSxLQUFLLElBQUksR0FBRyxHQUFHLENBQUMsRUFBRSxHQUFHLEdBQUcsSUFBSSxFQUFFLEVBQUUsR0FBRyxFQUFFO0FBQzdDO0FBQ0EsWUFBWSxNQUFNLEdBQUcsR0FBRyxDQUFDLEdBQUcsRUFBRSxHQUFHLEVBQUUsQ0FBQyxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUMsRUFBQztBQUN2RCxZQUFZLElBQUksR0FBRyxHQUFHLE1BQU0sQ0FBQyxJQUFJLENBQUMsR0FBRyxDQUFDLEdBQUcsQ0FBQyxDQUFDLElBQUksQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDMUQ7QUFDQSxZQUFZLEdBQUcsR0FBRyxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUMsRUFBQztBQUM1QjtBQUNBLFlBQVksTUFBTSxDQUFDLEdBQUcsR0FBRyxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsU0FBUyxFQUFFLENBQUMsQ0FBQztBQUMvQyxZQUFZLE1BQU0sRUFBRSxZQUFZLEVBQUUsQ0FBQyxFQUFFLEdBQUdBLDZCQUEyQixDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUMxRTtBQUNBLFlBQVksTUFBTSxLQUFLLEdBQUcsTUFBTSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUN6QztBQUNBLFlBQVksTUFBTSxHQUFHLEdBQUcsS0FBSyxDQUFDLFNBQVMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxLQUFLLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxHQUFHLElBQUksQ0FBQyxJQUFJLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDL0UsWUFBWSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM1QyxnQkFBZ0IsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDaEQsb0JBQW9CLENBQUMsQ0FBQyxTQUFTLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsS0FBSyxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsSUFBSSxDQUFDLEtBQUssQ0FBQyxHQUFHLENBQUMsR0FBRyxDQUFDLEVBQUUsR0FBRyxHQUFHLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ2hILGlCQUFpQjtBQUNqQixhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0E7QUFDQSxRQUFRLE1BQU0sRUFBRSxZQUFZLEVBQUUsQ0FBQyxFQUFFLEdBQUdBLDZCQUEyQixDQUFDLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUM7QUFDMUUsUUFBUSxJQUFJLENBQUMsQ0FBQyxHQUFHLE1BQU0sQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDLFNBQVMsRUFBRSxDQUFDO0FBQ3JEO0FBQ0E7QUFDQSxRQUFRLE9BQU8sSUFBSSxDQUFDLFVBQVUsQ0FBQztBQUMvQixLQUFLO0FBQ0w7O0FDcEVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTyxNQUFNLElBQUksU0FBUyxFQUFFLENBQUM7QUFDN0I7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJLFdBQVcsQ0FBQyxDQUFDLEVBQUUsVUFBVSxDQUFDLEVBQUUsRUFBRSxPQUFPLENBQUMsRUFBRSxFQUFFLENBQUMsQ0FBQyxDQUFDLEVBQUUsTUFBTSxDQUFDLFNBQVMsRUFBRSxJQUFJLENBQUMsSUFBSSxFQUFFO0FBQ2hGLFFBQVEsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsTUFBTSxFQUFFLElBQUksQ0FBQyxDQUFDO0FBQ2xDLFFBQVEsS0FBSyxDQUFDLGNBQWMsR0FBRyxDQUFDLFlBQVksRUFBRSxTQUFTLENBQUMsQ0FBQztBQUN6RCxRQUFRLEVBQUUsSUFBSSxDQUFDLEVBQUUsRUFBRSxJQUFJLENBQUMsRUFBRSxFQUFFLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUM7QUFDNUMsUUFBUSxJQUFJLENBQUMsU0FBUyxDQUFDLFlBQVksRUFBRSxJQUFJLENBQUMsR0FBRyxDQUFDLFVBQVUsRUFBRSxJQUFJLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDeEUsUUFBUSxJQUFJLENBQUMsU0FBUyxDQUFDLFNBQVMsRUFBRSxPQUFPLENBQUMsQ0FBQztBQUMzQyxRQUFRLElBQUksQ0FBQyxLQUFLLEdBQUcsQ0FBQyxDQUFDO0FBQ3ZCLFFBQVEsSUFBSSxDQUFDLENBQUMsR0FBRyxJQUFJLE1BQU0sQ0FBQyxJQUFJLENBQUMsRUFBRSxFQUFFLElBQUksQ0FBQyxFQUFFLEVBQUUsTUFBTSxJQUFJLENBQUMsV0FBVyxDQUFDLE1BQU0sQ0FBQyxDQUFDO0FBQzdFLFFBQVEsT0FBTyxJQUFJLENBQUM7QUFDcEIsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksSUFBSSxDQUFDLGVBQWUsQ0FBQyxJQUFJLEVBQUU7QUFDL0I7QUFDQSxRQUFRLE1BQU0sT0FBTyxHQUFHLElBQUksQ0FBQyxHQUFHLENBQUMsSUFBSSxDQUFDLFdBQVcsQ0FBQyxDQUFDO0FBQ25ELFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLEVBQUUsQ0FBQztBQUMxQixRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDMUIsUUFBUSxNQUFNLE1BQU0sR0FBRyxJQUFJLENBQUMsT0FBTyxDQUFDO0FBQ3BDLFFBQVEsTUFBTSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQztBQUN6QixRQUFRLElBQUksS0FBSyxDQUFDO0FBQ2xCLFFBQVEsSUFBSSxlQUFlLEVBQUU7QUFDN0IsWUFBWSxLQUFLLEdBQUcsZUFBZSxDQUFDO0FBQ3BDLFNBQVMsTUFBTTtBQUNmLFlBQVksS0FBSyxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUNyQyxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDeEMsZ0JBQWdCLE1BQU0sR0FBRyxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDckMsZ0JBQWdCLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ2hELG9CQUFvQixNQUFNLFFBQVEsR0FBRyxNQUFNLENBQUMsR0FBRyxFQUFFLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLEVBQUM7QUFDMUQsb0JBQW9CLEtBQUssQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxRQUFRLENBQUMsQ0FBQztBQUNwRCxvQkFBb0IsS0FBSyxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDO0FBQ3BELGlCQUFpQjtBQUNqQixhQUFhO0FBQ2I7QUFDQSxTQUFTO0FBQ1Q7QUFDQSxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsT0FBTyxDQUFDLENBQUM7QUFDNUM7QUFDQSxRQUFRLElBQUksQ0FBQyxNQUFNLEdBQUcsSUFBSSxNQUFNLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxPQUFPLENBQUMsQ0FBQztBQUNoRCxRQUFRLElBQUksQ0FBQyxNQUFNLEdBQUcsSUFBSSxNQUFNLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUMxQztBQUNBO0FBQ0EsUUFBUSxJQUFJLElBQUksR0FBRyxJQUFJLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDeEMsUUFBUSxNQUFNLEdBQUcsR0FBRyxJQUFJLENBQUM7QUFDekIsUUFBUSxNQUFNLFFBQVEsR0FBRyxFQUFFLENBQUM7QUFDNUIsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3BDLFlBQVksSUFBSSxPQUFPLEdBQUcsQ0FBQyxRQUFRLENBQUM7QUFDcEMsWUFBWSxJQUFJLE9BQU8sR0FBRyxRQUFRLENBQUM7QUFDbkMsWUFBWSxJQUFJLElBQUksR0FBRyxDQUFDLENBQUM7QUFDekIsWUFBWSxJQUFJLElBQUksR0FBRyxLQUFLLENBQUM7QUFDN0I7QUFDQSxZQUFZLElBQUksR0FBRyxHQUFHLENBQUMsQ0FBQztBQUN4QixZQUFZLE1BQU0sQ0FBQyxJQUFJLEVBQUU7QUFDekIsZ0JBQWdCLElBQUksSUFBSSxHQUFHLENBQUMsQ0FBQztBQUM3QixnQkFBZ0IsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM1QyxvQkFBb0IsSUFBSSxFQUFFLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDO0FBQ2pFLG9CQUFvQixJQUFJLENBQUMsS0FBSyxDQUFDLEVBQUUsRUFBRSxHQUFHLENBQUMsQ0FBQztBQUN4QyxvQkFBb0IsSUFBSSxDQUFDLENBQUMsQ0FBQyxHQUFHLEVBQUUsQ0FBQztBQUNqQyxvQkFBb0IsSUFBSSxJQUFJLEVBQUUsQ0FBQztBQUMvQixpQkFBaUI7QUFDakIsZ0JBQWdCLElBQUksS0FBSyxHQUFHLENBQUMsQ0FBQztBQUM5QixnQkFBZ0IsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM1QyxvQkFBb0IsSUFBSSxFQUFFLEdBQUcsQ0FBQyxJQUFJLEtBQUssQ0FBQyxJQUFJLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQyxDQUFDLEdBQUcsSUFBSSxDQUFDO0FBQy9ELG9CQUFvQixJQUFJLENBQUMsQ0FBQyxDQUFDLEdBQUcsRUFBRSxDQUFDO0FBQ2pDLG9CQUFvQixJQUFJLEVBQUUsR0FBRyxJQUFJLEVBQUU7QUFDbkMsd0JBQXdCLEtBQUssSUFBSSxFQUFFLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsQ0FBQztBQUNuRCxxQkFBcUI7QUFDckIsaUJBQWlCO0FBQ2pCLGdCQUFnQixJQUFJLEtBQUssR0FBRyxPQUFPLEVBQUU7QUFDckMsb0JBQW9CLE9BQU8sR0FBRyxJQUFJLENBQUM7QUFDbkMsb0JBQW9CLElBQUksR0FBRyxDQUFDLE9BQU8sS0FBSyxRQUFRLEtBQUssSUFBSSxHQUFHLENBQUMsS0FBSyxDQUFDLElBQUksR0FBRyxPQUFPLElBQUksQ0FBQyxDQUFDLENBQUM7QUFDeEYsaUJBQWlCLE1BQU07QUFDdkIsb0JBQW9CLE9BQU8sR0FBRyxJQUFJLENBQUM7QUFDbkMsb0JBQW9CLElBQUksR0FBRyxDQUFDLE9BQU8sS0FBSyxDQUFDLFFBQVEsS0FBSyxJQUFJLEdBQUcsQ0FBQyxLQUFLLENBQUMsSUFBSSxHQUFHLE9BQU8sSUFBSSxDQUFDLENBQUMsQ0FBQztBQUN6RixpQkFBaUI7QUFDakIsZ0JBQWdCLEVBQUUsR0FBRyxDQUFDO0FBQ3RCLGdCQUFnQixJQUFJLElBQUksQ0FBQyxHQUFHLENBQUMsS0FBSyxHQUFHLE9BQU8sQ0FBQyxHQUFHLEdBQUcsRUFBRSxJQUFJLEdBQUcsSUFBSSxDQUFDO0FBQ2pFLGdCQUFnQixJQUFJLEdBQUcsSUFBSSxRQUFRLEVBQUUsSUFBSSxHQUFHLElBQUksQ0FBQztBQUNqRCxhQUFhO0FBQ2I7QUFDQSxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDeEMsZ0JBQWdCLENBQUMsQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUMzQyxhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0E7QUFDQSxRQUFRLE1BQU0sSUFBSSxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsT0FBTyxFQUFDO0FBQzlDLFFBQVEsTUFBTSxFQUFFLEdBQUcsQ0FBQyxHQUFHLENBQUMsQ0FBQztBQUN6QixRQUFRLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDcEMsWUFBWSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3hDLGdCQUFnQixNQUFNLENBQUMsR0FBRyxJQUFJLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLElBQUksRUFBRSxFQUFFLE1BQU0sQ0FBQyxDQUFDO0FBQ2pGLGdCQUFnQixJQUFJLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDeEMsZ0JBQWdCLElBQUksQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUN4QyxhQUFhO0FBQ2IsU0FBUztBQUNULFFBQVEsSUFBSSxDQUFDLEVBQUUsR0FBRyxJQUFJLENBQUM7QUFDdkIsUUFBUSxPQUFPLElBQUksQ0FBQztBQUNwQixLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxTQUFTLENBQUMsVUFBVSxDQUFDLEdBQUcsRUFBRTtBQUM5QixRQUFRLElBQUksQ0FBQyxVQUFVLEVBQUUsQ0FBQztBQUMxQixRQUFRLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxVQUFVLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDN0MsWUFBWSxJQUFJLENBQUMsSUFBSSxFQUFFLENBQUM7QUFDeEIsU0FBUztBQUNULFFBQVEsT0FBTyxJQUFJLENBQUMsVUFBVSxDQUFDO0FBQy9CLEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJLEVBQUUsU0FBUyxDQUFDLFVBQVUsQ0FBQyxHQUFHLEVBQUU7QUFDaEMsUUFBUSxJQUFJLENBQUMsVUFBVSxFQUFFLENBQUM7QUFDMUIsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsVUFBVSxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQzdDLFlBQVksSUFBSSxDQUFDLElBQUksRUFBRSxDQUFDO0FBQ3hCLFlBQVksTUFBTSxJQUFJLENBQUMsVUFBVSxDQUFDO0FBQ2xDLFNBQVM7QUFDVCxRQUFRLE9BQU8sSUFBSSxDQUFDLFVBQVUsQ0FBQztBQUMvQixLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxJQUFJLEdBQUc7QUFDWCxRQUFRLE1BQU0sSUFBSSxHQUFHLEVBQUUsSUFBSSxDQUFDLEtBQUssQ0FBQztBQUNsQyxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDMUIsUUFBUSxNQUFNLEtBQUssR0FBRyxJQUFJLENBQUMsTUFBTSxDQUFDO0FBQ2xDLFFBQVEsTUFBTSxLQUFLLEdBQUcsSUFBSSxDQUFDLE1BQU0sQ0FBQztBQUNsQyxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDMUIsUUFBUSxNQUFNLE9BQU8sR0FBRyxJQUFJLENBQUMsUUFBUSxDQUFDO0FBQ3RDLFFBQVEsTUFBTSxHQUFHLEdBQUcsSUFBSSxDQUFDLEVBQUUsQ0FBQztBQUM1QixRQUFRLElBQUksQ0FBQyxHQUFHLElBQUksQ0FBQyxDQUFDLENBQUM7QUFDdkI7QUFDQTtBQUNBLFFBQVEsTUFBTSxJQUFJLEdBQUcsSUFBSSxHQUFHLEdBQUcsR0FBRyxDQUFDLEdBQUcsQ0FBQyxDQUFDO0FBQ3hDO0FBQ0E7QUFDQSxRQUFRLE1BQU0sRUFBRSxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsT0FBTyxFQUFDO0FBQzVDLFFBQVEsSUFBSSxJQUFJLEdBQUcsQ0FBQyxDQUFDO0FBQ3JCLFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUNwQyxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQzVDLGdCQUFnQixJQUFJLElBQUksR0FBRyxDQUFDLENBQUM7QUFDN0IsZ0JBQWdCLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxHQUFHLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDOUMsb0JBQW9CLE1BQU0sS0FBSyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDO0FBQ2hFLG9CQUFvQixJQUFJLElBQUksS0FBSyxHQUFHLEtBQUssQ0FBQztBQUMxQyxpQkFBaUI7QUFDakIsZ0JBQWdCLE1BQU0sRUFBRSxHQUFHLENBQUMsSUFBSSxDQUFDLEdBQUcsSUFBSSxDQUFDLENBQUM7QUFDMUMsZ0JBQWdCLEVBQUUsQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxFQUFFLENBQUMsQ0FBQztBQUN2QyxnQkFBZ0IsRUFBRSxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxDQUFDO0FBQ3ZDLGdCQUFnQixJQUFJLElBQUksQ0FBQyxHQUFHLEVBQUUsQ0FBQztBQUMvQixhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0E7QUFDQSxRQUFRLE1BQU0sQ0FBQyxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxFQUFDO0FBQ3JDLFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUNwQyxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQzVDLGdCQUFnQixNQUFNLEdBQUcsR0FBRyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLElBQUksRUFBRSxNQUFNLENBQUMsQ0FBQztBQUNwRSxnQkFBZ0IsQ0FBQyxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxDQUFDO0FBQ3ZDLGdCQUFnQixDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsR0FBRyxDQUFDLENBQUM7QUFDdkMsYUFBYTtBQUNiLFNBQVM7QUFDVDtBQUNBLFFBQVEsTUFBTSxJQUFJLEdBQUcsSUFBSSxNQUFNLENBQUMsQ0FBQyxFQUFFLEdBQUcsRUFBRSxPQUFPLENBQUMsQ0FBQztBQUNqRCxRQUFRLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDcEMsWUFBWSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3hDLGdCQUFnQixNQUFNLE9BQU8sR0FBRyxDQUFDLElBQUksSUFBSSxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLEdBQUcsRUFBRSxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDNUYsZ0JBQWdCLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxHQUFHLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDOUMsb0JBQW9CLElBQUksQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxJQUFJLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsR0FBRyxPQUFPLElBQUksQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3ZHLGlCQUFpQjtBQUNqQixhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0E7QUFDQSxRQUFRLElBQUksS0FBSyxHQUFHLElBQUksWUFBWSxDQUFDLEdBQUcsQ0FBQyxDQUFDO0FBQzFDLFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUNwQyxZQUFZLEtBQUssSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxHQUFHLEVBQUUsRUFBRSxDQUFDLEVBQUU7QUFDMUMsZ0JBQWdCLE1BQU0sR0FBRyxHQUFHLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDO0FBQzdDLGdCQUFnQixNQUFNLEdBQUcsR0FBRyxLQUFLLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUM5QyxnQkFBZ0IsTUFBTSxNQUFNLEdBQUcsS0FBSyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDakQ7QUFDQSxnQkFBZ0IsSUFBSSxPQUFPLEdBQUcsSUFBSSxDQUFDLElBQUksQ0FBQyxHQUFHLENBQUMsS0FBSyxJQUFJLENBQUMsSUFBSSxDQUFDLEdBQUcsQ0FBQyxHQUFHLE1BQU0sR0FBRyxFQUFFLEdBQUcsTUFBTSxHQUFHLEVBQUUsQ0FBQztBQUM1RixnQkFBZ0IsSUFBSSxPQUFPLEdBQUcsR0FBRyxFQUFFLE9BQU8sR0FBRyxHQUFHLENBQUM7QUFDakQsZ0JBQWdCLEtBQUssQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxPQUFPLENBQUMsQ0FBQztBQUMvQztBQUNBLGdCQUFnQixNQUFNLE1BQU0sR0FBRyxJQUFJLEdBQUcsR0FBRyxHQUFHLEVBQUUsR0FBRyxFQUFFLENBQUM7QUFDcEQsZ0JBQWdCLE1BQU0sTUFBTSxHQUFHLE1BQU0sR0FBRyxHQUFHLEdBQUcsT0FBTyxHQUFHLE9BQU8sR0FBRyxHQUFHLENBQUM7QUFDdEUsZ0JBQWdCLEtBQUssQ0FBQyxTQUFTLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxNQUFNLENBQUMsQ0FBQztBQUM5QztBQUNBLGdCQUFnQixDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxDQUFDLEtBQUssQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLEdBQUcsTUFBTSxDQUFDLENBQUM7QUFDMUQsZ0JBQWdCLEtBQUssQ0FBQyxDQUFDLENBQUMsSUFBSSxDQUFDLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQztBQUMxQyxhQUFhO0FBQ2IsU0FBUztBQUNUO0FBQ0EsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3BDLFlBQVksS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUN4QyxnQkFBZ0IsQ0FBQyxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLENBQUMsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLEtBQUssQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLEVBQUM7QUFDL0QsYUFBYTtBQUNiLFNBQVM7QUFDVDtBQUNBLFFBQVEsT0FBTyxJQUFJLENBQUMsQ0FBQyxDQUFDO0FBQ3RCLEtBQUs7QUFDTDs7QUM1T0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDZSxlQUFRLEVBQUUsQ0FBQyxFQUFFLEVBQUUsRUFBRSxRQUFRLEdBQUcsR0FBRyxFQUFFO0FBQ2hELElBQUksTUFBTSxPQUFPLEdBQUcsSUFBSSxDQUFDO0FBQ3pCLElBQUksTUFBTSxDQUFDLEdBQUcsRUFBRSxDQUFDLE1BQU0sQ0FBQztBQUN4QixJQUFJLElBQUksS0FBSyxHQUFHLElBQUksQ0FBQztBQUNyQixJQUFJLElBQUksR0FBRyxHQUFHLEtBQUssQ0FBQztBQUNwQixJQUFJLElBQUksQ0FBQyxHQUFHLEVBQUUsQ0FBQyxLQUFLLEVBQUUsQ0FBQztBQUN2QixJQUFJLElBQUksRUFBRSxHQUFHLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUNsQixJQUFJLElBQUksV0FBVyxHQUFHLEtBQUssQ0FBQztBQUM1QjtBQUNBLElBQUksT0FBTyxRQUFRLEVBQUUsSUFBSSxDQUFDLElBQUksQ0FBQyxXQUFXLEVBQUU7QUFDNUMsUUFBUSxXQUFXLEdBQUcsSUFBSSxDQUFDO0FBQzNCLFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLENBQUMsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUNwQyxZQUFZLENBQUMsQ0FBQyxDQUFDLENBQUMsSUFBSSxJQUFJLENBQUM7QUFDekIsWUFBWSxJQUFJLEdBQUcsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDM0IsWUFBWSxDQUFDLENBQUMsQ0FBQyxDQUFDLElBQUksSUFBSSxDQUFDO0FBQ3pCLFlBQVksSUFBSSxFQUFFLEdBQUcsQ0FBQyxHQUFHLEdBQUcsRUFBRSxJQUFJLElBQUksQ0FBQztBQUN2QyxZQUFZLElBQUksSUFBSSxDQUFDLEdBQUcsQ0FBQyxFQUFFLENBQUMsR0FBRyxPQUFPLEVBQUU7QUFDeEMsZ0JBQWdCLFdBQVcsR0FBRyxLQUFLLENBQUM7QUFDcEMsYUFBYTtBQUNiLFlBQVksQ0FBQyxDQUFDLENBQUMsQ0FBQyxJQUFJLEtBQUssR0FBRyxFQUFFLENBQUM7QUFDL0IsWUFBWSxFQUFFLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3RCLFNBQVM7QUFDVCxRQUFRLEtBQUssSUFBSSxHQUFHLElBQUksRUFBRSxHQUFHLElBQUksR0FBRyxHQUFHLENBQUM7QUFDeEMsUUFBUSxHQUFHLEdBQUcsRUFBRSxDQUFDO0FBQ2pCLEtBQUs7QUFDTCxJQUFJLE9BQU8sQ0FBQyxDQUFDO0FBQ2I7O0FDMUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTyxNQUFNLElBQUksU0FBUyxFQUFFLENBQUM7QUFDN0I7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSSxXQUFXLENBQUMsQ0FBQyxFQUFFLFdBQVcsQ0FBQyxFQUFFLEVBQUUsa0JBQWtCLENBQUMsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxDQUFDLENBQUMsRUFBRSxNQUFNLENBQUMsU0FBUyxFQUFFLElBQUksQ0FBQyxJQUFJLEVBQUU7QUFDdkcsUUFBUSxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsRUFBRSxNQUFNLEVBQUUsSUFBSSxFQUFDO0FBQ2pDLFFBQVEsS0FBSyxDQUFDLGNBQWMsR0FBRyxDQUFDLGFBQWEsRUFBRSxvQkFBb0IsRUFBRSxVQUFVLENBQUMsQ0FBQztBQUNqRixRQUFRLEVBQUUsSUFBSSxDQUFDLEVBQUUsRUFBRSxJQUFJLENBQUMsRUFBRSxFQUFFLEdBQUcsSUFBSSxDQUFDLENBQUMsQ0FBQyxLQUFLLENBQUM7QUFDNUMsUUFBUSxXQUFXLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsRUFBRSxHQUFHLENBQUMsRUFBRSxXQUFXLENBQUMsQ0FBQztBQUN6RCxRQUFRLElBQUksQ0FBQyxTQUFTLENBQUMsYUFBYSxFQUFFLFdBQVcsQ0FBQyxDQUFDO0FBQ25ELFFBQVEsSUFBSSxDQUFDLFNBQVMsQ0FBQyxvQkFBb0IsRUFBRSxJQUFJLENBQUMsR0FBRyxDQUFDLGtCQUFrQixFQUFFLFdBQVcsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQzVGLFFBQVEsSUFBSSxDQUFDLFNBQVMsQ0FBQyxVQUFVLEVBQUUsUUFBUSxDQUFDLENBQUM7QUFDN0MsUUFBUSxJQUFJLENBQUMsS0FBSyxHQUFHLENBQUMsQ0FBQztBQUN2QixRQUFRLElBQUksQ0FBQyxPQUFPLEdBQUcsQ0FBQyxDQUFDO0FBQ3pCLFFBQVEsSUFBSSxDQUFDLGlCQUFpQixHQUFHLENBQUMsQ0FBQztBQUNuQyxRQUFRLElBQUksQ0FBQyxtQkFBbUIsR0FBRyxDQUFDLENBQUM7QUFDckMsUUFBUSxJQUFJLENBQUMscUJBQXFCLEdBQUcsQ0FBQyxDQUFDO0FBQ3ZDLFFBQVEsSUFBSSxDQUFDLFNBQVMsR0FBRyxHQUFHLENBQUM7QUFDN0IsUUFBUSxJQUFJLENBQUMsY0FBYyxHQUFHLENBQUMsQ0FBQztBQUNoQyxRQUFRLElBQUksQ0FBQyxDQUFDLEdBQUcsSUFBSSxNQUFNLENBQUMsSUFBSSxDQUFDLEVBQUUsRUFBRSxJQUFJLENBQUMsRUFBRSxFQUFFLE1BQU0sSUFBSSxDQUFDLFdBQVcsQ0FBQyxNQUFNLENBQUMsQ0FBQztBQUM3RSxRQUFRLE9BQU8sSUFBSSxDQUFDO0FBQ3BCLEtBQUs7QUFDTDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksZUFBZSxDQUFDLE1BQU0sRUFBRSxRQUFRLEVBQUU7QUFDdEMsUUFBUSxNQUFNLEtBQUssR0FBRyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsQ0FBQyxLQUFLLENBQUMsSUFBSSxDQUFDLEdBQUcsQ0FBQyxHQUFHLElBQUksQ0FBQyxHQUFHLENBQUMsQ0FBQyxFQUFFLENBQUMsR0FBRyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3BFLFFBQVEsTUFBTSxFQUFFLEdBQUcsUUFBUSxDQUFDLENBQUMsRUFBRSxNQUFNLEdBQUcsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxDQUFDO0FBQ2hELFFBQVEsTUFBTSxFQUFFLEdBQUcsUUFBUSxDQUFDLENBQUMsRUFBRSxNQUFNLEdBQUcsQ0FBQyxFQUFFLEdBQUcsQ0FBQyxDQUFDO0FBQ2hEO0FBQ0EsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsRUFBRSxDQUFDLE1BQU0sRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ25ELFlBQVksTUFBTSxJQUFJLEdBQUcsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQy9CLFlBQVksRUFBRSxDQUFDLENBQUMsQ0FBQyxJQUFJLElBQUksR0FBRyxRQUFRLEdBQUcsQ0FBQyxHQUFHLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxJQUFJLEdBQUcsUUFBUSxDQUFDLEdBQUcsTUFBTSxDQUFDLENBQUMsQ0FBQztBQUNsRixTQUFTO0FBQ1Q7QUFDQSxRQUFRLE1BQU0sR0FBRyxHQUFHLENBQUMsQ0FBQyxLQUFLO0FBQzNCLFlBQVksTUFBTSxLQUFLLEdBQUcsUUFBUSxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxLQUFLLEVBQUUsQ0FBQyxDQUFDLENBQUMsR0FBRyxLQUFLLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQzNGLFlBQVksT0FBTyxJQUFJLENBQUMsSUFBSSxDQUFDLFdBQVcsQ0FBQyxLQUFLLENBQUMsR0FBRyxDQUFDLENBQUMsSUFBSSxDQUFDLEdBQUcsQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ2pFLFVBQVM7QUFDVDtBQUNBLFFBQVEsT0FBTyxNQUFNLENBQUMsR0FBRyxFQUFFLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDbkMsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FB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SyxDQUFDLE9BQU8sQ0FBQyxDQUFDLENBQUMsQ0FBQyxFQUFFLENBQUMsT0FBTyxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxDQUFDO0FBQ3BELFlBQVksSUFBSSxVQUFVLEdBQUcsT0FBTyxDQUFDLENBQUMsQ0FBQyxFQUFFLFVBQVUsR0FBRyxPQUFPLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDakUsU0FBUztBQUNULFFBQVEsSUFBSSxZQUFZLEdBQUcsQ0FBQyxRQUFRLENBQUM7QUFDckMsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsVUFBVSxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQzdDLFlBQVksT0FBTyxDQUFDLENBQUMsQ0FBQyxJQUFJLFVBQVUsQ0FBQztBQUNyQyxZQUFZLE9BQU8sQ0FBQyxDQUFDLENBQUMsSUFBSSxLQUFLLENBQUM7QUFDaEMsWUFBWSxPQUFPLENBQUMsQ0FBQyxDQUFDLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQyxDQUFDLEdBQUcsVUFBVSxHQUFHLE9BQU8sQ0FBQyxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQy9ELFlBQVksSUFBSSxZQUFZLEdBQUcsT0FBTyxDQUFDLENBQUMsQ0FBQyxFQUFFLFlBQVksR0FBRyxPQUFPLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDckUsU0FBUztBQUNULFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLFVBQVUsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM3QyxZQUFZLE9BQU8sQ0FBQyxDQUFDLENBQUMsSUFBSSxZQUFZLENBQUM7QUFDdkMsU0FBUztBQUNULFFBQVEsT0FBTztBQUNmLFlBQVksVUFBVSxFQUFFLFFBQVE7QUFDaEMsWUFBWSxTQUFTLEVBQUUsT0FBTztBQUM5QixTQUFTO0FBQ1QsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUksT0FBTyxDQUFDLGFBQWEsRUFBRSxHQUFHLEVBQUUsSUFBSSxF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O0FBQ0EsSUFBSSxLQUFLLENBQUMsQ0FBQyxFQUFFO0FBQ2IsUUFBUSxNQUFNLFNBQVMsR0FBRyxJQUFJLENBQUMsU0FBUyxDQUFDO0FBQ3pDLFFBQVEsTUFBTSxVQUFVLEdBQUcsSUFBSSxDQUFDLFVBQVUsQ0FBQztBQUMzQyxRQUFRLE1BQU0sUUFBUSxHQUFHLElBQUksQ0FBQyxRQUFRLENBQUM7QUFDdkMsUUFBUSxNQUFNLE9BQU8sR0FBRyxJQUFJLENBQUMsT0FBTyxDQUFDO0FBQ3JDLFFBQVEsTUFBTSxDQUFDLENBQUMsRUFBRSxHQUFHLENBQUMsR0FBRyxDQUFDLENBQUMsS0FBSyxDQUFDO0FBQ2pDLFFBQVEsTUFBTSxVQUFVLEdBQUcsUUFBUSxDQUFDLEtBQUssQ0FBQyxDQUFDLENBQUMsQ0FBQztBQUM3QyxRQUFRLE1BQU0sSUFBSSxHQUFHLElBQUksTUFBTSxDQUFDLENBQUMsRUFBRSxHQUFHLEVBQUUsQ0FBQyxDQUFDLENBQUM7QUFDM0MsUUFBUSxJQUFJLElBQUksR0FBRyxJQUFJLEtBQUssQ0FBQyxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDMUMsUUFBUSxJQUFJLElBQUksR0FBRyxJQUFJLEtBQUssQ0FBQyxHQUFHLENBQUMsQ0FBQyxJQUFJLENBQUMsQ0FBQyxDQUFDLENBQUM7QUFDMUMsUUFBUSxJQUFJLElBQUksR0FBRyxDQUFDLENBQUM7QUFDckIsUUFBUSxJQUFJLElBQUksR0FBRyxDQUFDLENBQUM7QUFDckIsUUFBUSxJQUFJLE1BQU0sR0FBRyxDQUFDLENBQUM7QUFDdkIsUUFBUSxJQUFJLElBQUksR0FBRyxDQUFDLENBQUM7QUFDckIsUUFBUSxNQUFNLGNBQWMsR0FBRyxDQUFDLEdBQUcsU0FBUyxHQUFHLFVBQVUsQ0FBQztBQUMxRDtBQUNBLFFBQVEsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLFVBQVUsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUM3QyxZQUFZLE1BQU0sQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFF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QUM7QUFDN0IsUUFBUSxNQUFNLElBQUksR0FBRyxJQUFJLENBQUMsSUFBSSxDQUFDO0FBQy9CLFFBQVEsTUFBTSxHQUFHLEdBQUcsSUFBSSxDQUFDLEdBQUcsQ0FBQztBQUM3QixRQUFRLE1BQU0sRUFBRSxHQUFHLElBQUksQ0FBQyxFQUFFLENBQUM7QUFDM0IsUUFBUSxLQUFLLElBQUksQ0FBQyxHQUFHLENBQUMsRUFBRSxDQUFDLEdBQUcsQ0FBQyxFQUFFLEVBQUUsQ0FBQyxFQUFFO0FBQ3BDLFlBQVksS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLEdBQUcsRUFBRSxFQUFFLENBQUMsRUFBRTtBQUMxQyxnQkFBZ0IsTUFBTSxRQUFRLEdBQUcsQ0FBQyxJQUFJLENBQUMsSUFBSSxDQUFDLEdBQUcsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxDQUFDLElBQUksSUFBSSxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxJQUFJLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLEVBQUUsR0FBRyxJQUFJLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLEVBQUUsRUFBRSxRQUFRLENBQUMsQ0FBQztBQUNqSyxnQkFBZ0IsSUFBSSxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUUsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDO0FBQy9DLGdCQUFnQixHQUFHLENBQUMsU0FBUyxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsS0FBSyxHQUFHLEdBQUcsQ0FBQyxLQUFLLENBQUMsQ0FBQyxFQUFFLENBQUMsQ0FBQyxHQUFHLEVBQUUsR0FBRyxJQUFJLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsR0FBRyxJQUFJLENBQUMsS0FBSyxDQUFDLENBQUMsRUFBRSxDQUFDLENBQUMsQ0FBQyxDQUFDO0FBQ3hHLGdC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
