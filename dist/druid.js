// https://renecutura.eu v0.6.0 Copyright 2022 Rene Cutura
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.druid = global.druid || {}));
})(this, (function (exports) { 'use strict';

/**
 * Computes the euclidean distance (<code>l<sub>2</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the euclidean distance between <code>a</code> and <code>b</code>.
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
 * @param {Number[]} summands - Array of values to sum up.
 * @returns {Number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
function neumair_sum (summands) {
    const n = summands.length;
    let sum = 0;
    let compensation = 0;

    for (let i = 0; i < n; ++i) {
        const summand = summands[i];
        const t = sum + summand;
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
 * Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the squared euclidean distance between <code>a</code> and <code>b</code>.
 */
function euclidean_squared (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    const s = new Float64Array(n);
    for (let i = 0; i < n; ++i) {
        const x = a[i];
        const y = b[i];
        const x_y = x - y;
        s[i] = x_y * x_y;
    }
    return neumair_sum(s);
}

/**
 * Computes the cosine distance (not similarity) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} The cosine distance between {@link a} and {@link b}.
 * 
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 * 
 * druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4;
 * 
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
 * Computes the manhattan distance (<code>l<sub>1</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the manhattan distance between <code>a</code> and <code>b</code>.
 */ 
function manhattan (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
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
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the chebyshev distance between {@link a} and {@link b}.
 */
function chebyshev (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res);
}

/**
 * Computes the canberra distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias canberra
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the canberra distance between <code>a</code> and <code>b</code>.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
function canberra(a, b) {
    if (a.length !== b.length) return undefined;
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += (Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i])));
    }
    return sum;
}

/**
 * Computes the jaccard distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias jaccard
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the jaccard distance between <code>a</code> and <code>b</code>.
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
 * Computes the hamming distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias hamming
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the hamming distance between <code>a</code> and <code>b</code>.
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
 * Computes the Sokal-Michener distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias sokal_michener
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the Sokal-Michener distance between <code>a</code> and <code>b</code>.  
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
 * Computes the yule distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias yule
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the yule distance between <code>a</code> and <code>b</code>.
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
 * Computes the k-nearest neighbors of each row of {@link A}.
 * @memberof module:matrix
 * @alias k_nearest_neigbhors
 * @param {Matrix} A - Either the data matrix, or a distance matrix.
 * @param {Number} k - The number of neighbors to compute.
 * @param {Function|"precomputed"} [metric=euclidean]
 * @returns {Array<Object>} -
 */
function k_nearest_neighbors (A, k, metric = euclidean) {
    const rows = A.shape[0];
    let D = metric == "precomputed" ? A : distance_matrix(A, metric);
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
 * Computes the distance matrix of datamatrix {@link A}.
 * @memberof module:matrix
 * @alias distance_matrix
 * @param {Matrix} A - Matrix.
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
 * Creates an Array containing {@link number} numbers from {@link start} to {@link end}.
 * If <code>{@link number} = null</null>.
 * @memberof module:matrix
 * @alias linspace
 * @param {Number} start - Start value.
 * @param {Number} end - End value.
 * @param {Number} [number = null] - Number of number between {@link start} and {@link end}.
 * @returns {Array} - An array with {@link number} entries, beginning at {@link start} ending at {@link end}.
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

/**
 * Computes the norm of a vector, by computing its distance to **0**.
 * @memberof module:matrix
 * @alias norm
 * @param {Matrix|Array<Number>|Float64Array} v - Vector.
 * @param {Function} [metric = euclidean] - Which metric should be used to compute the norm.
 * @returns {Number} - The norm of {@link v}.
 */
function norm (v, metric = euclidean) {
    let vector = null;
    if (v instanceof Matrix) {
        let [rows, cols] = v.shape;
        if (rows === 1) vector = v.row(0);
        else if (cols === 1) vector = v.col(0);
        else throw new Error("Matrix must be 1d!");
    } else {
        vector = v;
    }
    const n = vector.length;
    const zeros = new Float64Array(n);
    return metric(vector, zeros);
}

/**
 * Normalizes Vector {@link v}.
 * @memberof module:matrix
 * @alias normalize
 * @param {Array<Number>|Float64Array} v - Vector
 * @param {Function} metric 
 * @returns {Array<Number>|Float64Array} - The normalized vector with length 1.
 */
function normalize(v, metric = euclidean)  {
    const v_norm = norm(v, metric);
    return v.map(value => value / v_norm);
}

/**
 * Computes the QR Decomposition of the Matrix {@link A} using Gram-Schmidt process.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */
function qr_gramschmidt (A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, cols, "identity");
    const R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        let v = A.col(j);
        for (let i = 0; i < j; ++i) {
            const q = Q.col(i);
            const q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            R.set_entry(i, j, q_dot_v);
            v = v.map((v_, k) => v_ - q_dot_v * q[k]);
        }
        const v_norm = norm(v, euclidean);
        for (let k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j, j, v_norm);
    }
    return { R, Q };
}

/**
 * Computes the QR Decomposition of the Matrix {@link A} with householder transformations.
 * @memberof module:linear_algebra
 * @alias qr_householder
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */
function qr_householder (A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, rows, "I");
    const R = A.clone();

    for (let j = 0; j < cols; ++j) {
        const x = Matrix.from(R.col(j).slice(j));
        const x_norm = norm(x);
        const x0 = x.entry(0, 0);
        const rho = -Math.sign(x0);
        const u1 = x0 - rho * x_norm;
        const u = x.divide(u1).set_entry(0, 0, 1);
        const beta = (-rho * u1) / x_norm;

        const u_outer_u = u.outer(u);
        const R_block = R.get_block(j, 0);
        const new_R = R_block.sub(u_outer_u.dot(R_block).mult(beta));
        const Q_block = Q.get_block(0, j);
        const new_Q = Q_block.sub(Q_block.dot(u_outer_u).mult(beta));
        R.set_block(j, 0, new_R);
        Q.set_block(0, j, new_Q);
    }
    return { R, Q };
}

/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @memberof module:linear_algebra
 * @alias simultaneous_poweriteration
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Object} parameters - Object containing parameterization of the simultanious poweriteration method.
 * @param {Number} [parameters.max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [parameters.seed=1212] - The seed value or a randomizer used in the algorithm.
 * @param {Function} [parameters.qr=qr_gramschmidt] - The QR technique to use.
 * @param {Number} [parameters.tol=1e-8] - Allowed error for stopping criteria
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */
function simultaneous_poweriteration (A, k = 2, {seed = 1212, max_iterations = 100, qr = qr_gramschmidt, tol = 1e-8} = {}) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q, R } = qr(new Matrix(n, k, () => (randomizer.random - .5) * 2));
    while (max_iterations--) {
        const oldQ = Q.clone();
        const Z = A.dot(Q);
        const QR = qr(Z);
        Q = QR.Q;
        R = QR.R;
        const error = euclidean_squared(Q.values, oldQ.values);
        if (error < tol) {
            break;
        }
    }

    const eigenvalues = R.diag;
    const eigenvectors = Q.transpose().to2dArray;
    return { eigenvalues, eigenvectors };
}

/**
 * Computes the inner product between two arrays of the same length.
 * @memberof module:linear_algebra
 * @alias inner_product
 * @param {Array|Float64Array} a - Array a
 * @param {Array|Float64Array} b - Array b
 * @returns The inner product between {@link a} and {@link b}
 */
function inner_product (a, b) {
    const N = a.length;
    if (N != b.length) {
        throw new Error("Array a and b must have the same length!")
    }
    let sum = 0;
    for (let i = 0; i < N; ++i) {
        sum += a * b;
    }
    return sum;
}

/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
class Matrix {
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
     * @param {Number} row
     * @param {Array} values
     * @returns {Matrix}
     */
    set_row(row, values) {
        const cols = this._cols;
        if ((Array.isArray(values) || values instanceof Float64Array) && values.length === cols) {
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
     * @param {Function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {Function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
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
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
class DR {
    /**
     * Takes the default parameters and seals them, remembers the type of input {@link X}, and initializes the random number generator.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X, default_parameters, parameters) {
        this._parameters = Object.assign(Object.seal(default_parameters), parameters);
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("No valid type for X!");
        }
        [this._N, this._D] = this.X.shape;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._is_initialized = false;
        return this;
    }

    /**
     * Set and get parameters
     * @param {String} [name = null] - Name of the parameter. If not given then returns all parameters as an Object.
     * @param {any} [value = null] - Value of the parameter to set. If <code>name</code> is set and <code>value</code> is not given, returns the value of the respective parameter.
     * @returns {DR|any|Object} 
     * On setting a parameter, this function returns the DR object. 
     * If <code>name</code> is set and <code>value == null</code> then return actual parameter value.
     * If <code>name</code> is not given, then returns all parameters as an Object.
     * 
     * @example
     * '''
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
     * DR.parameter("d"); // returns 3,
     * DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
     * '''
     */
    parameter(name = null, value = null) {
        if (name === null) {
            return Object.assign({}, this._parameters);
        }
        if (!this._parameters.hasOwnProperty(name)) {
            throw new Error(`${name} is not a valid parameter!`);
        }
        if (value !== null) {
            this._parameters[name] = value;
            this._is_initialized = false;
            return this;
        } else {
            return this._parameters[name];
        }
    }

    para(name = null, value = null) {
        return this.parameter(name, value);
    }

    p(name = null, value = null) {
        return this.parameter(name, value);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} the projection.
     */
    transform() {
        this.check_init();
        return this.projection;
    }

    /**
     * Computes the projection.
     * @yields {Matrix|Number[][]} the intermediate steps of the projection.
     */
    *generator() {
        return this.transform();
    }

    /**
     * If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.
     * @returns {DR}
     */
    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
        return this;
    }

    /**
     * @returns {Matrix|Number[][]} the projection in the type of input <code>X</code>.
     */
    get projection() {
        if (this.hasOwnProperty("Y")) {
            this.check_init();
            return this._type === "matrix" ? this.Y : this.Y.to2dArray;
        } else {
            throw new Error("The dataset is not transformed yet!");
        }
    }

    /**
     * Computes the projection.
     * @param  {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<Matrix|Number[][]>} the dimensionality reduced dataset.
     */
    async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Matrix|Array} the dimensionality reduced dataset.
     */
    static transform(...args) {
        let dr = new this(...args);
        return dr.transform();
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise} a promise yielding the dimensionality reduced dataset.
     */
    static async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator} a generator yielding the intermediate steps of the dimensionality reduction method.
     */
    static *generator(...args) {
        const dr = new this(...args);
        const generator = dr.generator();
        for (const result of generator) {
            yield result;
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
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @returns {PCA}
     */
    constructor(X, parameters) {
        super(X, { d: 2, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
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
            return Matrix.from(A).dot(V).asArray;
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
        const { d, eig_args } = this._parameters;
        const X = this.X;
        const means = Matrix.from(X.meanCols);
        const X_cent = X.sub(means);
        const C = X_cent.transpose().dot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration(C, d, eig_args);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }

    static principal_components(X, parameters) {
        const dr = new this(X, parameters);
        return dr.principal_components();
    }
}

/**
 * @class
 * @alias MDS
 * @extends DR
 */
class MDS extends DR {
    /**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|Array}
     */
    transform() {
        const X = this.X;
        const rows = X.shape[0];
        const { d, metric, eig_args } = this._parameters;
        const A = metric === "precomputed" ? X : distance_matrix(X, metric);
        const ai_ = A.meanCols;
        const a_j = A.meanRows;
        const a__ = A.mean;

        this._d_X = A;
        const B = new Matrix(rows, rows, (i, j) => A.entry(i, j) - ai_[i] - a_j[j] + a__);

        const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
        this.Y = Matrix.from(V).transpose();

        return this.projection;
    }

    /**
     * @returns {Number} - the stress of the projection.
     */
    stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X;
        const d_Y = new Matrix();
        d_Y.shape = [
            N,
            N,
            (i, j) => {
                return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
            },
        ];
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
     * Isometric feature mapping (ISOMAP).
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(this._parameters.neighbors ?? Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
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
        const { d, metric, eig_args, neighbors } = this._parameters;
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix();
        D.shape = [rows, rows, (i, j) => (i <= j ? metric(X.row(i), X.row(j)) : D.entry(j, i))];
        const kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            const row = [];
            for (let j = 0; j < rows; ++j) {
                row.push({
                    index: j,
                    distance: D.entry(i, j),
                });
            }
            const H = new Heap(row, (d) => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, neighbors + 1));
        }

        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
        const G = new Matrix(rows, rows, (i, j) => {
            const other = kNearestNeighbors[i].find((n) => n.index === j);
            return other ? other.distance : Infinity;
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
        const A = new Matrix(rows, rows, (i, j) => {
            let val = G.entry(i, j);
            val = val === Infinity ? 0 : val;
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        });

        ai_ = ai_.map((v) => v / rows);
        a_j = a_j.map((v) => v / rows);
        a__ /= rows ** 2;
        const B = new Matrix(rows, rows, (i, j) => A.entry(i, j) - ai_[i] - a_j[j] + a__);

        // compute d eigenvectors
        const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
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
class FASTMAP extends DR {
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212 }, parameters);
        return this;
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {Function} dist
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */
    _choose_distant_objects(dist) {
        const X = this.X;
        const N = X.shape[0];
        let a_index = (this._randomizer.random_int % N) - 1;
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
        const { d, metric } = this._parameters;
        const Y = new Matrix(N, d, 0);
        let dist = (a, b) => metric(X.row(a), X.row(b));

        for (let _col = 0; _col < d; ++_col) {
            let old_dist = dist;
            // choose pivot objects
            const [a_index, b_index, d_ab] = this._choose_distant_objects(dist);
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
        // return embedding.
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
     * Linear Discriminant Analysis.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - The high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Array} parameters.labels - The labels / classes for each data point.
     * @param {number} [parameters.d = 2] - The dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
    constructor(X, parameters) {
        super(X, { labels: null, d: 2, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const [rows, cols] = X.shape;
        const { d, labels, eig_args } = this._parameters;
        if (labels === null || labels.length != rows) {
            throw new Error("LDA needs parameter label to every datapoint to work!");
        }
        const unique_labels = {};
        let label_id = 0;
        labels.forEach((l, i) => {
            if (l in unique_labels) {
                unique_labels[l].count++;
                unique_labels[l].rows.push(X.row(i));
            } else {
                unique_labels[l] = {
                    id: label_id++,
                    count: 1,
                    rows: [X.row(i)],
                };
            }
        });

        // create X_mean and vector means;
        const X_mean = X.mean;
        const V_mean = new Matrix(label_id, cols);
        for (const label in unique_labels) {
            const V = Matrix.from(unique_labels[label].rows);
            const v_mean = V.meanCols;
            for (let j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }
        }
        // scatter_between
        let S_b = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const m = new Matrix(cols, 1, (j) => v[j] - X_mean);
            const N = unique_labels[label].count;
            S_b = S_b.add(m.dot(m.transpose()).mult(N));
        }

        // scatter_within
        let S_w = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const m = new Matrix(cols, 1, (j) => v[j]);
            const R = unique_labels[label].rows;
            for (let i = 0, n = unique_labels[label].count; i < n; ++i) {
                const row_v = new Matrix(cols, 1, (j, _) => R[i][j] - m.entry(j, 0));
                S_w = S_w.add(row_v.dot(row_v.transpose()));
            }
        }

        let { eigenvectors: V } = simultaneous_poweriteration(S_w.inverse().dot(S_b), d, eig_args);
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
     * Locally Linear Embedding.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const rows = this._N;
        const cols = this._D;
        const { neighbors, d, eig_args, metric } = this._parameters;
        const nN = k_nearest_neighbors(X, neighbors, metric);
        const O = new Matrix(neighbors, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const nN_row = nN[row];
            const Z = new Matrix(neighbors, cols, (i, j) => X.entry(nN_row[i].j, j) - X.entry(row, j));
            const C = Z.dot(Z.T);
            if (neighbors > cols) {
                const C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < neighbors; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }
            // reconstruct;
            let w = Matrix.solve_CG(C, O, this._randomizer);
            w = w.divide(w.sum);
            for (let j = 0; j < neighbors; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        const I = new Matrix(rows, rows, "identity");
        const IW = I.sub(W);
        const M = IW.T.dot(IW);
        const { eigenvectors: V } = simultaneous_poweriteration(M.T.inverse(), d + 1, eig_args);
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
     * Local Tangent Space Alignment
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link LTSA} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        if (this._D <= this.parameter("d")) {
            throw new Error(`Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${this.parameter("d")})!`);
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const [rows, D] = X.shape;
        const { d, neighbors, metric, eig_args } = this._parameters;
        // 1.1 determine k nearest neighbors
        const nN = k_nearest_neighbors(X, neighbors, metric);
        // center matrix
        const O = new Matrix(D, D, "center");
        const B = new Matrix(rows, rows, 0);

        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            const I_i = [row, ...nN[row].map((n) => n.j)];
            let X_i = Matrix.from(I_i.map((n) => X.row(n)));
            // center X_i
            X_i = X_i.dot(O);
            // correlation matrix
            const C = X_i.dot(X_i.transpose());
            const { eigenvectors: g } = simultaneous_poweriteration(C, d, eig_args);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            const G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            const W_i = G_i_t.transpose()
                .dot(G_i_t)
                .add(1 / Math.sqrt(neighbors + 1));
            for (let i = 0; i < neighbors + 1; ++i) {
                for (let j = 0; j < neighbors + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0) + W_i.entry(i, j));
                }
            }
        }

        // 3. Aligning global coordinates
        const { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1, eig_args);
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
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.perplexity = 50] - perplexity.
     * @param {Number} [parameters.epsilon = 10] - learning parameter.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TSNE}
     */
    constructor(X, parameters) {
        super(X, { perplexity: 50, epsilon: 10, d: 2, metric: euclidean, seed: 1212 }, parameters);
        [this._N, this._D] = this.X.shape;
        this._iter = 0;
        this.Y = new Matrix(this._N, this.parameter("d"), () => this._randomizer.random);
        return this;
    }

    /**
     *
     * @returns {TSNE}
     */
    init() {
        // init
        const Htarget = Math.log(this.parameter("perplexity"));
        const N = this._N;
        const D = this._D;
        const {metric} = this._parameters;
        const X = this.X;
        let Delta;
        if (metric =="precomputed") {
            Delta = druid.Matrix.from(X);
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
        let prow = new Float64Array(N);
        const tol = 1e-4;
        const maxtries = 50;
        for (let i = 0; i < N; ++i) {
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let done = false;

            let num = 0;
            while (!done) {
                let psum = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = Math.exp(-Delta.entry(i, j) * beta);
                    if (i === j) pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                let Hhere = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = psum === 0 ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) {
                        Hhere -= pj * Math.log(pj);
                    }
                }
                if (Hhere > Htarget) {
                    betamin = beta;
                    beta = betamax === Infinity ? beta * 2 : (beta + betamax) / 2;
                } else {
                    betamax = beta;
                    beta = betamin === -Infinity ? beta / 2 : (beta + betamin) / 2;
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
     * @param {Number} [iterations=500] - Number of iterations.
     * @returns {Matrix|Number[][]} the projection.
     */
    transform(iterations = 500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     *
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Number[][]} - the projection.
     */
    *generator(iterations = 500) {
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
        const { d: dim, epsilon} = this._parameters;
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

                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * 0.8 : gainid + 0.2;
                if (newgain < 0.01) newgain = 0.01;
                gains.set_entry(i, d, newgain);

                const momval = iter < 250 ? 0.5 : 0.8;
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
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [parameters.local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [parameters.min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.
     * @param {Number} [parameters._spread = 1] - The effective scale of embedded points. (In combination with {@link parameters.min_dist})
     * @param {Number} [parameters._set_op_mix_ratio = 1] - Interpolate between union and intersection.
     * @param {Number} [parameters._repulsion_strength = 1]  - Weighting applied to negative samples.
     * @param {Number} [parameters._negative_sample_rate = 5] - The number of negative samples per positive sample.
     * @param {Number} [parameters._n_epochs = 350] - The number of training epochs.
     * @param {Number} [parameter._initial_alpha = 1] - The initial learning rate for the optimization.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {UMAP}
     */
    constructor(X, parameters) {
        super(X, { n_neighbors: 15, local_connectivity: 1, min_dist: 1, d: 2, metric: euclidean, seed: 1212, _spread: 1, _set_op_mix_ratio: 1, _repulsion_strength: 1, _negative_sample_rate: 5, _n_epochs: 350, _initial_alpha: 1 }, parameters);
        [this._N, this._D] = this.X.shape;
        /* let n_neighbors = Math.min(this._N - 1, parameters.n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(this.parameter("local_connectivity"), n_neighbors - 1)); */
        if (this.parameter("n_neighbors") > this._N) {
            throw new Error(`Parameter n_neighbors (=${this.parameter("n_neighbors")}) needs to be smaller than dataset size (N=${this._N})!`);
        }
        if (this.parameter("local_connectivity") > this.parameter("n_neighbors")) {
            throw new Error(`Parameter local_connectivity (=${this.parameter("local_connectivity")}) needs to be smaller than parameter n_neighbors (=${this.parameter("n_neighbors")})`);
        }
        this._iter = 0;
        const randomizer = this._randomizer;
        this.Y = new Matrix(this._N, this.parameter("d"), () => randomizer.random);
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
            yv[i] = xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread);
        }

        const err = (p) => {
            const error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map((e) => e * e)));
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
        const { local_connectivity, metric } = this._parameters;
        const target = Math.log2(k);
        const rhos = [];
        const sigmas = [];
        const X = this.X;
        const N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());

        const distances = [];
        if (metric === "precomputed") {
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
            const non_zero_dist = search_result.filter((d) => d.value > 0);
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
                    psum += d > 0 ? Math.exp(-(d / mid)) : 1;
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
            distances: distances,
            sigmas: sigmas,
            rhos: rhos,
        };
    }

    /**
     * @private
     * @param {Matrix} X
     * @param {Number} n_neighbors
     * @returns {Matrix}
     */
    _fuzzy_simplicial_set(X, n_neighbors) {
        const N = X.shape[0];
        const { metric, _set_op_mix_ratio } = this._parameters;
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
            .mult(_set_op_mix_ratio)
            .add(prod_matrix.mult(1 - _set_op_mix_ratio));
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
        const n_samples = weights.map((w) => n_epochs * (w / weights_max));
        for (let i = 0; i < result.length; ++i) if (n_samples[i] > 0) result[i] = Math.round(n_epochs / n_samples[i]);
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
        const [rows_n, cols_n] = graph.shape;
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
            rows: rows,
            cols: cols,
            data: data,
        };
    }

    /**
     * Computes all necessary
     * @returns {UMAP}
     */
    init() {
        const { _spread, min_dist, n_neighbors, _n_epochs, _negative_sample_rate } = this._parameters;
        const [a, b] = this._find_ab_params(_spread, min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, n_neighbors);
        const { rows, cols, data: weights } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        this._weights = weights;
        this._epochs_per_sample = this._make_epochs_per_sample(_n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map((d) => d * _negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        return this;
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
    transform(iterations = 350) {
        if (this.parameter("_n_epochs") != iterations) {
            this.parameter("_n_epochs", iterations);
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
    *generator(iterations = 350) {
        if (this.parameter("_n_epochs") != iterations) {
            this.parameter("_n_epochs", iterations);
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
        const randomizer = this._randomizer;
        const { _repulsion_strength, d: dim } = this._parameters;
        const { _alpha: alpha, _a: a, _b: b, _epochs_per_sample: epochs_per_sample, _epochs_per_negative_sample: epochs_per_negative_sample, _epoch_of_next_negative_sample: epoch_of_next_negative_sample, _epoch_of_next_sample: epoch_of_next_sample, _clip: clip } = this;
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
                    const k = randomizer.random_int % tail_length;
                    const other = tail_embedding.row(tail[k]);
                    const dist = euclidean_squared(current, other);
                    let grad_coeff = 0;
                    if (dist > 0) {
                        grad_coeff = (2 * _repulsion_strength * b) / ((0.01 + dist) * (a * Math.pow(dist, b) + 1));
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
                epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
            }
        }
        return head_embedding;
    }

    /**
     * @private
     * @returns {Matrix}
     */
    next() {
        const iter = ++this._iter;
        const Y = this.Y;
        const { _initial_alpha, _n_epochs } = this._parameters;
        this._alpha = _initial_alpha * (1 - iter / _n_epochs);
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }
}

/**
 * @class
 * @alias TriMap
 * @extends DR
 */
class TriMap extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.weight_adj = 500] - scaling factor.
     * @param {Number} [parameters.c = 5] - number of triplets multiplier.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.tol = 1e-8] -
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X, parameters) {
        super(X, { weight_adj: 500, c: 5, d: 2, metric: euclidean, tol: 1e-8, seed: 1212 }, parameters);
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
        const { d, metric, c } = this._parameters;
        this.n_inliers = 2 * c;
        this.n_outliers = 1 * c;
        this.n_random = 1 * c;
        this.Y = pca || new PCA(X, d).transform();
        this.knn = knn || new BallTree(X.to2dArray, metric);
        const { triplets, weights } = this._generate_triplets(this.n_inliers, this.n_outliers, this.n_random);
        this.triplets = triplets;
        this.weights = weights;
        this.lr = (1000 * N) / triplets.shape[0];
        this.C = Infinity;
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
        const { metric, weight_adj } = this._parameters;
        const X = this.X;
        const N = X.shape[0];
        const knn = this.knn;
        const n_extra = Math.min(n_inliers + 20, N);
        const nbrs = new Matrix(N, n_extra);
        const knn_distances = new Matrix(N, n_extra);
        for (let i = 0; i < N; ++i) {
            knn.search(X.row(i), n_extra + 1)
                .raw_data()
                .filter((d) => d.value != 0)
                .sort((a, b) => a.value - b.value)
                .forEach((d, j) => {
                    nbrs.set_entry(i, j, d.element.index);
                    knn_distances.set_entry(i, j, d.value);
                });
        }
        // scale parameter
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            sig[i] = Math.max((knn_distances.entry(i, 3) + knn_distances.entry(i, 4) + knn_distances.entry(i, 5) + knn_distances.entry(i, 6)) / 4, 1e-10);
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
            const { random_triplets, random_weights } = this._sample_random_triplets(X, n_random, sig);
            triplets = triplets.concat(random_triplets, "vertical");
            weights = Float64Array.from([...weights, ...random_weights]);
        }
        n_triplets = triplets.shape[0];
        let max_weight = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            if (isNaN(weights[i])) {
                weights[i] = 0;
            }
            if (max_weight < weights[i]) max_weight = weights[i];
        }
        let max_weight_2 = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight;
            weights[i] += 0.0001;
            weights[i] = Math.log(1 + weight_adj * weights[i]);
            if (max_weight_2 < weights[i]) max_weight_2 = weights[i];
        }
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight_2;
        }
        return {
            triplets: triplets,
            weights: weights,
        };
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
            return Math.exp(-(knn_distances.entry(i, j) ** 2 / sig[i] / sig[nbrs.entry(i, j)]));
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
            const sort_indices = this.__argsort(P.row(i).map((d) => -d));
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
        return A.map((d, i) => {
            return { d: d, i: i };
        })
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
        const interval = linspace(0, max_int - 1).filter((d) => rejects.indexOf(d) < 0);
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
        const metric = this.parameter("metric");
        const randomizer = this._randomizer;
        const N = X.shape[0];
        const random_triplets = new Matrix(N * n_random, 3);
        const random_weights = new Float64Array(N * n_random);
        for (let i = 0; i < N; ++i) {
            const n_i = i * n_random;
            const indices = [...linspace(0, i - 1), ...linspace(i + 1, N - 1)];
            for (let j = 0; j < n_random; ++j) {
                let [sim, out] = randomizer.choice(indices, 2);
                let p_sim = Math.exp(-(metric(X.row(i), X.row(sim)) ** 2 / (sig[i] * sig[sim])));
                if (p_sim < 1e-20) p_sim = 1e-20;
                let p_out = Math.exp(-(metric(X.row(i), X.row(out)) ** 2 / (sig[i] * sig[out])));
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
            random_triplets: random_triplets,
            random_weights: random_weights,
        };
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
        let y_ij = new Float64Array(dim);
        let y_ik = new Float64Array(dim);
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
                    d_ij += y_ij[d] ** 2;
                    d_ik += y_ik[d] ** 2;
                }
                // update y_ik and d_ik only
            } else {
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_kd = Y.entry(k, d);
                    y_ik[d] = Y_id - Y_kd;
                    d_ik += y_ik[d] ** 2;
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
        return { grad, loss, n_viol };
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
     * @param {Number} max_iteration
     * @yields {Matrix}
     * @returns {Matrix}
     */
    *generator(max_iteration = 800) {
        this.check_init();
        for (let iter = 0; iter < max_iteration; ++iter) {
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
        const gamma = iter > 150 ? 0.5 : 0.3;
        const old_C = this.C;
        const vel = this.vel;
        const Y = this.Y.add(vel.mult(gamma));
        const { grad, loss, n_viol } = this._grad(Y);
        this.C = loss;
        this.Y = this._update_embedding(Y, iter, grad);
        this.lr *= old_C > loss + this._parameters.tol ? 1.01 : 0.9;
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
        const gamma = iter > 150 ? 0.9 : 0.5; // moment parameter
        const min_gain = 0.01;
        const gain = this.gain;
        const vel = this.vel;
        const lr = this.lr;
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const new_gain = Math.sign(vel.entry(i, d)) != Math.sign(grad.entry(i, d)) ? gain.entry(i, d) + 0.2 : Math.max(gain.entry(i, d) * 0.8, min_gain);
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
     * Least Squares Projection.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.neighbors = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {Number} [parameters.control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     * @todo accept precomputed distance matrix.
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, control_points: undefined, d: 2, metric: euclidean, seed: 1212 }, parameters);
        this.parameter("neighbors", Math.min(parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        this.parameter("control_points", Math.min(parameters.control_points ?? Math.ceil(Math.sqrt(this._N)), this._N - 1));
        this._is_initialized = false;
        return this;
    }

    /**
     *
     * @param {DR} DR - method used for position control points.
     * @param {Object} DR_parameters - Object containing parameters for the DR method which projects the control points
     * @returns {LSP}
     */
    init(DR = MDS, DR_parameters = {}, KNN = BallTree) {
        if (this._is_initialized) return this;
        const X = this.X;
        const N = this._N;
        const K = this.parameter("neighbors");
        const d = this.parameter("d");
        const seed = this.parameter("seed");
        const metric = this.parameter("metric");
        DR_parameters = Object.assign({d, metric, seed }, DR_parameters);
        const nc = this.parameter("control_points");
        const control_points = new KMedoids(X, nc, null, metric).get_clusters().medoids;
        const C = new Matrix(nc, N, "zeros");
        control_points.forEach((c_i, i) => {
            C.set_entry(i, c_i, 1);
        });
        const Y_C = new DR(Matrix.from(control_points.map((c_i) => X.row(c_i))), DR_parameters).transform();

        const XA = X.to2dArray;
        const knn = new KNN(XA, metric);
        const L = new Matrix(N, N, "I");
        const alpha = -1 / K;
        XA.forEach((x_i, i) => {
            for (const { index: j } of knn.search(x_i, K).iterate()) {
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
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TopoMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { metric: euclidean, seed: 1212 }, parameters);
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
        const { metric} = this._parameters;
        this.Y = new Matrix(this._N, 2, 0);
        this._Emst = this._make_minimum_spanning_tree(metric);
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
        const Y = this.Y.to2dArray;
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
        const Y = this.Y.to2dArray;
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
     * SAMMON's Mapping
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {"PCA"|"MDS"|"random"} [parameters.init = "random"] - Either "PCA" or "MDS", with which SAMMON initialiates the projection. With "random" a random matrix gets used as starting point.
     * @param {Object} [parameters.init_parameters] - Parameters for the {@link init}-DR method.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { magic: 0.1, d: 2, metric: euclidean, seed: 1212, init_DR: "random", init_parameters: {} }, parameters);
        return this;
    }

    /**
     * initializes the projection.
     * @private
     */
    init() {
        const N = this.X.shape[0];
        const { d, metric, init_DR: init_DR, init_parameters: DR_parameters } = this._parameters;
        if (init_DR === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (["PCA", "MDS"].includes(init_DR)) {
            this.Y = Matrix.from(init_DR == "PCA" ? PCA.transform(this.X, DR_parameters) : MDS.transform(this.X, DR_parameters));
        } else {
            throw new Error('init_DR needs to be either "random" or a DR method!')
        }
        this.distance_matrix = metric == "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Matrix|Array} - The projection of {@link X}.
     */
    transform(max_iter = 200) {
        if (!this._is_initialized) this.init();
        for (let j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Generator} - A generator yielding the intermediate steps of the projection of {@link X}.
     */
    *generator(max_iter = 200) {
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
        const N = this.X.shape[0];
        const { d, metric } = this._parameters;
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
                    e1[k] += (delta[k] * dq) / dr;
                    e2[k] += (dq - (Math.pow(delta[k], 2) * (1 + dq / dY)) / dY) / dr;
                }
            }
            for (let k = 0; k < d; ++k) {
                const val = Y.entry(i, k) + ((MAGIC * e1[k]) / Math.abs(e2[k]) || 0);
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

class SQDMDS extends DR {
    /**
     * SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE and UMAP.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @param {Matrix|Number[][]} X
     * @param {Object} [parameters]
     * @param {Number} [parameters.d=2]
     * @param {Function} [parameters.metric = euclidean]
     * @param {Number} [parameters.decay_start = 0.1] - Percentage of iterations using exaggeration phase. If random init: it is recommended to start the decay later to give the time for the global config to adjust with big steps.
     * @param {Number} [parameters.decay_cte = 0.34] - Controls the decay of the learning parameter.
     * @param {Object} [parameters.init_DR]
     * @returns {SQDMDS}
     * @see {@link https://arxiv.org/pdf/2202.12087.pdf}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                d: 2,
                metric: euclidean,
                seed: 1212,
                decay_start: 0.1,
                decay_cte: 0.34, // 0.34
                init_DR: {type: "random"}
            },
            parameters
        );

        return this;
    }

    /**
     * @private
     */
    init() {
        const N = this._N;
        const d = this.parameter("d");

        // initialize helpers.
        this._add = this.__add(d);
        this._sub_div = this.__sub_div(d);
        this._minus = this.__minus(d);
        this._mult = this.__mult(d);
        this._LR_init = Math.max(2, 0.005 * N);
        this._LR = this._LR_init;
        this._offset = -Math.exp(-1 / this.parameter("decay_cte"));
        this._momentums = new Matrix(N, d, 0);
        this._grads = new Matrix(N, d, 0);
        this._indices = linspace(0, N - 1);
        // initialize projection.
        const R = this._randomizer;
        this.Y = new Matrix(N, d, () => R.random - 0.5);

        // preparing metric for optimization.
        const this_metric = this.parameter("metric");
        if (this_metric === "precomputed") {
            this._HD_metric = function (i, j, X) {
                return X.entry(i, j);
            };
            this._HD_metric_exaggeration = function (i, j, X) {
                return Math.pow(X.entry(i, j), 2);
            };
        } else {
            this._HD_metric = function (i, j, X) {
                return this_metric(X.row(i), X.row(j));
            };
            if (this_metric == euclidean) {
                this._HD_metric_exaggeration = function (i, j, X) {
                    return euclidean_squared(X.row(i), X.row(j));
                };
            } else {
                this._HD_metric_exaggeration = function (i, j, X) {
                    return Math.pow(this_metric(X.row(i), X.row(j)), 2);
                };
            }
        }
        return;
    }

    /**
     * Computes the projection.
     * @param {Number} [iterations=500] - Number of iterations.
     * @returns {Matrix|Number[][]} the projection.
     */
    transform(iterations = 500) {
        this.check_init();
        this._decay_start = Math.round(this.parameter("decay_start") * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
        }
        return this.projection;
    }

    /**
     * Computes the projection.
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Number[][]} the intermediate steps of the projection.
     */
    *generator(iterations = 500) {
        this.check_init();
        this._decay_start = Math.round(this.parameter("decay_start") * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * Performs an optimization step.
     * @private
     * @param {Number} i - Acutal iteration.
     * @param {Number} iterations - Number of iterations.
     */
    _step(i, iterations) {
        const decay_start = this._decay_start;
        if (i > decay_start) {
            const decay_cte = this.parameter("decay_cte");
            const offset = this._offset;
            const ratio = (i - decay_start) / (iterations - decay_start);
            this._LR = this._LR_init * (Math.exp(-(ratio * ratio) / decay_cte) + offset);
            this._distance_exaggeration = false;
        } else {
            this._distance_exaggeration = true;
        }
        this._nestrov_iteration(this._distance_exaggeration);
    }

    /**
     * Creates quartets of non overlapping indices.
     * @private
     * @returns {Number[][]}
     */
    __quartets() {
        const N = this._N;
        const max_N = N - (N % 4);
        const R = this._randomizer;
        const shuffled_indices = R.choice(this._indices, max_N);
        const result = [];
        for (let i = 0; i < max_N; i += 4) {
            result.push(Uint32Array.of(shuffled_indices[i], shuffled_indices[i + 1], shuffled_indices[i + 2], shuffled_indices[i + 3]));
        }
        return result;
    }

    /**
     * Computes and applies gradients, and updates momentum.
     * @private
     * @param {Boolean} distance_exaggeration
     */
    _nestrov_iteration(distance_exaggeration) {
        const momentums = this._momentums.mult(0.99, { inline: true });
        const LR = this._LR;
        const grads = this._fill_MDS_grads(this.Y.add(momentums), this._grads, distance_exaggeration);
        const [n, d] = momentums.shape;
        for (let i = 0; i < n; ++i) {
            const g_i = grads.row(i);
            const g_i_norm = norm(g_i);
            if (g_i_norm == 0) continue;
            const mul = LR / g_i_norm;
            const m_i = momentums.row(i);
            for (let j = 0; j < d; ++j) {
                m_i[j] -= mul * g_i[j];
            }
        } // momentums -= (LR / norm) * grads
        this.Y.add(momentums, { inline: true });
    }

    /**
     * Computes the gradients.
     * @param {Matrix} Y - The Projection.
     * @param {Matrix} grads - The gradients.
     * @param {Boolean} [exaggeration = false] - Whether or not to use early exaggeration.
     * @param {Boolean} [zero_grad = true] - Whether or not to reset the gradient in the beginning.
     * @returns {Matrix} the gradients.
     */
    _fill_MDS_grads(Y, grads, exaggeration = false, zero_grad = true) {
        if (zero_grad) {
            // compute new gradients
            grads.values.fill(0);
        }
        const add = this._add;
        const X = this.X;
        let HD_metric;
        if (exaggeration == true) {
            HD_metric = this._HD_metric_exaggeration;
        } else {
            HD_metric = this._HD_metric;
        }

        const D_quartet = new Float64Array(6);
        const quartets = this.__quartets();
        for (const [i, j, k, l] of quartets) {
            // compute quartet's HD distances.
            D_quartet[0] = HD_metric(i, j, X);
            D_quartet[1] = HD_metric(i, k, X);
            D_quartet[2] = HD_metric(i, l, X);
            D_quartet[3] = HD_metric(j, k, X);
            D_quartet[4] = HD_metric(j, l, X);
            D_quartet[5] = HD_metric(k, l, X);

            const D_quartet_sum = neumair_sum(D_quartet);

            if (D_quartet_sum > 0) {
                for (let i = 0; i < 6; ++i) {
                    D_quartet[i] /= D_quartet_sum;
                    D_quartet[i] += 1e-11;
                }
            }
            const [gi, gj, gk, gl] = this._compute_quartet_grads(Y, [i, j, k, l], D_quartet);

            // add is inline, row acces the matrix
            add(grads.row(i), gi);
            add(grads.row(j), gj);
            add(grads.row(k), gk);
            add(grads.row(l), gl);
        }
        return grads;
    }

    /**
     * Quartet gradients for a projection.
     * @private
     * @param {Matrix} Y - The acutal projection.
     * @param {Number[]} quartet - The indices of the quartet.
     * @param {Number[]} D_hd - The high-dimensional distances of the quartet.
     * @returns {Number[][]} the gradients for the quartet.
     */
    _compute_quartet_grads(Y, quartet, [p_ab, p_ac, p_ad, p_bc, p_bd, p_cd]) {
        const [a, b, c, d] = quartet.map((index) => Y.row(index));
        // LD distances, add a small number just in case
        const d_ab = euclidean(a, b) + 1e-12;
        const d_ac = euclidean(a, c) + 1e-12;
        const d_ad = euclidean(a, d) + 1e-12;
        const d_bc = euclidean(b, c) + 1e-12;
        const d_bd = euclidean(b, d) + 1e-12;
        const d_cd = euclidean(c, d) + 1e-12;
        const sum_LD_dist = neumair_sum([d_ab, d_ac, d_ad, d_bc, d_bd, d_cd]);

        // for each element of the sum: use the same gradient function and just permute the points given in input.
        const [gA1, gB1, gC1, gD1] = this._ABCD_grads(a, b, c, d, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd, p_ab, sum_LD_dist);
        const [gA2, gC2, gB2, gD2] = this._ABCD_grads(a, c, b, d, d_ac, d_ab, d_ad, d_bc, d_cd, d_bd, p_ac, sum_LD_dist);
        const [gA3, gD3, gC3, gB3] = this._ABCD_grads(a, d, c, b, d_ad, d_ac, d_ab, d_cd, d_bd, d_bc, p_ad, sum_LD_dist);
        const [gB4, gC4, gA4, gD4] = this._ABCD_grads(b, c, a, d, d_bc, d_ab, d_bd, d_ac, d_cd, d_ad, p_bc, sum_LD_dist);
        const [gB5, gD5, gA5, gC5] = this._ABCD_grads(b, d, a, c, d_bd, d_ab, d_bc, d_ad, d_cd, d_ac, p_bd, sum_LD_dist);
        const [gC6, gD6, gA6, gB6] = this._ABCD_grads(c, d, a, b, d_cd, d_ac, d_bc, d_ad, d_bd, d_ab, p_cd, sum_LD_dist);

        const add = this._add;
        const gA = add(gA1, gA2, gA3, gA4, gA5, gA6);
        const gB = add(gB1, gB2, gB3, gB4, gB5, gB6);
        const gC = add(gC1, gC2, gC3, gC4, gC5, gC6);
        const gD = add(gD1, gD2, gD3, gD4, gD5, gD6);

        return [gA, gB, gC, gD];
    }

    /**
     * Gradients for one element of the loss function's sum.
     * @private
     */
    _ABCD_grads(a, b, c, d, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd, p_ab, sum_LD_dist) {
        const ratio = d_ab / sum_LD_dist;
        const twice_ratio = 2 * ((p_ab - ratio) / sum_LD_dist);
        const minus = this._minus;
        const add = this._add;
        const mult = this._mult;
        const sub_div = this._sub_div;
        // no side effects because sub_div creates new arrays, and the inline functions work on this new created arrays.
        const gA = mult(minus(mult(add(sub_div(a, b, d_ab), sub_div(a, c, d_ac), sub_div(a, d, d_ad)), ratio), sub_div(a, b, d_ab)), twice_ratio);
        const gB = mult(minus(mult(add(sub_div(b, a, d_ab), sub_div(b, c, d_bc), sub_div(b, d, d_bd)), ratio), sub_div(b, a, d_ab)), twice_ratio);
        const gC = mult(add(sub_div(c, a, d_ac), sub_div(c, b, d_bc), sub_div(c, d, d_cd)), ratio * twice_ratio);
        const gD = mult(add(sub_div(d, a, d_ad), sub_div(d, b, d_bd), sub_div(d, c, d_cd)), ratio * twice_ratio);
        return [gA, gB, gC, gD];
    }

    /**
     * Inline!
     */
    __minus(d) {
        return (a, b) => {
            for (let i = 0; i < d; ++i) {
                a[i] -= b[i];
            }
            return a;
        };
    }

    /**
     * Inline!
     */
    __add(d) {
        return (...summands) => {
            const n = summands.length;
            const s1 = summands[0];
            for (let j = 1; j < n; ++j) {
                const summand = summands[j];
                for (let i = 0; i < d; ++i) {
                    s1[i] += summand[i];
                }
            }
            return s1;
        };
    }

    /**
     * Inline!
     */
    __mult(d) {
        return (a, v) => {
            for (let i = 0; i < d; ++i) {
                a[i] *= v;
            }
            return a;
        };
    }

    /**
     * Creates a new array <code>(x - y) / div</code>
     */
    __sub_div(d) {
        return (x, y, div) => {
            return Float64Array.from({ length: d }, (_, i) => (x[i] - y[i]) / div);
        };
    }
}

var version="0.6.0";

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
exports.SQDMDS = SQDMDS;
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
exports.inner_product = inner_product;
exports.jaccard = jaccard;
exports.k_nearest_neighbors = k_nearest_neighbors;
exports.kahan_sum = kahan_sum;
exports.linspace = linspace;
exports.manhattan = manhattan;
exports.max = max;
exports.min = min;
exports.neumair_sum = neumair_sum;
exports.norm = norm;
exports.normalize = normalize;
exports.powell = powell;
exports.qr = qr_gramschmidt;
exports.qr_householder = qr_householder;
exports.simultaneous_poweriteration = simultaneous_poweriteration;
exports.sokal_michener = sokal_michener;
exports.version = version;
exports.yule = yule;

Object.defineProperty(exports, '__esModule', { value: true });

}));
//# sourceMappingURL=druid.js.map
