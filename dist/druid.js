// https://renecutura.eu v0.5.1 Copyright 2022 Rene Cutura
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('tslib')) :
typeof define === 'function' && define.amd ? define(['exports', 'tslib'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.druid = global.druid || {}, global.tslib));
})(this, (function (exports, tslib) { 'use strict';

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
    var n = summands.length;
    var sum = 0;
    var compensation = 0;
    var y, t;
    for (var i = 0; i < n; ++i) {
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
    var n = summands.length;
    var sum = 0;
    var compensation = 0;
    for (var i = 0; i < n; ++i) {
        var summand = summands[i];
        var t = sum + summand;
        if (Math.abs(sum) >= Math.abs(summand)) {
            compensation += sum - t + summand;
        }
        else {
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var s = new Array(n);
    for (var i = 0; i < n; ++i) {
        var x = a[i];
        var y = b[i];
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
    if (a.length !== b.length)
        return undefined;
    var n = a.length;
    var sum = 0;
    var sum_a = 0;
    var sum_b = 0;
    for (var i = 0; i < n; ++i) {
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var sum = 0;
    for (var i = 0; i < n; ++i) {
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var res = [];
    for (var i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max.apply(Math, res);
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
function canberra (a, b) {
    if (a.length !== b.length)
        return undefined;
    var n = a.length;
    var sum = 0;
    for (var i = 0; i < n; ++i) {
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var num_non_zero = 0;
    var num_equal = 0;
    for (var i = 0; i < n; ++i) {
        var x = a[i] != 0;
        var y = b[i] != 0;
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var disagree = 0;
    for (var i = 0; i < n; ++i) {
        var x = a[i];
        var y = b[i];
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
function sokal_michener (a, b) {
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var num_not_equal = 0;
    for (var i = 0; i < n; ++i) {
        var x = a[i] != 0;
        var y = b[i] != 0;
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
    if (a.length != b.length)
        return undefined;
    var n = a.length;
    var num_true_true = 0;
    var num_true_false = 0;
    var num_false_true = 0;
    for (var i = 0; i < n; ++i) {
        var x = a[i] != 0;
        var y = b[i] != 0;
        num_true_true += x && y;
        num_true_false += x && !y;
        num_false_true += !x && x;
    }
    var num_false_false = n - num_true_true - num_true_false - num_false_true;
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
function k_nearest_neighbors (A, k, metric) {
    if (metric === void 0) { metric = euclidean; }
    var rows = A.shape[0];
    var D = metric == "precomputed" ? A : distance_matrix(A, metric);
    var nN = new Array(rows);
    var _loop_1 = function (row) {
        nN[row] = Array.from(D.row(row))
            .map(function (distance, col) {
            return {
                i: row,
                j: col,
                distance: distance
            };
        })
            .sort(function (a, b) { return a.distance - b.distance; })
            .slice(1, k + 1);
    };
    for (var row = 0; row < rows; ++row) {
        _loop_1(row);
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
function distance_matrix (A, metric) {
    if (metric === void 0) { metric = euclidean; }
    var n = A.shape[0];
    var D = new Matrix(n, n);
    for (var i = 0; i < n; ++i) {
        var A_i = A.row(i);
        for (var j = i + 1; j < n; ++j) {
            var dist = metric(A_i, A.row(j));
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
function linspace (start, end, number) {
    if (number === void 0) { number = null; }
    if (!number) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return number === 1 ? [start] : [];
    }
    var result = new Array(number);
    number -= 1;
    for (var i = number; i >= 0; --i) {
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
function norm (v, metric) {
    if (metric === void 0) { metric = euclidean; }
    var vector = null;
    if (v instanceof Matrix) {
        var _a = v.shape, rows = _a[0], cols = _a[1];
        if (rows === 1)
            vector = v.row(0);
        else if (cols === 1)
            vector = v.col(0);
        else
            throw new Error("Matrix must be 1d!");
    }
    else {
        vector = v;
    }
    var n = vector.length;
    var zeros = Float64Array.from({ length: n }, function () { return 0; });
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
function normalize (v, metric) {
    if (metric === void 0) { metric = euclidean; }
    var v_norm = norm(v, metric);
    return v.map(function (value) { return value / v_norm; });
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
    var _a = A.shape, rows = _a[0], cols = _a[1];
    var Q = new Matrix(rows, cols, "identity");
    var R = new Matrix(cols, cols, 0);
    var _loop_1 = function (j) {
        var v = A.col(j);
        var _loop_2 = function (i) {
            var q = Q.col(i);
            var q_dot_v = neumair_sum(q.map(function (q_, k) { return q_ * v[k]; }));
            R.set_entry(i, j, q_dot_v);
            v = v.map(function (v_, k) { return v_ - q_dot_v * q[k]; });
        };
        for (var i = 0; i < j; ++i) {
            _loop_2(i);
        }
        var v_norm = norm(v, euclidean);
        for (var k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j, j, v_norm);
    };
    for (var j = 0; j < cols; ++j) {
        _loop_1(j);
    }
    return { R: R, Q: Q };
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
    var _a = A.shape, rows = _a[0], cols = _a[1];
    var Q = new Matrix(rows, rows, "I");
    var R = A.clone();
    for (var j = 0; j < cols; ++j) {
        var x = Matrix.from(R.col(j).slice(j));
        var x_norm = norm(x);
        var x0 = x.entry(0, 0);
        var rho = -Math.sign(x0);
        var u1 = x0 - rho * x_norm;
        var u = x.divide(u1).set_entry(0, 0, 1);
        var beta = (-rho * u1) / x_norm;
        var u_outer_u = u.outer(u);
        var R_block = R.get_block(j, 0);
        var new_R = R_block.sub(u_outer_u.dot(R_block).mult(beta));
        var Q_block = Q.get_block(0, j);
        var new_Q = Q_block.sub(Q_block.dot(u_outer_u).mult(beta));
        R.set_block(j, 0, new_R);
        Q.set_block(0, j, new_Q);
    }
    return { R: R, Q: Q };
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
function simultaneous_poweriteration (A, k, _a) {
    if (k === void 0) { k = 2; }
    var _b = _a === void 0 ? {} : _a, _c = _b.seed, seed = _c === void 0 ? 1212 : _c, _d = _b.max_iterations, max_iterations = _d === void 0 ? 100 : _d, _e = _b.qr, qr = _e === void 0 ? qr_gramschmidt : _e, _f = _b.tol, tol = _f === void 0 ? 1e-8 : _f;
    var randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix))
        A = Matrix.from(A);
    var n = A.shape[0];
    var _g = qr(new Matrix(n, k, function () { return (randomizer.random - .5) * 2; })), Q = _g.Q, R = _g.R;
    while (max_iterations--) {
        var oldQ = Q.clone();
        var Z = A.dot(Q);
        var QR = qr(Z);
        Q = QR.Q;
        R = QR.R;
        var error = euclidean_squared(Q.values, oldQ.values);
        if (error < tol) {
            break;
        }
    }
    var eigenvalues = R.diag;
    var eigenvectors = Q.transpose().to2dArray;
    return { eigenvalues: eigenvalues, eigenvectors: eigenvectors };
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
    var N = a.length;
    if (N != b.length) {
        throw new Error("Array a and b must have the same length!");
    }
    var sum = 0;
    for (var i = 0; i < N; ++i) {
        sum += a * b;
    }
    return sum;
}

/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
var Matrix = /** @class */ (function () {
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
    function Matrix(rows, cols, value) {
        if (rows === void 0) { rows = null; }
        if (cols === void 0) { cols = null; }
        if (value === void 0) { value = null; }
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
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
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
                    for (var row = 0; row < rows; ++row) {
                        this._data[row * cols + row] = 1;
                    }
                    return this;
                }
                if (value === "center" && rows == cols) {
                    this._data = new Float64Array(rows * cols);
                    value = function (i, j) { return (i === j ? 1 : 0) - 1 / rows; };
                    for (var row = 0; row < rows; ++row) {
                        for (var col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = value(row, col);
                        }
                    }
                    return this;
                }
            }
            if (typeof value === "number") {
                this._data = new Float64Array(rows * cols);
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
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
    Matrix.from = function (A, type) {
        if (type === void 0) { type = "row"; }
        if (A instanceof Matrix) {
            return A.clone();
        }
        else if (Array.isArray(A) || A instanceof Float64Array) {
            var m = A.length;
            if (m === 0)
                throw new Error("Array is empty");
            // 1d
            if (!Array.isArray(A[0]) && !(A[0] instanceof Float64Array)) {
                if (type === "row") {
                    return new Matrix(1, m, function (_, j) { return A[j]; });
                }
                else if (type === "col") {
                    return new Matrix(m, 1, function (i) { return A[i]; });
                }
                else if (type === "diag") {
                    return new Matrix(m, m, function (i, j) { return (i == j ? A[i] : 0); });
                }
                else {
                    throw new Error("1d array has NaN entries");
                }
                // 2d
            }
            else if (Array.isArray(A[0]) || A[0] instanceof Float64Array) {
                var n = A[0].length;
                for (var row = 0; row < m; ++row) {
                    if (A[row].length !== n) {
                        throw new Error("various array lengths");
                    }
                }
                return new Matrix(m, n, function (i, j) { return A[i][j]; });
            }
        }
        else if (typeof A === "number") {
            return new Matrix(1, 1, A);
        }
        else {
            throw new Error("error");
        }
    };
    /**
     * Returns the {@link row}<sup>th</sup> row from the Matrix.
     * @param {Number} row
     * @returns {Float64Array}
     */
    Matrix.prototype.row = function (row) {
        var data = this.values;
        var cols = this._cols;
        return data.subarray(row * cols, (row + 1) * cols);
    };
    /**
     * Returns an generator yielding each row of the Matrix.
     * @yields {Float64Array}
     */
    Matrix.prototype.iterate_rows = function () {
        var cols, rows, data, row;
        return tslib.__generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    cols = this._cols;
                    rows = this._rows;
                    data = this.values;
                    row = 0;
                    _a.label = 1;
                case 1:
                    if (!(row < rows)) return [3 /*break*/, 4];
                    return [4 /*yield*/, data.subarray(row * cols, (row + 1) * cols)];
                case 2:
                    _a.sent();
                    _a.label = 3;
                case 3:
                    ++row;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/];
            }
        });
    };
    /**
     * Makes a {@link Matrix} object an iterable object.
     * @yields {Float64Array}
     */
    Matrix.prototype[Symbol.iterator] = function () {
        var _i, _a, row;
        return tslib.__generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    _i = 0, _a = this.iterate_rows();
                    _b.label = 1;
                case 1:
                    if (!(_i < _a.length)) return [3 /*break*/, 4];
                    row = _a[_i];
                    return [4 /*yield*/, row];
                case 2:
                    _b.sent();
                    _b.label = 3;
                case 3:
                    _i++;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/];
            }
        });
    };
    /**
     * Sets the entries of {@link row}<sup>th</sup> row from the Matrix to the entries from {@link values}.
     * @param {int} row
     * @param {Array} values
     * @returns {Matrix}
     */
    Matrix.prototype.set_row = function (row, values) {
        var cols = this._cols;
        if (Array.isArray(values) && values.length === cols) {
            var offset = row * cols;
            for (var col = 0; col < cols; ++col) {
                this.values[offset + col] = values[col];
            }
        }
        else if (values instanceof Matrix && values.shape[1] === cols && values.shape[0] === 1) {
            var offset = row * cols;
            for (var col = 0; col < cols; ++col) {
                this.values[offset + col] = values._data[col];
            }
        }
        return this;
    };
    /**
     * Returns the {@link col}<sup>th</sup> column from the Matrix.
     * @param {int} col
     * @returns {Array}
     */
    Matrix.prototype.col = function (col) {
        var result_col = new Float64Array(this._rows);
        for (var row = 0; row < this._rows; ++row) {
            result_col[row] = this.values[row * this._cols + col];
        }
        return result_col;
    };
    /**
     * Returns the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @returns {float64}
     */
    Matrix.prototype.entry = function (row, col) {
        return this.values[row * this._cols + col];
    };
    /**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given {@link value}.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */
    Matrix.prototype.set_entry = function (row, col, value) {
        this.values[row * this._cols + col] = value;
        return this;
    };
    /**
     * Returns a new transposed Matrix.
     * @returns {Matrix}
     */
    Matrix.prototype.transpose = function () {
        var _this = this;
        var B = new Matrix(this._cols, this._rows, function (row, col) { return _this.entry(col, row); });
        return B;
    };
    Object.defineProperty(Matrix.prototype, "T", {
        /**
         * Returns a new transposed Matrix. Short-form of {@function transpose}.
         * @returns {Matrix}
         */
        get: function () {
            return this.transpose();
        },
        enumerable: false,
        configurable: true
    });
    /**
     * Returns the inverse of the Matrix.
     * @returns {Matrix}
     */
    Matrix.prototype.inverse = function () {
        var _this = this;
        var rows = this._rows;
        var cols = this._cols;
        var B = new Matrix(rows, 2 * cols, function (i, j) {
            if (j >= cols) {
                return i === j - cols ? 1 : 0;
            }
            else {
                return _this.entry(i, j);
            }
        });
        var h = 0;
        var k = 0;
        while (h < rows && k < cols) {
            var i_max = 0;
            var max_val = -Infinity;
            for (var i = h; i < rows; ++i) {
                var val = Math.abs(B.entry(i, k));
                if (max_val < val) {
                    i_max = i;
                    max_val = val;
                }
            }
            if (B.entry(i_max, k) == 0) {
                k++;
            }
            else {
                // swap rows
                for (var j = 0; j < 2 * cols; ++j) {
                    var h_val = B.entry(h, j);
                    var i_val = B.entry(i_max, j);
                    B.set_entry(h, j, h_val);
                    B.set_entry(i_max, j, i_val);
                }
                for (var i = h + 1; i < rows; ++i) {
                    var f = B.entry(i, k) / B.entry(h, k);
                    B.set_entry(i, k, 0);
                    for (var j = k + 1; j < 2 * cols; ++j) {
                        B.set_entry(i, j, B.entry(i, j) - B.entry(h, j) * f);
                    }
                }
                h++;
                k++;
            }
        }
        for (var row = 0; row < rows; ++row) {
            var f = B.entry(row, row);
            for (var col = row; col < 2 * cols; ++col) {
                B.set_entry(row, col, B.entry(row, col) / f);
            }
        }
        for (var row = rows - 1; row >= 0; --row) {
            var B_row_row = B.entry(row, row);
            for (var i = 0; i < row; i++) {
                var B_i_row = B.entry(i, row);
                var f = B_i_row / B_row_row;
                for (var j = i; j < 2 * cols; ++j) {
                    var B_i_j = B.entry(i, j);
                    var B_row_j = B.entry(row, j);
                    B_i_j = B_i_j - B_row_j * f;
                    B.set_entry(i, j, B_i_j);
                }
            }
        }
        return new Matrix(rows, cols, function (i, j) { return B.entry(i, j + cols); });
    };
    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    Matrix.prototype.dot = function (B) {
        if (B instanceof Matrix) {
            var A_1 = this;
            if (A_1.shape[1] !== B.shape[0]) {
                throw new Error("A.dot(B): A is a ".concat(A_1.shape.join(" ⨯ "), "-Matrix, B is a ").concat(B.shape.join(" ⨯ "), "-Matrix: \n                A has ").concat(A_1.shape[1], " cols and B ").concat(B.shape[0], " rows. \n                Must be equal!"));
            }
            var I_1 = A_1.shape[1];
            var C = new Matrix(A_1.shape[0], B.shape[1], function (row, col) {
                var A_i = A_1.row(row);
                var B_i = B.col(col);
                var sum = 0;
                for (var i = 0; i < I_1; ++i) {
                    sum += A_i[i] * B_i[i];
                }
                return sum;
            });
            return C;
        }
        else if (Array.isArray(B) || B instanceof Float64Array) {
            var rows = this._rows;
            if (B.length !== rows) {
                throw new Error("A.dot(B): A has ".concat(rows, " cols and B has ").concat(B.length, " rows. Must be equal!"));
            }
            var C = new Array(rows);
            var _loop_1 = function (row) {
                C[row] = neumair_sum(this_1.row(row).map(function (e) { return e * B[row]; }));
            };
            var this_1 = this;
            for (var row = 0; row < rows; ++row) {
                _loop_1(row);
            }
            return C;
        }
        else {
            throw new Error("B must be Matrix or Array");
        }
    };
    /**
     * Computes the outer product from {@link this} and {@link B}.
     * @param {Matrix} B
     * @returns {Matrix}
     */
    Matrix.prototype.outer = function (B) {
        var A = this;
        var l = A._data.length;
        var r = B._data.length;
        if (l != r)
            return undefined;
        var C = new Matrix();
        C.shape = [
            l,
            l,
            function (i, j) {
                if (i <= j) {
                    return A._data[i] * B._data[j];
                }
                else {
                    return C.entry(j, i);
                }
            },
        ];
        return C;
    };
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
    Matrix.prototype.concat = function (B, type) {
        if (type === void 0) { type = "horizontal"; }
        var A = this;
        var _a = A.shape, rows_A = _a[0], cols_A = _a[1];
        var _b = B.shape, rows_B = _b[0], cols_B = _b[1];
        if (type == "horizontal") {
            if (rows_A != rows_B) {
                throw new Error("A.concat(B, \"horizontal\"): A and B need same number of rows, A has ".concat(rows_A, " rows, B has ").concat(rows_B, " rows."));
            }
            var X = new Matrix(rows_A, cols_A + cols_B, "zeros");
            X.set_block(0, 0, A);
            X.set_block(0, cols_A, B);
            return X;
        }
        else if (type == "vertical") {
            if (cols_A != cols_B) {
                throw new Error("A.concat(B, \"vertical\"): A and B need same number of columns, A has ".concat(cols_A, " columns, B has ").concat(cols_B, " columns."));
            }
            var X = new Matrix(rows_A + rows_B, cols_A, "zeros");
            X.set_block(0, 0, A);
            X.set_block(rows_A, 0, B);
            return X;
        }
        else if (type == "diag") {
            var X = new Matrix(rows_A + rows_B, cols_A + cols_B, "zeros");
            X.set_block(0, 0, A);
            X.set_block(rows_A, cols_A, B);
            return X;
        }
        else {
            throw new Error("type must be \"horizontal\" or \"vertical\", but type is ".concat(type, "!"));
        }
    };
    /**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row
     * @param {int} offset_col
     * @param {Matrix} B
     * @returns {Matrix}
     */
    Matrix.prototype.set_block = function (offset_row, offset_col, B) {
        var _a = B.shape, rows = _a[0], cols = _a[1];
        for (var row = 0; row < rows; ++row) {
            if (row > this._rows) {
                continue;
            }
            for (var col = 0; col < cols; ++col) {
                if (col > this._cols) {
                    continue;
                }
                this.set_entry(row + offset_row, col + offset_col, B.entry(row, col));
            }
        }
        return this;
    };
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
    Matrix.prototype.get_block = function (start_row, start_col, end_row, end_col) {
        if (end_row === void 0) { end_row = null; }
        if (end_col === void 0) { end_col = null; }
        var _a = this.shape, rows = _a[0], cols = _a[1];
        end_row = end_row !== null && end_row !== void 0 ? end_row : rows;
        end_col = end_col !== null && end_col !== void 0 ? end_col : cols;
        if (end_row <= start_row || end_col <= start_col) {
            throw new Error("\n                end_row must be greater than start_row, and \n                end_col must be greater than start_col, but\n                end_row = ".concat(end_row, ", start_row = ").concat(start_row, ", end_col = ").concat(end_col, ", and start_col = ").concat(start_col, "!"));
        }
        var X = new Matrix(end_row - start_row, end_col - start_col, "zeros");
        for (var row = start_row, new_row = 0; row < end_row; ++row, ++new_row) {
            for (var col = start_col, new_col = 0; col < end_col; ++col, ++new_col) {
                X.set_entry(new_row, new_col, this.entry(row, col));
            }
        }
        return X;
        //return new Matrix(end_row - start_row, end_col - start_col, (i, j) => this.entry(i + start_row, j + start_col));
    };
    /**
     * Returns a new array gathering entries defined by the indices given by argument.
     * @param {Array<Number>} row_indices - Array consists of indices of rows for gathering entries of this matrix
     * @param {Array<Number>} col_indices  - Array consists of indices of cols for gathering entries of this matrix
     * @returns {Matrix}
     */
    Matrix.prototype.gather = function (row_indices, col_indices) {
        var N = row_indices.length;
        var D = col_indices.length;
        var R = new Matrix(N, D);
        for (var i = 0; i < N; ++i) {
            var row_index = row_indices[i];
            for (var j = 0; j < N; ++j) {
                var col_index = col_indices[j];
                R.set_entry(i, j, this.entry(row_index, col_index));
            }
        }
        return R;
    };
    /**
     * Applies a function to each entry of the matrix.
     * @private
     * @param {function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    Matrix.prototype._apply_array = function (f, v) {
        var data = this.values;
        var _a = this.shape, rows = _a[0], cols = _a[1];
        for (var row = 0; row < rows; ++row) {
            var offset = row * cols;
            for (var col = 0; col < cols; ++col) {
                var i = offset + col;
                data[i] = f(data[i], v(row, col));
            }
        }
        return this;
    };
    Matrix.prototype._apply_rowwise_array = function (values, f) {
        return this._apply_array(f, function (_, j) { return values[j]; });
    };
    Matrix.prototype._apply_colwise_array = function (values, f) {
        var data = this.values;
        var _a = this.shape, rows = _a[0], cols = _a[1];
        for (var row = 0; row < rows; ++row) {
            var offset = row * cols;
            for (var col = 0; col < cols; ++col) {
                var i = offset + col;
                data[i] = f(data[i], values[row]);
            }
        }
        return this;
    };
    Matrix.prototype._apply = function (value, f) {
        var data = this.values;
        if (value instanceof Matrix) {
            var _a = value.shape, value_rows = _a[0], value_cols = _a[1];
            var _b = this.shape, rows = _b[0], cols = _b[1];
            if (value_rows === 1) {
                if (cols !== value_cols) {
                    throw new Error("cols !== value_cols");
                }
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(0, col));
                    }
                }
            }
            else if (value_cols === 1) {
                if (rows !== value_rows) {
                    throw new Error("rows !== value_rows");
                }
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, 0));
                    }
                }
            }
            else if (rows == value_rows && cols == value_cols) {
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, col));
                    }
                }
            }
            else {
                throw new Error("error");
            }
        }
        else if (Array.isArray(value)) {
            var rows = this._rows;
            var cols = this._cols;
            if (value.length === rows) {
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[row]);
                    }
                }
            }
            else if (value.length === cols) {
                for (var row = 0; row < rows; ++row) {
                    for (var col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[col]);
                    }
                }
            }
            else {
                throw new Error("error");
            }
        }
        else {
            for (var i = 0, n = this._rows * this._cols; i < n; ++i) {
                data[i] = f(data[i], value);
            }
        }
        return this;
    };
    /**
     * Clones the Matrix.
     * @returns {Matrix}
     */
    Matrix.prototype.clone = function () {
        var B = new Matrix();
        B._rows = this._rows;
        B._cols = this._cols;
        B._data = this.values.slice(0);
        return B;
    };
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
    Matrix.prototype.mult = function (value) {
        return this.clone()._apply(value, function (a, b) { return a * b; });
    };
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
    Matrix.prototype.divide = function (value) {
        return this.clone()._apply(value, function (a, b) { return a / b; });
    };
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
    Matrix.prototype.add = function (value) {
        return this.clone()._apply(value, function (a, b) { return a + b; });
    };
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
    Matrix.prototype.sub = function (value) {
        return this.clone()._apply(value, function (a, b) { return a - b; });
    };
    Object.defineProperty(Matrix.prototype, "shape", {
        /**
         * Returns the number of rows and columns of the Matrix.
         * @returns {Array} An Array in the form [rows, columns].
         */
        get: function () {
            return [this._rows, this._cols];
        },
        /**
         * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
         * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
         * @returns {Matrix}
         */
        set: function (_a) {
            var rows = _a[0], cols = _a[1], _b = _a[2], value = _b === void 0 ? function () { return 0; } : _b;
            this._rows = rows;
            this._cols = cols;
            this._data = new Float64Array(rows * cols);
            for (var row = 0; row < rows; ++row) {
                for (var col = 0; col < cols; ++col) {
                    this._data[row * cols + col] = value(row, col);
                }
            }
            return this;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "to2dArray", {
        /**
         * Returns the Matrix as a Array of Float64Arrays.
         * @returns {Array<Float64Array>}
         */
        get: function () {
            var result = [];
            for (var _i = 0, _a = this.iterate_rows(); _i < _a.length; _i++) {
                var row = _a[_i];
                result.push(row);
            }
            return result;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "asArray", {
        /**
         * Returns the Matrix as a Array of Arrays.
         * @returns {Array<Array>}
         */
        get: function () {
            var result = [];
            for (var _i = 0, _a = this.iterate_rows(); _i < _a.length; _i++) {
                var row = _a[_i];
                result.push(Array.from(row));
            }
            return result;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "diag", {
        /**
         * Returns the diagonal of the Matrix.
         * @returns {Float64Array}
         */
        get: function () {
            var rows = this._rows;
            var cols = this._cols;
            var min_row_col = Math.min(rows, cols);
            var result = new Float64Array(min_row_col);
            for (var i = 0; i < min_row_col; ++i) {
                result[i] = this.entry(i, i);
            }
            return result;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "mean", {
        /**
         * Returns the mean of all entries of the Matrix.
         * @returns {Number}
         */
        get: function () {
            var sum = this.sum;
            var n = this._rows * this._cols;
            return sum / n;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "sum", {
        /**
         * Returns the sum oof all entries of the Matrix.
         * @returns {Number}
         */
        get: function () {
            var data = this.values;
            return neumair_sum(data);
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "values", {
        /**
         * Returns the sum oof all entries of the Matrix.
         * @returns {Float64Array}
         */
        get: function () {
            var data = this._data;
            return data;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "meanRows", {
        /**
         * Returns the mean of each row of the matrix.
         * @returns {Float64Array}
         */
        get: function () {
            var data = this.values;
            var rows = this._rows;
            var cols = this._cols;
            var result = Float64Array.from({ length: rows });
            for (var row = 0; row < rows; ++row) {
                result[row] = 0;
                for (var col = 0; col < cols; ++col) {
                    result[row] += data[row * cols + col];
                }
                result[row] /= cols;
            }
            return result;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Matrix.prototype, "meanCols", {
        /** Returns the mean of each column of the matrix.
         * @returns {Float64Array}
         */
        get: function () {
            var data = this.values;
            var rows = this._rows;
            var cols = this._cols;
            var result = Float64Array.from({ length: cols });
            for (var col = 0; col < cols; ++col) {
                result[col] = 0;
                for (var row = 0; row < rows; ++row) {
                    result[col] += data[row * cols + col];
                }
                result[col] /= rows;
            }
            return result;
        },
        enumerable: false,
        configurable: true
    });
    /**
     * Solves the equation {@link A}x = {@link b} using the conjugate gradient method. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer} [randomizer=null]
     * @param {Number} [tol=1e-3]
     * @returns {Matrix}
     */
    Matrix.solve_CG = function (A, b, randomizer, tol) {
        if (tol === void 0) { tol = 1e-3; }
        if (randomizer === null) {
            randomizer = new Randomizer();
        }
        var rows = A.shape[0];
        var cols = b.shape[1];
        var result = new Matrix(rows, 0);
        for (var i = 0; i < cols; ++i) {
            var b_i = Matrix.from(b.col(i)).T;
            var x = new Matrix(rows, 1, function () { return randomizer.random; });
            var r = b_i.sub(A.dot(x));
            var d = r.clone();
            do {
                var z = A.dot(d);
                var alpha = r.T.dot(r).entry(0, 0) / d.T.dot(z).entry(0, 0);
                x = x.add(d.mult(alpha));
                var r_next = r.sub(z.mult(alpha));
                var beta = r_next.T.dot(r_next).entry(0, 0) / r.T.dot(r).entry(0, 0);
                d = r_next.add(d.mult(beta));
                r = r_next;
            } while (Math.abs(r.mean) > tol);
            result = result.concat(x, "horizontal");
        }
        return result;
    };
    /**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    Matrix.solve = function (A, b) {
        var _a = "L" in A && "U" in A ? A : Matrix.LU(A), L = _a.L, U = _a.U;
        var rows = L.shape[0];
        var x = b.clone();
        // forward
        for (var row = 0; row < rows; ++row) {
            for (var col = 0; col < row - 1; ++col) {
                x.set_entry(0, row, x.entry(0, row) - L.entry(row, col) * x.entry(1, col));
            }
            x.set_entry(0, row, x.entry(0, row) / L.entry(row, row));
        }
        // backward
        for (var row = rows - 1; row >= 0; --row) {
            for (var col = rows - 1; col > row; --col) {
                x.set_entry(0, row, x.entry(0, row) - U.entry(row, col) * x.entry(0, col));
            }
            x.set_entry(0, row, x.entry(0, row) / U.entry(row, row));
        }
        return x;
    };
    /**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */
    Matrix.LU = function (A) {
        var rows = A.shape[0];
        var L = new Matrix(rows, rows, "zeros");
        var U = new Matrix(rows, rows, "identity");
        for (var j = 0; j < rows; ++j) {
            for (var i = j; i < rows; ++i) {
                var sum = 0;
                for (var k = 0; k < j; ++k) {
                    sum += L.entry(i, k) * U.entry(k, j);
                }
                L.set_entry(i, j, A.entry(i, j) - sum);
            }
            for (var i = j; i < rows; ++i) {
                if (L.entry(j, j) === 0) {
                    return undefined;
                }
                var sum = 0;
                for (var k = 0; k < j; ++k) {
                    sum += L.entry(j, k) * U.entry(k, i);
                }
                U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
            }
        }
        return { L: L, U: U };
    };
    /**
     * Computes the determinante of {@link A}, by using the LU decomposition of {@link A}.
     * @param {Matrix} A
     * @returns {Number} det - Returns the determinate of the Matrix {@link A}.
     */
    Matrix.det = function (A) {
        var rows = A.shape[0];
        var _a = Matrix.LU(A), L = _a.L, U = _a.U;
        var L_diag = L.diag;
        var U_diag = U.diag;
        var det = L_diag[0] * U_diag[0];
        for (var row = 1; row < rows; ++row) {
            det *= L_diag[row] * U_diag[row];
        }
        return det;
    };
    /**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M
     * @param {int} [k=2]
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */
    Matrix.SVD = function (M, k) {
        if (k === void 0) { k = 2; }
        var MT = M.T;
        var MtM = MT.dot(M);
        var MMt = M.dot(MT);
        var _a = simultaneous_poweriteration(MtM, k), V = _a.eigenvectors, Sigma = _a.eigenvalues;
        var U = simultaneous_poweriteration(MMt, k).eigenvectors;
        return { U: U, Sigma: Sigma.map(function (sigma) { return Math.sqrt(sigma); }), V: V };
        //Algorithm 1a: Householder reduction to bidiagonal form:
        /* const [m, n] = A.shape;
        let U = new Matrix(m, n, (i, j) => i == j ? 1 : 0);
        console.log(U.to2dArray)
        let V = new Matrix(n, m, (i, j) => i == j ? 1 : 0);
        console.log(V.to2dArray)
        let B = Matrix.bidiagonal(A.clone(), U, V);
        console.log(U,V,B)
        return { U: U, "Sigma": B, V: V }; */
    };
    return Matrix;
}());

/**
 * @class
 * @memberof module:utils
 * @alias Randomizer
 */
var Randomizer = /** @class */ (function () {
    /**
     * Mersenne Twister random number generator.
     * @constructor
     * @param {Number} [_seed=new Date().getTime()] - The seed for the random number generator. If <code>_seed == null</code> then the actual time gets used as seed.
     * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
     */
    function Randomizer(_seed) {
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
    Object.defineProperty(Randomizer.prototype, "seed", {
        /**
         * Returns the seed of the random number generator.
         * @returns {Number} - The seed.
         */
        get: function () {
            return this._seed;
        },
        set: function (_seed) {
            this._seed = _seed;
            var mt = this._mt;
            mt[0] = _seed >>> 0;
            for (this._mti = 1; this._mti < this._N; this._mti += 1) {
                var mti = this._mti;
                var s = mt[mti - 1] ^ (mt[mti - 1] >>> 30);
                mt[mti] = ((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253 + mti;
                mt[mti] >>>= 0;
            }
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Randomizer.prototype, "random", {
        /**
         * Returns a float between 0 and 1.
         * @returns {Number} - A random number between [0, 1]
         */
        get: function () {
            return this.random_int * (1.0 / 4294967296.0);
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Randomizer.prototype, "random_int", {
        /**
         * Returns an integer between 0 and MAX_INTEGER.
         * @returns {Integer} - A random integer.
         */
        get: function () {
            var y, mag01 = new Array(0x0, this._MATRIX_A);
            if (this._mti >= this._N) {
                var kk = void 0;
                /* if (this._mti == this._N + 1) {
                    this.seed = 5489;
                } */
                var N_M = this._N - this._M;
                var M_N = this._M - this._N;
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
        },
        enumerable: false,
        configurable: true
    });
    /**
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Number[]|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @returns {Number[]} - A random selection form {@link A} of {@link n} samples.
     */
    Randomizer.prototype.choice = function (A, n) {
        if (A instanceof Matrix) {
            var rows = A.shape[0];
            if (n > rows) {
                throw new Error("n bigger than A!");
            }
            var sample = new Array(n);
            var index_list = linspace(0, rows - 1);
            for (var i = 0, l = index_list.length; i < n; ++i, --l) {
                var random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map(function (d) { return A.row(d); });
        }
        else if (Array.isArray(A) || A instanceof Float64Array) {
            var rows = A.length;
            if (n > rows) {
                throw new Error("n bigger than A!");
            }
            var sample = new Array(n);
            var index_list = linspace(0, rows - 1);
            for (var i = 0, l = index_list.length; i < n; ++i, --l) {
                var random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map(function (d) { return A[d]; });
        }
    };
    /**
     * @static
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Number[]|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @param {Number|Randomizer} seed - The seed for the random number generator.
     * @returns {Number[]} - A random selection form {@link A} of {@link n} samples.
     */
    Randomizer.choice = function (A, n, seed) {
        if (seed === void 0) { seed = 1212; }
        var R = seed instanceof Randomizer ? seed : new Randomizer(seed);
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
    };
    return Randomizer;
}());

/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias max
 * @param {Array} values
 * @returns {Number}
 */
function max (values) {
    var max;
    for (var _i = 0, values_1 = values; _i < values_1.length; _i++) {
        var value = values_1[_i];
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
    var min;
    for (var _i = 0, values_1 = values; _i < values_1.length; _i++) {
        var value = values_1[_i];
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
var Heap = /** @class */ (function () {
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
    function Heap(elements, accessor, comparator) {
        if (elements === void 0) { elements = null; }
        if (accessor === void 0) { accessor = function (d) { return d; }; }
        if (comparator === void 0) { comparator = "min"; }
        if (elements) {
            return Heap.heapify(elements, accessor, comparator);
        }
        else {
            this._accessor = accessor;
            this._container = [];
            if (comparator == "min") {
                this._comparator = function (a, b) { return a < b; };
            }
            else if (comparator == "max") {
                this._comparator = function (a, b) { return a > b; };
            }
            else {
                this._comparator = comparator;
            }
            return this;
        }
    }
    /**
     * Creates a Heap from an Array
     * @param {Array|Set} elements - Contains the elements for the Heap.
     * @param {Function=} [accessor = (d) => d] - Function returns the value of the element.
     * @param {(String=|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     */
    Heap.heapify = function (elements, accessor, comparator) {
        if (accessor === void 0) { accessor = function (d) { return d; }; }
        if (comparator === void 0) { comparator = "min"; }
        var heap = new Heap(null, accessor, comparator);
        var container = heap._container;
        for (var _i = 0, elements_1 = elements; _i < elements_1.length; _i++) {
            var e = elements_1[_i];
            container.push({
                "element": e,
                "value": accessor(e)
            });
        }
        for (var i = Math.floor((elements.length / 2) - 1); i >= 0; --i) {
            heap._heapify_down(i);
        }
        return heap;
    };
    /**
     * Swaps elements of container array.
     * @private
     * @param {Number} index_a
     * @param {Number} index_b
     */
    Heap.prototype._swap = function (index_a, index_b) {
        var _a;
        var container = this._container;
        _a = [container[index_a], container[index_b]], container[index_b] = _a[0], container[index_a] = _a[1];
        return;
    };
    /**
     * @private
     */
    Heap.prototype._heapify_up = function () {
        var container = this._container;
        var index = container.length - 1;
        while (index > 0) {
            var parentIndex = Math.floor((index - 1) / 2);
            if (!this._comparator(container[index].value, container[parentIndex].value)) {
                break;
            }
            else {
                this._swap(parentIndex, index);
                index = parentIndex;
            }
        }
    };
    /**
     * Pushes the element to the heap.
     * @param {} element
     * @returns {Heap}
     */
    Heap.prototype.push = function (element) {
        var value = this._accessor(element);
        //const node = new Node(element, value);
        var node = { "element": element, "value": value };
        this._container.push(node);
        this._heapify_up();
        return this;
    };
    /**
     * @private
     * @param {Number} [start_index = 0]
     */
    Heap.prototype._heapify_down = function (start_index) {
        if (start_index === void 0) { start_index = 0; }
        var container = this._container;
        var comparator = this._comparator;
        var length = container.length;
        var left = 2 * start_index + 1;
        var right = 2 * start_index + 2;
        var index = start_index;
        if (index > length)
            throw "index higher than length";
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
    };
    /**
     * Removes and returns the top entry of the heap.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    Heap.prototype.pop = function () {
        var container = this._container;
        if (container.length === 0) {
            return null;
        }
        else if (container.length === 1) {
            return container.pop();
        }
        this._swap(0, container.length - 1);
        var item = container.pop();
        this._heapify_down();
        return item;
    };
    Object.defineProperty(Heap.prototype, "first", {
        /**
         * Returns the top entry of the heap without removing it.
         * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
         */
        get: function () {
            return this._container.length > 0 ? this._container[0] : null;
        },
        enumerable: false,
        configurable: true
    });
    /**
     * Yields the raw data
     * @yields {Object} Object consists of the element and its value (computed by {@link accessor}).
     */
    Heap.prototype.iterate = function () {
        var i, n;
        return tslib.__generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    i = 0, n = this._container.length;
                    _a.label = 1;
                case 1:
                    if (!(i < n)) return [3 /*break*/, 4];
                    return [4 /*yield*/, this._container[i].element];
                case 2:
                    _a.sent();
                    _a.label = 3;
                case 3:
                    ++i;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/];
            }
        });
    };
    /**
     * Returns the heap as ordered array.
     * @returns {Array} Array consisting the elements ordered by {@link comparator}.
     */
    Heap.prototype.toArray = function () {
        var _this = this;
        return this.data()
            .sort(function (a, b) { return _this._comparator(a, b) ? -1 : 0; });
    };
    /**
     * Returns elements of container array.
     * @returns {Array} Array consisting the elements.
     */
    Heap.prototype.data = function () {
        return this._container
            .map(function (d) { return d.element; });
    };
    /**
     * Returns the container array.
     * @returns {Array} The container array.
     */
    Heap.prototype.raw_data = function () {
        return this._container;
    };
    Object.defineProperty(Heap.prototype, "length", {
        /**
         * The size of the heap.
         * @returns {Number}
         */
        get: function () {
            return this._container.length;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Heap.prototype, "empty", {
        /**
         * Returns false if the the heap has entries, true if the heap has no entries.
         * @returns {Boolean}
         */
        get: function () {
            return this.length === 0;
        },
        enumerable: false,
        configurable: true
    });
    return Heap;
}());

/**
 * @class
 * @alias DisjointSet
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */
var DisjointSet = /** @class */ (function () {
    /**
     * @constructor
     * @alias DisjointSet
     * @memberof module:datastructure
     * @param {Array=} elements
     * @returns {DisjointSet}
     */
    function DisjointSet(elements) {
        if (elements === void 0) { elements = null; }
        this._list = new Set();
        if (elements) {
            for (var _i = 0, elements_1 = elements; _i < elements_1.length; _i++) {
                var e = elements_1[_i];
                this.make_set(e);
            }
        }
        return this;
    }
    DisjointSet.prototype.make_set = function (x) {
        var list = this._list;
        if (!list.has(x)) {
            list.add(x);
            x.__disjoint_set = {};
            x.__disjoint_set.parent = x;
            x.__disjoint_set.children = new Set([x]);
            x.__disjoint_set.size = 1;
        }
        return this;
    };
    DisjointSet.prototype.find = function (x) {
        var _a;
        var list = this._list;
        if (list.has(x)) {
            if (x.__disjoint_set.parent !== x) {
                (_a = x.__disjoint_set.children).add.apply(_a, x);
                x.__disjoint_set.parent = this.find(x.__disjoint_set.parent);
                return x.__disjoint_set.parent;
            }
            else {
                return x;
            }
        }
        else {
            return null;
        }
    };
    DisjointSet.prototype.union = function (x, y) {
        var _a;
        var node_x = this.find(x);
        var node_y = this.find(y);
        if (node_x === node_y)
            return this;
        if (node_x.__disjoint_set.size < node_y.__disjoint_set.size)
            _a = [node_y, node_x], node_x = _a[0], node_y = _a[1];
        node_y.__disjoint_set.parent = node_x;
        // keep track of children?
        node_y.__disjoint_set.children.forEach(node_x.__disjoint_set.children.add, node_x.__disjoint_set.children);
        node_x.__disjoint_set.size += node_y.__disjoint_set.size;
        return this;
    };
    return DisjointSet;
}());

/**
 * @class
 * @alias BallTree
 */
var BallTree = /** @class */ (function () {
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
    function BallTree(elements, metric) {
        if (elements === void 0) { elements = null; }
        if (metric === void 0) { metric = euclidean; }
        this._Node = /** @class */ (function () {
            function _Node(pivot, child1, child2, radius) {
                if (child1 === void 0) { child1 = null; }
                if (child2 === void 0) { child2 = null; }
                if (radius === void 0) { radius = null; }
                this.pivot = pivot;
                this.child1 = child1;
                this.child2 = child2;
                this.radius = radius;
            }
            return _Node;
        }());
        this._Leaf = /** @class */ (function () {
            function _Leaf(points) {
                this.points = points;
            }
            return _Leaf;
        }());
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
    BallTree.prototype.add = function (elements) {
        elements = elements.map(function (element, index) {
            return { index: index, element: element };
        });
        this._root = this._construct(elements);
        return this;
    };
    /**
     * @private
     * @param {Array<*>} elements
     * @returns {Node} root of balltree.
     */
    BallTree.prototype._construct = function (elements) {
        var _this = this;
        if (elements.length === 1) {
            return new this._Leaf(elements);
        }
        else {
            var c_1 = this._greatest_spread(elements);
            var sorted_elements = elements.sort(function (a, b) { return a.element[c_1] - b.element[c_1]; });
            var n = sorted_elements.length;
            var p_index = Math.floor(n / 2);
            var p_1 = elements[p_index];
            var L = sorted_elements.slice(0, p_index);
            var R = sorted_elements.slice(p_index, n);
            var radius = Math.max.apply(Math, elements.map(function (d) { return _this._metric(p_1.element, d.element); }));
            var B = void 0;
            if (L.length > 0 && R.length > 0) {
                B = new this._Node(p_1, this._construct(L), this._construct(R), radius);
            }
            else {
                B = new this._Leaf(elements);
            }
            return B;
        }
    };
    /**
     * @private
     * @param {Node} B
     * @returns {Number}
     */
    BallTree.prototype._greatest_spread = function (B) {
        var d = B[0].element.length;
        var start = new Array(d);
        for (var i = 0; i < d; ++i) {
            start[i] = [Infinity, -Infinity];
        }
        var spread = B.reduce(function (acc, current) {
            for (var i = 0; i < d; ++i) {
                acc[i][0] = Math.min(acc[i][0], current.element[i]);
                acc[i][1] = Math.max(acc[i][1], current.element[i]);
            }
            return acc;
        }, start);
        spread = spread.map(function (d) { return d[1] - d[0]; });
        var c = 0;
        for (var i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c;
    };
    /**
     *
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    BallTree.prototype.search = function (t, k) {
        var _this = this;
        if (k === void 0) { k = 5; }
        return this._search(t, k, new Heap(null, function (d) { return _this._metric(d.element, t); }, "max"), this._root);
    };
    /**
     * @private
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @param {Heap} Q - Heap consists of the currently found {@link k} nearest neighbors.
     * @param {Node|Leaf} B
     */
    BallTree.prototype._search = function (t, k, Q, B) {
        // B is Node
        if (Q.length >= k && B.pivot && B.radius && this._metric(t, B.pivot.element) - B.radius >= Q.first.value) {
            return Q;
        }
        if (B.child1)
            this._search(t, k, Q, B.child1);
        if (B.child2)
            this._search(t, k, Q, B.child2);
        // B is leaf
        if (B.points) {
            for (var i = 0, n = B.points.length; i < n; ++i) {
                var p = B.points[i];
                if (k > Q.length) {
                    Q.push(p);
                }
                else {
                    Q.push(p);
                    Q.pop();
                }
            }
        }
        return Q;
    };
    return BallTree;
}());

/**
 * @class
 * @alias KNN
 */
var KNN = /** @class */ (function () {
    /**
     * Generates a KNN list with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias KNN
     * @param {Array=} elements - Elements which should be added to the KNN list
     * @param {Function|"precomputed"} [metric = euclidean] metric is either precomputed or a function to use: (a, b) => distance
     * @returns {KNN}
     */
    function KNN(elements, metric) {
        if (elements === void 0) { elements = null; }
        if (metric === void 0) { metric = euclidean; }
        this._metric = metric;
        this._elements = elements instanceof Matrix ? elements : Matrix.from(elements);
        var N = this._elements.shape[0];
        if (metric === "precomputed") {
            this._D = this._elements.clone();
        }
        else {
            this._D = distance_matrix(this._elements, metric);
        }
        this.KNN = [];
        for (var row = 0; row < N; ++row) {
            var distances = this._D.row(row);
            var H = new Heap(null, function (d) { return d.value; }, "min");
            for (var j = 0; j < N; ++j) {
                H.push({
                    value: distances[j],
                    index: j
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
    KNN.prototype.search = function (t, k) {
        if (k === void 0) { k = 5; }
        var metric = this._metric;
        var KNN = this.KNN;
        var H;
        if (Array.isArray(t)) {
            if (this._metric == "precomputed") {
                throw "Search by query element is only possible when not using a precomputed distance matrix!";
            }
            var elements = this._elements;
            var N = KNN.length;
            var nearest_element_index = null;
            var nearest_dist = Infinity;
            for (var i = 0; i < N; ++i) {
                var element = elements.row(i);
                var dist = metric(t, element);
                if (dist < nearest_dist) {
                    nearest_element_index = i;
                    nearest_dist = dist;
                }
            }
            H = KNN[nearest_element_index];
        }
        else if (Number.isInteger(t)) {
            H = KNN[t];
        }
        var result = [];
        for (var i = 0; i < k; ++i) {
            result.push(H.pop());
        }
        result.forEach(function (res) { return H.push(res.element); });
        return result;
    };
    return KNN;
}());

/**
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
var DR = /** @class */ (function () {
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
    function DR(X, default_parameters, parameters) {
        var _a;
        this._parameters = Object.assign(Object.seal(default_parameters), parameters);
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        }
        else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        }
        else {
            throw new Error("No valid type for X!");
        }
        _a = this.X.shape, this._N = _a[0], this._D = _a[1];
        this._randomizer = new Randomizer(this._parameters.seed);
        this._is_initialized = false;
        return this;
    }
    /**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {any} [value = null] - value of the parameter to set.
     * @returns {DR|any} - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.
     * @example
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
     * DR.parameter("d"); // returns 3,
     * DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
     */
    DR.prototype.parameter = function (name, value) {
        if (value === void 0) { value = null; }
        if (!this._parameters.hasOwnProperty(name)) {
            throw new Error("".concat(name, " is not a valid parameter!"));
        }
        if (value) {
            this._parameters[name] = value;
            this._is_initialized = false;
            return this;
        }
        else {
            return this._parameters[name];
        }
    };
    DR.prototype.para = function (name, value) {
        if (value === void 0) { value = null; }
        return this.parameter(name, value);
    };
    DR.prototype.p = function (name, value) {
        if (value === void 0) { value = null; }
        return this.parameter(name, value);
    };
    /**
     * Computes the projection.
     * @returns {Matrix} - Returns the projection.
     */
    DR.prototype.transform = function () {
        this.check_init();
        return this.projection;
    };
    /**
     * Computes the projection.
     * @returns {Generator} - A generator yielding the intermediate steps of the dimensionality reduction method.
     */
    DR.prototype.generator = function () {
        return tslib.__generator(this, function (_a) {
            return [2 /*return*/, this.transform()];
        });
    };
    /**
     * If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.
     * @returns {DR}
     */
    DR.prototype.check_init = function () {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
        return this;
    };
    Object.defineProperty(DR.prototype, "projection", {
        /**
         * @returns {Matrix|Array} Returns the projection.
         */
        get: function () {
            if (this.hasOwnProperty("Y")) {
                this.check_init();
                return this._type === "matrix" ? this.Y : this.Y.to2dArray;
            }
            else {
                throw new Error("The dataset is not transformed yet!");
            }
        },
        enumerable: false,
        configurable: true
    });
    /**
     *
     * @param  {...any} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise} - A promise yielding the dimensionality reduced dataset.
     */
    DR.prototype.transform_async = function () {
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        return tslib.__awaiter(this, void 0, void 0, function () {
            return tslib.__generator(this, function (_a) {
                return [2 /*return*/, this.transform.apply(this, args)];
            });
        });
    };
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Matrix|Array} - The dimensionality reduced dataset.
     */
    DR.transform = function () {
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        var dr = new (this.bind.apply(this, tslib.__spreadArray([void 0], args, false)))();
        return dr.transform();
    };
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise} - A promise yielding the dimensionality reduced dataset.
     */
    DR.transform_async = function () {
        var args = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        return tslib.__awaiter(this, void 0, void 0, function () {
            return tslib.__generator(this, function (_a) {
                return [2 /*return*/, this.transform.apply(this, args)];
            });
        });
    };
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator} - A generator yielding the intermediate steps of the dimensionality reduction method.
     */
    DR.generator = function () {
        var _i, dr, generator, _a, generator_1, result;
        var args = [];
        for (_i = 0; _i < arguments.length; _i++) {
            args[_i] = arguments[_i];
        }
        return tslib.__generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    dr = new (this.bind.apply(this, tslib.__spreadArray([void 0], args, false)))();
                    generator = dr.generator();
                    _a = 0, generator_1 = generator;
                    _b.label = 1;
                case 1:
                    if (!(_a < generator_1.length)) return [3 /*break*/, 4];
                    result = generator_1[_a];
                    return [4 /*yield*/, result];
                case 2:
                    _b.sent();
                    _b.label = 3;
                case 3:
                    _a++;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/];
            }
        });
    };
    return DR;
}());

/**
 * @class
 * @alias PCA
 * @augments DR
 */
var PCA = /** @class */ (function (_super) {
    tslib.__extends(PCA, _super);
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
    function PCA(X, parameters) {
        var _this = _super.call(this, X, { d: 2, seed: 1212, eig_args: {} }, parameters) || this;
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        return _this;
    }
    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}. If parameter {@link A} is given, then project {@link A} with the principal components of {@link X}.
     * @param {null|Matrix|Array} [A = null] - If given, the data to project.
     * @returns {Matrix|Array} - The projected data.
     */
    PCA.prototype.transform = function (A) {
        if (A === void 0) { A = null; }
        var V = this.principal_components();
        if (A == null) {
            var X = this.X;
            this.Y = X.dot(V);
            return this.projection;
        }
        else if (Array.isArray(A)) {
            return Matrix.from(A).dot(V).asArray;
        }
        else if (A instanceof Matrix) {
            return A.dot(V);
        }
        else {
            throw new Error("No valid type for A!");
        }
    };
    /**
     * Computes the {@link d} principal components of Matrix {@link X}.
     * @returns {Matrix}
     */
    PCA.prototype.principal_components = function () {
        if (this.V) {
            return this.V;
        }
        var _a = this._parameters, d = _a.d, eig_args = _a.eig_args;
        var X = this.X;
        var means = Matrix.from(X.meanCols);
        var X_cent = X.sub(means);
        var C = X_cent.transpose().dot(X_cent);
        var V = simultaneous_poweriteration(C, d, eig_args).eigenvectors;
        this.V = Matrix.from(V).transpose();
        return this.V;
    };
    PCA.principal_components = function (X, parameters) {
        var dr = new this(X, parameters);
        return dr.principal_components();
    };
    return PCA;
}(DR));

/**
 * @class
 * @alias MDS
 * @extends DR
 */
var MDS = /** @class */ (function (_super) {
    tslib.__extends(MDS, _super);
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
    function MDS(X, parameters) {
        var _this = _super.call(this, X, { d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters) || this;
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        return _this;
    }
    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|Array}
     */
    MDS.prototype.transform = function () {
        var X = this.X;
        var rows = X.shape[0];
        var _a = this._parameters, d = _a.d, metric = _a.metric, eig_args = _a.eig_args;
        var A = metric === "precomputed" ? X : distance_matrix(X, metric);
        var ai_ = A.meanCols;
        var a_j = A.meanRows;
        var a__ = A.mean;
        this._d_X = A;
        var B = new Matrix(rows, rows, function (i, j) { return A.entry(i, j) - ai_[i] - a_j[j] + a__; });
        var V = simultaneous_poweriteration(B, d, eig_args).eigenvectors;
        this.Y = Matrix.from(V).transpose();
        return this.projection;
    };
    /**
     * @returns {Number} - the stress of the projection.
     */
    MDS.prototype.stress = function () {
        var N = this.X.shape[0];
        var Y = this.Y;
        var d_X = this._d_X;
        var d_Y = new Matrix();
        d_Y.shape = [
            N,
            N,
            function (i, j) {
                return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
            },
        ];
        var top_sum = 0;
        var bottom_sum = 0;
        for (var i = 0; i < N; ++i) {
            for (var j = i + 1; j < N; ++j) {
                top_sum += Math.pow(d_X.entry(i, j) - d_Y.entry(i, j), 2);
                bottom_sum += Math.pow(d_X.entry(i, j), 2);
            }
        }
        return Math.sqrt(top_sum / bottom_sum);
    };
    return MDS;
}(DR));

/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */
var ISOMAP = /** @class */ (function (_super) {
    tslib.__extends(ISOMAP, _super);
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
    function ISOMAP(X, parameters) {
        var _a;
        var _this = _super.call(this, X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters) || this;
        _this.parameter("neighbors", Math.min((_a = _this._parameters.neighbors) !== null && _a !== void 0 ? _a : Math.max(Math.floor(_this.X.shape[0] / 10), 2), _this._N - 1));
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        return _this;
    }
    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    ISOMAP.prototype.transform = function () {
        this.check_init();
        var X = this.X;
        var rows = this._N;
        var _a = this._parameters, d = _a.d, metric = _a.metric, eig_args = _a.eig_args, neighbors = _a.neighbors;
        // TODO: make knn extern and parameter for constructor or transform?
        var D = new Matrix();
        D.shape = [rows, rows, function (i, j) { return (i <= j ? metric(X.row(i), X.row(j)) : D.entry(j, i)); }];
        var kNearestNeighbors = [];
        for (var i = 0; i < rows; ++i) {
            var row = [];
            for (var j = 0; j < rows; ++j) {
                row.push({
                    index: j,
                    distance: D.entry(i, j)
                });
            }
            var H = new Heap(row, function (d) { return d.distance; }, "min");
            kNearestNeighbors.push(H.toArray().slice(1, neighbors + 1));
        }
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
        var G = new Matrix(rows, rows, function (i, j) {
            var other = kNearestNeighbors[i].find(function (n) { return n.index === j; });
            return other ? other.distance : Infinity;
        });
        for (var i = 0; i < rows; ++i) {
            for (var j = 0; j < rows; ++j) {
                for (var k = 0; k < rows; ++k) {
                    G.set_entry(i, j, Math.min(G.entry(i, j), G.entry(i, k) + G.entry(k, j)));
                }
            }
        }
        var ai_ = new Float64Array(rows);
        var a_j = new Float64Array(rows);
        var a__ = 0;
        var A = new Matrix(rows, rows, function (i, j) {
            var val = G.entry(i, j);
            val = val === Infinity ? 0 : val;
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        });
        ai_ = ai_.map(function (v) { return v / rows; });
        a_j = a_j.map(function (v) { return v / rows; });
        a__ /= Math.pow(rows, 2);
        var B = new Matrix(rows, rows, function (i, j) { return A.entry(i, j) - ai_[i] - a_j[j] + a__; });
        // compute d eigenvectors
        var V = simultaneous_poweriteration(B, d, eig_args).eigenvectors;
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.projection;
    };
    return ISOMAP;
}(DR));

/**
 * @class
 * @alias FASTMAP
 * @extends DR
 */
var FASTMAP = /** @class */ (function (_super) {
    tslib.__extends(FASTMAP, _super);
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
    function FASTMAP(X, parameters) {
        var _this = _super.call(this, X, { d: 2, metric: euclidean, seed: 1212 }, parameters) || this;
        return _this;
    }
    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {Function} dist
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */
    FASTMAP.prototype._choose_distant_objects = function (dist) {
        var X = this.X;
        var N = X.shape[0];
        var a_index = (this._randomizer.random_int % N) - 1;
        var b_index = null;
        var max_dist = -Infinity;
        for (var i = 0; i < N; ++i) {
            var d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        max_dist = -Infinity;
        for (var i = 0; i < N; ++i) {
            var d_bi = dist(b_index, i);
            if (d_bi > max_dist) {
                max_dist = d_bi;
                a_index = i;
            }
        }
        return [a_index, b_index, max_dist];
    };
    /**
     * Computes the projection.
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */
    FASTMAP.prototype.transform = function () {
        var X = this.X;
        var N = X.shape[0];
        var _a = this._parameters, d = _a.d, metric = _a.metric;
        var Y = new Matrix(N, d, 0);
        var dist = function (a, b) { return metric(X.row(a), X.row(b)); };
        var _loop_1 = function (_col) {
            var old_dist = dist;
            // choose pivot objects
            var _b = this_1._choose_distant_objects(dist), a_index = _b[0], b_index = _b[1], d_ab = _b[2];
            if (d_ab !== 0) {
                // project the objects on the line (O_a, O_b)
                for (var i = 0; i < N; ++i) {
                    var d_ai = dist(a_index, i);
                    var d_bi = dist(b_index, i);
                    var y_i = (Math.pow(d_ai, 2) + Math.pow(d_ab, 2) - Math.pow(d_bi, 2)) / (2 * d_ab);
                    Y.set_entry(i, _col, y_i);
                }
                // consider the projections of the objects on a
                // hyperplane perpendicluar to the line (a, b);
                // the distance function D'() between two
                // projections is given by Eq.4
                dist = function (a, b) { return Math.sqrt(Math.pow(old_dist(a, b), 2) - Math.pow((Y.entry(a, _col) - Y.entry(b, _col)), 2)); };
            }
        };
        var this_1 = this;
        for (var _col = 0; _col < d; ++_col) {
            _loop_1(_col);
        }
        // return embedding.
        this.Y = Y;
        return this.projection;
    };
    return FASTMAP;
}(DR));

/**
 * @class
 * @alias LDA
 * @extends DR
 */
var LDA = /** @class */ (function (_super) {
    tslib.__extends(LDA, _super);
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
    function LDA(X, parameters) {
        var _this = _super.call(this, X, { labels: null, d: 2, seed: 1212, eig_args: {} }, parameters) || this;
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        return _this;
    }
    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    LDA.prototype.transform = function () {
        var X = this.X;
        var _a = X.shape, rows = _a[0], cols = _a[1];
        var _b = this._parameters, d = _b.d, labels = _b.labels, eig_args = _b.eig_args;
        if (labels === null || labels.length != rows) {
            throw new Error("LDA needs parameter label to every datapoint to work!");
        }
        var unique_labels = {};
        var label_id = 0;
        labels.forEach(function (l, i) {
            if (l in unique_labels) {
                unique_labels[l].count++;
                unique_labels[l].rows.push(X.row(i));
            }
            else {
                unique_labels[l] = {
                    id: label_id++,
                    count: 1,
                    rows: [X.row(i)]
                };
            }
        });
        // create X_mean and vector means;
        var X_mean = X.mean;
        var V_mean = new Matrix(label_id, cols);
        for (var label in unique_labels) {
            var V_1 = Matrix.from(unique_labels[label].rows);
            var v_mean = V_1.meanCols;
            for (var j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }
        }
        // scatter_between
        var S_b = new Matrix(cols, cols);
        var _loop_1 = function (label) {
            var v = V_mean.row(unique_labels[label].id);
            var m = new Matrix(cols, 1, function (j) { return v[j] - X_mean; });
            var N = unique_labels[label].count;
            S_b = S_b.add(m.dot(m.transpose()).mult(N));
        };
        for (var label in unique_labels) {
            _loop_1(label);
        }
        // scatter_within
        var S_w = new Matrix(cols, cols);
        var _loop_2 = function (label) {
            var v = V_mean.row(unique_labels[label].id);
            var m = new Matrix(cols, 1, function (j) { return v[j]; });
            var R = unique_labels[label].rows;
            var _loop_3 = function (i, n) {
                var row_v = new Matrix(cols, 1, function (j, _) { return R[i][j] - m.entry(j, 0); });
                S_w = S_w.add(row_v.dot(row_v.transpose()));
            };
            for (var i = 0, n = unique_labels[label].count; i < n; ++i) {
                _loop_3(i);
            }
        };
        for (var label in unique_labels) {
            _loop_2(label);
        }
        var V = simultaneous_poweriteration(S_w.inverse().dot(S_b), d, eig_args).eigenvectors;
        V = Matrix.from(V).transpose();
        this.Y = X.dot(V);
        // return embedding
        return this.projection;
    };
    return LDA;
}(DR));

/**
 * @class
 * @alias LLE
 * @extends DR
 */
var LLE = /** @class */ (function (_super) {
    tslib.__extends(LLE, _super);
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
    function LLE(X, parameters) {
        var _a;
        var _this = _super.call(this, X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters) || this;
        _this.parameter("neighbors", Math.min((_a = parameters.neighbors) !== null && _a !== void 0 ? _a : Math.max(Math.floor(_this._N / 10), 2), _this._N - 1));
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        return _this;
    }
    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    LLE.prototype.transform = function () {
        var X = this.X;
        var rows = this._N;
        var cols = this._D;
        var _a = this._parameters, neighbors = _a.neighbors, d = _a.d, eig_args = _a.eig_args, metric = _a.metric;
        var nN = k_nearest_neighbors(X, neighbors, metric);
        var O = new Matrix(neighbors, 1, 1);
        var W = new Matrix(rows, rows);
        var _loop_1 = function (row) {
            var nN_row = nN[row];
            var Z = new Matrix(neighbors, cols, function (i, j) { return X.entry(nN_row[i].j, j) - X.entry(row, j); });
            var C = Z.dot(Z.T);
            if (neighbors > cols) {
                var C_trace = neumair_sum(C.diag) / 1000;
                for (var j = 0; j < neighbors; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }
            // reconstruct;
            var w = Matrix.solve_CG(C, O, this_1._randomizer);
            w = w.divide(w.sum);
            for (var j = 0; j < neighbors; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        };
        var this_1 = this;
        for (var row = 0; row < rows; ++row) {
            _loop_1(row);
        }
        // comp embedding
        var I = new Matrix(rows, rows, "identity");
        var IW = I.sub(W);
        var M = IW.T.dot(IW);
        var V = simultaneous_poweriteration(M.T.inverse(), d + 1, eig_args).eigenvectors;
        this.Y = Matrix.from(V.slice(1, 1 + d)).T;
        // return embedding
        return this.projection;
    };
    return LLE;
}(DR));

/**
 * @class
 * @alias LTSA
 * @extends DR
 */
var LTSA = /** @class */ (function (_super) {
    tslib.__extends(LTSA, _super);
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
    function LTSA(X, parameters) {
        var _a;
        var _this = _super.call(this, X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters) || this;
        _this.parameter("neighbors", Math.min((_a = parameters.neighbors) !== null && _a !== void 0 ? _a : Math.max(Math.floor(_this._N / 10), 2), _this._N - 1));
        if (!_this._parameters.eig_args.hasOwnProperty("seed")) {
            _this._parameters.eig_args.seed = _this._randomizer;
        }
        if (_this._D <= _this.parameter("d")) {
            throw new Error("Dimensionality of X (D = ".concat(_this._D, ") must be greater than the required dimensionality of the result (d = ").concat(_this.parameter("d"), ")!"));
        }
        return _this;
    }
    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    LTSA.prototype.transform = function () {
        var X = this.X;
        var _a = X.shape, rows = _a[0], D = _a[1];
        var _b = this._parameters, d = _b.d, neighbors = _b.neighbors, metric = _b.metric, eig_args = _b.eig_args;
        // 1.1 determine k nearest neighbors
        var nN = k_nearest_neighbors(X, neighbors, metric);
        // center matrix
        var O = new Matrix(D, D, "center");
        var B = new Matrix(rows, rows, 0);
        for (var row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            var I_i = tslib.__spreadArray([row], nN[row].map(function (n) { return n.j; }), true);
            var X_i = Matrix.from(I_i.map(function (n) { return X.row(n); }));
            // center X_i
            X_i = X_i.dot(O);
            // correlation matrix
            var C = X_i.dot(X_i.transpose());
            var g = simultaneous_poweriteration(C, d, eig_args).eigenvectors;
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            var G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            var W_i = G_i_t.transpose()
                .dot(G_i_t)
                .add(1 / Math.sqrt(neighbors + 1));
            for (var i = 0; i < neighbors + 1; ++i) {
                for (var j = 0; j < neighbors + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0) + W_i.entry(i, j));
                }
            }
        }
        // 3. Aligning global coordinates
        var Y = simultaneous_poweriteration(B, d + 1, eig_args).eigenvectors;
        this.Y = Matrix.from(Y.slice(1)).transpose();
        // return embedding
        return this.projection;
    };
    return LTSA;
}(DR));

/**
 * @class
 * @alias TSNE
 * @extends DR
 */
var TSNE = /** @class */ (function (_super) {
    tslib.__extends(TSNE, _super);
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
    function TSNE(X, parameters) {
        var _a;
        var _this = _super.call(this, X, { perplexity: 50, epsilon: 10, d: 2, metric: euclidean, seed: 1212 }, parameters) || this;
        _a = _this.X.shape, _this._N = _a[0], _this._D = _a[1];
        _this._iter = 0;
        _this.Y = new Matrix(_this._N, _this.parameter("d"), function () { return _this._randomizer.random; });
        return _this;
    }
    /**
     *
     * @param {Matrix} distance_matrix - accepts a precomputed distance matrix
     * @returns {TSNE}
     */
    TSNE.prototype.init = function () {
        // init
        var Htarget = Math.log(this.parameter("perplexity"));
        var N = this._N;
        var D = this._D;
        var metric = this._parameters.metric;
        var X = this.X;
        var Delta;
        if (metric == "precomputed") {
            Delta = druid.Matrix.from(X);
        }
        else {
            Delta = new Matrix(N, N);
            for (var i = 0; i < N; ++i) {
                var X_i = X.row(i);
                for (var j = i + 1; j < N; ++j) {
                    var distance = metric(X_i, X.row(j));
                    Delta.set_entry(i, j, distance);
                    Delta.set_entry(j, i, distance);
                }
            }
        }
        var P = new Matrix(N, N, "zeros");
        this._ystep = new Matrix(N, D, "zeros");
        this._gains = new Matrix(N, D, 1);
        // search for fitting sigma
        var prow = new Float64Array(N);
        var tol = 1e-4;
        var maxtries = 50;
        for (var i = 0; i < N; ++i) {
            var betamin = -Infinity;
            var betamax = Infinity;
            var beta = 1;
            var done = false;
            var num = 0;
            while (!done) {
                var psum = 0;
                for (var j = 0; j < N; ++j) {
                    var pj = Math.exp(-Delta.entry(i, j) * beta);
                    if (i === j)
                        pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                var Hhere = 0;
                for (var j = 0; j < N; ++j) {
                    var pj = psum === 0 ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) {
                        Hhere -= pj * Math.log(pj);
                    }
                }
                if (Hhere > Htarget) {
                    betamin = beta;
                    beta = betamax === Infinity ? beta * 2 : (beta + betamax) / 2;
                }
                else {
                    betamax = beta;
                    beta = betamin === -Infinity ? beta / 2 : (beta + betamin) / 2;
                }
                ++num;
                if (Math.abs(Hhere - Htarget) < tol)
                    done = true;
                if (num >= maxtries)
                    done = true;
            }
            for (var j = 0; j < N; ++j) {
                P.set_entry(i, j, prow[j]);
            }
        }
        //compute probabilities
        var Pout = new Matrix(N, N, "zeros");
        var N2 = N * 2;
        for (var i = 0; i < N; ++i) {
            for (var j = i; j < N; ++j) {
                var p = Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100);
                Pout.set_entry(i, j, p);
                Pout.set_entry(j, i, p);
            }
        }
        this._P = Pout;
        return this;
    };
    /**
     *
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */
    TSNE.prototype.transform = function (iterations) {
        if (iterations === void 0) { iterations = 500; }
        this.check_init();
        for (var i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    };
    /**
     *
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */
    TSNE.prototype.generator = function (iterations) {
        var i;
        if (iterations === void 0) { iterations = 500; }
        return tslib.__generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    this.check_init();
                    i = 0;
                    _a.label = 1;
                case 1:
                    if (!(i < iterations)) return [3 /*break*/, 4];
                    this.next();
                    return [4 /*yield*/, this.projection];
                case 2:
                    _a.sent();
                    _a.label = 3;
                case 3:
                    ++i;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/, this.projection];
            }
        });
    };
    /**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */
    TSNE.prototype.next = function () {
        var iter = ++this._iter;
        var P = this._P;
        var ystep = this._ystep;
        var gains = this._gains;
        var N = this._N;
        var _a = this._parameters, dim = _a.d, epsilon = _a.epsilon;
        var Y = this.Y;
        //calc cost gradient;
        var pmul = iter < 100 ? 4 : 1;
        // compute Q dist (unnormalized)
        var Qu = new Matrix(N, N, "zeros");
        var qsum = 0;
        for (var i = 0; i < N; ++i) {
            for (var j = i + 1; j < N; ++j) {
                var dsum = 0;
                for (var d = 0; d < dim; ++d) {
                    var dhere = Y.entry(i, d) - Y.entry(j, d);
                    dsum += dhere * dhere;
                }
                var qu = 1 / (1 + dsum);
                Qu.set_entry(i, j, qu);
                Qu.set_entry(j, i, qu);
                qsum += 2 * qu;
            }
        }
        // normalize Q dist
        var Q = new Matrix(N, N, 0);
        for (var i = 0; i < N; ++i) {
            for (var j = i + 1; j < N; ++j) {
                var val = Math.max(Qu.entry(i, j) / qsum, 1e-100);
                Q.set_entry(i, j, val);
                Q.set_entry(j, i, val);
            }
        }
        var grad = new Matrix(N, dim, "zeros");
        for (var i = 0; i < N; ++i) {
            for (var j = 0; j < N; ++j) {
                var premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (var d = 0; d < dim; ++d) {
                    grad.set_entry(i, d, grad.entry(i, d) + premult * (Y.entry(i, d) - Y.entry(j, d)));
                }
            }
        }
        // perform gradient step
        var ymean = new Float64Array(dim);
        for (var i = 0; i < N; ++i) {
            for (var d = 0; d < dim; ++d) {
                var gid = grad.entry(i, d);
                var sid = ystep.entry(i, d);
                var gainid = gains.entry(i, d);
                var newgain = Math.sign(gid) === Math.sign(sid) ? gainid * 0.8 : gainid + 0.2;
                if (newgain < 0.01)
                    newgain = 0.01;
                gains.set_entry(i, d, newgain);
                var momval = iter < 250 ? 0.5 : 0.8;
                var newsid = momval * sid - epsilon * newgain * gid;
                ystep.set_entry(i, d, newsid);
                Y.set_entry(i, d, Y.entry(i, d) + newsid);
                ymean[d] += Y.entry(i, d);
            }
        }
        for (var i = 0; i < N; ++i) {
            for (var d = 0; d < 2; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - ymean[d] / N);
            }
        }
        return this.Y;
    };
    return TSNE;
}(DR));

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
function powell (f, x0, max_iter) {
    if (max_iter === void 0) { max_iter = 300; }
    var epsilon = 1e-2;
    var n = x0.length;
    var alpha = 1e-3;
    var pfx = 10000;
    var x = x0.slice();
    var fx = f(x);
    var convergence = false;
    while (max_iter-- >= 0 && !convergence) {
        convergence = true;
        for (var i = 0; i < n; ++i) {
            x[i] += 1e-6;
            var fxi = f(x);
            x[i] -= 1e-6;
            var dx = (fxi - fx) / 1e-6;
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
var UMAP = /** @class */ (function (_super) {
    tslib.__extends(UMAP, _super);
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
    function UMAP(X, parameters) {
        var _c;
        var _this = _super.call(this, X, { n_neighbors: 15, local_connectivity: 1, min_dist: 1, d: 2, metric: euclidean, seed: 1212, _spread: 1, _set_op_mix_ratio: 1, _repulsion_strength: 1, _negative_sample_rate: 5, _n_epochs: 350, _initial_alpha: 1 }, parameters) || this;
        _c = _this.X.shape, _this._N = _c[0], _this._D = _c[1];
        /* let n_neighbors = Math.min(this._N - 1, parameters.n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(this.parameter("local_connectivity"), n_neighbors - 1)); */
        if (_this.parameter("n_neighbors") > _this._N) {
            throw new Error("Parameter n_neighbors (=".concat(_this.parameter("n_neighbors"), ") needs to be smaller than dataset size (N=").concat(_this._N, ")!"));
        }
        if (_this.parameter("local_connectivity") > _this.parameter("n_neighbors")) {
            throw new Error("Parameter local_connectivity (=".concat(_this.parameter("local_connectivity"), ") needs to be smaller than parameter n_neighbors (=").concat(_this.parameter("n_neighbors"), ")"));
        }
        _this._iter = 0;
        var randomizer = _this._randomizer;
        _this.Y = new Matrix(_this._N, _this.parameter("d"), function () { return randomizer.random; });
        return _this;
    }
    /**
     * @private
     * @param {Number} spread
     * @param {Number} min_dist
     * @returns {Array}
     */
    UMAP.prototype._find_ab_params = function (spread, min_dist) {
        var curve = function (x, a, b) { return 1 / (1 + a * Math.pow(x, 2 * b)); };
        var xv = linspace(0, spread * 3, 300);
        var yv = linspace(0, spread * 3, 300);
        for (var i = 0, n = xv.length; i < n; ++i) {
            var xv_i = xv[i];
            yv[i] = xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread);
        }
        var err = function (p) {
            var error = linspace(1, 300).map(function (_, i) { return yv[i] - curve(xv[i], p[0], p[1]); });
            return Math.sqrt(neumair_sum(error.map(function (e) { return e * e; })));
        };
        return powell(err, [1, 1]);
    };
    /**
     * @private
     * @param {Array<Array>} distances
     * @param {Array<Number>} sigmas
     * @param {Array<Number>} rhos
     * @returns {Array}
     */
    UMAP.prototype._compute_membership_strengths = function (distances, sigmas, rhos) {
        for (var i = 0, n = distances.length; i < n; ++i) {
            for (var j = 0, m = distances[i].length; j < m; ++j) {
                var v = distances[i][j].value - rhos[i];
                distances[i][j].value = v > 0 ? Math.exp(-v / sigmas[i]) : 1;
            }
        }
        return distances;
    };
    /**
     * @private
     * @param {KNN|BallTree} knn
     * @param {Number} k
     * @returns {Object}
     */
    UMAP.prototype._smooth_knn_dist = function (knn, k) {
        var _c, _d, _e;
        var SMOOTH_K_TOLERANCE = 1e-5;
        var MIN_K_DIST_SCALE = 1e-3;
        var n_iter = 64;
        var _f = this._parameters, local_connectivity = _f.local_connectivity, metric = _f.metric;
        var target = Math.log2(k);
        var rhos = [];
        var sigmas = [];
        var X = this.X;
        var N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());
        var distances = [];
        if (metric === "precomputed") {
            for (var i = 0; i < N; ++i) {
                distances.push(knn.search(i, k).reverse());
            }
        }
        else {
            for (var _i = 0, X_1 = X; _i < X_1.length; _i++) {
                var x_i = X_1[_i];
                distances.push(knn.search(x_i, k).raw_data().reverse());
            }
        }
        for (var i = 0; i < N; ++i) {
            var lo = 0;
            var hi = Infinity;
            var mid = 1;
            var search_result = distances[i];
            var non_zero_dist = search_result.filter(function (d) { return d.value > 0; });
            var non_zero_dist_length = non_zero_dist.length;
            if (non_zero_dist_length >= local_connectivity) {
                var index = Math.floor(local_connectivity);
                var interpolation = local_connectivity - index;
                if (index > 0) {
                    rhos.push(non_zero_dist[index - 1]);
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rhos[i].value += interpolation * (non_zero_dist[index].value - non_zero_dist[index - 1]);
                    }
                }
                else {
                    rhos[i].value = interpolation * non_zero_dist[0].value;
                }
            }
            else if (non_zero_dist_length > 0) {
                rhos[i] = non_zero_dist[non_zero_dist_length - 1].value;
            }
            for (var x = 0; x < n_iter; ++x) {
                var psum = 0;
                for (var j = 0; j < k; ++j) {
                    var d = search_result[j].value - rhos[i];
                    psum += d > 0 ? Math.exp(-(d / mid)) : 1;
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    _c = [mid, (lo + hi) / 2], hi = _c[0], mid = _c[1];
                }
                else {
                    if (hi === Infinity) {
                        _d = [mid, mid * 2], lo = _d[0], mid = _d[1];
                    }
                    else {
                        _e = [mid, (lo + hi) / 2], lo = _e[0], mid = _e[1];
                    }
                }
            }
            sigmas[i] = mid;
            var mean_ithd = search_result.reduce(function (a, b) { return a + b.value; }, 0) / search_result.length;
            //let mean_d = null;
            if (rhos[i] > 0) {
                if (sigmas[i] < MIN_K_DIST_SCALE * mean_ithd) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_ithd;
                }
            }
            else {
                var mean_d = distances.reduce(function (acc, res) { return acc + res.reduce(function (a, b) { return a + b.value; }, 0) / res.length; });
                if (sigmas[i] > MIN_K_DIST_SCALE * mean_d) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_d;
                }
            }
        }
        return {
            distances: distances,
            sigmas: sigmas,
            rhos: rhos
        };
    };
    /**
     * @private
     * @param {Matrix} X
     * @param {Number} n_neighbors
     * @returns {Matrix}
     */
    UMAP.prototype._fuzzy_simplicial_set = function (X, n_neighbors) {
        var N = X.shape[0];
        var _c = this._parameters, metric = _c.metric, _set_op_mix_ratio = _c._set_op_mix_ratio;
        var knn = metric === "precomputed" ? new KNN(X, "precomputed") : new BallTree(X.to2dArray, metric);
        var _d = this._smooth_knn_dist(knn, n_neighbors), distances = _d.distances, sigmas = _d.sigmas, rhos = _d.rhos;
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        var result = new Matrix(N, N, "zeros");
        for (var i = 0; i < N; ++i) {
            var distances_i = distances[i];
            for (var j = 0; j < distances_i.length; ++j) {
                result.set_entry(i, distances_i[j].element.index, distances_i[j].value);
            }
        }
        var transposed_result = result.T;
        var prod_matrix = result.mult(transposed_result);
        return result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(_set_op_mix_ratio)
            .add(prod_matrix.mult(1 - _set_op_mix_ratio));
    };
    /**
     * @private
     * @param {Number} n_epochs
     * @returns {Array}
     */
    UMAP.prototype._make_epochs_per_sample = function (n_epochs) {
        var weights = this._weights;
        var result = new Float32Array(weights.length).fill(-1);
        var weights_max = max(weights);
        var n_samples = weights.map(function (w) { return n_epochs * (w / weights_max); });
        for (var i = 0; i < result.length; ++i)
            if (n_samples[i] > 0)
                result[i] = Math.round(n_epochs / n_samples[i]);
        return result;
    };
    /**
     * @private
     * @param {Matrix} graph
     * @returns {Object}
     */
    UMAP.prototype._tocoo = function (graph) {
        var rows = [];
        var cols = [];
        var data = [];
        var _c = graph.shape, rows_n = _c[0], cols_n = _c[1];
        for (var row = 0; row < rows_n; ++row) {
            for (var col = 0; col < cols_n; ++col) {
                var entry = graph.entry(row, col);
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
            data: data
        };
    };
    /**
     * Computes all necessary
     * @returns {UMAP}
     */
    UMAP.prototype.init = function () {
        var _c = this._parameters, _spread = _c._spread, min_dist = _c.min_dist, n_neighbors = _c.n_neighbors, _n_epochs = _c._n_epochs, _negative_sample_rate = _c._negative_sample_rate;
        var _d = this._find_ab_params(_spread, min_dist), a = _d[0], b = _d[1];
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, n_neighbors);
        var _e = this._tocoo(this._graph), rows = _e.rows, cols = _e.cols, weights = _e.data;
        this._head = rows;
        this._tail = cols;
        this._weights = weights;
        this._epochs_per_sample = this._make_epochs_per_sample(_n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map(function (d) { return d * _negative_sample_rate; });
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        return this;
    };
    UMAP.prototype.graph = function () {
        this.check_init();
        return { cols: this._head, rows: this._tail, weights: this._weights };
    };
    /**
     *
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    UMAP.prototype.transform = function (iterations) {
        if (iterations === void 0) { iterations = 350; }
        if (this.parameter("_n_epochs") != iterations) {
            this.parameter("_n_epochs", iterations);
            this.init();
        }
        this.check_init();
        for (var i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    };
    /**
     *
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    UMAP.prototype.generator = function (iterations) {
        var i;
        if (iterations === void 0) { iterations = 350; }
        return tslib.__generator(this, function (_c) {
            switch (_c.label) {
                case 0:
                    if (this.parameter("_n_epochs") != iterations) {
                        this.parameter("_n_epochs", iterations);
                        this.init();
                    }
                    this.check_init();
                    i = 0;
                    _c.label = 1;
                case 1:
                    if (!(i < iterations)) return [3 /*break*/, 4];
                    this.next();
                    return [4 /*yield*/, this.projection];
                case 2:
                    _c.sent();
                    _c.label = 3;
                case 3:
                    ++i;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/, this.projection];
            }
        });
    };
    /**
     * @private
     * @param {Number} x
     * @returns {Number}
     */
    UMAP.prototype._clip = function (x) {
        if (x > 4)
            return 4;
        if (x < -4)
            return -4;
        return x;
    };
    /**
     * performs the optimization step.
     * @private
     * @param {Matrix} head_embedding
     * @param {Matrix} tail_embedding
     * @param {Matrix} head
     * @param {Matrix} tail
     * @returns {Matrix}
     */
    UMAP.prototype._optimize_layout = function (head_embedding, tail_embedding, head, tail) {
        var randomizer = this._randomizer;
        var _c = this._parameters, _repulsion_strength = _c._repulsion_strength, dim = _c.d;
        var _d = this, alpha = _d._alpha, a = _d._a, b = _d._b, epochs_per_sample = _d._epochs_per_sample, epochs_per_negative_sample = _d._epochs_per_negative_sample, epoch_of_next_negative_sample = _d._epoch_of_next_negative_sample, epoch_of_next_sample = _d._epoch_of_next_sample, clip = _d._clip;
        var tail_length = tail.length;
        for (var i = 0, n = epochs_per_sample.length; i < n; ++i) {
            if (epoch_of_next_sample[i] <= this._iter) {
                var j = head[i];
                var k = tail[i];
                var current = head_embedding.row(j);
                var other = tail_embedding.row(k);
                var dist = euclidean_squared(current, other);
                var grad_coeff = 0;
                if (dist > 0) {
                    grad_coeff = (-2 * a * b * Math.pow(dist, b - 1)) / (a * Math.pow(dist, b) + 1);
                }
                for (var d = 0; d < dim; ++d) {
                    var grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                    var c = current[d] + grad_d;
                    var o = other[d] - grad_d;
                    current[d] = c;
                    other[d] = o;
                    head_embedding.set_entry(j, d, c);
                    tail_embedding.set_entry(k, d, o);
                }
                epoch_of_next_sample[i] += epochs_per_sample[i];
                var n_neg_samples = (this._iter - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i];
                for (var p = 0; p < n_neg_samples; ++p) {
                    var k_1 = randomizer.random_int % tail_length;
                    var other_1 = tail_embedding.row(tail[k_1]);
                    var dist_1 = euclidean_squared(current, other_1);
                    var grad_coeff_1 = 0;
                    if (dist_1 > 0) {
                        grad_coeff_1 = (2 * _repulsion_strength * b) / ((0.01 + dist_1) * (a * Math.pow(dist_1, b) + 1));
                    }
                    else if (j === k_1) {
                        continue;
                    }
                    for (var d = 0; d < dim; ++d) {
                        var grad_d = clip(grad_coeff_1 * (current[d] - other_1[d])) * alpha;
                        var c = current[d] + grad_d;
                        var o = other_1[d] - grad_d;
                        current[d] = c;
                        other_1[d] = o;
                        head_embedding.set_entry(j, d, c);
                        tail_embedding.set_entry(tail[k_1], d, o);
                    }
                }
                epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
            }
        }
        return head_embedding;
    };
    /**
     * @private
     * @returns {Matrix}
     */
    UMAP.prototype.next = function () {
        var iter = ++this._iter;
        var Y = this.Y;
        var _c = this._parameters, _initial_alpha = _c._initial_alpha, _n_epochs = _c._n_epochs;
        this._alpha = _initial_alpha * (1 - iter / _n_epochs);
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);
        return this.Y;
    };
    return UMAP;
}(DR));

/**
 * @class
 * @alias TriMap
 * @extends DR
 */
var TriMap = /** @class */ (function (_super) {
    tslib.__extends(TriMap, _super);
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
    function TriMap(X, parameters) {
        var _this = _super.call(this, X, { weight_adj: 500, c: 5, d: 2, metric: euclidean, tol: 1e-8, seed: 1212 }, parameters) || this;
        return _this;
    }
    /**
     *
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used).
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used).
     */
    TriMap.prototype.init = function (pca, knn) {
        if (pca === void 0) { pca = null; }
        if (knn === void 0) { knn = null; }
        var X = this.X;
        var N = X.shape[0];
        var _a = this._parameters, d = _a.d, metric = _a.metric, c = _a.c;
        this.n_inliers = 2 * c;
        this.n_outliers = 1 * c;
        this.n_random = 1 * c;
        this.Y = pca || new PCA(X, d).transform();
        this.knn = knn || new BallTree(X.to2dArray, metric);
        var _b = this._generate_triplets(this.n_inliers, this.n_outliers, this.n_random), triplets = _b.triplets, weights = _b.weights;
        this.triplets = triplets;
        this.weights = weights;
        this.lr = (1000 * N) / triplets.shape[0];
        this.C = Infinity;
        this.vel = new Matrix(N, d, 0);
        this.gain = new Matrix(N, d, 1);
        return this;
    };
    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {Number} n_inliers
     * @param {Number} n_outliers
     * @param {Number} n_random
     */
    TriMap.prototype._generate_triplets = function (n_inliers, n_outliers, n_random) {
        var _a = this._parameters, metric = _a.metric, weight_adj = _a.weight_adj;
        var X = this.X;
        var N = X.shape[0];
        var knn = this.knn;
        var n_extra = Math.min(n_inliers + 20, N);
        var nbrs = new Matrix(N, n_extra);
        var knn_distances = new Matrix(N, n_extra);
        var _loop_1 = function (i) {
            knn.search(X.row(i), n_extra + 1)
                .raw_data()
                .filter(function (d) { return d.value != 0; })
                .sort(function (a, b) { return a.value - b.value; })
                .forEach(function (d, j) {
                nbrs.set_entry(i, j, d.element.index);
                knn_distances.set_entry(i, j, d.value);
            });
        };
        for (var i = 0; i < N; ++i) {
            _loop_1(i);
        }
        // scale parameter
        var sig = new Float64Array(N);
        for (var i = 0; i < N; ++i) {
            sig[i] = Math.max((knn_distances.entry(i, 3) + knn_distances.entry(i, 4) + knn_distances.entry(i, 5) + knn_distances.entry(i, 6)) / 4, 1e-10);
        }
        var P = this._find_p(knn_distances, sig, nbrs);
        var triplets = this._sample_knn_triplets(P, nbrs, n_inliers, n_outliers);
        var n_triplets = triplets.shape[0];
        var outlier_distances = new Float64Array(n_triplets);
        for (var i = 0; i < n_triplets; ++i) {
            var j = triplets.entry(i, 0);
            var k = triplets.entry(i, 2);
            outlier_distances[i] = metric(X.row(j), X.row(k));
        }
        var weights = this._find_weights(triplets, P, nbrs, outlier_distances, sig);
        if (n_random > 0) {
            var _b = this._sample_random_triplets(X, n_random, sig), random_triplets = _b.random_triplets, random_weights = _b.random_weights;
            triplets = triplets.concat(random_triplets, "vertical");
            weights = Float64Array.from(tslib.__spreadArray(tslib.__spreadArray([], weights, true), random_weights, true));
        }
        n_triplets = triplets.shape[0];
        var max_weight = -Infinity;
        for (var i = 0; i < n_triplets; ++i) {
            if (isNaN(weights[i])) {
                weights[i] = 0;
            }
            if (max_weight < weights[i])
                max_weight = weights[i];
        }
        var max_weight_2 = -Infinity;
        for (var i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight;
            weights[i] += 0.0001;
            weights[i] = Math.log(1 + weight_adj * weights[i]);
            if (max_weight_2 < weights[i])
                max_weight_2 = weights[i];
        }
        for (var i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight_2;
        }
        return {
            triplets: triplets,
            weights: weights
        };
    };
    /**
     * Calculates the similarity matrix P
     * @private
     * @param {Matrix} knn_distances - matrix of pairwise knn distances
     * @param {Float64Array} sig - scaling factor for the distances
     * @param {Matrix} nbrs - nearest neighbors
     * @returns {Matrix} pairwise similarity matrix
     */
    TriMap.prototype._find_p = function (knn_distances, sig, nbrs) {
        var _a = knn_distances.shape, N = _a[0], n_neighbors = _a[1];
        return new Matrix(N, n_neighbors, function (i, j) {
            return Math.exp(-(Math.pow(knn_distances.entry(i, j), 2) / sig[i] / sig[nbrs.entry(i, j)]));
        });
    };
    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix {@link P}. Row i corresponds to the i-th point.
     * @param {Number} n_inliers - Number of inlier points.
     * @param {Number} n_outliers - Number of outlier points.
     *
     */
    TriMap.prototype._sample_knn_triplets = function (P, nbrs, n_inliers, n_outliers) {
        var N = nbrs.shape[0];
        var triplets = new Matrix(N * n_inliers * n_outliers, 3);
        for (var i = 0; i < N; ++i) {
            var n_i = i * n_inliers * n_outliers;
            var sort_indices = this.__argsort(P.row(i).map(function (d) { return -d; }));
            for (var j = 0; j < n_inliers; ++j) {
                var n_j = j * n_outliers;
                var sim = nbrs.entry(i, sort_indices[j]);
                var samples = this._rejection_sample(n_outliers, N, sort_indices.slice(0, j + 1));
                for (var k = 0; k < n_outliers; ++k) {
                    var index = n_i + n_j + k;
                    var out = samples[k];
                    triplets.set_entry(index, 0, i);
                    triplets.set_entry(index, 1, sim);
                    triplets.set_entry(index, 2, out);
                }
            }
        }
        return triplets;
    };
    /**
     * Should do the same as np.argsort()
     * @private
     * @param {Array} A
     */
    TriMap.prototype.__argsort = function (A) {
        return A.map(function (d, i) {
            return { d: d, i: i };
        })
            .sort(function (a, b) { return a.d - b.d; })
            .map(function (d) { return d.i; });
    };
    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are in the {@link rejects}.
     * @private
     * @param {*} n_samples
     * @param {*} max_int
     * @param {*} rejects
     */
    TriMap.prototype._rejection_sample = function (n_samples, max_int, rejects) {
        var randomizer = this._randomizer;
        var interval = linspace(0, max_int - 1).filter(function (d) { return rejects.indexOf(d) < 0; });
        return randomizer.choice(interval, Math.min(n_samples, interval.length - 2));
    };
    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - scaling factor for the distances.
     */
    TriMap.prototype._find_weights = function (triplets, P, nbrs, outlier_distances, sig) {
        var n_triplets = triplets.shape[0];
        var weights = new Float64Array(n_triplets);
        for (var t = 0; t < n_triplets; ++t) {
            var i = triplets.entry(t, 0);
            var sim = nbrs.row(i).indexOf(triplets.entry(t, 1));
            var p_sim = P.entry(i, sim);
            var p_out = Math.exp(-(Math.pow(outlier_distances[t], 2) / (sig[i] * sig[triplets.entry(t, 2)])));
            if (p_out < 1e-20)
                p_out = 1e-20;
            weights[t] = p_sim / p_out;
        }
        return weights;
    };
    /**
     * Sample uniformly ranom triplets
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {Number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    TriMap.prototype._sample_random_triplets = function (X, n_random, sig) {
        var _a, _b;
        var metric = this.parameter("metric");
        var randomizer = this._randomizer;
        var N = X.shape[0];
        var random_triplets = new Matrix(N * n_random, 3);
        var random_weights = new Float64Array(N * n_random);
        for (var i = 0; i < N; ++i) {
            var n_i = i * n_random;
            var indices = tslib.__spreadArray(tslib.__spreadArray([], linspace(0, i - 1), true), linspace(i + 1, N - 1), true);
            for (var j = 0; j < n_random; ++j) {
                var _c = randomizer.choice(indices, 2), sim = _c[0], out = _c[1];
                var p_sim = Math.exp(-(Math.pow(metric(X.row(i), X.row(sim)), 2) / (sig[i] * sig[sim])));
                if (p_sim < 1e-20)
                    p_sim = 1e-20;
                var p_out = Math.exp(-(Math.pow(metric(X.row(i), X.row(out)), 2) / (sig[i] * sig[out])));
                if (p_out < 1e-20)
                    p_out = 1e-20;
                if (p_sim < p_out) {
                    _a = [out, sim], sim = _a[0], out = _a[1];
                    _b = [p_out, p_sim], p_sim = _b[0], p_out = _b[1];
                }
                var index = n_i + j;
                random_triplets.set_entry(index, 0, i);
                random_triplets.set_entry(index, 1, sim);
                random_triplets.set_entry(index, 2, out);
                random_weights[index] = p_sim / p_out;
            }
        }
        return {
            random_triplets: random_triplets,
            random_weights: random_weights
        };
    };
    /**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */
    TriMap.prototype._grad = function (Y) {
        var n_inliers = this.n_inliers;
        var n_outliers = this.n_outliers;
        var triplets = this.triplets;
        var weights = this.weights;
        var _a = Y.shape, N = _a[0], dim = _a[1];
        var n_triplets = triplets.shape[0];
        var grad = new Matrix(N, dim, 0);
        var y_ij = new Float64Array(dim);
        var y_ik = new Float64Array(dim);
        var d_ij = 1;
        var d_ik = 1;
        var n_viol = 0;
        var loss = 0;
        var n_knn_triplets = N * n_inliers * n_outliers;
        for (var t = 0; t < n_triplets; ++t) {
            var _b = triplets.row(t), i = _b[0], j = _b[1], k = _b[2];
            // update y_ij, y_ik, d_ij, d_ik
            if (t % n_outliers == 0 || t >= n_knn_triplets) {
                d_ij = 1;
                d_ik = 1;
                for (var d = 0; d < dim; ++d) {
                    var Y_id = Y.entry(i, d);
                    var Y_jd = Y.entry(j, d);
                    var Y_kd = Y.entry(k, d);
                    y_ij[d] = Y_id - Y_jd;
                    y_ik[d] = Y_id - Y_kd;
                    d_ij += Math.pow(y_ij[d], 2);
                    d_ik += Math.pow(y_ik[d], 2);
                }
                // update y_ik and d_ik only
            }
            else {
                d_ik = 1;
                for (var d = 0; d < dim; ++d) {
                    var Y_id = Y.entry(i, d);
                    var Y_kd = Y.entry(k, d);
                    y_ik[d] = Y_id - Y_kd;
                    d_ik += Math.pow(y_ik[d], 2);
                }
            }
            if (d_ij > d_ik)
                ++n_viol;
            loss += weights[t] / (1 + d_ik / d_ij);
            var w = Math.pow((weights[t] / (d_ij + d_ik)), 2);
            for (var d = 0; d < dim; ++d) {
                var gs = y_ij[d] * d_ik * w;
                var go = y_ik[d] * d_ij * w;
                grad.set_entry(i, d, grad.entry(i, d) + gs - go);
                grad.set_entry(j, d, grad.entry(j, d) - gs);
                grad.set_entry(k, d, grad.entry(k, d) + go);
            }
        }
        return { grad: grad, loss: loss, n_viol: n_viol };
    };
    /**
     *
     * @param {Number} max_iteration
     */
    TriMap.prototype.transform = function (max_iteration) {
        if (max_iteration === void 0) { max_iteration = 400; }
        this.check_init();
        for (var iter = 0; iter < max_iteration; ++iter) {
            this._next(iter);
        }
        return this.projection;
    };
    /**
     * @param {Number} max_iteration
     * @yields {Matrix}
     * @returns {Matrix}
     */
    TriMap.prototype.generator = function (max_iteration) {
        var iter;
        if (max_iteration === void 0) { max_iteration = 800; }
        return tslib.__generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    this.check_init();
                    iter = 0;
                    _a.label = 1;
                case 1:
                    if (!(iter < max_iteration)) return [3 /*break*/, 4];
                    this._next(iter);
                    return [4 /*yield*/, this.projection];
                case 2:
                    _a.sent();
                    _a.label = 3;
                case 3:
                    ++iter;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/, this.projection];
            }
        });
    };
    /**
     * Does the iteration step.
     * @private
     * @param {Number} iter
     */
    TriMap.prototype._next = function (iter) {
        var gamma = iter > 150 ? 0.5 : 0.3;
        var old_C = this.C;
        var vel = this.vel;
        var Y = this.Y.add(vel.mult(gamma));
        var _a = this._grad(Y), grad = _a.grad, loss = _a.loss; _a.n_viol;
        this.C = loss;
        this.Y = this._update_embedding(Y, iter, grad);
        this.lr *= old_C > loss + this._parameters.tol ? 1.01 : 0.9;
        return this.Y;
    };
    /**
     * Updates the embedding.
     * @private
     * @param {Matrix} Y
     * @param {Number} iter
     * @param {Matrix} grad
     */
    TriMap.prototype._update_embedding = function (Y, iter, grad) {
        var _a = Y.shape, N = _a[0], dim = _a[1];
        var gamma = iter > 150 ? 0.9 : 0.5; // moment parameter
        var min_gain = 0.01;
        var gain = this.gain;
        var vel = this.vel;
        var lr = this.lr;
        for (var i = 0; i < N; ++i) {
            for (var d = 0; d < dim; ++d) {
                var new_gain = Math.sign(vel.entry(i, d)) != Math.sign(grad.entry(i, d)) ? gain.entry(i, d) + 0.2 : Math.max(gain.entry(i, d) * 0.8, min_gain);
                gain.set_entry(i, d, new_gain);
                vel.set_entry(i, d, gamma * vel.entry(i, d) - lr * gain.entry(i, d) * grad.entry(i, d));
                Y.set_entry(i, d, Y.entry(i, d) + vel.entry(i, d));
            }
        }
        return Y;
    };
    return TriMap;
}(DR));

/**
 * @class
 * @alias Hierarchical_Clustering
 */
var Hierarchical_Clustering = /** @class */ (function () {
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
    function Hierarchical_Clustering(matrix, linkage, metric) {
        if (linkage === void 0) { linkage = "complete"; }
        if (metric === void 0) { metric = euclidean; }
        this._id = 0;
        this._matrix = matrix instanceof Matrix ? matrix : Matrix.from(matrix);
        this._metric = metric;
        this._linkage = linkage;
        if (metric === "precomputed" && this._matrix.shape[0] !== this._matrix.shape[1]) {
            throw new Error("If metric is 'precomputed', then matrix has to be square!");
        }
        this.init();
        this.root = this["do"]();
        return this;
    }
    /**
     *
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    Hierarchical_Clustering.prototype.get_clusters = function (value, type) {
        if (type === void 0) { type = "distance"; }
        var clusters = [];
        var accessor;
        switch (type) {
            case "distance":
                accessor = function (d) { return d.dist; };
                break;
            case "depth":
                accessor = function (d) { return d.depth; };
                break;
            default:
                throw new Error("invalid type");
        }
        this._traverse(this.root, accessor, value, clusters);
        return clusters;
    };
    /**
     * @private
     * @param {} node
     * @param {*} f
     * @param {*} value
     * @param {*} result
     */
    Hierarchical_Clustering.prototype._traverse = function (node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        }
        else {
            this._traverse(node.left, f, value, result);
            this._traverse(node.right, f, value, result);
        }
    };
    /**
     * computes the tree.
     */
    Hierarchical_Clustering.prototype.init = function () {
        var metric = this._metric;
        var A = this._matrix;
        var n = (this._n = A.shape[0]);
        var d_min = (this._d_min = new Float64Array(n));
        var distance_matrix;
        if (metric !== "precomputed") {
            distance_matrix = new Matrix(n, n, 0); //new Array(n);
            for (var i = 0; i < n; ++i) {
                d_min[i] = 0;
                //distance_matrix[i] = new Float64Array(n);
                for (var j = 0; j < n; ++j) {
                    distance_matrix.set_entry(i, j, i === j ? Infinity : metric(A.row(i), A.row(j)));
                    if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        else {
            distance_matrix = this._matrix.clone();
            for (var i = 0; i < n; ++i) {
                for (var j = 0; j < n; ++j) {
                    if (i === j) {
                        distance_matrix.set_entry(i, j, Infinity);
                    }
                    else if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        this._distance_matrix = distance_matrix;
        var clusters = (this._clusters = new Array(n));
        var c_size = (this._c_size = new Uint16Array(n));
        for (var i = 0; i < n; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        return this;
    };
    /**
     * computes the tree.
     */
    Hierarchical_Clustering.prototype["do"] = function () {
        var n = this._n;
        var d_min = this._d_min;
        var D = this._distance_matrix;
        var clusters = this._clusters;
        var c_size = this._c_size;
        var linkage = this._linkage;
        var root = null;
        for (var p = 0, p_max = n - 1; p < p_max; ++p) {
            var c1 = 0;
            for (var i = 0; i < n; ++i) {
                var D_i_min = D.entry(i, d_min[i]);
                for (var j = i + 1; j < n; ++j) {
                    if (D_i_min > D.entry(i, j)) {
                        d_min[i] = j;
                        D_i_min = D.entry(i, d_min[i]);
                    }
                }
            }
            for (var i = 0; i < n; ++i) {
                if (D.entry(i, d_min[i]) < D.entry(c1, d_min[c1])) {
                    c1 = i;
                }
            }
            var c2 = d_min[c1];
            var c1_cluster = clusters[c1][0];
            var c2_cluster = clusters[c2][0];
            var c1_cluster_indices = c1_cluster.isLeaf ? [c1_cluster.index] : c1_cluster.index;
            var c2_cluster_indices = c2_cluster.isLeaf ? [c2_cluster.index] : c2_cluster.index;
            var indices = c1_cluster_indices.concat(c2_cluster_indices);
            var new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D.entry(c1, c2), null, indices);
            c1_cluster.parent = new_cluster;
            c2_cluster.parent = new_cluster;
            clusters[c1].unshift(new_cluster);
            c_size[c1] += c_size[c2];
            for (var j = 0; j < n; ++j) {
                var D_c1_j = D.entry(c1, j);
                var D_c2_j = D.entry(c2, j);
                var value = void 0;
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
            for (var i = 0; i < n; ++i) {
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
    };
    return Hierarchical_Clustering;
}());
var Cluster = /** @class */ (function () {
    function Cluster(id, left, right, dist, centroid, index, size, depth) {
        this.id = id;
        this.left = left;
        this.right = right;
        this.dist = dist;
        this.index = index;
        this.size = size !== null && size !== void 0 ? size : left.size + right.size;
        this.depth = depth !== null && depth !== void 0 ? depth : 1 + Math.max(left.depth, right.depth);
        this.centroid = centroid !== null && centroid !== void 0 ? centroid : this._calculate_centroid(left, right);
        this.parent = null;
        return this;
    }
    Cluster.prototype._calculate_centroid = function (left, right) {
        var l_size = left.size;
        var r_size = right.size;
        var l_centroid = left.centroid;
        var r_centroid = right.centroid;
        var size = this.size;
        var n = left.centroid.length;
        var new_centroid = new Float64Array(n);
        for (var i = 0; i < n; ++i) {
            new_centroid[i] = (l_size * l_centroid[i] + r_size * r_centroid[i]) / size;
        }
        return new_centroid;
    };
    Object.defineProperty(Cluster.prototype, "isLeaf", {
        get: function () {
            return this.depth === 0;
        },
        enumerable: false,
        configurable: true
    });
    Cluster.prototype.leaves = function () {
        if (this.isLeaf)
            return [this];
        var left = this.left;
        var right = this.right;
        return (left.isLeaf ? [left] : left.leaves()).concat(right.isLeaf ? [right] : right.leaves());
    };
    Cluster.prototype.descendants = function () {
        if (this.isLeaf)
            return [this];
        var left_descendants = this.left.descendants();
        var right_descendants = this.right.descendants();
        return left_descendants.concat(right_descendants).concat([this]);
    };
    return Cluster;
}());

/**
 * @class
 * @alias KMeans
 */
var KMeans = /** @class */ (function () {
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
    function KMeans(matrix, K, metric, seed, init) {
        if (metric === void 0) { metric = euclidean; }
        if (seed === void 0) { seed = 1987; }
        if (init === void 0) { init = true; }
        this._metric = metric;
        this._matrix = matrix;
        this._K = K;
        var _a = matrix.shape, N = _a[0], D = _a[1];
        this._N = N;
        this._D = D;
        if (K > N)
            K = N;
        this._randomizer = new Randomizer(seed);
        this._clusters = new Array(N).fill(undefined);
        this._cluster_centroids = this._get_random_centroids(K);
        if (init)
            this.init(K, this._cluster_centroids);
        return this;
    }
    /**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    KMeans.prototype.get_clusters = function () {
        var K = this._K;
        var clusters = this._clusters;
        var result = new Array(K).fill().map(function () { return new Array(); });
        clusters.forEach(function (c, i) { return result[c].push(i); });
        return result;
    };
    /**
     * @private
     * @param {Array} points
     * @param {Array} candidates
     */
    KMeans.prototype._furthest_point = function (points, candidates) {
        var A = this._matrix;
        var metric = this._metric;
        var i = points.length;
        var H = Heap.heapify(candidates, function (d) {
            var Ad = A.row(d);
            var sum = 0;
            for (var j = 0; j < i; ++j) {
                sum += metric(Ad, points[j]);
            }
            return sum;
        }, "max");
        return H.pop().element;
    };
    KMeans.prototype._get_random_centroids = function (K) {
        var N = this._N;
        var randomizer = this._randomizer;
        var A = this._matrix;
        var cluster_centroids = new Array(K).fill();
        var indices = linspace(0, N - 1);
        var random_point = randomizer.random_int % (N - 1);
        cluster_centroids[0] = A.row(random_point);
        var init_points = [random_point];
        var sample_size = Math.floor((N - K) / K); // / K
        for (var i = 1; i < K; ++i) {
            // sampling + kmeans++ improvement?
            var sample = randomizer.choice(indices.filter(function (d) { return init_points.indexOf(d) == -1; }), sample_size);
            var furthest_point = this._furthest_point(cluster_centroids.slice(0, i), sample);
            init_points.push(furthest_point);
            cluster_centroids[i] = A.row(furthest_point);
        }
        return cluster_centroids;
    };
    KMeans.prototype._iteration = function (cluster_centroids) {
        var K = cluster_centroids.length;
        var N = this._N;
        var D = this._D;
        var A = this._matrix;
        var metric = this._metric;
        var clusters = this._clusters;
        var clusters_changed = false;
        // find nearest cluster centroid.
        for (var i = 0; i < N; ++i) {
            var Ai = A.row(i);
            var min_dist = Infinity;
            var min_cluster = null;
            for (var j = 0; j < K; ++j) {
                var d = metric(cluster_centroids[j], Ai);
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
        for (var i = 0; i < K; ++i) {
            var centroid = cluster_centroids[i];
            for (var j = 0; j < D; ++j) {
                centroid[j] = 0;
            }
        }
        // compute centroid
        this._compute_centroid(cluster_centroids);
        return {
            "clusters_changed": clusters_changed,
            "cluster_centroids": cluster_centroids
        };
    };
    KMeans.prototype._compute_centroid = function (cluster_centroids) {
        var K = cluster_centroids.length;
        var N = this._N;
        var D = this._D;
        var A = this._matrix;
        var clusters = this._clusters;
        var cluster_counter = new Array(K).fill(0);
        for (var i = 0; i < N; ++i) {
            var Ai = A.row(i);
            var ci = clusters[i];
            cluster_counter[ci]++;
            var centroid = cluster_centroids[ci];
            for (var j = 0; j < D; ++j) {
                centroid[j] += Ai[j];
            }
        }
        var _loop_1 = function (i) {
            var n = cluster_counter[i];
            cluster_centroids[i] = cluster_centroids[i].map(function (c) { return c / n; });
        };
        for (var i = 0; i < K; ++i) {
            _loop_1(i);
        }
    };
    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    KMeans.prototype.init = function (K, cluster_centroids) {
        if (!K)
            K = this._K;
        if (!cluster_centroids)
            cluster_centroids = this._get_random_centroids(K);
        var clusters_changed = false;
        do {
            var iteration_result = this._iteration(cluster_centroids);
            cluster_centroids = iteration_result.cluster_centroids;
            clusters_changed = iteration_result.clusters_changed;
        } while (clusters_changed);
    };
    return KMeans;
}());

/**
 * @class
 * @alias KMedoids
 */
var KMedoids = /** @class */ (function () {
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
    function KMedoids(matrix, K, max_iter, metric, seed) {
        if (max_iter === void 0) { max_iter = null; }
        if (metric === void 0) { metric = euclidean; }
        if (seed === void 0) { seed = 1212; }
        this._metric = metric;
        this._matrix = matrix;
        this._A = this._matrix.to2dArray;
        this._K = K;
        var _a = matrix.shape, N = _a[0], D = _a[1];
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
        if (K > N)
            K = N;
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
    KMedoids.prototype.get_clusters = function () {
        var _this = this;
        var K = this._K;
        var A = this._A;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        var result = new Array(K).fill().map(function () { return new Array(); });
        A.forEach(function (x_j, j) {
            result[_this._nearest_medoid(x_j, j).index_nearest].push(j);
        });
        result.medoids = this._cluster_medoids;
        return result;
    };
    KMedoids.prototype.generator = function () {
        return tslib.__asyncGenerator(this, arguments, function generator_1() {
            var max_iter, finish, i;
            return tslib.__generator(this, function (_a) {
                switch (_a.label) {
                    case 0:
                        max_iter = this._max_iter;
                        return [4 /*yield*/, tslib.__await(this.get_clusters())];
                    case 1: return [4 /*yield*/, _a.sent()];
                    case 2:
                        _a.sent();
                        finish = false;
                        i = 0;
                        _a.label = 3;
                    case 3:
                        finish = this._iteration();
                        return [4 /*yield*/, tslib.__await(this.get_clusters())];
                    case 4: return [4 /*yield*/, _a.sent()];
                    case 5:
                        _a.sent();
                        _a.label = 6;
                    case 6:
                        if (!finish && ++i < max_iter) return [3 /*break*/, 3];
                        _a.label = 7;
                    case 7: return [2 /*return*/];
                }
            });
        });
    };
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
    KMedoids.prototype._iteration = function () {
        var _this = this;
        var A = this._A;
        var K = this._K;
        var medoids = this._cluster_medoids;
        var cache = A.map(function (x_o, o) { return _this._nearest_medoid(x_o, o); });
        // empty best candidates array
        var DeltaTD = new Array(K).fill(0);
        var xs = new Array(K).fill(null);
        A.forEach(function (x_j, j) {
            if (medoids.findIndex(function (m) { return m === j; }) < 0) {
                var d_j = cache[j].distance_nearest; // distance to current medoid
                var deltaTD_1 = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach(function (x_o, o) {
                    if (j === o)
                        return;
                    var d_oj = _this._get_distance(o, j, x_o, x_j); // distance to new medoid
                    var _a = cache[o], n = _a["index_nearest"], d_n = _a["distance_nearest"], d_s = _a["distance_second"]; // cached
                    deltaTD_1[n] += Math.min(d_oj, d_s) - d_n; // loss change for x_o
                    // Reassignment check
                    if (d_oj < d_n) {
                        // update loss change
                        for (var i = 0; i < K; ++i) {
                            if (i !== n)
                                deltaTD_1[i] += d_oj - d_n;
                        }
                    }
                });
                // remember best swap for i;
                deltaTD_1
                    .map(function (d, i) { return [d, i]; })
                    .filter(function (_a) {
                    var d = _a[0], i = _a[1];
                    return d < DeltaTD[i];
                })
                    .forEach(function (_a) {
                    var d = _a[0], i = _a[1];
                    if (d < DeltaTD[i]) {
                        DeltaTD[i] = d;
                        xs[i] = j;
                    }
                });
            }
        });
        // stop if no improvements were found
        if (min(DeltaTD) >= 0)
            return true;
        var _loop_1 = function () {
            // swap roles of medoid m_i and non_medoid xs_i
            var i = DeltaTD
                .map(function (d, i) { return [d, i]; })
                .sort(function (_a, _b) {
                var a = _a[0];
                var b = _b[0];
                return a - b;
            })[0][1];
            if (medoids.filter(function (m) { return m == xs[i]; }).length == 0) {
                medoids[i] = xs[i];
            }
            // disable the swap just performed
            DeltaTD[i] = 0;
            // recompute TD for remaining swap candidates
            DeltaTD
                .map(function (d_j, j) { return [d_j, j]; })
                .filter(function (_a) {
                var d_j = _a[0];
                return d_j < 0;
            })
                .forEach(function (_a) {
                _a[0]; var j = _a[1];
                var x_j = A[j];
                var sum = 0;
                A.forEach(function (x_o, o) {
                    if (medoids.findIndex(function (m) { return m != j && m == o; }) >= 0)
                        return;
                    if (i == j)
                        return;
                    if (cache[o].index_nearest === medoids[j])
                        sum += (Math.min(_this._get_distance(o, j, x_o, x_j), cache[o].distance_second) - cache[o].distance_nearest);
                    else {
                        sum += (Math.min(_this._get_distance(o, j, x_o, x_j) - cache[o].distance_nearest, 0));
                    }
                });
                DeltaTD[j] = sum;
            });
        };
        // execute all improvements
        while (min(DeltaTD) < 0) {
            _loop_1();
        }
        this._cluster_medoids = medoids;
        return false;
    };
    KMedoids.prototype._get_distance = function (i, j, x_i, x_j) {
        if (x_i === void 0) { x_i = null; }
        if (x_j === void 0) { x_j = null; }
        if (i === j)
            return 0;
        var D = this._distance_matrix;
        var A = this._A;
        var metric = this._metric;
        var d_ij = D.entry(i, j);
        if (d_ij === 0) {
            d_ij = metric(x_i || A[i], x_j || A[j]);
            D.set_entry(i, j, d_ij);
            D.set_entry(j, i, d_ij);
        }
        return d_ij;
    };
    KMedoids.prototype._nearest_medoid = function (x_j, j) {
        var _this = this;
        var medoids = this._cluster_medoids;
        var A = this._A;
        var _a = medoids
            .map(function (m, i) {
            var x_m = A[m];
            return [_this._get_distance(j, m, x_j, x_m), i];
        })
            .sort(function (m1, m2) { return m1[0] - m2[0]; }), nearest = _a[0], second = _a[1];
        return {
            "distance_nearest": nearest[0],
            "index_nearest": nearest[1],
            "distance_second": second[0],
            "index_second": second[1]
        };
    };
    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    KMedoids.prototype.init = function (K, cluster_medoids) {
        if (!K)
            K = this._K;
        if (!cluster_medoids)
            cluster_medoids = this._get_random_medoids(K);
        var max_iter = this._max_iter;
        var finish = false;
        var i = 0;
        do {
            finish = this._iteration();
        } while (!finish && ++i < max_iter);
        return this;
    };
    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @param {number} K - number of clusters
     *
     */
    KMedoids.prototype._get_random_medoids = function (K) {
        var _this = this;
        var N = this._N;
        var A = this._A;
        var indices = linspace(0, N - 1);
        var randomizer = this._randomizer;
        var n = Math.min(N, 10 + Math.ceil(Math.sqrt(N)));
        var TD = new Array(n).fill(Infinity);
        var medoids = [];
        // first medoid
        var TD0 = Infinity;
        var S = randomizer.choice(indices, n);
        for (var j = 0; j < n; ++j) {
            var S_j = S[j];
            var x_j = A[S_j];
            for (var o = 0; o < n; ++o) {
                if (o === j)
                    continue;
                var x_o = A[S[o]];
                TD[j] += this._get_distance(j, o, x_j, x_o);
            }
            if (TD[j] < TD0) {
                TD0 = TD[j]; // smallest distance sum
                medoids.push(S_j);
            }
        }
        // other medoids
        for (var i = 1; i < K; ++i) {
            var DeltaTD = Infinity;
            S = randomizer.choice(indices.filter(function (index) { return medoids.findIndex(function (d) { return d === index; }) < 0; }), n);
            for (var j = 0; j < n; ++j) {
                var deltaTD = 0;
                var S_j = S[j];
                var x_j = A[S_j];
                var _loop_2 = function (o) {
                    if (o === j)
                        return "continue";
                    var S_o = S[o];
                    var x_o = A[S_o];
                    var delta = this_1._get_distance(S_j, S_o, x_j, x_o) - min(medoids.map(function (m) { return _this._get_distance(S_o, m, x_o); }));
                    if (delta < 0) {
                        deltaTD = deltaTD + delta;
                    }
                };
                var this_1 = this;
                for (var o = 0; o < n; ++o) {
                    _loop_2(o);
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
    };
    return KMedoids;
}());

/**
 * @class
 * @alias OPTICS
 */
var OPTICS = /** @class */ (function () {
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
    function OPTICS(matrix, epsilon, min_points, metric) {
        if (metric === void 0) { metric = euclidean; }
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
    OPTICS.prototype.init = function () {
        var ordered_list = this._ordered_list;
        var matrix = this._matrix;
        var N = matrix.shape[0];
        var DB = this._DB;
        var clusters = this._clusters;
        var cluster_index = this._cluster_index = 0;
        for (var i = 0; i < N; ++i) {
            DB[i] = {
                "element": matrix.row(i),
                "index": i,
                "reachability_distance": undefined,
                "processed": false
            };
        }
        for (var _i = 0, DB_1 = DB; _i < DB_1.length; _i++) {
            var p = DB_1[_i];
            if (p.processed)
                continue;
            p.neighbors = this._get_neighbors(p);
            p.processed = true;
            clusters.push([p.index]);
            cluster_index = clusters.length - 1;
            ordered_list.push(p);
            if (this._core_distance(p) != undefined) {
                var seeds = new Heap(null, function (d) { return d.reachability_distance; }, "min");
                this._update(p, seeds);
                this._expand_cluster(seeds, clusters[cluster_index]);
            }
        }
        return this;
    };
    /**
     *
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Array} An array consisting of the {@link epsilon}-neighborhood of {@link p}.
     */
    OPTICS.prototype._get_neighbors = function (p) {
        if ("neighbors" in p)
            return p.neighbors;
        var DB = this._DB;
        var metric = this._metric;
        var epsilon = this._epsilon;
        var neighbors = [];
        for (var _i = 0, DB_2 = DB; _i < DB_2.length; _i++) {
            var q = DB_2[_i];
            if (q.index == p.index)
                continue;
            if (metric(p.element, q.element) < epsilon) {
                neighbors.push(q);
            }
        }
        return neighbors;
    };
    /**
     *
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Number} The distance to the {@link min_points}-th nearest point of {@link p}, or undefined if the {@link epsilon}-neighborhood has fewer elements than {@link min_points}.
     */
    OPTICS.prototype._core_distance = function (p) {
        var min_points = this._min_points;
        var metric = this._metric;
        if (p.neighbors && p.neighbors.length <= min_points) {
            return undefined;
        }
        return metric(p.element, p.neighbors[min_points].element);
    };
    /**
     * Updates the reachability distance of the points.
     * @private
     * @param {Object} p
     * @param {Heap} seeds
     */
    OPTICS.prototype._update = function (p, seeds) {
        var metric = this._metric;
        var core_distance = this._core_distance(p);
        var neighbors = this._get_neighbors(p); //p.neighbors;
        var _loop_1 = function (q) {
            if (q.processed)
                return "continue";
            var new_reachability_distance = Math.max(core_distance, metric(p.element, q.element));
            //if (q.reachability_distance == undefined) { // q is not in seeds
            if (seeds.raw_data().findIndex(function (d) { return d.element == q; }) < 0) {
                q.reachability_distance = new_reachability_distance;
                seeds.push(q);
            }
            else { // q is in seeds
                if (new_reachability_distance < q.reachability_distance) {
                    q.reachability_distance = new_reachability_distance;
                    seeds = Heap.heapify(seeds.data(), function (d) { return d.reachability_distance; }, "min"); // seeds change key =/
                }
            }
        };
        for (var _i = 0, neighbors_1 = neighbors; _i < neighbors_1.length; _i++) {
            var q = neighbors_1[_i];
            _loop_1(q);
        }
    };
    /**
     * Expands the {@link cluster} with points in {@link seeds}.
     * @private
     * @param {Heap} seeds
     * @param {Array} cluster
     */
    OPTICS.prototype._expand_cluster = function (seeds, cluster) {
        var ordered_list = this._ordered_list;
        while (!seeds.empty) {
            var q = seeds.pop().element;
            q.neighbors = this._get_neighbors(q);
            q.processed = true;
            cluster.push(q.index);
            ordered_list.push(q);
            if (this._core_distance(q) != undefined) {
                this._update(q, seeds);
                this._expand_cluster(seeds, cluster);
            }
        }
    };
    /**
     * Returns an array of clusters.
     * @returns {Array<Array>} Array of clusters with the indices of the rows in given {@link matrix}.
     */
    OPTICS.prototype.get_clusters = function () {
        var clusters = [];
        var outliers = [];
        var min_points = this._min_points;
        for (var _i = 0, _a = this._clusters; _i < _a.length; _i++) {
            var cluster = _a[_i];
            if (cluster.length < min_points) {
                outliers.push.apply(outliers, cluster);
            }
            else {
                clusters.push(cluster);
            }
        }
        clusters.push(outliers);
        return clusters;
    };
    /**
     * @returns {Array} Returns an array, where the ith entry defines the cluster affirmation of the ith point of {@link matrix}. (-1 stands for outlier)
     */
    OPTICS.prototype.get_cluster_affirmation = function () {
        var N = this._matrix.shape[0];
        var result = new Array(N).fill();
        var clusters = this.get_clusters();
        for (var i = 0, n = clusters.length; i < n; ++i) {
            var cluster = clusters[i];
            for (var _i = 0, cluster_1 = cluster; _i < cluster_1.length; _i++) {
                var index = cluster_1[_i];
                result[index] = (i < n - 1) ? i : -1;
            }
        }
        return result;
    };
    return OPTICS;
}());

/**
 * @class
 * @alias LSP
 * @extends DR
 */
var LSP = /** @class */ (function (_super) {
    tslib.__extends(LSP, _super);
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
    function LSP(X, parameters) {
        var _a, _c;
        var _this = _super.call(this, X, { neighbors: undefined, control_points: undefined, d: 2, metric: euclidean, seed: 1212 }, parameters) || this;
        _this.parameter("neighbors", Math.min((_a = parameters.neighbors) !== null && _a !== void 0 ? _a : Math.max(Math.floor(_this._N / 10), 2), _this._N - 1));
        _this.parameter("control_points", Math.min((_c = parameters.control_points) !== null && _c !== void 0 ? _c : Math.ceil(Math.sqrt(_this._N)), _this._N - 1));
        _this._is_initialized = false;
        return _this;
    }
    /**
     *
     * @param {DR} DR - method used for position control points.
     * @param {Object} DR_parameters - Object containing parameters for the DR method which projects the control points
     * @returns {LSP}
     */
    LSP.prototype.init = function (DR, DR_parameters, KNN) {
        if (DR === void 0) { DR = MDS; }
        if (DR_parameters === void 0) { DR_parameters = {}; }
        if (KNN === void 0) { KNN = BallTree; }
        if (this._is_initialized)
            return this;
        var X = this.X;
        var N = this._N;
        var K = this.parameter("neighbors");
        var d = this.parameter("d");
        var seed = this.parameter("seed");
        var metric = this.parameter("metric");
        DR_parameters = Object.assign({ d: d, metric: metric, seed: seed }, DR_parameters);
        var nc = this.parameter("control_points");
        var control_points = new KMedoids(X, nc, null, metric).get_clusters().medoids;
        var C = new Matrix(nc, N, "zeros");
        control_points.forEach(function (c_i, i) {
            C.set_entry(i, c_i, 1);
        });
        var Y_C = new DR(Matrix.from(control_points.map(function (c_i) { return X.row(c_i); })), DR_parameters).transform();
        var XA = X.to2dArray;
        var knn = new KNN(XA, metric);
        var L = new Matrix(N, N, "I");
        var alpha = -1 / K;
        XA.forEach(function (x_i, i) {
            for (var _i = 0, _a = knn.search(x_i, K).iterate(); _i < _a.length; _i++) {
                var j = _a[_i].index;
                if (i === j)
                    continue;
                L.set_entry(i, j, alpha);
            }
        });
        var A = L.concat(C, "vertical");
        var z = new Matrix(N, d, "zeros");
        var b = z.concat(Y_C, "vertical");
        this._A = A;
        this._b = b;
        this._is_initialized = true;
        return this;
    };
    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    LSP.prototype.transform = function () {
        this.check_init();
        var A = this._A;
        var AT = A.T;
        var b = this._b;
        var ATA = AT.dot(A);
        var ATb = AT.dot(b);
        this.Y = Matrix.solve_CG(ATA, ATb, this._randomizer);
        return this.projection;
    };
    return LSP;
}(DR));

/**
 * @class
 * @alias TopoMap
 * @memberof module:dimensionality_reduction
 * @extends DR
 */
var TopoMap = /** @class */ (function (_super) {
    tslib.__extends(TopoMap, _super);
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
    function TopoMap(X, parameters) {
        var _a;
        var _this = _super.call(this, X, { metric: euclidean, seed: 1212 }, parameters) || this;
        _a = _this.X.shape, _this._N = _a[0], _this._D = _a[1];
        _this._distance_matrix = new Matrix(_this._N, _this._N, 0);
        return _this;
    }
    /**
     * @private
     */
    TopoMap.prototype.__lazy_distance_matrix = function (i, j, metric) {
        var D = this._distance_matrix;
        var X = this.X;
        var D_ij = D.entry(i, j);
        if (D_ij === 0) {
            var dist = metric(X.row(i), X.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
            return dist;
        }
        return D_ij;
    };
    /**
     * Computes the minimum spanning tree, using a given metric
     * @private
     * @param {Function} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    TopoMap.prototype._make_minimum_spanning_tree = function (metric) {
        if (metric === void 0) { metric = euclidean; }
        var N = this._N;
        var X = tslib.__spreadArray([], this.X, true);
        var disjoint_set = new DisjointSet(X);
        var F = [];
        var E = [];
        for (var i = 0; i < N; ++i) {
            for (var j = i + 1; j < N; ++j) {
                E.push([i, j, this.__lazy_distance_matrix(i, j, metric)]);
            }
        }
        E = E.sort(function (a, b) { return a[2] - b[2]; });
        for (var _i = 0, E_1 = E; _i < E_1.length; _i++) {
            var _a = E_1[_i], u = _a[0], v = _a[1], w = _a[2];
            var set_u = disjoint_set.find(X[u]);
            var set_v = disjoint_set.find(X[v]);
            if (set_u !== set_v) {
                F.push([u, v, w]);
                disjoint_set.union(set_u, set_v);
            }
        }
        return F.sort(function (a, b) { return a[2] - b[2]; });
    };
    /**
     * initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    TopoMap.prototype.init = function () {
        var metric = this._parameters.metric;
        this.Y = new Matrix(this._N, 2, 0);
        this._Emst = this._make_minimum_spanning_tree(metric);
        this._is_initialized = true;
        return this;
    };
    /**
     * Returns true if Point C is left of line AB.
     * @private
     * @param {Array} PointA - Point A of line AB
     * @param {Array} PointB - Point B of line AB
     * @param {Array} PointC - Point C
     * @returns {Boolean}
     */
    TopoMap.prototype.__hull_cross = function (_a, _b, _c) {
        var ax = _a[0], ay = _a[1];
        var bx = _b[0], by = _b[1];
        var sx = _c[0], sy = _c[1];
        return (bx - ax) * (sy - ay) - (by - ay) * (sx - ax) <= 0;
    };
    /**
     * Computes the convex hull of the set of Points S
     * @private
     * @param {Array} S - Set of Points.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     * @returns {Array} convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     */
    TopoMap.prototype.__hull = function (S) {
        var points = S.sort(function (_a, _b) {
            var x1 = _a[0], y1 = _a[1];
            var x2 = _b[0], y2 = _b[1];
            return y1 - y2 || x1 - x2;
        });
        var N = points.length;
        if (N <= 2)
            return points;
        var lower = [];
        for (var i = 0; i < N; ++i) {
            while (lower.length >= 2 && this.__hull_cross(lower[lower.length - 2], lower[lower.length - 1], points[i])) {
                lower.pop();
            }
            lower.push(points[i]);
        }
        var upper = [];
        for (var i = N - 1; i >= 0; --i) {
            while (upper.length >= 2 && this.__hull_cross(upper[upper.length - 2], upper[upper.length - 1], points[i])) {
                upper.pop();
            }
            upper.push(points[i]);
        }
        upper.pop();
        lower.pop();
        return lower.concat(upper);
    };
    /**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     * @private
     * @param {Array} PointA
     * @param {Array} PointB
     * @return {Object} Object containing the sinus- and cosinus-values for a rotation.
     */
    TopoMap.prototype.__findAngle = function (_a, _b) {
        var p1x = _a[0], p1y = _a[1];
        var p2x = _b[0], p2y = _b[1];
        var n = euclidean([p1x, p1y], [p2x, p2y]);
        if (n === 0)
            return {
                sin: 0,
                cos: 1
            };
        var vec = [(p2x - p1x) / n, (p2y - p1y) / n];
        var cos = vec[0];
        var sin = Math.sqrt(1 - cos * cos);
        sin = vec[1] >= 0 ? -sin : sin;
        return {
            sin: sin,
            cos: cos
        };
    };
    /**
     * @private
     * @param {Array} hull
     * @param {Array} p
     * @param {Bool} topEdge
     */
    TopoMap.prototype.__align_hull = function (hull, p, topEdge) {
        var v = -1;
        var d2;
        for (var i = 0; i < hull.length; ++i) {
            var d = euclidean(hull[i], p);
            if (v === -1) {
                d2 = d;
                v = i;
            }
            else {
                if (d2 > d) {
                    d2 = d;
                    v = i;
                }
            }
        }
        var v1;
        var v2;
        if (topEdge) {
            v1 = hull[v];
            v2 = hull[(v + 1) % hull.length];
        }
        else {
            if (v == 0)
                v = hull.length - 1;
            v1 = hull[v];
            v2 = hull[(v - 1) % hull.length];
        }
        var transformation = {
            tx: -hull[v][0],
            ty: -hull[v][1]
        };
        if (hull.length >= 2) {
            var _a = this.__findAngle(v1, v2), sin = _a.sin, cos = _a.cos;
            transformation.sin = sin;
            transformation.cos = cos;
        }
        else {
            transformation.sin = 0;
            transformation.cos = 1;
        }
        return transformation;
    };
    /**
     * @private
     * @param {Array} Point - The point which should get transformed.
     * @param {Object} Transformation - contains the values for translation and rotation.
     */
    TopoMap.prototype.__transform = function (_a, _b) {
        var px = _a[0], py = _a[1];
        var tx = _b.tx, ty = _b.ty, sin = _b.sin, cos = _b.cos;
        var x = px + tx;
        var y = py + ty;
        var xx = x * cos - y * sin;
        var yy = x * sin + y * cos;
        return [xx, yy];
    };
    /**
     * Calls {@link __transform} for each point in Set C
     * @private
     * @param {Array} C - Set of points.
     * @param {Object} t - Transform object.
     * @param {Number} yOffset - value to offset set C.
     */
    TopoMap.prototype.__transform_component = function (C, t, yOffset) {
        var N = C.length;
        for (var i = 0; i < N; ++i) {
            var c = C[i];
            var _a = this.__transform(c, t), cx = _a[0], cy = _a[1];
            c[0] = cx;
            c[1] = cy + yOffset;
        }
    };
    /**
     * @private
     * @param {Array} u - point u
     * @param {Array} v - point v
     * @param {Number} w - edge weight w
     */
    TopoMap.prototype.__align_components = function (u, v, w) {
        var points_u = tslib.__spreadArray([], u.__disjoint_set.children, true);
        var points_v = tslib.__spreadArray([], v.__disjoint_set.children, true);
        var hull_u = this.__hull(points_u);
        var hull_v = this.__hull(points_v);
        var t_u = this.__align_hull(hull_u, u, false);
        var t_v = this.__align_hull(hull_v, v, true);
        this.__transform_component(points_u, t_u, 0);
        this.__transform_component(points_v, t_v, w);
    };
    /**
     * Transforms the inputdata {@link X} to dimensionality 2.
     */
    TopoMap.prototype.transform = function () {
        if (!this._is_initialized)
            this.init();
        var Emst = this._Emst;
        var Y = this.Y.to2dArray;
        var components = new DisjointSet(Y.map(function (y, i) {
            y.i = i;
            return y;
        }));
        for (var _i = 0, Emst_1 = Emst; _i < Emst_1.length; _i++) {
            var _a = Emst_1[_i], u = _a[0], v = _a[1], w = _a[2];
            var component_u = components.find(Y[u]);
            var component_v = components.find(Y[v]);
            if (component_u === component_v)
                continue;
            this.__align_components(component_u, component_v, w);
            components.union(component_u, component_v);
        }
        return this.projection;
    };
    TopoMap.prototype.generator = function () {
        var Emst, Y, components, _i, Emst_2, _a, u, v, w, component_u, component_v;
        return tslib.__generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    if (!this._is_initialized)
                        this.init();
                    Emst = this._Emst;
                    Y = this.Y.to2dArray;
                    components = new DisjointSet(Y.map(function (y, i) {
                        y.i = i;
                        return y;
                    }));
                    _i = 0, Emst_2 = Emst;
                    _b.label = 1;
                case 1:
                    if (!(_i < Emst_2.length)) return [3 /*break*/, 4];
                    _a = Emst_2[_i], u = _a[0], v = _a[1], w = _a[2];
                    component_u = components.find(Y[u]);
                    component_v = components.find(Y[v]);
                    if (component_u === component_v)
                        return [3 /*break*/, 3];
                    this.__align_components(component_u, component_v, w);
                    components.union(component_u, component_v);
                    return [4 /*yield*/, this.projection];
                case 2:
                    _b.sent();
                    _b.label = 3;
                case 3:
                    _i++;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/, this.projection];
            }
        });
    };
    return TopoMap;
}(DR));

/**
 * @class
 * @alias SAMMON
 * @extends DR
 */
var SAMMON = /** @class */ (function (_super) {
    tslib.__extends(SAMMON, _super);
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
    function SAMMON(X, parameters) {
        var _this = _super.call(this, X, { magic: 0.1, d: 2, metric: euclidean, seed: 1212, init_DR: "random", init_parameters: {} }, parameters) || this;
        return _this;
    }
    /**
     * initializes the projection.
     * @private
     */
    SAMMON.prototype.init = function () {
        var N = this.X.shape[0];
        var _a = this._parameters, d = _a.d, metric = _a.metric, init_DR = _a.init_DR, DR_parameters = _a.init_parameters;
        if (init_DR === "random") {
            var randomizer_1 = this._randomizer;
            this.Y = new Matrix(N, d, function () { return randomizer_1.random; });
        }
        else if (["PCA", "MDS"].includes(init_DR)) {
            this.Y = Matrix.from(init_DR == "PCA" ? PCA.transform(this.X, DR_parameters) : MDS.transform(this.X, DR_parameters));
        }
        else {
            throw new Error('init_DR needs to be either "random" or a DR method!');
        }
        this.distance_matrix = metric == "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        return this;
    };
    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Matrix|Array} - The projection of {@link X}.
     */
    SAMMON.prototype.transform = function (max_iter) {
        if (max_iter === void 0) { max_iter = 200; }
        if (!this._is_initialized)
            this.init();
        for (var j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    };
    /**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Generator} - A generator yielding the intermediate steps of the projection of {@link X}.
     */
    SAMMON.prototype.generator = function (max_iter) {
        var j;
        if (max_iter === void 0) { max_iter = 200; }
        return tslib.__generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    if (!this._is_initialized)
                        this.init();
                    j = 0;
                    _a.label = 1;
                case 1:
                    if (!(j < max_iter)) return [3 /*break*/, 4];
                    this._step();
                    return [4 /*yield*/, this.projection];
                case 2:
                    _a.sent();
                    _a.label = 3;
                case 3:
                    ++j;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/, this.projection];
            }
        });
    };
    SAMMON.prototype._step = function () {
        var MAGIC = this.parameter("magic");
        var D = this.distance_matrix;
        var N = this.X.shape[0];
        var _a = this._parameters, d = _a.d, metric = _a.metric;
        var Y = this.Y;
        var G = new Matrix(N, d, 0);
        var sum = new Float64Array(d);
        for (var i = 0; i < N; ++i) {
            var e1 = new Float64Array(d);
            var e2 = new Float64Array(d);
            var Yi = Y.row(i);
            for (var j = 0; j < N; ++j) {
                if (i === j)
                    continue;
                var Yj = Y.row(j);
                var delta = new Float64Array(d);
                for (var k = 0; k < d; ++k) {
                    delta[k] = Yi[k] - Yj[k];
                }
                var dY = metric(Yi, Yj);
                var dX = D.entry(i, j);
                var dq = dX - dY;
                var dr = Math.max(dX * dY, 1e-2);
                for (var k = 0; k < d; ++k) {
                    e1[k] += (delta[k] * dq) / dr;
                    e2[k] += (dq - (Math.pow(delta[k], 2) * (1 + dq / dY)) / dY) / dr;
                }
            }
            for (var k = 0; k < d; ++k) {
                var val = Y.entry(i, k) + ((MAGIC * e1[k]) / Math.abs(e2[k]) || 0);
                G.set_entry(i, k, val);
                sum[k] += val;
            }
        }
        for (var k = 0; k < d; ++k) {
            sum[k] /= N;
        }
        for (var i = 0; i < N; ++i) {
            for (var k = 0; k < d; ++k) {
                Y.set_entry(i, k, G.entry(i, k) - sum[k]);
            }
        }
        return Y;
    };
    return SAMMON;
}(DR));

var version="0.5.1";

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
