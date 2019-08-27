// https://renecutura.eu v0.0.2 Copyright 2019 Rene Cutura
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = global || self, factory(global.druid = global.druid || {}));
}(this, function (exports) { 'use strict';

/*import { neumair_sum } from "../numerical/index";

export default function(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        let x = a[i];
        let y = b[i]
        s[i] = ((x - y) * (x - y))
    }
    return Math.sqrt(neumair_sum(s))
}*/

function euclidean(a, b) {
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
function kahan_sum(summands) {
    let n = summands.length;
    let sum = 0;
    let compensation = 0;
    let y, t;

    for (let i = 0; i < n; ++i) {
        y = summands[i] - compensation;
        t = sum + y;
        compensation = (t - sum) - y;
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
function neumair_sum(summands) {
    let n = summands.length;
    let sum = 0;
    let compensation = 0;

    for (let i = 0; i < n; ++i) {
        let summand = summands[i];
        let t = sum + summand;
        if (Math.abs(sum) >= Math.abs(summand)) {
            compensation += (sum - t) + summand;
        } else {
            compensation += (summand - t) + sum;
        }
        sum = t;
    }
    return sum + compensation;
}

/**
 * @module numerical
 */

function euclidean_squared(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        let x = a[i];
        let y = b[i];
        s[i] = ((x - y) * (x - y));
    }
    return neumair_sum(s);
}

/**
 * Computes the Cosine distance between vector {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Array} a 
 * @param {Array} b 
 * @returns {float64} The Cosine distance between vector {@link a} and {@link b}.  
 */
function cosine(a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += (a[i] * b[i]);
        sum_a += (a[i] * a[i]);
        sum_b += (b[i] * b[i]);
    }
    return sum / ((Math.sqrt(sum_a) * Math.sqrt(sum_b)));
}

function manhattan(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum
}

/**
 * Computes the Chebyshev distance between vector {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Array} a 
 * @param {Array} b 
 * @returns {float64} the Chebyshev distance between vector {@link a} and {@link b}.  
 */
function chebyshev(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res)
}

/**
 * @module metrics
 */

function k_nearest_neighbors(A, k, distance_matrix = null, metric = euclidean) {
    let n = A.length;
    let D = distance_matrix || dmatrix(A, metric);
    for (let i = 0; i < n; ++i) {
        D[i] = D[i].map((d,j) => {
            return {
                i: i, j: j, distance: D[i][j]
            }
        }).sort((a, b) => a.distance - b.distance)
        .slice(1, k + 1);
    }
    return D
}

function dmatrix(A, metric = euclidean) {
    let distance = metric;
    if (distance === undefined) return undefined;
    let n = A.length;
    let D = new Array(n);
    for (let i = 0; i < n; ++i) {
        D[i] = new Float64Array(n);
    }
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = distance(A[i], A[j]);
        }
    }
    return D;
}

function linspace(start, end, number = null) {
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
    return result
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
 * Computes the QR Decomposition of the Matrix {@link A}.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A 
 */
function qr(A) {
    let [rows, cols] = A.shape;
    let Q = new Matrix(rows, cols, "identity");
    let R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        let v = A.col(j);
        for (let i = 0; i < j; ++i) {
            let q = Q.col(i);
            let q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            R.set_entry(i,j, q_dot_v);
            v = v.map((v_, k) => v_ - q_dot_v * q[k]);
        }
        let v_norm = norm(v, euclidean);
        for (let k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j,j, v_norm);
    }
    /*let R = matrix.clone()

    for (let j = 0; j < cols; ++j) {
        let z = new Matrix()
        z.shape = [rows - j, 1, (i,_) => R.entry(i+j, j)]
        let norm_z = norm(z, euclidean);
        let rho = z.entry(0, 0) < 0 ? 1 : -1;
        let u1 = z.entry(0, 0) - rho * norm_z;
        let u = z.divide(u1);
        u.set_entry(0, 0, 1);
        let beta = -rho * u1 / norm_z;

        let row_j = rows - j
        let u_outer_u = u.dot(u.transpose())
        let R_j__ = new Matrix()
        R_j__.shape = [row_j, row_j, (r, c) => R.entry(row_j + r, row_j + c)]
        let R_sub = u_outer_u.dot(R_j__).mult(-beta)
        for (let r = 0; r < row_j; ++r) {
            for (let c = 0; c < row_j; ++c) {
                R.set_entry(row_j + r, row_j + c, R.entry(row_j + r, row_j + c) + R_sub.entry(r, c))
            }
        }
        let Q__j_ = new Matrix();
        Q__j_.shape = [rows, row_j, (r, c) => Q.entry(r, row_j + c)]
        let Q_sub = Q__j_.dot(u_outer_u).mult(-beta)
        for (let r = 0; r < rows; ++r) {
            for ( let c = 0; c < row_j; ++c) {
                Q.set_entry(r, row_j + c, Q.entry(r, row_j + c) + Q_sub.entry(r, c))
            }
        }

        console.log(j, Q, R)
    }*/

    return { R: R, Q: Q };
}

function qr_householder(A) {
    let [rows, cols] = A.shape;
    let Q = new Matrix(rows, rows, "identity");
    let R = A.clone();

    for (let j = 0; j < cols; ++j) {
        let x = R.get_block(j, j).col(0);
        let x_norm = norm(x);
        let x0 = x[0];
        let rho = -Math.sign(x0);
        let u1 = x0 - rho * x_norm;
        let u = x.map(e => e / u1);
        u[0] = 1;
        let beta = -rho * u1 / x_norm;

        //let u_outer_u = new Matrix(cols - j, cols - j, (i, j) => u[i] * u[j]);
        u = Matrix.from(u, "col");
        let R_j_0 = R.get_block(j, 0);
        R.set_block(j, 0, R_j_0.sub(u.dot(u.T.dot(R_j_0)).mult(beta)));
        let Q_0_j = Q.get_block(0, j);
        Q.set_block(0, j, Q_0_j.sub(Q_0_j.dot(u).dot(u.T))); 
    }

    // repair R // numerical unstability?
    for (let i = 1; i < rows; ++i) {
        for (let j = 0; j < i; ++j) {
            R.set_entry(i, j, 0);
        }
    }
    return { R: R, Q: Q };
}

function simultaneous_poweriteration(A, k = 2, max_iterations = 100, seed = 19870307) {
    let randomizer = new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0];
    let { Q: Q, R: R } = qr(new Matrix(n, k, () => randomizer.random));
    while(max_iterations--) {
        let oldR = R.clone();
        let Z = A.dot(Q);
        let QR = qr(Z);
        [ Q, R ] = [ QR.Q, QR.R ]; 
        if (neumair_sum(R.sub(oldR).diag) / n < 1e-12) {
            max_iterations = 0;
        }        
    }

    let eigenvalues = R.diag;
    let eigenvectors = Q.transpose().to2dArray;//.map((d,i) => d.map(dd => dd * eigenvalues[i]))
    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors
    };
}

// crout algorithm
// https://en.wikipedia.org/wiki/Crout_matrix_decomposition
function lu(A) {
    let rows = A.shape[0];
    let L = new Matrix(rows, rows, "zeros");
    let U = new Matrix(rows, rows, "identity");
    let sum;

    for (let j = 0; j < rows; ++j) {
        for (let i = j; i < rows; ++i) {
            sum = 0;
            for (let k = 0; k < j; ++k) {
                sum += L.entry(i, k) * U.entry(k, j);
            }
            /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(i, k) * U.entry(k, j)))
            console.log("lsum1", sum)
            sum = []
            for (let k = 0; k < j; ++k) {
                sum.push(L.entry(i, k) * U.entry(k, j))
            }
            sum = neumair_sum(sum)
            console.log("lsum2", sum)*/
            L.set_entry(i, j, A.entry(i, j) - sum);
        }
        for (let i = j; i < rows; ++i) {
            if (L.entry(j, j) === 0) {
                return undefined;
            }
            sum = 0;
            for (let k = 0; k < j; ++k) {
                sum += L.entry(j, k) * U.entry(k, i);
            }
            /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(j, k) * U.entry(k, i)))
            console.log("usum1", sum)
            sum = []
            for (let k = 0; k < j; ++k) {
                sum.push(L.entry(j, k) * U.entry(k, i))
            }
            sum = neumair_sum("usum2", sum)
            console.log(sum)*/
            U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
        }
    }

    return { L: L, U: U };
}

// doolittle algorithm
/*export default function(M) {
    let [rows, cols] = M.shape;
    let P = new Matrix();
    P.shape = [rows + 1, 1, i => i];
    let A = M.clone();
    let L = new Matrix();
    let U = new Matrix();
    let I = new Matrix();
    I.shape = [rows, rows, (i, j) => i === j ? 1 : 0];

    for (let i = 0; i < rows; ++i) {
        let max_A = 0;
        let i_max = i;
        let abs_A;
        for (let k = i; k < rows; ++k) {
            abs_A = Math.abs(A.entry(k, i))
            if (abs_A > max_A) {
                [ max_A, i_max ] = [ abs_A, k ];
            }

            if (max_A < 1e-12) return undefined;

            if (i_max !== i) {
                let p = P.row(i);
                P.set_row(i, P.row(i_max));
                P.set_row(i_max, p);

                let A_row = A.row(i);
                A.set_row(i, A.row(i_max));
                A.set_row(i_max, A_row);

                P.set_entry(rows + 1, 0, P.entry(rows + 1, 0) + 1)
            }
        }

        for (let j = i + 1; j < rows; ++j) {
            A.set_entry(j, i,  A.entry(j, i) / A.entry(i, i));
            for (let k = i + 1; k < rows; ++k) {
                A.set_entry(j, k, A.entry(j, k) - A.entry(j, i) * A.entry(i, k));
            }
        }
    }

    L.shape = [rows, rows, (i, j) => {
        if ( i > j ) return A.entry(i, j);
        if ( i === j ) return A.entry(i, j) + 1;
        return 0
    }]

    U = A.sub(L)

    return {L: L, U: U, P: P};   
}*/

/**
 * Computes the eigenvector of {@link X} with an accelerated stochastic power iteration algorithm.
 * @memberof module:linear_algebra 
 * @alias svrg
 * @see {@link https://arxiv.org/abs/1707.02670}
 * @param {Matrix} data - the data matrix
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} epoch - number of epochs
 * @param {number} m - epoch length
 * @param {number} s - mini-batch size
 * @param {number} seed - seed for the random number generator
 */
function svrg(data, x, beta, epoch=20, m=10, s=1, seed) {
    let [n, d] = data.shape;
    const randomizer = new Randomizer(seed);
    x = new Matrix(d, 1, () => randomizer.random);
    x = x.divide(norm(x));
    let x0 = x.clone();
    let A = data.T.dot(data).divide(n);
    let x_tilde = x.clone();
    
    for (let t = 0; t < epoch; ++t) {
        const gx = A.dot(x_tilde);
        for (let i = 0; i < m; ++i) {
            const ang = x.T.dot(x_tilde).entry(0, 0);
            const sample = Matrix.from(Randomizer.choice(data, s));
            const sampleT_dot_sample = sample.T.dot(sample);
            const x_tmp = x.clone();
            const XTXx = sampleT_dot_sample
                    .dot(x.divide(s));
            const XTXx_tilde = sampleT_dot_sample
                    .dot(x_tilde.mult(ang / s));
            x = XTXx.sub(XTXx_tilde)
                    .add(gx.mult(ang).sub(x0.mult(beta)));
            x0 = x_tmp;
            const x_norm = norm(x);
            x = x.divide(x_norm);
            x0 = x0.divide(x_norm);        
        }  
        x_tilde = x.clone();
    }
    return x;

}

/**
 * 
 * @param {Matrix} data - the data matrix
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} max_iter - maximum number of iterations
 * @param {number} seed - seed for the random number generator
 */
function poweriteration_m(data, x0, beta, max_iter=20, seed) {
    let randomizer = new Randomizer(seed);
    let [ n, d ] = data.shape;
    let A = data.T.dot(data).divide(n);
    if (x0 === null) x0 = new Matrix(d, 1, () => randomizer.random);
    x0 = x0.divide(norm(x0));
    let x = x0.clone();
    for (let i = 0; i < max_iter; ++i) {
        let x_tmp = x.clone();
        x = A.dot(x).sub(x0.mult(beta));
        x0 = x_tmp;
        let z = norm(x);
        x = x.divide(z);
        x0 = x0.divide(z);
    }
    return x;
}

/**
 * @typedef {Eigenpair} Eigenpair
 * @property {Array} Eigenvalues - Array of Eigenvalues
 * @property {Array[]} Eigenvectors - Array of Eigenvectors 
 */


/**
 * Computes the {@link n} biggest Eigenpair of the Matrix {@link data}.
 * @memberof module:linear_algebra
 * @alias poweriteration_n
 * @param {Matrix} data - the data matrix
 * @param {int} n - Number of Eigenvalues / Eigenvectors
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} max_iter - maximum number of iterations
 * @param {number} seed - seed for the random number generator
 * @returns {Eigenpair} The {@link n} Eigenpairs.
 */
function poweriteration_n(data, n, x0, beta, max_iter=100, seed) {
    const randomizer = new Randomizer(seed);
    const N = data.shape[0];
    //let b = new Matrix(N, n, () => randomizer.random);
    let b = [];

    if (x0 == null) {
        x0 = new Array(n);//new Matrix(N, n, () => randomizer.random)
        
        for (let i = 0; i < n; ++i) {
            x0[i] = new Float64Array(N);
            b[i] = new Float64Array(N);
            for (let j = 0; j < N; ++j) {
                const value = randomizer.random;
                x0[i][j] = value;
                b[i][j] = value;
            }
            let x0_i_norm = norm(x0[i]);
            x0[i] = x0[i].map(x => x / x0_i_norm);
        }
        //x0 = Matrix.from(x0).T;
        //b = Matrix.from(b).T;
    }
    //x0 = x0.divide(norm(x0));
    for (let k = 0; k < n; ++k) {
        let bk = b[k];
        for (let s = 0; s < max_iter; ++s) {
            // Orthogonalize vector
            for (let l = 0; l < k; ++l) {
                const row = b[l];
                const d = neumair_sum((new Float64Array(N)).map((_, i) => bk[i] * row[i]));
                for (let i = 0; i < N; ++i) {
                    bk[i] = bk[i] - (d * row[i]);
                }
            }
            let tmp = data.dot(bk);
            const tmp_norm = norm(tmp);
            x0[k] = tmp.map(t => t / tmp_norm);
            if (neumair_sum((new Float64Array(N)).map((_, i) => tmp[i] * bk[i])) > (1 - 1e-12)) {
                break;
            }
            [bk, tmp] = [tmp, bk];
        }
    }

    return {
        "eigenvalues": b,
        "eigenvectors": x0
    }
}

/**
 * @module linear_algebra
 */

/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
class Matrix{
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
    constructor(rows=null, cols=null, value=null) {
        this._rows = rows;
        this._cols = cols;
        this._data = null;
        if (rows && cols) {
            if (!value) {
                this._data = new Float64Array(rows * cols);
                return this;
            }
            if (typeof(value) === "function") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value(row, col);
                    }
                }
                return this;
            }
            if (typeof(value) === "string") {
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
                    value = (i, j) => (i === j ? 1 : 0) - (1 / rows);
                    for (let row = 0; row < rows; ++row) {
                        for (let col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = value(row, col);
                        }
                    }
                    return this;
                }
            }
            if (typeof(value) === "number") {
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
     * @param {string} [type = "row"] - If {@link A} is a Array or Float64Array, then type defines if it is a row- or a column vector. 
     * @returns {Matrix}
     * 
     * @example
     * let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
     */
    static from(A, type="row") {
        if (A instanceof Matrix) {
            return A.clone();
        } else if (Array.isArray(A) || A instanceof Float64Array) {
            let m = A.length;
            if (m === 0) throw "Array is empty";
            // 1d
            if (!Array.isArray(A[0]) && !(A[0] instanceof Float64Array)) {
                if (type === "row") {  
                    return new Matrix(1, m, (_, j) => A[j]);
                } else if (type === "col") {
                    return new Matrix(m, 1, (i) => A[i]);
                } else {
                    throw "1d array has NaN entries"
                }
            // 2d
            } else if (Array.isArray(A[0]) || A[0] instanceof Float64Array) {
                let n = A[0].length;
                for (let row = 0; row < m; ++row) {
                    if (A[row].length !== n) throw "various array lengths";
                }
                return new Matrix(m, n, (i, j) => A[i][j])
            }
        } else if (typeof(A) === "number") {
            return new Matrix(1, 1, A);
        } else {
            throw "error"
        }
    }

    /**
     * Returns the {@link row}th row from the Matrix.
     * @param {int} row 
     * @returns {Array}
     */
    row(row) {
        let result_row = new Array(this._cols);
        for (let col = 0; col < this._cols; ++col) {
            result_row[col] = this._data[row * this._cols + col];
        }
        return result_row;
    }

    /**
     * Sets the entries of {@link row}th row from the Matrix to the entries from {@link values}.
     * @param {int} row 
     * @param {Array} values 
     * @returns {Matrix}
     */
    set_row(row, values) {
        let cols = this._cols;
        if (Array.isArray(values) && values.length === cols) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this._data[offset + col] = values[col];
            }
        } else if (values instanceof Matrix && values.shape[1] === cols && values.shape[0] === 1) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this._data[offset + col] = values._data[col];
            }
        }
        return this;
    }

    /**
     * Returns the {@link col}th column from the Matrix.
     * @param {int} col 
     * @returns {Array}
     */
    col(col) {
        let result_col = new Array(this._rows);
        for (let row = 0; row < this._rows; ++row) {
            result_col[row] = this._data[row * this._cols + col];
        }
        return result_col;
    }

    /**
     * Returns the {@link col}th entry from the {@link row}th row of the Matrix.
     * @param {int} row 
     * @param {int} col 
     * @returns {float64}
     */
    entry(row, col) {
        return this._data[row * this._cols + col];
    }

    /**
     * Sets the {@link col}th entry from the {@link row}th row of the Matrix to the given {@link value}.
     * @param {int} row 
     * @param {int} col 
     * @param {float64} value
     * @returns {Matrix}
     */
    set_entry(row, col, value) {
        this._data[row * this._cols + col] = value;
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
        let B = new Matrix(rows, 2 * cols, (i,j) => {
            if (j >= cols) {
                return (i === (j - cols)) ? 1 : 0;
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
                let val = Math.abs(B.entry(i,k));
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
                    let B_i_j = B.entry(i,j);
                    let B_row_j = B.entry(row, j);
                    B_i_j = B_i_j - B_row_j * f;
                    B.set_entry(i, j, B_i_j);
                }
            }
        }

        return new Matrix(rows, cols, (i,j) => B.entry(i, j + cols));
    }

    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    dot(B) {
        if (B instanceof Matrix) {
            let A = this;
            if (A.shape[1] !== B.shape[0]) return undefined;
            let I = A.shape[1];
            let C = new Matrix(A.shape[0], B.shape[1], (row, col) => {
                let A_i = A.row(row);
                let B_i = B.col(col);
                for (let i = 0; i < I; ++i) {
                    A_i[i] *= B_i[i];
                }
                return neumair_sum(A_i);
            });
            return C;
        } else if (Array.isArray(B) || (B instanceof Float64Array)) {
            let rows = this._rows;
            if (B.length !== rows) return undefined;
            let C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map(e => e * B[row]));
            }
            return C;
        } else {
            return undefined;
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
        C.shape = [l, l, (i, j) => {
            if (i <= j) {
                return A._data[i] * B._data[j];
            } else {
                return C.entry(j, i);
            }
        }];
        return C;
    }

    /**
     * Transforms A to bidiagonal form. (With Householder transformations)
     * @param {Matrix} A 
     */
    static bidiagonal(A) {
        A = A.clone();
        let [ rows, cols ] = A.shape;
        for (let k = 0; k < cols; ++k) {
            let x = A.col(k).slice(k);
            let x_norm = norm(x);
            let u = A.col(k).slice(k);
            u[0] = x[0] + (x[0] < 0 ? -1 : 1) * x_norm;
            let u_norm = norm(u);
            u = u.map(u_i => u_i /= u_norm);
            u = Matrix.from(u, "col");
            let A_ = A.get_block(k, k);
            A_ = A_.sub(u.mult(2).outer(u).dot(A_));
            A.set_block(k, k, A_); /** @todo repair numerical unstability? */
            // repair numerical unstability?
            for (let r = k + 1; r < rows; ++r) {
                A.set_entry(r, k, 0);
            }
            if (k <= cols - 2) {
                let y = A.row(k).slice(k + 1);
                let y_norm = norm(y);
                let v = A.row(k).slice(k + 1);
                v[0] = y[0] + (y[0] < 0 ? -1 : 1) * y_norm;
                let v_norm = norm(v);
                v = v.map(v_i => v_i /= v_norm);
                v = Matrix.from(v, "col");
                let _A = A.get_block(k, k + 1);
                _A = _A.sub(_A.dot(v.mult(2).outer(v))); 
                A.set_block(k, k + 1, _A);
                // repair numerical unstability?
                for (let c = k + 2; c < cols; ++c) {
                    A.set_entry(k, c, 0);
                }
            }
        }
        return A;
    }

    /**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row 
     * @param {int} offset_col 
     * @param {Matrix} B 
     * @returns {Matrix}
     */
    set_block(offset_row, offset_col, B) {
        let [ rows, cols ] = B.shape;
        for (let row = 0; row < rows; ++row) {
            if (row > this._rows) continue;
            for (let col = 0; col < cols; ++col) {
                if (col > this._cols) continue;
                this.set_entry(row + offset_row, col + offset_col, B.entry(row, col));
            }
        }
        return this;
    }

    /**
     * Extracts the entries of A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row 
     * @param {int} offset_col 
     * @returns {Matrix}
     */
    get_block(offset_rows, offset_cols) {
        let [ rows, cols ] = this.shape;
        return new Matrix(rows - offset_rows, cols - offset_cols, (i, j) => this.entry(i + offset_rows, j + offset_cols));
    }

    /**
     * Applies a function to each entry of the matrix.
     * @param {function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    _apply_array(f, v) {
        const data = this._data;
        const [ rows, cols ] = this.shape;
        for (let row = 0; row < rows; ++row) {
            const o = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = o + col;
                const d = data[i];
                data[i] = f(d, v(row, col));
            }
        }
        return this; 
    }

    _apply_rowwise_array(values, f) {
        return this._apply_array(f, (i, j) => values[j]);
    }

    _apply_colwise_array(values, f) {
        const data = this._data;
        const [ rows, cols ] = this.shape;
        for (let row = 0; row < rows; ++row) {
            const o = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = o + col;
                const d = data[i];
                data[i] = f(d, values[row]);
            }
        }
        return this; 
    }

    /*_apply_pointwise_number(value, f) {

    }*/

    _apply(value, f) {
        let data = this._data;
        if (value instanceof Matrix) {
            let [ value_rows, value_cols ] = value.shape;
            let [ rows, cols ] = this.shape;
            if (value_rows === 1) {
                if (cols !== value_cols) return undefined;
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(0, col));
                    }
                }
            } else if (value_cols === 1) {
                if (rows !== value_rows) return undefined;
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
            } else return undefined;
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
                return undefined;
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
        B._data = this._data.slice(0);
        return B;
    }

    mult(value) {
        return this.clone()._apply(value, (a,b) => a * b)
    }

    divide(value) {
        return this.clone()._apply(value, (a,b) => a / b)
    }

    add(value) {
        return this.clone()._apply(value, (a,b) => a + b)
    }

    sub(value) {
        return this.clone()._apply(value, (a,b) => a - b)
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
    set shape([ rows, cols, value = () => 0 ]) {
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
     * Returns the Matrix as a two-dimensional Array.
     * @returns {Array}
     */
    get to2dArray() {
        const rows = this._rows;
        const cols = this._cols;
        let result = new Array(rows);
        for (let row = 0; row < rows; ++row) {
            let result_col = new Array(cols);
            for (let col = 0; col < cols; ++col) {
                result_col[col] = this.entry(row, col);
            }
            result[row] = result_col;
        }
        return result;
    }

    /**
     * Returns the diagonal of the Matrix.
     * @returns {Array}
     */
    get diag() {
        const rows = this._rows;
        const cols = this._cols;
        const min_row_col = Math.min(rows, cols);
        let result = new Array(min_row_col);
        for (let i = 0; i < min_row_col; ++i) {
            result[i] = this.entry(i,i);
        }
        return result;
    }

    /**
     * Returns the mean of all entries of the Matrix.
     * @returns {float64}
     */
    get mean() {
        const data = this._data;
        const n = this._rows * this._cols;
        let sum = 0;
        for (let i = 0; i < n; ++i) {
            sum += data[i];
        }
        return sum / n;
    }

    /**
     * Returns the mean of each row of the matrix.
     * @returns {Array}
     */
    get meanRows() {
        const data = this._data;
        const rows = this._rows;
        const cols = this._cols;
        let result = [];
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
     * @returns {Array}
     */
    get meanCols() {
        const data = this._data;
        const rows = this._rows;
        const cols = this._cols;
        let result = [];
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
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A, b) {
        let rows = A.shape[0];
        let { L: L, U: U } = Matrix.LU(A);//lu(A);
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
        let rows = A.shape[0];
        let L = new Matrix(rows, rows, "zeros");
        let U = new Matrix(rows, rows, "identity");
        let sum;

        for (let j = 0; j < rows; ++j) {
            for (let i = j; i < rows; ++i) {
                sum = 0;
                for (let k = 0; k < j; ++k) {
                    sum += L.entry(i, k) * U.entry(k, j);
                }
                /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(i, k) * U.entry(k, j)))
                console.log("lsum1", sum)
                sum = []
                for (let k = 0; k < j; ++k) {
                    sum.push(L.entry(i, k) * U.entry(k, j))
                }
                sum = neumair_sum(sum)
                console.log("lsum2", sum)*/
                L.set_entry(i, j, A.entry(i, j) - sum);
            }
            for (let i = j; i < rows; ++i) {
                if (L.entry(j, j) === 0) {
                    return undefined;
                }
                sum = 0;
                for (let k = 0; k < j; ++k) {
                    sum += L.entry(j, k) * U.entry(k, i);
                }
                /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(j, k) * U.entry(k, i)))
                console.log("usum1", sum)
                sum = []
                for (let k = 0; k < j; ++k) {
                    sum.push(L.entry(j, k) * U.entry(k, i))
                }
                sum = neumair_sum("usum2", sum)
                console.log(sum)*/
                U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
            }
        }

        return { L: L, U: U };
    }

    /**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M 
     * @param {int} [k=2] 
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */
    static SVD(M, k=2) {
        let MtM = M.transpose().dot(M);
        let MMt = M.dot(M.transpose());
        let { eigenvectors: V, eigenvalues: Sigma } = simultaneous_poweriteration(MtM, k);
        let { eigenvectors: U } = simultaneous_poweriteration(MMt, k);

        return { U: U, Sigma: Sigma.map(sigma => Math.sqrt(sigma)), V: V };
        /*const [ rows, cols ] = matrix.shape;
        // 1.
        let U = new Matrix();
        U.shape = [rows, cols, (i, j) => i === j ? 1 : 0];
        let Sigma = matrix.clone();
        let Vt = new Matrix();
        Vt.shape = [cols, cols, (i, j) => i === j ? 1 : 0];
        let I = new Matrix();
        I.shape = [cols, cols, (i, j) => i === j ? 1 : 0];

        function householder_vector(x) {
            let dot_1on = norm(x.col(0).slice(1));
            let v = x.clone();
            v.set_entry(0, 0, 1);
            let beta = 0;
            if (dot_1on > 1e-63) {
                let x_0 = x.entry(0, 0);
                let x_norm = Math.sqrt(x_0 ** 2 + dot_1on);
                let v_0;
                if (x_0 <= 0) {
                    v_0 = x_0 - x_norm;
                } else {
                    v_0 = -dot_1on / (x_0 + x_norm);
                }
                v.set_entry(0, 0, v_0)
                beta = 2 * v_0 ** 2 / (dot_1on + v_0 ** 2)
                v.divide(v_0)
            }
            return {"v": v, "beta": beta}
        }

        function householder_matrix(size, col, v, beta) {
            let vvt = v.dot(v.transpose())
            let V = new Matrix()
            V.shape = [size, size, (r, c) => {
                if (r < col || c < col) {
                    return r === c ? 1 : 0;
                }
                return (r === c ? 1 : 0) - (beta * vvt.entry(r - col, c - col));
            }]

            console.log(V, vvt)
            return V;
        }

        for (let col = 0; col < cols; ++col) {
            let x = Matrix.from(Sigma.col(col), "col")
            let { v: u, beta: beta } = householder_vector(x);
            let Q = householder_matrix(cols, col, u, beta);
            console.log("u", u, Q)
            U = U.dot(Q)
            Sigma = Q.dot(Sigma);

            if (col < cols - 2) {
                let x = Matrix.from(Sigma.row(col).slice(1), "col")
                let { v: v, beta: beta } = householder_vector(x);
                Q = householder_matrix(cols, col + 1, v, beta);
                console.log("v", v, Q)
                Vt = Q.dot(Vt)
                Sigma = Sigma.dot(Q);
            }

            /*let u = Sigma.col(col);
            let u_norm = norm(u);
            let sign_u = u[col] >= 0 ? 1 : -1;
            u[col] = sign_u * (Math.abs(u[col]) + u_norm);
            u_norm = norm(u);
            //u = u.map(u_i => u_i / u_norm);
            for (let c = 0; c < cols; ++c) {
                u[c] = u[c] / u_norm;
            }
            u = Matrix.from(u, "col")
            let H = I.sub(u.dot(u.transpose()).mult(2))
            console.log(u.to2dArray, H.to2dArray)
            Sigma = H.dot(Sigma);
            U = U.dot(H.transpose())

            if (col < cols - 1) {

            }*/

            /*let u = Sigma.col(col).slice(col);
            let u_norm = norm(u);
            let sign_u = u[0] >= 0 ? 1 : -1;
            u[0] = sign_u * (Math.abs(u[0]) + u_norm);//+= sign_u * u_norm;
            u_norm = norm(u)//Math.sqrt(2 * u_norm * (u_norm + u[0]))//norm(u);
            u = u.map(u_i => u_i / u_norm);
            u = Matrix.from(u, "col")
            let A_k_m = new Matrix()
            A_k_m.shape = [rows - col, cols - col, (i, j) => Sigma.entry(i + col, j + col)];
            let U_k = u.dot(u.transpose().dot(A_k_m)).mult(2)
            for (let r = 0; r < rows - col; ++r) {
                for (let c = 0; c < cols - col; ++c) {
                    let val = Sigma.entry(col + r, col + c) - U_k.entry(r, c);
                    val = Math.abs(val) < 1e-15 ? 0 : val;
                    Sigma.set_entry(col + r, col + c, val)
                }
            }
            console.log(Sigma)



            /*if (col <= cols - 2) {
                let col_1 = col + 1
                let v = Sigma.row(col).slice(col_1);
                let v_norm = norm(v);
                let sign_v = v[0] >= 0 ? 1 : -1;
                v[0] += sign_v * (Math.abs(v[0]) + v_norm);//sign_v * v_norm;
                v_norm = norm(v);
                v = v.map(v_i => v_i / v_norm);
                v = Matrix.from(v, "col")
                let A_k1_n = new Matrix()
                A_k1_n.shape = [rows - col, cols - col_1, (i, j) => Sigma.entry(col + i, col_1 + j)];
                let V_k = A_k1_n.dot(v.dot(v.transpose())).mult(2)
                for (let r = 0; r < rows - col; ++r) {
                    for (let c = 0; c < cols - col_1; ++c) {
                        let val = Sigma.entry(col + r, col_1 + c) - V_k.entry(r, c)
                        console.log("V", [col + r, col_1 + c], 2 * Sigma.entry(col + r, col_1 + c) == V_k.entry(r, c))
                        val = Math.abs(val) < 1e-15 ? 0 : val;
                        Sigma.set_entry(col + r, col_1 + c, val)
                    }
                }
            }

            */
            /*console.log(Sigma)
        }*/
    }

    
}

// https://github.com/jasondavies/science.js/blob/master/src/stats/hcluster.js

class Hierarchical_Clustering {
    constructor(matrix, linkage="single", metric=euclidean) {
        this._id = 0;
        this._matrix = matrix;
        this._metric = metric;
        this._linkage = linkage;
        this.init();
        this.root = this.do();
        return this;
    }

    get_clusters(value, type="distance") {
        let clusters = [];
        let accessor;
        switch (type) {
            case "distance":
                accessor = d => d.dist;
                break;
            case "depth":
                accessor = d => d.depth;
                break;
            default:
                throw "invalid type";
        }
        this._traverse(this.root, accessor, value, clusters);
        return clusters
    }

    _traverse(node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        } else {
            this._traverse(node.left, f, value, result);
            this._traverse(node.right, f, value, result);
        }
    }

    init() {
        const metric = this._metric;
        const A = this._matrix;
        const n = this._n = A.shape[0];
        const d_min = this._d_min = new Float64Array(n);
        const distance_matrix = this._distance_matrix = new Array(n);
        for (let i = 0; i < n; ++i) {
            d_min[i] = 0;
            distance_matrix[i] = new Float64Array(n);
            for (let j = 0; j < n; ++j) {
                distance_matrix[i][j] = i === j ? Infinity : metric(A.row(i), A.row(j));
                if (distance_matrix[i][d_min[i]] > distance_matrix[i][j]) {
                    d_min[i] = j;
                }
            }
        }

        const clusters = this._clusters = new Array(n);
        const c_size = this._c_size = new Uint16Array(n);
        for (let i = 0; i < n; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        return this;
    }

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
                if (D[i][d_min[i]] < D[c1][d_min[c1]]) {
                    c1 = i;
                }
            }
            let c2 = d_min[c1];

            let c1_cluster = clusters[c1][0];
            let c2_cluster = clusters[c2][0];

            let new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D[c1][c2]);
            clusters[c1].unshift(new_cluster);
            c_size[c1] += c_size[c2];

            for (let j = 0; j < n; ++j) {
                switch(linkage) {
                    case "single":
                        if (D[c1][j] > D[c2][j]) {
                            D[j][c1] = D[c1][j] = D[c2][j];
                        }
                        break;
                    case "complete":
                        if (D[c1][j] < D[c2][j]) {
                            D[j][c1] = D[c1][j] = D[c2][j];
                        }
                        break;
                    case "average":
                        D[j][c1] = D[c1][j] = (c_size[c1] * D[c1][j] + c_size[c2] * D[c2][j]) / (c_size[c1] + c_size[j]);
                        break;
                }
            }

            D[c1][c1] = Infinity;

            for (let i = 0; i < n; ++i) {
                D[i][c2] = D[c2][i] = Infinity;
            }
            for (let j = 0; j < n; ++j) {
                if (d_min[j] === c2) {
                    d_min[j] = c1;
                }
                if (D[c1][j] < D[c1][d_min[c1]]) {
                    d_min[c1] = j;
                }
            }

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
        this.size = size != null ? size : left.size + right.size;
        this.depth = depth != null ? depth : 1 + Math.max(left.depth, right.depth);
        this.centroid = centroid != null ? centroid : this._calculate_centroid(left, right);
        /*if (centroid == null) {
            const n = left.centroid.length;
            const new_centroid = this.centroid = new Float64Array(n);
            for (let i = 0; i < n; ++i) {
                new_centroid[i] = (left.size * left.centroid[i] + right.size * right.centroid[i]) / this.size;
            }
        } else {
            this.centroid = centroid;
        }
        console.log(left, right, this)*/
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
        if (this.isLeaf) return [this.index];
        const left = this.left;
        const right = this.right;
        return (left.isLeaf ? [left.index] : left.leaves())
            .concat(right.isLeaf ? [right.index] : right.leaves())
    }
}

// https://github.com/jdfekete/reorder.js/blob/master/reorder.v1.js
class Reorder{
    constructor(A) {
        this._A = A;
        this._optimal_leaf_order = null;
    }

    get available() {
        return [
            "optimal_leaf_order",
            "spectral_order",
        ]
    }

    reorder(type="optimal_leaf_order", metric=euclidean) {
        let result = null;
        switch (type) {
            case "optimal_leaf_order":
                this._optimal_leaf_order = new Optimal_Leaf_Order(this._A, metric);
                result = this._optimal_leaf_order.ordering;
                break;

            case "spectral_order":
                this._spectral_order = new Spectral_Order(this._A, metric);
                result = this._spectral_order.ordering;
                break;
            case "barycenter_order":
                break;
        }
        return result;
    }
}

//class Barycenter_Order{    
//}

class Spectral_Order{
    constructor(A, metric=euclidean) {
        this._A = A;
        this._metric = metric;
        this._N = A.shape[0];
        const fiedler_vector = this._fiedler_vector(A);
        this._ordering = linspace(0, this._N - 1).sort((a, b) => fiedler_vector[a] - fiedler_vector[b]);
        return this;
    }

    get ordering() {
        return this._ordering;
    }

    _fiedler_vector(B) {
        const g = this._gershgorin_bound(B);
        const N = B.shape[0];
        const B_hat = new Matrix(N, N, (i, j) => i === j ? g - B.entry(i, j) : -B.entry(i, j));
        const eig = simultaneous_poweriteration(B_hat, 2);
        return eig.eigenvectors[1]
    }

    _gershgorin_bound(B) {
        let max = 0;
        let N = B.shape[0];
        for (let i = 0; i < N; ++i) {
            let t = B.entry(i, i);
            for (let j = 0; j < N; ++j) {
                if (i !== j) {
                    t += Math.abs(B.entry(i, j));
                }
            }
            max = max > t ? max : t;
        }
        return max;
    }
}

class Optimal_Leaf_Order{
    constructor(A, metric=euclidean) {
        this._A = A;
        const N = A.shape[0];
        const hclust = this._hclust = new Hierarchical_Clustering(A, "complete", metric);
        const distance_matrix = this._distance_matrix = new Array(N);
        for (let i = 0; i < N; ++i) {
            distance_matrix[i] = new Float64Array(N);
            for (let j = 0; j < N; ++j) {
                distance_matrix[i][j] = i === j ? Infinity : metric(A.row(i), A.row(j));
            }
        }
        this._order_map = new Map();
        let min = Infinity;
        this._optimal_order = null;
        let left = hclust.root.left.leaves();
        let right = hclust.root.right.leaves();

        for (let i = 0, n = left.length; i < n; ++i) {
            for (let j = 0, m = right.length; j < m; ++j) {
                let so = this.order(hclust.root, left[i], right[j]);
                if (so[0] < min) {
                    min = so[0];
                    this._optimal_order = so[1];
                }
            }
        }
        return this;
    }

    get ordering() {
        return this._optimal_order;
    }

    order(v, i, j) {
        const order_map = this._order_map;
        const key = `k${v.id}-${i}-${j}`; // ugly key
        /*if (key in order_map) 
            return order_map[key];
        return (order_map[key] = this._order(v, i, j))*/
        /*if (order_map.has(v)) {
            const v_map = order_map.get(v)
            if (v_map.has(`${i},${j}`)) {
                return v_map.get(`${i},${j}`)
            } else {
                let value = this._order(v, i, j);
                v_map.set(`${i},${j}`, value);
                return value;
            }
        } else {
            let value = this._order(v, i, j);
            const v_map = new Map();
            v_map.set(`${i},${j}`, value);
            order_map.set(v, v_map);
            return value;
        }*/
        if (order_map.has(key)) {
            return order_map.get(key);
        } else {
            let value = this._order(v, i, j);
            order_map.set(key, value);
            return value;
        }
    }

    _order(v, i, j) {
        if (v.isLeaf) {
            return [0, [v.index]];
        }
        const D = this._distance_matrix;
        let l = v.left;
        let r = v.right;
        let L = l ? l.leaves() : [];
        let R = r ? r.leaves() : [];
        let w;
        let x;
        if (L.indexOf(i) !== -1 && R.indexOf(j) !== -1) {
            w = l; 
            x = r;
        } else if (R.indexOf(i) !== -1 && L.indexOf(j) !== -1) {
            w = r;
            x = l;
        } else {
            throw "Node is not common ancestor of i and j";
        }

        let Wl = w.left ? w.left.leaves() : [];
        let Wr = w.right ? w.right.leaves() : [];
        let Ks = Wr.indexOf(i) != -1 ? Wl : Wr;
        if (Ks.length === 0) { 
            Ks = [i];
        }

        let Xl = x.left ? x.left.leaves() : [];
        let Xr = x.right ? x.right.leaves() : [];
        let Ls = Xr.indexOf(j) != -1 ? Xl : Xr;
        if (Ls.length === 0) {
            Ls = [j];
        }

        let min = Infinity;
        let optimal_order = [];
        for (let k = 0, Ks_length = Ks.length; k < Ks_length; ++k) {
            let w_min = this.order(w, i, Ks[k]);
            for (let m = 0, Ls_length = Ls.length; m < Ls_length; ++m) {
                let x_min = this.order(x, Ls[m], j);
                let dist = w_min[0] + D[Ks[k]][Ls[m]] + x_min[0];
                if (dist < min) {
                    min = dist;
                    optimal_order = w_min[1].concat(x_min[1]);
                }
            }
        }

        return [min, optimal_order];
    }
}

/**
 * @module matrix
 */

class Randomizer {
    // https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
    /**
     * Mersenne Twister
     * @param {*} _seed 
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
            mt[mti] = (((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253) + mti;
            mt[mti] >>>= 0;
        }
    }

    get seed() {
        return this._seed;
    }

    get random() {
        return this.random_int * (1.0 / 4294967296.0)
    }

    get random_int() {
        let y, mag01 = new Array(0x0, this._MATRIX_A);
        if (this._mti >= this._N) {
            let kk;

            if (this._mti == this._N + 1) {
                this.seed = 5489;
            }

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

        y = this._mt[this._mti += 1];
        y ^= (y >>> 11);
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= (y >>> 18);

        return y >>> 0;
    }

    choice(A, n) {
        if (A instanceof Matrix) {
            let [rows, cols] = A.shape;
            if (n > rows) throw "n bigger than A!";
            let sample = new Array(n);
            let index_list = linspace(0, rows - 1);
            for (let i = 0, l = index_list.length; i < n; ++i, --l) {
                let random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map(d => A.row(d));
        } else if (Array.isArray(A) || A instanceof Float64Array) {
            let rows = A.length;
            if (n > rows) {
                throw "n bigger than A!";
            }
            let sample = new Array(n);
            let index_list = linspace(0, rows - 1);
            for (let i = 0, l = index_list.length; i < n; ++i, --l) {
                let random_index = this.random_int % l;
                sample[i] = index_list.splice(random_index, 1)[0];
            }
            return sample.map(d => A[d]);
        }
    }

    static choice(A, n, seed=19870307) {
        let [rows, cols] = A.shape;
        if (n > rows) throw "n bigger than A!"
        let rand = new Randomizer(seed);
        let sample = new Array(n);
        let index_list = linspace(0, rows - 1);
        /*let index_list = new Array(rows);
        for (let i = 0; i < rows; ++i) {
            index_list[i] = i;
        }*/
        //let result = new Matrix(n, cols);
        for (let i = 0, l = index_list.length; i < n; ++i, --l) {
            let random_index = rand.random_int % l;
            sample[i] = index_list.splice(random_index, 1)[0];
            //random_index = index_list.splice(random_index, 1)[0];
            //result.set_row(i, A.row(random_index))
        }
        //return result;
        //return new Matrix(n, cols, (row, col) => A.entry(sample[row], col))
        return sample.map(d => A.row(d));
    }
}

/*export class Heap {
    constructor(arr = null, accessor = (d) => d, comparator = "min") {
        this.root = null;
        this.accessor = accessor;

        if (comparator == "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator == "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d));
        }
    }

    push(element) {
        const value = this.accessor(element);
        const newNode = new Node(element, value);
        if (!this.root || this._comparator(value, this.root.value)) {
            newNode.next = this.root;
            this.root = newNode;
        } else {
            let pointer = this.root;
            while (pointer.next && !this._comparator(value, pointer.next.value)) {
                pointer = pointer.next;
            }
            newNode.next = pointer.next;
            pointer.next = newNode;
        }
        return this;
    }

    pop() {
        if (!this.root) {
            return null;
        }
        const root = this.root;
        this.root = this.root.next;
        return root;
    }

    get first() {
        return this.root;
    }

    * iterate() {
        let pointer = this.root;
        while (pointer) {
            yield pointer.element;
            pointer = pointer.next;
        }
    }

    toArray() {
        let res = [];
        let pointer = this.root;
        while (pointer) {
            res.push(pointer.element)
            pointer = pointer.next;
        }
        return res;
    }

    get length() {
        let len = 0;
        let pointer = this.root;
        while (pointer) {
            len += 1;
            pointer = pointer.next;
        }
        return len;
    }

    get empty() {
        return this.root === null;
    }
}

class Node {
    constructor(element, value) {
        this.element = element;
        this.value = value;
        this.next = null;
    }
}*/

class Heap {
    constructor(arr = null, accessor = (d) => d, comparator = "min") {
        //this.root = null;
        this._accessor = accessor;
        this._container = [];

        if (comparator == "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator == "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        //console.log(arr)
        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d));
        }

    }

    /* constructor(arr = null, accessor = (d) => d, comparator = "min") {
        this._accessor = accessor;

        if (comparator === "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator === "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        console.log(arr)
        this._container = []
        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d)) 
        }
    } */

    _get_left_child_index(index) {
        return (index * 2) + 1;
    }

    _get_right_child_index(index) {
        return (index * 2) + 2;
    }

    _get_parent_index(index) {
        return Math.floor((index - 1) / 2);
    }

    _has_parent(index) {
        return this._get_parent_index(index) >= 0;
    }

    _has_left_child(index) {
        return this._get_left_child_index(index) < this._container.length;
    }

    _has_right_child(index) {
        return this._get_right_child_index(index) < this._container.length;
    }

    _left_child(index) {
        return this._container[this._get_left_child_index(index)];
    }

    _right_child(index) {
        return this._container[this._get_right_child_index(index)];
    }

    _parent(index) {
        return this._container[this._get_parent_index(index)];
    }

    _swap(index_a, index_b) {
        //[this._container[index_b], this._container[index_a]] = [this._container[index_a], this._container[index_b]];
        let tmp = this._container[index_b];
        this._container[index_b] = this._container[index_a];
        this._container[index_a] = tmp;
    }

    push(element) {
        const value = this._accessor(element);
        const node = new Node(element, value);
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    _heapify_up(start_index) {
        let index = start_index || this._container.length - 1;
        while (this._has_parent(index) && !this._comparator(this._parent(index).value, this._container[index].value)) {
            this._swap(index, this._get_parent_index(index));
            index = this._get_parent_index(index);
        }
    }

    pop() {
        if (this._container.length === 0) {
            return null;
        }
        if (this._container.length === 1) {
            return this._container.pop().element;
        }
        
        const item = this._container[0];

        this._container[0] = this._container.pop();
        this._heapify_down();

        return item;
    }

    _heapify_down(start_index=0) {
        let index = start_index;
        let next_index = null;

        while (this._has_left_child(index)) {
            if (this._has_right_child(index) && this._comparator(this._right_child(index).value, this._left_child(index).value)) {
                next_index = this._get_right_child_index(index);
            } else {
                next_index = this._get_left_child_index(index);
            }

            if (this._comparator(this._comparator(this._container[index].value, this._container[next_index].value))) {
                break;
            }

            this._swap(index, next_index);
            index = next_index;
        }
    }

    // peek
    get first() {
        return this._container.length > 0 ? this._container[0] : null;
    }

    * iterate() {
        for (let i = 0, n = this._container.length; i < n; ++i) {
            yield this._container[i].element;
        }
    }

    toArray() {
        const comparator = this._comparator;
        const accessor = this._accessor;
        let container = this._container;//.map(d => d.element);
        return container.sort((a, b) => comparator(accessor(a.element), accessor(b.element)) ? -1 : 1)
        //return this._container.sort((a, b) => comparator(a.value, b.value))//.map(d => d.element);
    }

    get length() {
        return this._container.length;
    }

    get empty() {
        return this.length === 0;
    }
}

class Node {
    constructor(element, value) {
        this.element = element;
        this.value = value;
    }
}

class HNSW {
    /**
     * 
     * @param {*} metric metric to use: (a, b) => distance
     * @param {*} heuristic use heuristics or naive selection
     * @param {*} m max number of connections
     * @param {*} ef size of candidate list
     * @param {*} m0 max number of connections for ground layer 
     * @see {@link https://arxiv.org/abs/1603.09320}
     */
    constructor(metric = euclidean, heuristic = true, m = 5, ef = 200, m0 = null, mL = null) {
        this._metric = metric;
        this._select = heuristic ? this._select_heuristic : this._select_simple;
        this._m = m;
        this._ef = ef;
        this._m0 = m0 || 2 * m;
        this._graph = [];
        this._ep = null;
        this._L = null;
        this._mL = mL === null ? 1 / Math.log2(m) : mL;
        this.search = this.search;
    }

    addOne(element) {
        this.add([element]);
    }

    add(...elements) {
        const m = this._m;
        const ef = this._ef;
        const m0 = this._m0;
        //const metric = this._metric;
        const mL = this._mL;
        let graph = this._graph;
        for (const element of elements) {
            let ep = this._ep ? Array.from(this._ep): null;
            let W = [];
            let L = this._L;
            let l = Math.floor(-Math.log(Math.random() * mL));
            let min_L_l = Math.min(L, l);
            if (L) {
                for (let l_c = graph.length - 1; l_c > min_L_l; --l_c) {
                    ep = this._search_layer(element, ep, 1, l_c);
                }
                for (let l_c = min_L_l; l_c >= 0; --l_c) {
                    let layer_c = graph[l_c];
                    layer_c.points.push(element);
                    W = this._search_layer(element, ep, ef, l_c);
                    let neighbors = this._select(element, W, m, l_c);
                    neighbors.forEach(p => {
                        if (p !== element) {
                            //let distance = metric(p, element);
                            layer_c.edges.push({
                                idx1: p, 
                                idx2: element, 
                                ///distance: distance
                            });
                            layer_c.edges.push({
                                idx1: element, 
                                idx2: p, 
                                //distance: distance
                            });
                        }
                    });
                    let max = (l_c === 0 ? m0 : m);
                    for (let e of neighbors) {
                        let e_conn = layer_c.edges
                            .filter(edge => edge.idx1 === e)
                            .map(edge => edge.idx2);
                        if (e_conn.length > max) {
                            let neighborhood = this._select(e, e_conn, max, l_c);
                            layer_c.edges = layer_c.edges
                                .filter(edge => edge.idx1 !== e);
                            neighborhood.forEach(neighbor => {
                                if (e !== neighbor) {
                                    //let distance = metric(e, neighbor);
                                    layer_c.edges.push({
                                        idx1: e, 
                                        idx2: neighbor, 
                                        //distance: distance
                                    });
                                }
                            });
                        }
                    }
                    ep = W;
                }
            }
            if (graph.length < l || l > L) {
                for (let i = l, n = graph.length; i >= n; --i) {
                    let new_layer = {
                        l_c: i, 
                        points: [element], 
                        edges: new Array()
                    };
                    graph.push(new_layer);
                    if (i === l) {
                        this._ep = [element];
                        this._L = l;
                    }
                }
                graph = graph.sort((a, b) => a.l_c - b.l_c);
            }
        }
        return this;
    }

    _select_heuristic(q, candidates, M, l_c, extend_candidates = true, keep_pruned_connections = true) {
        if (l_c > this._graph.length - 1) return candidates
        const metric = this._metric;
        const layer = this._graph[l_c];
        let R = [];
        let W_set = new Set(candidates);
        if (extend_candidates) {
            for (let c of candidates) {
                for (let {idx2: c_adj} of layer.edges.filter(edge => edge.idx1 === c)) {
                    W_set.add(c_adj);
                }
            }
        }
        let W = new Heap(Array.from(W_set), d => metric(d, q), "min");
        let W_d = new Heap(null, d => metric(d, q), "min");
        while (W.first && R.length < M) {
            let e = W.pop();
            let random_r = Math.floor(Math.random() * R.length);
            if (R.length === 0 || e.value < metric(R[random_r], q)) {
                R.push(e.element);
            } else {
                W_d.push(e.element);
            }
        }
        if (keep_pruned_connections) {
            while (W_d.first && R.length < M) {
                R.push(W_d.pop().element);
            }
        }
        return R
    }

    _select_simple(q, C, M) {
        const metric = this._metric;
        let res = C.sort((a,b) => metric(a, q) - metric(b, q)).slice(0,M);
        return res
    }

    _search_layer(q, ep, ef, l_c) {
        const metric = this._metric;
        const layer = this._graph.find(l => l.l_c === l_c);
        let v = new Set(ep);
        let C = new Heap(ep, d => metric(d, q), "min");
        let W = new Heap(ep, d => metric(d, q), "max");
        while (C.length > 0) {
            let c = C.pop();
            let f = W.first;
            if (c.value > f.value) {
                break;
            }
            for (let {idx2: e} of layer.edges.filter(e => e.idx1 === c.element)) {
                if (!v.has(e)) {
                    v.add(e);
                    f = W.first.element;
                    if (metric(e, q) < metric(f, q) || W.length < ef) {
                        C.push(e);
                        W.push(e);
                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }
        return W.toArray().reverse().slice(0, ef);
    }

    search(q, K, ef = null) {
        ef = ef || 1;
        let ep = this._ep;
        let L = this._L;
        for (let l_c = L; l_c > 0; --l_c) {
            ep = this._search_layer(q, ep, ef, l_c);
        }
        ep = this._search_layer(q, ep, K, 0);
        return ep;
    }

    * search_iter(q, K, ef = null) {
        ef = ef || 1;
        let ep = this._ep ? Array.from(this._ep): null;
        let L = this._L;
        yield {l_c: L, ep: [q]};
        for (let l_c = L; l_c > 0; --l_c) {
            yield {l_c: l_c, ep: ep};
            ep = this._search_layer(q, ep, ef, l_c);
            yield {l_c: l_c, ep: ep};
        }
        yield {l_c: 0, ep: ep};
        ep = this._search_layer(q, ep, K, 0);
        yield {l_c: 0, ep: ep};
    }
}

/**
 * @memberof module:knn
 */
class BallTree {
    /**
     * Generates a BallTree with given {@link elements}.
     * @param {Array} elements - Elements which should be added to the BallTree
     * @param {function} metric metric to use: (a, b) => distance
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
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

    add(elements) {
        elements = elements.map((element, index) => {
            return {index: index, element: element}
        });
        this._root = this._construct(elements);
        return this;
    }

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
            //if (this._root === null) this._root = B;
            return B;
        } /*else {
            return null;
        }*/
    }

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
        return c
    }

    search(t, k=5) {
        return this._search(t, k, new Heap(null, d => this._metric(d.element, t), "max"), this._root);
    }

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
 * @memberof module:knn
 */
class NNDescent{
    /**
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     * @see {@link https://github.com/lmcinnes/pynndescent/blob/master/pynndescent/pynndescent_.py}
     */
    constructor(elements=null, metric=euclidean, seed=19870307) {
        this._metric = metric;
        const rand = this._randomizer = new Randomizer(seed);
        if (elements) {
            const k = this._k = 5;
            const B = this._knn_list = new Heap(rand.choice(elements, k), d => d, "min");
            let c = true;
            while (c) {
                let R = this._reverse(B);
                c = false;
            }

        }
        

        return this;   
    }

    _reverse(B) {
        return B;
    }

    /**
     * 
     * @param {Array} elements 
     */
    add(elements) {
        return this;
    }

    search(x, k=5) {
        return new Heap(x, k)
    }
}

/**
 * @module knn
 */

class PCA{
    constructor(X, d=2) {
        this.X = X;
        this.d = d;
    }

    transform() {
        let X = this.X;
        let D = X.shape[1];
        let O = new Matrix(D, D, "center");
        let X_cent = X.dot(O);

        let C = X_cent.transpose().dot(X_cent);
        let { eigenvectors: V } = simultaneous_poweriteration(C, this.d);
        console.log(V);
        V = Matrix.from(V).transpose();
        this.Y = X.dot(V);
        return this.Y
    }

    get projection() {
        return this.Y
    }

    

}

class MDS{
    constructor(X, d=2, metric=euclidean) {
        this.X = X;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        const X = this.X;
        //let sum_reduce = (a,b) => a + b
        const rows = X.shape[0];
        const metric = this._metric;
        let ai_ = [];
        let a_j = [];
        for (let i = 0; i < rows; ++i) {
            ai_.push(0);
            a_j.push(0);
        }
        let a__ = 0;
        const A = new Matrix();
        A.shape = [rows, rows, (i,j) => {
            let val = 0;
            if (i < j) {
                val = metric(X.row(i), X.row(j));
            } else if (i > j) {
                val = A.entry(j,i);
            }
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        }];
        this._d_X = A;
        ai_ = ai_.map(v => v / rows);
        a_j = a_j.map(v => v / rows);
        a__ /= (rows ** 2);
        const B = new Matrix(rows, rows, (i, j) => (A.entry(i, j) - ai_[i] - a_j[j] + a__));
        //B.shape = [rows, rows, (i,j) => (A.entry(i,j) - (A.row(i).reduce(sum_reduce) / rows) - (A.col(j).reduce(sum_reduce) / rows) + a__)]
                
        const { eigenvectors: V } = simultaneous_poweriteration(B, this.d);
        this.Y = Matrix.from(V).transpose();
        
        return this.Y
    }

    get projection() {
        return this.Y
    }

    get stress() {
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

class ISOMAP{
    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this.k = neighbors || Math.floor(this.X.shape[0] / 10);
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let rows = X.shape[0];
        // make knn extern and parameter for constructor or transform?
        let D = new Matrix();
        D.shape = [rows, rows, (i,j) => i <= j ? this._metric(X.row(i), X.row(j)) : D.entry(j,i)];
        let kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            let row = D.row(i).map((d,i) => { 
                return {
                    "index": i,
                    "distance": d
                }
            });
            let H = new Heap(row, d => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, this.k + 1));
        }
        
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        let G = new Matrix(rows, rows, (i,j) => {
            let other = kNearestNeighbors[i].find(n => n.index === j);
            return other ? other.distance : Infinity
        });

        for (let i = 0; i < rows; ++i) {
            for (let j = 0; j < rows; ++j) {
                for (let k = 0; k < rows; ++k) {
                    G.set_entry(i, j, Math.min(G.entry(i, j), G.entry(i, k) + G.entry(k, j)));
                }
            }
        }
        
        let ai_ = [];
        let a_j = [];
        for (let i = 0; i < rows; ++i) {
            ai_.push(0);
            a_j.push(0);
        }
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
        let B = new Matrix(rows, rows, (i,j) => (A.entry(i,j) - ai_[i] - a_j[j] + a__));
             
        // compute d eigenvectors
        let { eigenvectors: V } = simultaneous_poweriteration(B, this.d);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.Y
    }

    get projection() {
        return this.Y
    }

    

}

class FASTMAP{
    constructor(X, d=2, metric=euclidean) {
        this.X = X;
        this.d = d;
        this._metric = metric;
        this._col = -1;
        this.randomizer = new Randomizer(1212);
    }

    _choose_distant_objects(dist) {
        let X = this.X;
        let N = X.shape[0];
        let a_index = this.randomizer.random_int % N - 1;
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            let d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            let d_bi = dist(b_index, i);
            if (d_bi > max_dist) {
                max_dist = d_bi;
                a_index = i;
            }
        }
        return [a_index, b_index, max_dist];
    }

    transform() {
        let X = this.X;
        let [ rows, D ] = X.shape;
        let Y = new Matrix(rows, this.d);
        //let PA = [[], []];
        let dist = (a,b) => this._metric(X.row(a), X.row(b));
        let old_dist = dist;

        while(this._col < this.d - 1) {
            this._col += 1;
            let col = this._col;
            // choose pivot objects
            let [a_index, b_index, d_ab] = this._choose_distant_objects(dist);
            // record id of pivot objects
            //PA[0].push(a_index);
            //PA[1].push(b_index);
            if (d_ab === 0) {
                // because all inter-object distances are zeros
                for (let i = 0; i < rows; ++i) {
                    Y.set_entry(i, col, 0);
                }
            } else {
                // project the objects on the line (O_a, O_b)
                for (let i = 0; i < rows; ++i) {
                    let d_ai = dist(a_index, i);
                    let d_bi = dist(b_index, i);
                    let y_i = (d_ai ** 2 + d_ab ** 2 - d_bi ** 2) / (2 * d_ab);
                    Y.set_entry(i, col, y_i);
                }
                // consider the projections of the objects on a
                // hyperplane perpendicluar to the line (a, b);
                // the distance function D'() between two 
                // projections is given by Eq.4
                dist = (a,b) => Math.sqrt((old_dist(a,b) ** 2) - ((Y.entry(a, col) - Y.entry(b, col)) ** 2));
            }
        }
        // return embedding
        this.Y = Y;
        return this.Y;
    }

    get projection() {
        return this.Y
    }

    

}

class LDA{
    constructor(X, labels, d=2, metric=euclidean) {
        this.X = X;
        this._labels = labels;
        this.d = d;
        this._metric = metric;
    }

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

        let { eigenvectors: V } = simultaneous_poweriteration(S_w.inverse().dot(S_b), this.d);
        V = Matrix.from(V).transpose();
        this.Y = X.dot(V);

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}

class LLE{
    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [ rows, cols ] = X.shape;
        let k = this._k;
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        let O = new Matrix(k, 1, 1);
        let W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            let Z = new Matrix(k, cols, (i, j) => X.entry(nN[row][i].j, j) - X.entry(row, j));
            let C = Z.dot(Z.transpose());
            if ( k > cols ) {
                let C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < k; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }

            // reconstruct;
            let w = Matrix.solve(C, O);
            let w_sum = neumair_sum(w.col(0));
            w = w.divide(w_sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, nN[row][j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        let I = new Matrix(rows, rows, "identity");
        let IW = I.sub(W);
        let M = IW.transpose().dot(IW);
        let { eigenvectors: V } = simultaneous_poweriteration(M.transpose().inverse(), d + 1);
        
        this.Y = Matrix.from(V.slice(1, 1 + d)).transpose();

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}

class MLLE{
    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [ rows, cols ] = X.shape;
        let k = this._k;
        // 1.1 Determine a neighborset
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        let O = new Matrix(k, 1, 1);
        let W = new Matrix(rows, k);
        let Phi = new Matrix(rows, rows);

        let V = new Array(rows);
        let Lambda = new Array(rows);
        let P = new Array(rows);

        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map(n => n.j);
            let x_i = Matrix.from(X.row(row), "row");
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            X_i = X_i.sub(x_i);
            //X_i = X_i.dot(new Matrix(X_i._cols, X_i._cols, "center"))
            let C_i = X_i.dot(X_i.transpose()); // k by k

            let gamma = neumair_sum(C_i.diag) / 1000;
            for (let j = 0; j < k; ++j) {
                C_i.set_entry(j, j, C_i.entry(j, j) + gamma);
            }
            
            let { eigenvalues: Lambdas, eigenvectors: v } = simultaneous_poweriteration(C_i, k);
            V[row] = v; // k by k, rows are eigenvectors, big to small
            Lambda[row] = Lambdas; // 1 by k, cols are eigenvalues, big to small
            P.push(neumair_sum(Lambdas.slice(d + 1)) / neumair_sum(Lambdas.slice(0, d)));

            // reconstruct;
            let w = Matrix.solve(C_i, O); // k by 1
            let w_sum = neumair_sum(w.col(0));
            w = w.divide(w_sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, j, w.entry(j, 0));
            }
        }
        // find regularized weights // median
        let theta = P.sort((ρ_i, ρ_j) => ρ_i - ρ_j)[Math.ceil(rows / 2)];
        
        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map(n => n.j);
            let Lambdas = Lambda[row]; // 1 by k
            let s_i = Lambdas.map((Lambda, l) => {
                    return {
                        "l": l,
                        "ratio": neumair_sum(Lambdas.slice(k - l + 1)) / neumair_sum(Lambdas.slice(0, k - l)),
                        "Lambda": Lambda
                    }
                });
            //console.log(s_i)
            s_i = s_i
                .filter(s => s.ratio < theta && s.l <= k - d)
                .map(s => s.l).pop() || d;
            let V_i = V[row]; // k by k
            V_i = V_i.slice(k - s_i); // s_i by k
            let alpha_i = (1 / Math.sqrt(s_i)) * norm(V_i[0].map((_, j) => neumair_sum(V_i.map(r => r[j]))));
            V_i = Matrix.from(V_i); // s_i by k
            
            //https://github.com/scikit-learn/scikit-learn/blob/7b136e9/sklearn/manifold/locally_linear.py#L703

            let h = new Matrix(s_i, 1, alpha_i);
            let ones = new Matrix(k, 1, 1);
            h = h.sub(V_i.dot(ones));
            let h_norm = norm(h.col(0));
            h = h_norm < 1e-12 ? h.mult(0) : h.divide(h_norm);
            V_i = V_i.T;
            ones = new Matrix(s_i, 1, 1);
            let w_i = Matrix.from(W.row(row), "col");
            
            /*let H_i = new Matrix(s_i, s_i, "identity");
            H_i = H_i.sub(h.mult(2).outer(h));
            let W_i = V_i.sub(V_i.dot(h).dot(h.T).mult(2)).add(w_i.mult(1 - alpha_i))
            */
            let W_i = V_i.sub(V_i.dot(h).dot(h.T).mult(2)).add(w_i.mult(1 - alpha_i).dot(ones.T));
            
            W_i = W_i.dot(W_i.T);
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < s_i; ++j) {
                    Phi.set_entry(I_i[i], I_i[j], Phi.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }
        //let { eigenvectors: Y } = simultaneous_poweriteration(Phi.inverse(), d + 1);
        //this.Y = Matrix.from(Y.slice(1)).transpose()

        let { eigenvectors: Y } = simultaneous_poweriteration(Phi, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}

// https://epubs.siam.org/doi/abs/10.1137/S1064827502419154
class LTSA{

    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [ rows, D ] = X.shape;
        let k = this._k;
        // 1.1 determine k nearest neighbors
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        // center matrix
        let O = new Matrix(D, D, "center");
        let B = new Matrix(rows, rows, 0);
        
        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            let I_i = [row, ...nN[row].map(n => n.j)];
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            // center X_i
            X_i = X_i.dot(O);
            // correlation matrix
            let C = X_i.dot(X_i.transpose());
            let { eigenvectors: g } = simultaneous_poweriteration(C, d);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            let G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            let W_i = G_i_t.transpose().dot(G_i_t).add(1 / Math.sqrt(k + 1));
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < k + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }

        // 3. Aligning global coordinates
        let { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}

/*import { simultaneous_poweriteration} from "../linear_algebra/index";
import { k_nearest_neighbors } from "../matrix/index";
import { neumair_sum } from "../numerical/index";
import { norm } from "../matrix/index";
import { linspace } from "../matrix/index";*/

class TSNE{
    constructor(X, perplexity, epsilon, d=2, metric=euclidean, seed=1212) {
        this._X = X;
        this._d = d;
        [ this._N, this._D ] = X.shape;
        this._perplexity = perplexity;
        this._epsilon = epsilon;
        this._metric = metric;
        this._iter = 0;
        this.randomizer = new Randomizer(seed);
        this._Y = new Matrix(this._N, this._d, () => this.randomizer.random);
    }

    init(distance_matrix=null) {
        // init
        let Htarget = Math.log(this._perplexity);
        let D = distance_matrix || new Matrix(this._N, this._N, (i, j) => this._metric(this._X.row(i), this._X.row(j)));
        let P = new Matrix(this._N, this._N, "zeros");

        this._ystep = new Matrix(this._N, this._D, "zeros").to2dArray;
        this._gains = new Matrix(this._N, this._D, 1).to2dArray;

        // search for fitting sigma
        let prow = new Array(this._N).fill(0);
        for (let i = 0, N = this._N; i < N; ++i) {
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let done = false;
            let maxtries = 50;
            let tol = 1e-4;

            let num = 0;
            while(!done) {
                let psum = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = Math.exp(-D.entry(i, j) * beta);
                    if (i === j) pj = 0;
                    prow[j] = pj;
                    psum += pj;
                }
                let Hhere = 0;
                for (let j = 0; j < N; ++j) {
                    let pj = (psum === 0) ? 0 : prow[j] / psum;
                    prow[j] = pj;
                    if (pj > 1e-7) Hhere -= pj * Math.log(pj);
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
        let Pout = new Matrix(this._N, this._N, "zeros");
        let N2 = this._N * 2;
        for (let i = 0, N = this._N; i < N; ++i) {
            for (let j = 0; j < N; ++j) {
                Pout.set_entry(i, j, Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100));
            }
        }
        this._P = Pout;
        return this
    }

    set perplexity(value) {
        this._perplexity = value;
    }

    get perplexity() {
        return this._perplexity;
    }

    set epsilon(value) {
        this._epsilon = value;
    }

    get epsilon() {
        return this._epsilon;
    }

    transform(iterations=1000) {
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this._Y;
    }

    * transform_iter() {
        while (true) {
            this.next();
            yield this._Y;
        }
    }

    // perform optimization
    next() {
        let iter = ++this._iter;
        let P = this._P;
        let ystep = this._ystep;
        let gains = this._gains;
        let Y = this._Y;
        let N = this._N;
        let epsilon = this._epsilon;
        let dim = this._d;

        //calc cost gradient;
        let pmul = iter < 100 ? 4 : 1;
        
        // compute Q dist (unnormalized)
        let Qu = new Matrix(N, N, "zeros");
        let qsum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dsum = 0;
                for (let d = 0; d < dim; ++d) {
                    let dhere = Y.entry(i, d) - Y.entry(j, d);
                    dsum += dhere * dhere;
                }
                let qu = 1 / (1 + dsum);
                Qu.set_entry(i, j, qu);
                Qu.set_entry(j, i, qu);
                qsum += 2 * qu;
            }
        }

        // normalize Q dist
        let Q = new Matrix(N, N, (i, j) => Math.max(Qu.entry(i, j) / qsum, 1e-100));

        let cost = 0;
        let grad = [];
        for (let i = 0; i < N; ++i) {
            let gsum = new Array(dim).fill(0);
            for (let j = 0; j < N; ++j) {
                cost += -P.entry(i, j) * Math.log(Q.entry(i, j));
                let premult = 4 * (pmul * P.entry(i, j) - Q.entry(i, j)) * Qu.entry(i, j);
                for (let d = 0; d < dim; ++d) {
                    gsum[d] += premult * (Y.entry(i, d) - Y.entry(j, d));
                }
            }
            grad.push(gsum);
        }

        // perform gradient step
        let ymean = new Array(dim).fill(0);
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                let gid = grad[i][d];
                let sid = ystep[i][d];
                let gainid = gains[i][d];
                
                let newgain = Math.sign(gid) === Math.sign(sid) ? gainid * .8 : gainid + .2;
                if (newgain < .01) newgain = .01;
                gains[i][d] = newgain;

                let momval = iter < 250 ? .5 : .8;
                let newsid = momval * sid - epsilon * newgain * grad[i][d];
                ystep[i][d] = newsid;

                Y.set_entry(i, d, Y.entry(i, d) + newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < 2; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - ymean[d] / N);
            }
        }

        return this._Y;
    }

    get projection() {
        return this._Y;
    }
}

// http://optimization-js.github.io/optimization-js/optimization.js.html#line438
function powell(f, x0, max_iter=300) {
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
        alpha *= (pfx >= fx ? 1.05 : 0.4);
        pfx = fx;
    }
    return x;
}

class UMAP{
    constructor(X, local_connectivity, min_dist, d=2, metric=euclidean, seed=1212) {
        this._X = X;
        this._d = d;
        [ this._N, this._D ] = X.shape;
        this._local_connectivity = local_connectivity;
        this._min_dist = min_dist;
        this._metric = metric;
        this._iter = 0;
        this._n_neighbors = 11;
        this._spread = 1;
        this._set_op_mix_ratio = 1;
        this._repulsion_strength = 1;
        this._negative_sample_rate = 5;
        this._n_epochs = 200;
        this._initial_alpha = 1;
        this._randomizer = new Randomizer(seed);
        this._Y = new Matrix(this._N, this._d, () => this._randomizer.random);
    }

    _find_ab_params(spread, min_dist) {
        function curve(x, a, b) {
            return 1 / (1 + a * Math.pow(x, 2 * b));
        }
      
        var xv = linspace(0, spread * 3, 300);
        var yv = linspace(0, spread * 3, 300);
        
        for ( var i = 0, n = xv.length; i < n; ++i ) {
            if (xv[i] < min_dist) {
                yv[i] = 1;
            } else {
                yv[i] = Math.exp(-(xv[i] - min_dist) / spread);
            }
        }
      
        function err(p) {
            var error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map(e => e * e)));
        }
      
        var [ a, b ] = powell(err, [1,1]);
        return [ a, b ]
    }

    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            for (let j = 0, m = distances[i].length; j < m; ++j) {
                let v = distances[i][j].value - rhos[i];
                let value = 1;
                if (v > 0) {
                    value = Math.exp(-v / sigmas[i]);
                }
                distances[i][j].value = value;
            }
        }
        return distances;
    }

    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = this._local_connectivity;
        const bandwidth = 1;
        const target = Math.log2(k) * bandwidth;
        const rhos = [];
        const sigmas = [];
        const X = this._X;

        let distances = [];
        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            let x_i = X.row(i);
            distances.push(knn.search(x_i, Math.max(local_connectivity, k)).toArray().reverse());
        }

        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            let search_result = distances[i];
            rhos.push(search_result[0].value);

            let lo = 0;
            let hi = Infinity;
            let mid = 1;

            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    let d = search_result[j].value - rhos[i];
                    psum += (d > 0 ? Math.exp(-(d / mid)) : 1);
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    //[hi, mid] = [mid, (lo + hi) / 2];
                    hi = mid;
                    mid = (lo + hi) / 2; // PROBLEM mit hi?????
                } else {
                    lo = mid;
                    if (hi === Infinity) {
                        mid *= 2;
                    } else {
                        mid = (lo + hi) / 2;
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
        return {distances: distances, sigmas: sigmas, rhos: rhos}
    }

    _fuzzy_simplicial_set(X, n_neighbors) {
        const knn = new BallTree(X.to2dArray, euclidean);
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        let result = new Matrix(X.shape[0], X.shape[0], "zeros");
        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            for (let j = 0; j < n_neighbors; ++j) {
                result.set_entry(i, distances[i][j].element.index, distances[i][j].value);
            }
        }
        const transposed_result = result.T;
        const prod_matrix = result.mult(transposed_result);
        result = result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(this._set_op_mix_ratio)
            .add(prod_matrix.mult(1 - this._set_op_mix_ratio));
        return result;
    }

    _make_epochs_per_sample(graph, n_epochs) {
        const { data: weights } = this._tocoo(graph);
        let result = new Array(weights.length).fill(-1);
        const weights_max = Math.max(...weights);
        const n_samples = weights.map(w => n_epochs * (w / weights_max));
        result = result.map((d, i) => (n_samples[i] > 0) ? Math.round(n_epochs / n_samples[i]) : d);
        return result;
    }

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
        return {rows: rows, cols: cols, data: data};
    }

    init() {
        const [ a, b ] = this._find_ab_params(this._spread, this._min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this._X, this._n_neighbors);
        this._epochs_per_sample = this._make_epochs_per_sample(this._graph, this._n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map(d => d * this._negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        const { rows, cols } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        return this
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

    transform(iterations = 1000) {
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this._Y;
    }

    * transform_iter() {
        this._iter = 0;
        while (this._iter < this._n_epochs) {
            this.next();
            yield this._Y;
        }
        return this._Y;
    }

    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

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
                const dist = euclidean(current, other);//this._metric(current, other);
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
                    const dist = euclidean(current, other);//this._metric(current, other);
                    let grad_coeff = 0;
                    if (dist > 0) {
                        grad_coeff = (2 * repulsion_strength * b) / ((.01 + dist) * (a * Math.pow(dist, b) + 1));
                    } else if (j == k) {
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

    next() {
        let iter = ++this._iter;
        let Y = this._Y;

        this._alpha = (this._initial_alpha * (1 - iter / this._n_epochs));
        this._Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this._Y;
    }

    get projection() {
        return this._Y;
    }
}

//import { Matrix } from "../matrix/index";

/**
 * 
 */
class OAP {
    constructor(X, depth_field_lag, step_size, depth_weight, d = 2, metric = euclidean, seed = 1212) {
        this._X = X;
        this._d = d;
        [this._N, this._D] = X.shape;
        this._depth_field_lag = depth_field_lag;
        this._step_size = step_size;
        this._depth_weight = depth_weight;
        this._J = 3;
        this._max_iter = 1;
        this._metric = metric;
        this._seed = seed;
        this._randomizer = new Randomizer(seed);
    }

    _data_depth(technique = "chebyshev") {
        const X = this._X;
        const N = this._N;
        const h = new Float32Array(N);
        let deepest_point = 0;
        if (technique === "mdb") {
            h.fill(1);

            /*
            // Modified Band Depth 
            // https://www.tandfonline.com/doi/pdf/10.1198/jasa.2009.0108?casa_token=E1Uzntgs-5AAAAAA:Eo8mUpJDhpLQ5RHBkCB3Mdz0tbGM3Q0v78bwyCIAv7-peLGwfG3TcXLqShIaYuJLEqKc7GvaKlgvUg 
            const randomizer = this._randomizer;
            const h = new Float32Array(this._N);
            const J = this._J;
            const N = this._N;
            const D = this._D;
            const X = this._X;

            const one_div_by_n_choose_j = 1;
            for (let row = 0; row < N; ++row) {
                const x = X.row(row);
                const B_min = new Float32Array(D).fill(Infinity);
                const B_max = new Float32Array(D).fill(-Infinity);
                let r = Math.floor(randomizer.random * N);
                for (let i = 0; i < J; ++i) {
                    const x_j = X.row(r);
                    for (let d = 0; d < D; ++d) {
                        const x_jd = x_j[d]
                        B_min[d] = Math.min(B_min[d], x_jd);
                        B_max[d] = Math.max(B_max[d], x_jd);
                    }
                    r += Math.floor(randomizer.random * (N - 1));
                    r = r % N;
                }
                for (let d = 0; d < D; ++d) {
                    const x_d = x[d];
                    if (x_d >= B_min[d] && x_d <= B_max[d]) {
                        ++h[row]
                    }
                }
            }
            this._h = h;*/
        } else if (technique === "chebyshev") {
            // L∞ Depth
            // https://arxiv.org/pdf/1506.01332.pdf
            for (let i = 0; i < N; ++i) {
                let x = X.row(i);
                let sum = 0;
                for (let j = 0; j < N; ++j) {
                    if (i !== j) {
                        sum += chebyshev(x, X.row(j));
                    }
                }
                h[i] = 1 / (1 + sum / N);
                if (h[deepest_point] < h[i]) {
                    deepest_point = i;
                }
            }
        }
        this._h = h;
        this._deepest_point = deepest_point;

    }

    init() {
        this._iter = 0;
        // init with MDS
        const init_MDS = new MDS(this._X, this._d, this._metric);
        //console.log(init_MDS)
        this._Y = init_MDS.transform();

        // try häääh?
        this._X_distances = init_MDS._d_X;
        /*let max = -Infinity
        init_MDS._d_X._data.forEach(dx => max = Math.max(dx, max));
        this._X_distances = init_MDS._d_X.divide(max);*/
        // end try hääääh?
        
        // compute order statistics
        this._data_depth();
        this._M = this._monotonic_field(this._Y);
        //
        return this;
    }

    set depth_field_lag(value) {
        this._depth_field_lag = value;
    }

    get depth_field_lag() {
        return this._depth_field_lag;
    }

    set step_size(value) {
        this._step_size = value;
    }

    get step_size() {
        return this._step_size;
    }

    set depth_weight(value) {
        this._depth_weight = value;
    }

    get depth_weight() {
        return this._depth_weight;
    }

    transform(iterations = this._max_iter) {
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this._Y;
    }

    * transform_iter() {
        while (true) {
            this.next();
            yield this._Y;
        }
    }

    _monotonic_field(Y) {
        const h = this._h;
        const Y_ = this._Y_;
        const nn = new BallTree();
        nn.add(Y.to2dArray);

        const N = 5;
        let M = (x) => {
            let neighbors = nn.search(x, N).toArray();
            let d_sum = 0;//neighbors.reduce((a, b) => a + b.value, 0);
            let m = 0;
            for (let i = 0; i < N; ++i) {
                d_sum += neighbors[i].value;
                m += h[neighbors[i].element.index] * neighbors[i].value;
            }
            //console.log(m, d_sum)
            m /= d_sum;
            return m;
        };
        return M;
    }

    next() {
        const iter = ++this._iter;
        const l = this._depth_field_lag;
        const step_size = this._step_size;
        const w = this._depth_weight;
        const N = this._N;
        const dim = this._d;
        const d_X = this._X_distances;
        const h = this._h;
        let Y = this._Y;

        if ((iter % l) === 1) {
            // compute monotonic field
            this._Y_ = this._Y.clone();
            this._M = this._monotonic_field(Y);
        }
        const M = this._M;
        // perform gradient step

        // MDS stress step
        /*for (let i = 0; i < N; ++i) {
            const d_x = d_X.row(i);
            const y_i = Y.row(i)
            const delta_mds_stress = new Float32Array(dim);
            for (let j = 0; j < N; ++j) {
                if (i !== j) {
                    const y_j = Y.row(j)
                    const d_y = metric(y_i, y_j);
                    const d_x_j = d_x[j] === 0 ? 1e-2 : d_x[j]
                    const mult = 1 - (d_x_j / d_y)
                    for (let d = 0; d < dim; ++d) {
                        delta_mds_stress[d] += (mult * (y_i[d] - y_j[d]));
                    }
                }
            }
            for (let d = 0; d < dim; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - step_size * delta_mds_stress[d] / N)
            }
        }*/
        
        // MDS stress step
        const d_Y = new Matrix();
        d_Y.shape = [N, N, (i, j) => {
            return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
        }];
        const ratio = new Matrix();//d_X.divide(d_Y).mult(-1);
        ratio.shape = [N, N, (i, j) => {
            if (i === j) return 1e-8
            return i < j ? -d_X.entry(i, j) / d_Y.entry(i, j) : ratio.entry(j, i);
        }];
        for (let i = 0; i < N; ++i) {
            ratio.set_entry(i, i, ratio.entry(i, i) - neumair_sum(ratio.row(i)));
        }
        const mds_Y = ratio.dot(Y).divide(N);

        // Data depth step
        const diff_Y = new Matrix(N, dim, (i, j) => mds_Y.entry(i, j) - Y.entry(i, j));

        for (let i = 0; i < N; ++i) {
            const m = M(Y.row(i));
            const dm = M(mds_Y.row(i));
            const h_i = h[i];
            for (let d = 0; d < dim; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) + step_size * (diff_Y.entry(i, d) + w * 2 * (m - h_i) * dm));
            }
        }

        this._Y = Y;

        return this._Y;
    }

    get projection() {
        return this._Y;
    }
}

exports.Randomizer = Randomizer;
exports.kahan_sum = kahan_sum;
exports.neumair_sum = neumair_sum;
exports.euclidean = euclidean;
exports.euclidean_squared = euclidean_squared;
exports.cosine = cosine;
exports.manhattan = manhattan;
exports.chebyshev = chebyshev;
exports.k_nearest_neighbors = k_nearest_neighbors;
exports.distance_matrix = dmatrix;
exports.linspace = linspace;
exports.norm = norm;
exports.Matrix = Matrix;
exports.Reorder = Reorder;
exports.HNSW = HNSW;
exports.BallTree = BallTree;
exports.NNDescent = NNDescent;
exports.Heap = Heap;
exports.qr = qr;
exports.qr_householder = qr_householder;
exports.simultaneous_poweriteration = simultaneous_poweriteration;
exports.lu = lu;
exports.svrg = svrg;
exports.poweriteration_m = poweriteration_m;
exports.poweriteration_n = poweriteration_n;
exports.PCA = PCA;
exports.MDS = MDS;
exports.ISOMAP = ISOMAP;
exports.FASTMAP = FASTMAP;
exports.LDA = LDA;
exports.LLE = LLE;
exports.MLLE = MLLE;
exports.LTSA = LTSA;
exports.TSNE = TSNE;
exports.UMAP = UMAP;
exports.OAP = OAP;
exports.powell = powell;
exports.Hierarchical_Clustering = Hierarchical_Clustering;

Object.defineProperty(exports, '__esModule', { value: true });

}));
