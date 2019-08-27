// https://null.org v0.0.1 Copyright 2019 abc
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = global || self, factory(global.sci = global.sci || {}));
}(this, function (exports) { 'use strict';

function create_array(rows, columns, value_function) {
    if (rows < 0 || columns < 0) {
        return undefined;
    } else if (rows === 0 && columns === 0) {
        return value_function(0,0);
    } else if (columns === 0) {
        let A = new Array(rows);
        for (let i = 0; i < rows; ++i) {
            A[i] = value_function(i);
        }
        return A;
    } else {    
        let A = new Array(rows);
        for (let i = 0; i < rows; ++i) {
            A[i] = new Array(columns);
            for (let j = 0; j < columns; ++j) {
                A[i][j] = value_function(i,j);
            }
        }
        return A;
    }
}

function zeros(rows = 0, cols = 0) {
    /*if (n < 0 || m < 0) {
        return undefined
    } else if (n === 0 && m === 0) {
        return 0
    } else if (m === 0) {
        let A = new Array(n)
        for (let i = 0; i < n; ++i) {
            A[i] = 0
        }
        return A
    } else {    
        let A = new Array(n)
        for (let i = 0; i < n; ++i) {
            A[i] = new Array(m)
            for (let j = 0; j < m; ++j) {
                A[i][j] = 0
            }
        }
        return A
    }*/
    return create_array(rows, cols, () => 0)
}

// kahan summation algorithm
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm

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

// neumair summation algorithm
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements

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

function euclidean(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    /*let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += ((a[i] - b[i]) * (a[i] - b[i]))
    }*/
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        s[i] = ((a[i] - b[i]) * (a[i] - b[i]));
    }
    //return Math.sqrt(sum)
    return Math.sqrt(neumair_sum(s))
}

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

function chebyshev(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(res)
}

const euclidean$1 = euclidean;

function dmatrix(A, metric = euclidean$1) {
    let distance = metric;
    if (distance === undefined) return undefined
    let n = A.length;
    let D = zeros(n,n);
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = distance(A[i], A[j]);
        }
    }
    return D
}

const euclidean$2 = euclidean;

function k_nearest_neighbors(A, k, distance_matrix = null, metric = euclidean$2) {
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

function random_array(n, m) {
    return create_array(n, m, () => Math.random())
    /*if (n < 0 || m < 0) {
        return undefined;
    } else if (n === 0 && m === 0) {
        return Math.random();
    } else if (m === 0) {
        let A = new Array(n);
        for (let i = 0; i < n; ++i) {
            A[i] = Math.random();
        }
        return A;
    } else {    
        let A = new Array(n);
        for (let i = 0; i < n; ++i) {
            A[i] = new Array(m);
            for (let j = 0; j < m; ++j) {
                A[i][j] = Math.random();
            }
        }
        return A;
    }*/
}

function identity(n) {
    return create_array(n, n, (i,j) => (i === j) ? 1 : 0);
}

function dot(A, B) {
    let A_rows = A.length;
    let A_cols = Array.isArray(A[0]) ? A[0].length : 0;
    let B_rows = B.length;
    let B_cols = Array.isArray(B[0]) ? B[0].length : 0;
    let result = null;
    if (A_cols && B_cols) {
        if (A_cols !== B_rows) return undefined;
        result = zeros(A_rows,B_rows);
        for (let i = 0; i < A_rows; ++i) {
            for (let j = 0; j < B_cols; ++j) {
                for (let k = 0; k < A_cols; ++k) {
                    result[i][j] += (A[i][k] * B[k][j]);
                }
            }
        }
    } else if (A_cols && (B_cols == 0)) {
        result = zeros(A_rows,0);
        for (let i = 0; i < A_rows; ++i) {
            for (let k = 0; k < A_cols; ++k) {
                result[k] += (A[i][k] * B[k]);
            }
        }
    } else if ((A_cols == 0) && B_cols) {
        result = zeros(A_rows,0);
        for (let i = 0; i < A_rows; ++i) {
            for (let j = 0; j < B_rows; ++j) {
                result[j] = (A[i] * B[j][i]);
            }
        }
    } else if ((A_cols == 0) && (B_cols == 0)) {
        if (A_rows !== B_rows) return undefined;
        result = 0;
        for (let i = 0; i < A_rows; ++i) {
            result += (A[i] * B[i]);
        }
    } else {
        return undefined;
    }

    return result;
}

function mult(A, x) {
    if (!Array.isArray(A)) return A * x
    let n = A.length;
    let result;
    if (Array.isArray(A[0])) {
        let m = A[0].length;
        result = [];
        for (let i = 0; i < n; ++i) {
            let result_i = [];
            for (let j = 0; j < m; ++j) {
                result_i[j] = A[i][j] * x;    
            }
            result.push(result_i);
        }
    } else {
        result = [];
        for (let i = 0; i < n; ++i) {
            result.push(A[i] * x);
        }
    }
    return result;
}

function divide(A, x) {
    return mult(A, 1/x);
}

function norm(v, metric = euclidean) {
    let n = v.length;
    return metric(v, zeros(n,0));
}

function transpose(A) {
    let n = A.length;
    let result;
    if (Array.isArray(A[0])) {
        let m = A[0].length;
        result = new Array(m);
        for (let i = 0; i < m; ++i) {
            let result_i = new Array(n);
            for (let j = 0; j < n; ++j) {
                result_i[j] = A[j][i];
            }
            result[i] = result_i;
        }
    } else {
        result = new Array(n);
        for ( let i = 0; i < n; ++i) {
            result[i] = [A[i]];
        }
    }
    return result;
}

function add(A,x) {
    let n = A.length;
    let A_isArray = Array.isArray(A);
    let x_isArray = Array.isArray(x);
    let result;
    if (A_isArray && x_isArray) {
        let m = x.length;
        if (n !== m) return undefined;
        result = new Array(n);
        for (let i = 0; i < n; ++i) {
            result[i] = A[i] + x[i];
        }
        return result
    }
    if (Array.isArray(A[0])) {
        let m = A[0].length;
        result = new Array(n);
        for (let i = 0; i < n; ++i) {
            result[i] = new Array(m);
            for (let j = 0; j < m; ++j) {
                result[i][j] = A[i][j] + x;
            }
        }
    } else {
        result = new Array(n);
        for (let i = 0; i < n; ++i) {
            result[i] = A[i] + x;
        }
    }
    return result;
}

function sub(A,x) {
    return add(A, mult(x, -1));
}

function poweriteration(A, max_iterations = 100, metric = euclidean) {
    let n = A.length;
    let r = random_array(n,0);
    r = divide(r, norm(r, metric));
    while (max_iterations--) {
        let r_next = dot(A, r);
        r = divide(r_next, norm(r_next, metric));
    }
    let u = dot(dot(transpose(A), r), r);
    let l = dot(r, r);
    return {
        "eigenvector" : r, 
        "eigenvalue": u / l
    };
}

function lanczos(A, k = 2, max_iterations = 100, metric = euclidean) {
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0];
    let v = new Matrix(n, 1, () => Math.random());
    v = v.divide(norm(v._data, metric));
    
    let w_ = A.dot(v);
    let a = w_.dot(v);
    let w = w_.sub(a.dot(v));

    while (max_iterations--) {
        let b = norm(w._data, euclidean);
    }
    return {
        "eigenvector" : r, 
        "eigenvalue": u / l
    };
}

/**
 * Computes largest Eigenvalue / Eigenvector with 
 * Accelerated Stochastic Power Iteration
 * http://proceedings.mlr.press/v84/xu18a/xu18a-supp.pdf
 * @param {Matrix} A - The respective matrix
 * @param {Number} max_iterations - number of maximal iterations
 * @param {Number} batch_size - defines batchsize
 * @param {Number} beta - learning parameter
 * @param {Metric} metric - for computing the norm 
 */
function sapi(A, max_iterations = 100, batch_size = 10, beta = .05, metric = euclidean) {
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0];
    let r = new Matrix(n, 1, () => Math.random());
    let r_last = new Matrix(n, 1, 0);

    while (max_iterations--) {
        let A_r = A.dot(r);
        let beta_r_last = r_last.mult(beta);
        let r_next = A_r.sub(beta_r_last);
        let r_next_norm = norm(r_next._data, metric);
        r_last = r.divide(r_next_norm);
        r = r_next.divide(r_next_norm);
    }

    let u = r.transpose().dot(A).dot(r);
    let l = r.transpose().dot(r);
    let lambda = u.divide(l).entry(0,0);
    return {
        "eigenvector": r.transpose().to2dArray[0],
        "eigenvalue": lambda
    }

}

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

class Matrix{
    /*constructor(A = null) {
        this._data = null;
        this._rows = null;
        this._cols = null;
        if (A) {
            this.entries = A;
        }
        return this;
    }*/
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

    /*set entries(A) {
        if (!Array.isArray(A)) throw "not an array"
        if (A.length == 0) throw "empty array"
        if (!Array.isArray(A[0])) throw "not a 2d Array"
        if (A[0].length == 0) throw "empty array"
        let rows = A.length;
        let cols = A[0].length;
        this._rows = rows;
        this._cols = cols;
        let data = this._data = new Float64Array(rows * cols);
        for (let row = 0; row < rows; ++row) {
            let A_row = A[row]
            for (let col = 0; col < cols; ++col) {
                data[row * cols + col] = A_row[col];
            }
        }
        return this;
    }*/

    static from(A, type) {
        //let B;
        if (A instanceof Matrix) {
            return A.clone();
        } else if (Array.isArray(A)) {
            //B = new Matrix()
            let m = A.length;
            if (m === 0) return B;
            // 1d
            if (!Array.isArray(A[0])) {
                if (type === "row" || type === null) {  
                    //B.shape = [1, m, (_, j) => A[j]];
                    //return B;
                    return new Matrix(1, m, (_, j) => A[j]);
                } else if (type === "col") {
                    //B.shape = [m, 1, (i, _) => A[i]];
                    //return B;
                    return new Matrix(m, 1, (i, _) => A[i]);
                } else {
                    throw "1d array has NaN entries"
                }
            // 2d
            } else if (Array.isArray(A[0])) {
                let n = A[0].length;
                for (let row = 0; row < m; ++row) {
                    //let A_row = ;
                    if (A[row].length !== n) throw "various array lengths";
                    
                    /*for (let col = 0; col < n; ++col) {
                        if (!Number.isNaN(A_row[col])) {
                            B.set_entry(row, col, A_row[col]);
                        } else {
                            throw "2d array has NaN entries";
                        }
                    }*/
                }
                //B.shape = [m, n, (i, j) => A[i][j]];
                //return B;
                return new Matrix(m, n, (i, j) => A[i][j])
            }
        } else if (typeof(A) === "number") {
            return new Matrix(1, 1, A);
        } else {
            throw "error"
        }
    }

    row(row) {
        let result_row = new Array(this._cols);
        for (let col = 0; col < this._cols; ++col) {
            result_row[col] = this._data[row * this._cols + col];
        }
        return result_row;
    }

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

    col(col) {
        let result_col = new Array(this._rows);
        for (let row = 0; row < this._rows; ++row) {
            result_col[row] = this._data[row * this._cols + col];
        }
        return result_col;
    }

    entry(row, col) {
        return this._data[row * this._cols + col];
    }

    set_entry(row, col, value) {
        this._data[row * this._cols + col] = value;
    }

    transpose() {
        let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
        return B;
    }

    get T() {
        return this.transpose();
    }

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

    dot(B) {
        if (B instanceof Matrix) {
            let A = this;
            if (A.shape[1] !== B.shape[0]) return undefined;
            let I = A.shape[1];
            let C = new Matrix(A.shape[0], B.shape[1], (row, col) => {
                let A_i = A.row(row);
                let B_i = B.col(col);
                for (let i = 0; i < I; ++i) {
                    A_i[i] = A_i[i] * B_i[i];
                }
                return neumair_sum(A_i);
            });
            return C;
        } else if (Array.isArray(B)) {
            let rows = this._rows;
            if (B.length !== rows) return undefined;
            let C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map(e => e * B[row]));//.reduce((a,b) => a + b);
            }
            return C;
        } else {
            return undefined;
        }
    }

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

    householder(o) {
        let v = this.col(o).slice(o);
        let v_1 = v[0];
        let v_norm = norm(v);
        let α = -(v_1 < 0 ? -1 : 1) * v_norm;
        v[0] -= α;
        v = Matrix.from(v, "col");
        let β = 1 / (v_norm * (Math.abs(v_1) + v_norm));

        let H = v.mult(β).outer(v);
        let A = this.get_block(o, o);
        A = A.sub(H.dot(A));
        console.log(H, A);
        this.set_block(o, o, A);

        return this;
    }

    householder2(o) {
        /*let y = this.row(o).slice(o + 1);
        let y1 = y[0];
        let y_norm = norm(y);
        let alpha = (y1 < 0) ? y_norm : -y_norm;
        let omega = y;
        omega[0] -= alpha;
        omega = Matrix.from(omega, "col");
        let beta = 1 / (y_norm * (Math.abs(y1) + y_norm))

        let O = new Matrix(this._rows, this._cols, "identity")
        O.set_block(o + 1, o + 1, omega.mult(beta).outer(omega));
        console.log(O)
        return this.sub(O.dot(this));*/

        let v = this.row(o).slice(o);
        let n = v.length;
        let v_1 = v[1];
        let v_norm = norm(v);
        let α = -(v_1 < 0 ? -1 : 1) * v_norm;
        v[1] -= α;
        v = Matrix.from(v, "col");
        let β = 1 / (v_norm * (Math.abs(v_1) + v_norm));

        //let I = new Matrix(n + 1, n + 1, "identity");
        let H = v.mult(β).outer(v);
        //console.log(I, v,β,  H)
        //H = I.set_block(1, 1, H);

        let A = this.get_block(o, o);
        A = A.sub(H.dot(A));
        console.log(H, A);
        this.set_block(o, o, A);

        return this;
    }

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
            A.set_block(k, k, A_);
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

    get_block(offset_rows, offset_cols) {
        let [ rows, cols ] = this.shape;
        return new Matrix(rows - offset_rows, cols - offset_cols, (i, j) => this.entry(i + offset_rows, j + offset_cols));
    }

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

    get shape() {
        return [this._rows, this._cols];
    }

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

    get mean() {
        const data = this._data;
        const n = this._rows * this._cols;
        let sum = 0;
        for (let i = 0; i < n; ++i) {
            sum += data[i];
        }
        return sum / n;
    }

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

    static solve(A, b) {
        let [ rows, cols ] = A.shape;
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

    static SVD(M, k=2) {
        let MtM = M.transpose().dot(M);
        let MMt = M.dot(M.transpose());
        let { eigenvectors: V, eigenvalues: Σ } = simultaneous_poweriteration(MtM, k);
        let { eigenvectors: U } = simultaneous_poweriteration(MMt, k);

        return { U: U, Σ: Σ.map(ς => Math.sqrt(ς)), V: V };
        /*const [ rows, cols ] = matrix.shape;
        // 1.
        let U = new Matrix();
        U.shape = [rows, cols, (i, j) => i === j ? 1 : 0];
        let Σ = matrix.clone();
        let Vt = new Matrix();
        Vt.shape = [cols, cols, (i, j) => i === j ? 1 : 0];
        let I = new Matrix();
        I.shape = [cols, cols, (i, j) => i === j ? 1 : 0];

        function householder_vector(x) {
            let dot_1on = norm(x.col(0).slice(1));
            let v = x.clone();
            v.set_entry(0, 0, 1);
            let β = 0;
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
                β = 2 * v_0 ** 2 / (dot_1on + v_0 ** 2)
                v.divide(v_0)
            }
            return {"v": v, "β": β}
        }

        function householder_matrix(size, col, v, β) {
            let vvt = v.dot(v.transpose())
            let V = new Matrix()
            V.shape = [size, size, (r, c) => {
                if (r < col || c < col) {
                    return r === c ? 1 : 0;
                }
                return (r === c ? 1 : 0) - (β * vvt.entry(r - col, c - col));
            }]

            console.log(V, vvt)
            return V;
        }

        for (let col = 0; col < cols; ++col) {
            let x = Matrix.from(Σ.col(col), "col")
            let { v: u, β: β } = householder_vector(x);
            let Q = householder_matrix(cols, col, u, β);
            console.log("u", u, Q)
            U = U.dot(Q)
            Σ = Q.dot(Σ);

            if (col < cols - 2) {
                let x = Matrix.from(Σ.row(col).slice(1), "col")
                let { v: v, β: β } = householder_vector(x);
                Q = householder_matrix(cols, col + 1, v, β);
                console.log("v", v, Q)
                Vt = Q.dot(Vt)
                Σ = Σ.dot(Q);
            }

            /*let u = Σ.col(col);
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
            Σ = H.dot(Σ);
            U = U.dot(H.transpose())

            if (col < cols - 1) {

            }*/

            /*let u = Σ.col(col).slice(col);
            let u_norm = norm(u);
            let sign_u = u[0] >= 0 ? 1 : -1;
            u[0] = sign_u * (Math.abs(u[0]) + u_norm);//+= sign_u * u_norm;
            u_norm = norm(u)//Math.sqrt(2 * u_norm * (u_norm + u[0]))//norm(u);
            u = u.map(u_i => u_i / u_norm);
            u = Matrix.from(u, "col")
            let A_k_m = new Matrix()
            A_k_m.shape = [rows - col, cols - col, (i, j) => Σ.entry(i + col, j + col)];
            let U_k = u.dot(u.transpose().dot(A_k_m)).mult(2)
            for (let r = 0; r < rows - col; ++r) {
                for (let c = 0; c < cols - col; ++c) {
                    let val = Σ.entry(col + r, col + c) - U_k.entry(r, c);
                    val = Math.abs(val) < 1e-15 ? 0 : val;
                    Σ.set_entry(col + r, col + c, val)
                }
            }
            console.log(Σ)



            /*if (col <= cols - 2) {
                let col_1 = col + 1
                let v = Σ.row(col).slice(col_1);
                let v_norm = norm(v);
                let sign_v = v[0] >= 0 ? 1 : -1;
                v[0] += sign_v * (Math.abs(v[0]) + v_norm);//sign_v * v_norm;
                v_norm = norm(v);
                v = v.map(v_i => v_i / v_norm);
                v = Matrix.from(v, "col")
                let A_k1_n = new Matrix()
                A_k1_n.shape = [rows - col, cols - col_1, (i, j) => Σ.entry(col + i, col_1 + j)];
                let V_k = A_k1_n.dot(v.dot(v.transpose())).mult(2)
                for (let r = 0; r < rows - col; ++r) {
                    for (let c = 0; c < cols - col_1; ++c) {
                        let val = Σ.entry(col + r, col_1 + c) - V_k.entry(r, c)
                        console.log("V", [col + r, col_1 + c], 2 * Σ.entry(col + r, col_1 + c) == V_k.entry(r, c))
                        val = Math.abs(val) < 1e-15 ? 0 : val;
                        Σ.set_entry(col + r, col_1 + c, val)
                    }
                }
            }

            */
            /*console.log(Σ)
        }*/
    }
}

class Vector extends Array{
    constructor() {
        super();
        return this;
    }

    get norm() {
       return this; 
    }
}

class Heap {
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

        console.log(arr)
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

    *iterate() {
        let pointer = this.root;
        /*do {
            yield pointer.element;
        } while (pointer = pointer.next);*/
        while (pointer) {
            yield pointer.element;
            pointer = pointer.next;
        }
    }

    toArray() {
        let res = [];
        let pointer = this.root;
        while (pointer) {
            res.push(pointer.element);
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
}

class HNSW {
    /**
     * 
     * @param {*} metric metric to use: (a, b) => distance
     * @param {*} heuristic use heuristics or naive selection
     * @param {*} m max number of connections
     * @param {*} ef size of candidate list
     * @param {*} m0 max number of connections for ground layer 
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
        let X = this.X;
        let rows = X.shape[0];
        let ai_ = [];
        let a_j = [];
        for (let i = 0; i < rows; ++i) {
            ai_.push(0);
            a_j.push(0);
        }
        let a__ = 0;
        let A = new Matrix();
        A.shape = [rows, rows, (i,j) => {
            let val = 0;
            if (i < j) {
                val = this._metric(X.row(i), X.row(j));
            } else if (i > j) {
                val = A.entry(j,i);
            }
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        }];
        ai_ = ai_.map(v => v / rows);
        a_j = a_j.map(v => v / rows);
        a__ /= (rows ** 2);
        let B = new Matrix(rows, rows, (i, j) => (A.entry(i, j) - ai_[i] - a_j[j] + a__));
        //B.shape = [rows, rows, (i,j) => (A.entry(i,j) - (A.row(i).reduce(sum_reduce) / rows) - (A.col(j).reduce(sum_reduce) / rows) + a__)]
                
        let { eigenvectors: V } = simultaneous_poweriteration(B, this.d);
        this.Y = Matrix.from(V).transpose();
        
        return this.Y
    }

    get projection() {
        return this.Y
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
        let Φ = new Matrix(rows, rows);

        let V = new Array(rows);
        let Λ = new Array(rows);
        let P = new Array(rows);

        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map(n => n.j);
            let x_i = Matrix.from(X.row(row), "row");
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            X_i = X_i.sub(x_i);
            //X_i = X_i.dot(new Matrix(X_i._cols, X_i._cols, "center"))
            let C_i = X_i.dot(X_i.transpose()); // k by k

            let γ = neumair_sum(C_i.diag) / 1000;
            for (let j = 0; j < k; ++j) {
                C_i.set_entry(j, j, C_i.entry(j, j) + γ);
            }
            
            let { eigenvalues: λs, eigenvectors: v } = simultaneous_poweriteration(C_i, k);
            V[row] = v; // k by k, rows are eigenvectors, big to small
            Λ[row] = λs; // 1 by k, cols are eigenvalues, big to small
            P.push(neumair_sum(λs.slice(d + 1)) / neumair_sum(λs.slice(0, d)));

            // reconstruct;
            let w = Matrix.solve(C_i, O); // k by 1
            let w_sum = neumair_sum(w.col(0));
            w = w.divide(w_sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, j, w.entry(j, 0));
            }
        }
        // find regularized weights // median
        let η = P.sort((ρ_i, ρ_j) => ρ_i - ρ_j)[Math.ceil(rows / 2)];
        
        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map(n => n.j);
            let λs = Λ[row]; // 1 by k
            let s_i = λs.map((λ, l) => {
                    return {
                        "l": l,
                        "ratio": neumair_sum(λs.slice(k - l + 1)) / neumair_sum(λs.slice(0, k - l)),
                        "λ": λ
                    }
                });
            //console.log(s_i)
            s_i = s_i
                .filter(s => s.ratio < η && s.l <= k - d)
                .map(s => s.l).pop() || d;
            let V_i = V[row]; // k by k
            V_i = V_i.slice(k - s_i); // s_i by k
            let α_i = (1 / Math.sqrt(s_i)) * norm(V_i[0].map((_, j) => neumair_sum(V_i.map(r => r[j]))));
            V_i = Matrix.from(V_i); // s_i by k
            
            //https://github.com/scikit-learn/scikit-learn/blob/7b136e9/sklearn/manifold/locally_linear.py#L703

            let h = new Matrix(s_i, 1, α_i);
            let ones = new Matrix(k, 1, 1);
            h = h.sub(V_i.dot(ones));
            let h_norm = norm(h.col(0));
            h = h_norm < 1e-12 ? h.mult(0) : h.divide(h_norm);
            V_i = V_i.T;
            ones = new Matrix(s_i, 1, 1);
            let w_i = Matrix.from(W.row(row), "col");
            
            /*let H_i = new Matrix(s_i, s_i, "identity");
            H_i = H_i.sub(h.mult(2).outer(h));
            let W_i = V_i.sub(V_i.dot(h).dot(h.T).mult(2)).add(w_i.mult(1 - α_i))
            */
            let W_i = V_i.sub(V_i.dot(h).dot(h.T).mult(2)).add(w_i.mult(1 - α_i).dot(ones.T));
            
            W_i = W_i.dot(W_i.T);
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < s_i; ++j) {
                    Φ.set_entry(I_i[i], I_i[j], Φ.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }
        //let { eigenvectors: Y } = simultaneous_poweriteration(Φ.inverse(), d + 1);
        //this.Y = Matrix.from(Y.slice(1)).transpose()

        let { eigenvectors: Y } = simultaneous_poweriteration(Φ, d + 1);
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
        this._Y = new Matrix(this._N, this._N, () => this.randomizer.random);
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

class UMAP{
    constructor(X, local_connectivity, min_dist, d=2, metric=euclidean) {
        this._X = X;
        this._d = d;
        [ this._N, this._D ] = X.shape;
        this._local_connectivity = local_connectivity;
        this._min_dist = min_dist;
        this._metric = metric;
        this._iter = 0;
    }

    init(distance_matrix=null, seed=1212) {
        // init
        this.randomizer = new Randomizer(seed);
        this._Y = new Matrix(this._N, this._N, () => this.randomizer.random);

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

exports.k_nearest_neighbors = k_nearest_neighbors;
exports.distance_matrix = dmatrix;
exports.zeros = zeros;
exports.linspace = linspace;
exports.random_array = random_array;
exports.create_array = create_array;
exports.identity = identity;
exports.dot = dot;
exports.mult = mult;
exports.divide = divide;
exports.norm = norm;
exports.transpose = transpose;
exports.sub = sub;
exports.add = add;
exports.Matrix = Matrix;
exports.Vector = Vector;
exports.HNSW = HNSW;
exports.Heap = Heap;
exports.euclidean = euclidean;
exports.cosine = cosine;
exports.manhattan = manhattan;
exports.chebyshev = chebyshev;
exports.poweriteration = poweriteration;
exports.lanczos = lanczos;
exports.sapi = sapi;
exports.qr = qr;
exports.simultaneous_poweriteration = simultaneous_poweriteration;
exports.lu = lu;
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
exports.Randomizer = Randomizer;
exports.kahan_sum = kahan_sum;
exports.neumair_sum = neumair_sum;

Object.defineProperty(exports, '__esModule', { value: true });

}));
