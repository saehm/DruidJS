var version$1 = "0.8.0";
var pkg = {
	version: version$1};

/**
 * Computes the Bray-Curtis distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Bray-Curtis distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity}
 */
function bray_curtis(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    let sum_abs_diff = 0;
    let sum_ab = 0;
    for (let i = 0; i < a.length; ++i) {
        sum_abs_diff += Math.abs(a[i] - b[i]);
        sum_ab += a[i] + b[i];
    }
    return sum_abs_diff / sum_ab;
}

/**
 * Computes the canberra distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The canberra distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
function canberra(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i]));
    }
    return sum;
}

/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The chebyshev distance between `a` and `b`.
 */
function chebyshev(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    const res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res);
}

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
function cosine(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
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
 * Computes the squared euclidean distance (l_2^2) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The squared euclidean distance between `a` and `b`.

 */
function euclidean_squared(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        const a_b = a[i] - b[i];
        sum += a_b * a_b;
    }
    return sum;
}

/**
 * Computes the euclidean distance (`l_2`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The euclidean distance between `a` and `b`.
 */
function euclidean(a, b) {
    return Math.sqrt(euclidean_squared(a, b));
}

/**
 * Computes the Goodman-Kruskal gamma coefficient for ordinal association.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First categorical/ordinal variable
 * @param {number[] | Float64Array} b - Second categorical/ordinal variable
 * @returns {number} The Goodman-Kruskal gamma coefficient between `a` and `b` (-1 to 1).
 * @see {@link https://en.wikipedia.org/wiki/Goodman_and_Kruskal%27s_gamma}
 */
function goodman_kruskal(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    if (n < 2) return 0;

    let concordant = 0;
    let discordant = 0;
    let tie_a = 0;
    let tie_b = 0;

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const a_diff = a[i] - a[j];
            const b_diff = b[i] - b[j];
            const a_tied = a_diff === 0;
            const b_tied = b_diff === 0;

            if (a_tied && b_tied) ; else if (a_tied) {
                tie_a++;
            } else if (b_tied) {
                tie_b++;
            } else if (a_diff * b_diff > 0) {
                concordant++;
            } else {
                discordant++;
            }
        }
    }

    const denominator = concordant + discordant + tie_a + tie_b;
    if (denominator === 0) return 0;

    const numerator = concordant + discordant;
    if (numerator === 0) return 0;

    return (concordant - discordant) / numerator;
}

/**
 * Computes the hamming distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The hamming distance between `a` and `b`.
 */
function hamming(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let disagree = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i];
        const y = b[i];
        disagree += x !== y ? 1 : 0;
    }
    return disagree / n;
}

/**
 * Computes the Haversine distance between two points on a sphere of unit length 1. Multiply the result with the radius of the sphere. (For instance Earth's radius is 6371km)
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - Point [lat1, lon1] in radians
 * @param {number[] | Float64Array} b - Point [lat2, lon2] in radians
 * @returns {number} The Haversine distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Haversine_formula}
 */
function haversine(a, b) {
    if (a.length !== 2 || b.length !== 2)
        throw new Error("Haversine distance requires exactly 2 coordinates [lat, lon] for each point!");
    const lat1 = a[0];
    const lon1 = a[1];
    const lat2 = b[0];
    const lon2 = b[1];

    const dlat = lat2 - lat1;
    const dlon = lon2 - lon1;

    const sin_dlat2 = Math.sin(dlat / 2);
    const sin_dlon2 = Math.sin(dlon / 2);

    const x = sin_dlat2 * sin_dlat2 + Math.cos(lat1) * Math.cos(lat2) * sin_dlon2 * sin_dlon2;
    const c = 2 * Math.atan2(Math.sqrt(x), Math.sqrt(1 - x));

    return c;
}

/**
 * Computes the jaccard distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The jaccard distance between `a` and `b`.
 */
function jaccard(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let num_non_zero = 0;
    let num_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] !== 0;
        const y = b[i] !== 0;
        num_non_zero += x || y ? 1 : 0;
        num_equal += x && y ? 1 : 0;
    }
    return (num_non_zero - num_equal) / num_non_zero;
}

/**
 * Computes the manhattan distance (`l_1`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The manhattan distance between `a` and `b`.
 */
function manhattan(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum;
}

/**
 * Computes the Sokal-Michener distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Sokal-Michener distance between `a` and `b`.

 */
function sokal_michener(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let num_not_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] !== 0;
        const y = b[i] !== 0;
        num_not_equal += x !== y ? 1 : 0;
    }
    return (2 * num_not_equal) / (n + num_not_equal);
}

/**
 * Computes the 1D Wasserstein distance (Earth Mover's Distance) between two distributions.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First distribution (histogram or probability mass)
 * @param {number[] | Float64Array} b - Second distribution (histogram or probability mass)
 * @returns {number} The Wasserstein/EMD distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Wasserstein_metric}
 */
function wasserstein(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sumA = 0;
    let sumB = 0;
    for (let i = 0; i < n; i++) {
        sumA += a[i];
        sumB += b[i];
    }

    // Fallback if sums are 0
    if (sumA === 0 && sumB === 0) return 0;
    if (sumA === 0 || sumB === 0) return Infinity;

    let distance = 0;
    let cumA = 0;
    let cumB = 0;
    for (let i = 0; i < n; i++) {
        cumA += a[i] / sumA;
        cumB += b[i] / sumB;
        distance += Math.abs(cumA - cumB);
    }
    return distance;
}

/**
 * Computes the yule distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The yule distance between `a` and `b`.
 */
function yule(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let num_true_true = 0;
    let num_true_false = 0;
    let num_false_true = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] !== 0;
        const y = b[i] !== 0;
        num_true_true += x && y ? 1 : 0;
        num_true_false += x && !y ? 1 : 0;
        num_false_true += !x && y ? 1 : 0;
    }
    const num_false_false = n - num_true_true - num_true_false - num_false_true;
    return num_true_false === 0 || num_false_true === 0
        ? 0
        : (2 * num_true_false * num_false_true) / (num_true_true * num_false_false + num_true_false * num_false_true);
}

/** @import { Metric } from "../metrics/index.js" */

/**
 * @param {Matrix | Float64Array[] | number[][]} A
 * @returns {A is Matrix}
 */
function isMatrix(A) {
    return A instanceof Matrix;
}

/**
 * Computes the distance matrix of datamatrix `A`.
 *
 * @category Matrix
 * @param {Matrix | Float64Array[] | number[][]} A - Matrix.
 * @param {Metric} [metric=euclidean] - The diistance metric. Default is `euclidean`
 * @returns {Matrix} The distance matrix of `A`.
 */
function distance_matrix(A, metric = euclidean) {
    /** @type {number} */
    const n = isMatrix(A) ? A.shape[0] : A.length;
    const D = new Matrix(n, n);
    for (let i = 0; i < n; ++i) {
        const A_i = isMatrix(A) ? A.row(i) : A[i];
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, isMatrix(A) ? A.row(j) : A[j]);
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}

//@ts-check


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
function k_nearest_neighbors(A, k, metric = euclidean) {
    A = A instanceof Matrix ? A : Matrix.from(A);
    const rows = A.shape[0];
    const D = metric === "precomputed" ? A : distance_matrix(A, metric);
    /** @type {{ i: number; j: number; distance: number }[][]} */
    const nN = [];
    for (let row = 0; row < rows; ++row) {
        const res = Array.from(D.row(row))
            .map((distance, col) => {
                return {
                    i: row,
                    j: col,
                    distance: distance,
                };
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1);
        nN.push(res);
    }
    return nN;
}

/**
 * Creates an Array containing `number` numbers from `start` to `end`. If `number = null`.
 *
 * @category Matrix
 * @param {number} start - Start value.
 * @param {number} end - End value.
 * @param {number} [number] - Number of number between `start` and `end`.
 * @returns {number[]} An array with `number` entries, beginning at `start` ending at `end`.
 */
function linspace(start, end, number) {
    if (number === undefined || number === null) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return number === 1 ? [start] : [];
    }
    const result = new Array(number);
    number -= 1;
    for (let i = number; i >= 0; --i) {
        result[i] = (i * end + (number - i) * start) / number;
    }
    return result;
}

/**
 * Computes the inner product between two arrays of the same length.
 *
 * @category Linear Algebra
 * @param {number[] | Float64Array} a - Array a.
 * @param {number[] | Float64Array} b - Array b.
 * @returns The inner product between `a` and `b`.
 */
function inner_product(a, b) {
    const N = a.length;
    if (N !== b.length) {
        throw new Error("Array a and b must have the same length!");
    }
    let sum = 0;
    for (let i = 0; i < N; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

/**
 * Numerical stable summation with the Kahan summation algorithm.
 *
 * @category Numerical
 * @param {number[] | Float64Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
 */
function kahan_sum(summands) {
    const n = summands.length;
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
 *
 * @category Numerical
 * @param {number[] | Float64Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
function neumair_sum(summands) {
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
 * Computes the QR Decomposition of the Matrix `A` using Gram-Schmidt process.
 *
 * @category Linear Algebra
 * @param {Matrix} A
 * @returns {{ R: Matrix; Q: Matrix }}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */
function qr(A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, cols, "identity");
    const R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        const v = A.col(j);
        for (let i = 0; i < j; ++i) {
            const q = Q.col(i);
            const q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            for (let k = 0; k < rows; ++k) {
                v[k] -= q_dot_v * q[k];
            }
            R.set_entry(i, j, q_dot_v);
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
 * Computes the QR Decomposition of the Matrix `A` with householder transformations.
 *
 * @category Linear Algebra
 * @param {Matrix} A
 * @returns {{ R: Matrix; Q: Matrix }}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */
function qr_householder(A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, rows, "I");
    const R = A.clone();

    for (let j = 0; j < cols; ++j) {
        const x = Matrix.from_vector(R.col(j).slice(j), "row");
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
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
function max(values) {
    let max = -Infinity;
    for (const value of values) {
        if (value !== null && max < value) {
            max = value;
        }
    }
    return max;
}

/**
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
function min(values) {
    let min = Infinity;
    for (const value of values) {
        if (value !== null && min > value) {
            min = value;
        }
    }
    return min;
}

/**
 * @category Utils
 * @class
 */
class Randomizer {
    _N = 624;
    _M = 397;
    _MATRIX_A = 0x9908b0df;
    _UPPER_MASK = 0x80000000;
    _LOWER_MASK = 0x7fffffff;

    /** @type {number[]} */
    _mt;
    /** @type {number} */
    _mti;
    /** @type {number} */
    _seed;

    /**
     * Mersenne Twister random number generator.
     *
     * @param {number} [_seed=new Date().getTime()] - The seed for the random number generator. If `_seed == null` then
     *   the actual time gets used as seed. Default is `new Date().getTime()`
     * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
     */
    constructor(_seed) {
        this._mt = new Array(this._N);
        this._mti = this._N + 1;
        this._seed = _seed ?? Date.now();
        this.seed = this._seed;
    }

    /** @type {number} seed */
    set seed(_seed) {
        this._seed = _seed;
        const mt = this._mt;

        mt[0] = _seed >>> 0;
        for (this._mti = 1; this._mti < this._N; this._mti += 1) {
            const mti = this._mti;
            const s = mt[mti - 1] ^ (mt[mti - 1] >>> 30);
            mt[mti] = ((((s & 0xffff0000) >>> 16) * 1812433253) << 16) + (s & 0x0000ffff) * 1812433253 + mti;
            mt[mti] >>>= 0;
        }
    }

    /**
     * Returns the seed of the random number generator.
     *
     * @returns {number} - The seed.
     */
    get seed() {
        return this._seed;
    }

    /**
     * Returns a float between 0 and 1.
     *
     * @returns {number} - A random number between [0, 1]
     */
    get random() {
        return this.random_int * (1.0 / 4294967296.0);
    }

    /**
     * Returns an integer between 0 and MAX_INTEGER.
     *
     * @returns {number} - A random integer.
     */
    get random_int() {
        let y,
            mag01 = [0x0, this._MATRIX_A];
        if (this._mti >= this._N) {
            let kk;

            /* if (this._mti == this._N + 1) {
                this.seed = 5489;
            } */

            const N_M = this._N - this._M;
            const M_N = this._M - this._N;

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
        this._mti += 1;
        y = this._mt[this._mti];
        y ^= y >>> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >>> 18;

        return y >>> 0;
    }

    gauss_random() {
        let x, y, r;
        if (this._val != null) {
            x = this._val;
            this._val = null;
            return x;
        } else
            do {
                x = 2 * this.random - 1;
                y = 2 * this.random - 1;
                r = x * x + y * y;
            } while (!r || r > 1);
        const c = Math.sqrt((-2 * Math.log(r)) / r);
        this._val = y * c; // cache this for next function call for efficiency
        return x * c;
    }

    /**
     * @template T Returns samples from an input Matrix or Array.
     * @param {T[]} A - The input Matrix or Array.
     * @param {number} n - The number of samples.
     * @returns {T[]} A random selection form `A` of `n` samples.
     */
    choice(A, n) {
        if (!Array.isArray(A)) throw new Error("A must be an Array!");
        // if (A instanceof Matrix) {
        //     let rows = A.shape[0];
        //     if (n > rows) {
        //         throw new Error("n bigger than A!");
        //     }
        //     /** @type {number[]} */
        //     let sample = new Array(n);
        //     let index_list = linspace(0, rows - 1);
        //     for (let i = 0, l = index_list.length; i < n; ++i, --l) {
        //         let random_index = this.random_int % l;
        //         sample[i] = index_list.splice(random_index, 1)[0];
        //     }
        //     return sample.map((d) => A.row(d));
        // } else if (Array.isArray(A) || A instanceof Float64Array) {
        const rows = A.length;
        if (n > rows) {
            throw new Error("n bigger than A!");
        }
        const sample = new Array(n);
        const index_list = linspace(0, rows - 1);
        for (let i = 0, l = index_list.length; i < n; ++i, --l) {
            const random_index = this.random_int % l;
            sample[i] = index_list.splice(random_index, 1)[0];
        }
        return sample.map((d) => A[d]);
        //} else {
        //throw new Error("A must be of type Matrix or Float64Array or number[]!");
        // }
    }

    /**
     * @template T Returns samples from an input Matrix or Array.
     * @param {T[]} A - The input Matrix or Array.
     * @param {number} n - The number of samples.
     * @param {number} seed - The seed for the random number generator.
     * @returns {T[]} - A random selection form `A` of `n` samples.
     */
    static choice(A, n, seed = 1212) {
        const R = new Randomizer(seed);
        return R.choice(A, n);
    }
}

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
function simultaneous_poweriteration(
    A,
    k = 2,
    { seed = 1212, max_iterations = 100, qr: qr$1 = qr, tol = 1e-8 } = {},
) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q, R } = qr$1(new Matrix(n, k, () => (randomizer.random - 0.5) * 2));
    while (max_iterations--) {
        const oldQ = Q;
        const Z = A.dot(Q);
        const QR = qr$1(Z);
        Q = QR.Q;
        R = QR.R;
        const error = euclidean_squared(Q.values, oldQ.values);
        if (error < tol) {
            break;
        }
    }

    const eigenvalues = R.diag();
    const eigenvectors = Q.transpose().to2dArray();
    return { eigenvalues, eigenvectors };
}

/** @typedef {(i: number, j: number) => number} Accessor */

/**
 * @class
 * @category Matrix
 */
class Matrix {
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
    constructor(rows, cols, value = 0) {
        /** @type {number} */ this._rows = rows;
        /** @type {number} */ this._cols = cols;
        /** @type {Float64Array} */ this._data;

        if (rows && cols) {
            if (!value) {
                this._data = new Float64Array(rows * cols);
            }
            if (typeof value === "function") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value(row, col);
                    }
                }
            }
            if (typeof value === "string") {
                if (value === "zeros") {
                    this._data = new Float64Array(rows * cols);
                    for (let row = 0; row < rows; ++row) {
                        for (let col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = 0;
                        }
                    }
                }
                if (value === "identity" || value === "I") {
                    this._data = new Float64Array(rows * cols);
                    for (let row = 0; row < rows; ++row) {
                        this._data[row * cols + row] = 1;
                    }
                }
                if (value === "center" && rows === cols) {
                    this._data = new Float64Array(rows * cols);
                    value = (i, j) => (i === j ? 1 : 0) - 1 / rows;
                    for (let row = 0; row < rows; ++row) {
                        for (let col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = value(row, col);
                        }
                    }
                }
            }
            if (typeof value === "number") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value;
                    }
                }
            }
            if (Array.isArray(value)) {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value[row][col];
                    }
                }
            }
        }
    }

    /**
     * Creates a Matrix out of `A`.
     * @param {Matrix | Float64Array[] | number[][]} A - The matrix, array, or number, which should converted to a Matrix.
     * @returns {Matrix}
     * @example
     * let A = Matrix.from([ [1, 0], [0, 1], ]); //creates a two by two identity matrix.
     */
    static from(A) {
        if (A instanceof Matrix) {
            return A.clone();
        }
        if (Matrix.is2dArray(A)) {
            const m = A.length;
            const n = A[0].length;
            for (let row = 0; row < m; ++row) {
                if (A[row].length !== n) {
                    throw new Error("various array lengths");
                }
            }
            return new Matrix(m, n, (i, j) => A[i][j]);
        }
        throw new Error("error");
    }

    /**
     * Creates a Matrix with the diagonal being the values of `v`.
     *
     * @example let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     *
     * @param {number[] | Float64Array} v
     * @returns {Matrix}
     */
    static from_diag(v) {
        const N = v.length;
        return new Matrix(N, N, (i, j) => (i === j ? v[i] : 0));
    }

    /**
     * Creates a Matrix with the diagonal being the values of `v`.
     *
     * @example let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     *
     * @param {number[] | Float64Array} v
     * @param {"col" | "row"} type
     * @returns {Matrix}
     */
    static from_vector(v, type) {
        const N = v.length;
        if (type === "col") {
            return new Matrix(N, 1, (i, _) => v[i]);
        } else {
            return new Matrix(1, N, (_, j) => v[j]);
        }
    }

    /**
     * Returns the `row`<sup>th</sup> row from the Matrix.
     *
     * @param {number} row
     * @returns {Float64Array}
     */
    row(row) {
        const data = this.values;
        const cols = this._cols;
        return data.subarray(row * cols, (row + 1) * cols);
    }

    /**
     * Returns an generator yielding each row of the Matrix.
     *
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
     * Makes a `Matrix` object an iterable object.
     *
     * @yields {Float64Array}
     */
    *[Symbol.iterator]() {
        for (const row of this.iterate_rows()) {
            yield row;
        }
    }

    /**
     * Sets the entries of `row`<sup>th</sup> row from the Matrix to the entries from `values`.
     *
     * @param {number} row
     * @param {number[]} values
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
            throw new Error("Values not valid! Needs to be either an Array, a Float64Array, or a fitting Matrix!");
        }
        return this;
    }

    /**
     * Swaps the rows `row1` and `row2` of the Matrix.
     *
     * @param {number} row1
     * @param {number} row2
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
        return this;
    }

    /**
     * Returns the col<sup>th</sup> column from the Matrix.
     *
     * @param {number} col
     * @returns {Float64Array}
     */
    col(col) {
        const result_col = new Float64Array(this._rows);
        for (let row = 0; row < this._rows; ++row) {
            result_col[row] = this.values[row * this._cols + col];
        }
        return result_col;
    }

    /**
     * Returns the `col`<sup>th</sup> entry from the `row`<sup>th</sup> row of the Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @returns {number}
     */
    entry(row, col) {
        return this.values[row * this._cols + col];
    }

    /**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given
     * {@link value}.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    set_entry(row, col, value) {
        this.values[row * this._cols + col] = value;
        return this;
    }

    /**
     * Adds a given {@link value} to the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the
     * Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    add_entry(row, col, value) {
        this.values[row * this._cols + col] += value;
        return this;
    }

    /**
     * Subtracts a given {@link value} from the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the
     * Matrix.
     *
     * @param {number} row
     * @param {number} col
     * @param {number} value
     * @returns {Matrix}
     */
    sub_entry(row, col, value) {
        this.values[row * this._cols + col] -= value;
        return this;
    }

    /**
     * Returns a new transposed Matrix.
     *
     * @returns {Matrix}
     */
    transpose() {
        const B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
        return B;
    }

    /**
     * Returns a new transposed Matrix. Short-form of `transpose`.
     *
     * @returns {Matrix}
     */
    get T() {
        return this.transpose();
    }

    /**
     * Returns the inverse of the Matrix.
     *
     * @returns {Matrix}
     */
    inverse() {
        const rows = this._rows;
        const cols = this._cols;
        const A = this.clone();
        const B = new Matrix(rows, cols, "I");

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
                throw new Error("Cannot compute inverse of Matrix, determinant is zero");
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
                            A_row[s] -= f * A_col[s];
                        }
                        for (let s = 0; s < cols; ++s) {
                            B_row[s] -= f * B_col[s];
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
     * Returns the dot product. If `B` is an Array or Float64Array then an Array gets returned. If `B` is a Matrix then
     * a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    dot(B) {
        if (B instanceof Matrix) {
            const [rows_A, cols_A] = this.shape;
            const [rows_B, cols_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${this.shape.join(" ⨯ ")}-Matrix, B is a ${B.shape.join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows.
                Must be equal!`);
            }
            const C = new Matrix(rows_A, cols_B, 0);
            const A_val = this.values;
            const B_val = B.values;
            const C_val = C.values;

            for (let i = 0; i < rows_A; ++i) {
                const i_cols_A = i * cols_A;
                const i_cols_B = i * cols_B;
                for (let k = 0; k < cols_A; ++k) {
                    const aik = A_val[i_cols_A + k];
                    if (aik === 0) continue;
                    const k_cols_B = k * cols_B;
                    for (let j = 0; j < cols_B; ++j) {
                        C_val[i_cols_B + j] += aik * B_val[k_cols_B + j];
                    }
                }
            }
            return C;
        } else if (Matrix.isArray(B)) {
            // TODO: create Matrix directly
            const rows = this._rows;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            const C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map((e) => e * B[row]));
            }
            return Matrix.from(C);
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Transposes the current matrix and returns the dot product with `B`. If `B` is an Array or Float64Array then an
     * Array gets returned. If `B` is a Matrix then a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    transDot(B) {
        if (B instanceof Matrix) {
            const [cols_A, rows_A] = this.shape; // transpose matrix
            const [rows_B, cols_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${[rows_A, cols_A].join(" ⨯ ")}-Matrix, B is a ${B.shape.join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows, which must be equal!`);
            }
            // let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
            // this.values[row * this._cols + col];
            const C = new Matrix(rows_A, cols_B, 0);
            const A_val = this.values; // A is rows_B x rows_A (transposed)
            const B_val = B.values;
            const C_val = C.values;

            for (let k = 0; k < cols_A; ++k) {
                // cols_A is rows_B
                const k_rows_A = k * rows_A;
                const k_cols_B = k * cols_B;
                for (let i = 0; i < rows_A; ++i) {
                    const aki = A_val[k_rows_A + i];
                    if (aki === 0) continue;
                    for (let j = 0; j < cols_B; ++j) {
                        C_val[i * cols_B + j] += aki * B_val[k_cols_B + j];
                    }
                }
            }
            return C;
        } else if (Matrix.isArray(B)) {
            // TODO: create Matrix directly
            const rows = this._cols;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            const C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.col(row).map((e) => e * B[row]));
            }
            return Matrix.from(C);
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Returns the dot product with the transposed version of `B`. If `B` is an Array or Float64Array then an Array gets
     * returned. If `B` is a Matrix then a Matrix gets returned.
     *
     * @param {Matrix | number[] | Float64Array} B The right side
     * @returns {Matrix}
     */
    dotTrans(B) {
        if (B instanceof Matrix) {
            const [rows_A, cols_A] = this.shape;
            const [cols_B, rows_B] = B.shape;
            if (cols_A !== rows_B) {
                throw new Error(`A.dot(B): A is a ${this.shape.join(" ⨯ ")}-Matrix, B is a ${[rows_B, cols_B].join(" ⨯ ")}-Matrix:
                A has ${cols_A} cols and B ${rows_B} rows, which must be equal!`);
            }
            const C = new Matrix(rows_A, cols_B, (row, col) => {
                const A_i = this.row(row);
                const B_i = B.row(col);
                let sum = 0;
                for (let i = 0; i < cols_A; ++i) {
                    sum += A_i[i] * B_i[i];
                }
                return sum;
            });
            return C;
        } else if (Matrix.isArray(B)) {
            // TODO: create Matrix directly
            const rows = this._rows;
            if (B.length !== rows) {
                throw new Error(`A.dot(B): A has ${rows} cols and B has ${B.length} rows. Must be equal!`);
            }
            const C = new Array(rows);
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map((e) => e * B[row]));
            }
            return Matrix.from(C);
        } else {
            throw new Error(`B must be Matrix or Array`);
        }
    }

    /**
     * Computes the outer product from `this` and `B`.
     *
     * @param {Matrix} B
     * @returns {Matrix}
     */
    outer(B) {
        const l = this._data.length;
        const r = B._data.length;
        if (l !== r) throw new Error("Matrix A and B needs to be of the same length!");
        const C = new Matrix(
            l,
            l,
            /** @type {Accessor} */ (i, j) => {
                if (i <= j) {
                    return this._data[i] * B._data[j];
                } else {
                    return this.entry(j, i);
                }
            },
        );

        return C;
    }

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
    concat(B, type = "horizontal") {
        const [rows_A, cols_A] = this.shape;
        const [rows_B, cols_B] = B.shape;
        if (type === "horizontal") {
            if (rows_A !== rows_B) {
                throw new Error(
                    `A.concat(B, "horizontal"): A and B need same number of rows, A has ${rows_A} rows, B has ${rows_B} rows.`,
                );
            }
            const X = new Matrix(rows_A, cols_A + cols_B, "zeros");
            X.set_block(0, 0, this);
            X.set_block(0, cols_A, B);
            return X;
        } else if (type === "vertical") {
            if (cols_A !== cols_B) {
                throw new Error(
                    `A.concat(B, "vertical"): A and B need same number of columns, A has ${cols_A} columns, B has ${cols_B} columns.`,
                );
            }
            const X = new Matrix(rows_A + rows_B, cols_A, "zeros");
            X.set_block(0, 0, this);
            X.set_block(rows_A, 0, B);
            return X;
        } else if (type === "diag") {
            const X = new Matrix(rows_A + rows_B, cols_A + cols_B, "zeros");
            X.set_block(0, 0, this);
            X.set_block(rows_A, cols_A, B);
            return X;
        } else {
            throw new Error(`type must be "horizontal" or "vertical", but type is ${type}!`);
        }
    }

    /**
     * Writes the entries of B in A at an offset position given by `offset_row` and `offset_col`.
     *
     * @param {number} offset_row
     * @param {number} offset_col
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
    get_block(start_row, start_col, end_row, end_col) {
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
    }

    /**
     * Returns a new array gathering entries defined by the indices given by argument.
     *
     * @param {number[]} row_indices - Array consists of indices of rows for gathering entries of this matrix
     * @param {number[]} col_indices - Array consists of indices of cols for gathering entries of this matrix
     * @returns {Matrix}
     */
    gather(row_indices, col_indices) {
        const N = row_indices.length;
        const D = col_indices.length;

        const R = new Matrix(N, D);
        for (let i = 0; i < N; ++i) {
            const row_index = row_indices[i];
            for (let j = 0; j < D; ++j) {
                const col_index = col_indices[j];
                R.set_entry(i, j, this.entry(row_index, col_index));
            }
        }

        return R;
    }

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

    /**
     * @param {number[] | Float64Array} values
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
    _apply_rowwise_array(values, f) {
        return this._apply_array(f, (_, j) => values[j]);
    }
    /**
     * @param {number[] | Float64Array} values
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
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

    /**
     * @param {Matrix | number[] | Float64Array | number} value
     * @param {(d: number, v: number) => number} f
     * @returns {Matrix}
     */
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
                    const v = values[row];
                    for (let col = 0; col < cols; ++col, ++i) {
                        data[i] = f(data[i], v);
                    }
                }
            } else if (rows === value_rows && cols === value_cols) {
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
        } else {
            // scalar value
            for (let i = 0, n = rows * cols; i < n; ++i) {
                data[i] = f(data[i], value);
            }
        }
        return this;
    }

    /**
     * Clones the Matrix.
     *
     * @returns {Matrix}
     */
    clone() {
        const B = new Matrix(this._rows, this._cols);
        //B._rows = this._rows;
        //B._cols = this._cols;
        if (this._data) {
            B._data = this._data.slice(0);
        }
        return B;
    }

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
    mult(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a * b);
    }

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
    divide(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a / b);
    }

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
    add(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a + b);
    }

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
    sub(value, { inline = false } = {}) {
        const A = inline ? this : this.clone();
        return A._apply(value, (a, b) => a - b);
    }

    /**
     * Returns the number of rows and columns of the Matrix.
     *
     * @returns {number[]} An Array in the form [rows, columns].
     */
    get shape() {
        return [this._rows, this._cols];
    }

    /**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     *
     * @param {[number, number, Accessor]} parameter - Takes an Array in the form [rows, cols, value], where rows and
     *   cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row
     *   and col) which has to return a value for the colth entry of the rowth row.
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
    }

    /**
     * Returns the Matrix as a Array of Float64Arrays.
     *
     * @returns {Float64Array[]}
     */
    to2dArray() {
        const result = [];
        for (const row of this.iterate_rows()) {
            result.push(row);
        }
        return result;
    }

    /**
     * Returns the Matrix as a Array of Arrays.
     *
     * @returns {number[][]}
     */
    asArray() {
        const result = [];
        for (const row of this.iterate_rows()) {
            result.push(Array.from(row));
        }
        return result;
    }

    /**
     * Returns the diagonal of the Matrix.
     *
     * @returns {Float64Array}
     */
    diag() {
        const rows = this._rows;
        const cols = this._cols;
        const min_row_col = Math.min(rows, cols);
        const result = new Float64Array(min_row_col);
        for (let i = 0; i < min_row_col; ++i) {
            result[i] = this.entry(i, i);
        }
        return result;
    }

    /**
     * Returns the mean of all entries of the Matrix.
     *
     * @returns {number}
     */
    mean() {
        const sum = this.sum();
        const n = this._rows * this._cols;
        return sum / n;
    }

    /**
     * Returns the sum oof all entries of the Matrix.
     *
     * @returns {number}
     */
    sum() {
        const data = this.values;
        return neumair_sum(data);
    }

    /**
     * Returns the entries of the Matrix.
     *
     * @returns {Float64Array}
     */
    get values() {
        const data = this._data;
        return data;
    }

    /**
     * Returns the mean of each row of the matrix.
     *
     * @returns {Float64Array}
     */
    meanRows() {
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

    /**
     * Returns the mean of each column of the matrix.
     *
     * @returns {Float64Array}
     */
    meanCols() {
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
     * Solves the equation `Ax = b` using the conjugate gradient method. Returns the result `x`.
     *
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer | null} [randomizer]
     * @param {number} [tol=1e-3] Default is `1e-3`
     * @returns {Matrix}
     */
    static solve_CG(A, b, randomizer, tol = 1e-3) {
        if (!randomizer) {
            randomizer = new Randomizer();
        }
        const rows = A.shape[0];
        const cols = b.shape[1];
        let result = new Matrix(rows, 0);
        for (let i = 0; i < cols; ++i) {
            const b_i = Matrix.from_vector(b.col(i), "col");
            let x = new Matrix(rows, 1, () => randomizer.random);
            let r = b_i.sub(A.dot(x));
            let d = r.clone();
            let iter = 0;
            const max_iter = rows * 10; // Prevent infinite loops
            do {
                const z = A.dot(d);
                const alpha = r.transDot(r).entry(0, 0) / d.transDot(z).entry(0, 0);
                x = x.add(d.mult(alpha));
                const r_next = r.sub(z.mult(alpha));
                const beta = r_next.transDot(r_next).entry(0, 0) / r.transDot(r).entry(0, 0);
                d = r_next.add(d.mult(beta));
                r = r_next;
                iter++;
            } while (Math.abs(r.mean()) > tol && iter < max_iter);
            result = result.concat(x, "horizontal");
        }
        return result;
    }

    /**
     * Solves the equation `Ax = b`. Returns the result `x`.
     *
     * @param {Matrix | { L: Matrix; U: Matrix }} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A, b) {
        const { L, U } = "L" in A && "U" in A ? A : Matrix.LU(A);
        const rows = L.shape[0];
        const x = b.clone();

        // forward
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < row; ++col) {
                x.sub_entry(0, row, L.entry(row, col) * x.entry(0, col));
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
     * `LU` decomposition of the Matrix `A`. Creates two matrices, so that the dot product `LU` equals `A`.
     *
     * @param {Matrix} A
     * @returns {{ L: Matrix; U: Matrix }} The left triangle matrix `L` and the upper triangle matrix `U`.
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
                    throw new Error("L's diagonal not supposed to be 0!");
                }
                let sum = 0;
                for (let k = 0; k < j; ++k) {
                    sum += L.entry(j, k) * U.entry(k, i);
                }
                U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
            }
        }

        return { L, U };
    }

    /**
     * Computes the determinante of `A`, by using the `LU` decomposition of `A`.
     *
     * @param {Matrix} A
     * @returns {number} The determinate of the Matrix `A`.
     */
    static det(A) {
        const [rows, cols] = A.shape;

        if (rows === 2 && cols === 2) {
            return A.entry(0, 0) * A.entry(1, 1) - A.entry(0, 1) * A.entry(1, 0);
        }
        if (rows === 3 && cols === 3) {
            const a = A.entry(0, 0);
            const b = A.entry(0, 1);
            const c = A.entry(0, 2);
            const d = A.entry(1, 0);
            const e = A.entry(1, 1);
            const f = A.entry(1, 2);
            const g = A.entry(2, 0);
            const h = A.entry(2, 1);
            const i = A.entry(2, 2);
            return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        }

        const { L, U } = Matrix.LU(A);
        const L_diag = L.diag();
        const U_diag = U.diag();
        let det = L_diag[0] * U_diag[0];
        for (let row = 1; row < rows; ++row) {
            det *= L_diag[row] * U_diag[row];
        }
        return det;
    }

    /**
     * Computes the `k` components of the SVD decomposition of the matrix `M`.
     *
     * @param {Matrix} M
     * @param {number} [k=2] Default is `2`
     * @returns {{ U: Float64Array[]; Sigma: Float64Array; V: Float64Array[] }}
     */
    static SVD(M, k = 2) {
        const MtM = M.transDot(M);
        const MMt = M.dotTrans(M);
        const { eigenvectors: V, eigenvalues: Sigma } = simultaneous_poweriteration(MtM, k);
        const { eigenvectors: U } = simultaneous_poweriteration(MMt, k);
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

    /**
     * @param {unknown} A
     * @returns {A is unknown[]|number[]|Float64Array|Float32Array}
     */
    static isArray(A) {
        return Array.isArray(A) || A instanceof Float64Array || A instanceof Float32Array;
    }

    /**
     * @param {any[]} A
     * @returns {A is number[][]|Float64Array[]}
     */
    static is2dArray(A) {
        if (!Array.isArray(A) || A.length === 0) {
            return false;
        }
        const n = A[0].length;
        for (let i = 0; i < A.length; ++i) {
            if (!Array.isArray(A[i]) && !(A[i] instanceof Float64Array)) {
                return false;
            }
            if (A[i].length !== n) {
                return false;
            }
        }
        return true;
    }
}

/** @import { Metric } from "../metrics/index.js" */

/**
 * Computes the norm of a vector, by computing its distance to **0**.
 *
 * @category Matrix
 * @param {Matrix | number[] | Float64Array} v - Vector.
 * @param {Metric} [metric=euclidean] - Which metric should be used to compute the norm. Default is `euclidean`
 * @returns {number} - The norm of `v`.
 */
function norm(v, metric = euclidean) {
    let vector = null;
    if (v instanceof Matrix) {
        const [rows, cols] = v.shape;
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

/** @import { Metric } from "../metrics/index.js" */

/**
 * Normalizes Vector `v`.
 *
 * @category Matrix
 * @param {number[] | Float64Array} v - Vector
 * @param {Metric} metric
 * @returns {number[] | Float64Array} - The normalized vector with length 1.
 */
function normalize(v, metric = euclidean) {
    const v_norm = norm(v, metric);
    return v.map((value) => value / v_norm);
}

/** @import {InputType} from "../index.js" */

/**
 * Base class for all clustering algorithms.
 * @template Para
 */
class Clustering {
    /** @type {InputType} */
    _points;
    /** @type {Para} */
    _parameters;
    /** @type {Matrix} */
    _matrix;
    /** @type {number} */
    _N;
    /** @type {number} */
    _D;

    /**
     * Compute the respective Clustering with given parameters
     * @param {InputType} points
     * @param {Para} parameters
     */
    constructor(points, parameters) {
        this._points = points;
        this._parameters = parameters;

        this._matrix = points instanceof Matrix ? points : Matrix.from(points);
        const [N, D] = this._matrix.shape;
        this._N = N;
        this._D = D;
    }

    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[][]} An array with the indices of the clusters.
     */
    get_clusters(...args) {
        throw new Error("The function get_clusters must be implemented!");
    }

    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[]} An array with the clusters id's for each point.
     */
    get_cluster_list(...args) {
        throw new Error("The function get_cluster_list must be implemented!");
    }
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
class CURE extends Clustering {
    /** @type {number} */
    _K;
    /** @type {number} */
    _num_representatives;
    /** @type {number} */
    _shrink_factor;
    /**
     * @private
     * @type {CURECluster[]}
     */
    _clusters = [];
    /** @type {number[]} */
    _cluster_ids = [];

    /**
     * @param {InputType} points
     * @param {Partial<ParametersCURE>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersCURE} */ (
                Object.assign(
                    { K: 2, num_representatives: 5, shrink_factor: 0.5, metric: euclidean, seed: 1212 },
                    parameters,
                )
            ),
        );

        this._K = this._parameters.K ?? 2;
        this._num_representatives = this._parameters.num_representatives ?? 5;
        this._shrink_factor = this._parameters.shrink_factor ?? 0.5;

        // Initialize clusters
        this._initialize_clusters();
        // Run CURE algorithm
        this._cure();
    }

    /**
     * Initialize each point as its own cluster
     * @private
     */
    _initialize_clusters() {
        const N = this._N;
        //const D = this._D;
        this._clusters = [];

        for (let i = 0; i < N; ++i) {
            const point = this._matrix.row(i);
            const centroid = new Float64Array(point);
            // For single point, representative is the point itself
            const representatives = [new Float64Array(point)];

            this._clusters.push(new CURECluster([i], centroid, representatives));
        }
    }

    /**
     * Compute distance between two clusters using representative points
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {number}
     */
    _cluster_distance(cluster1, cluster2) {
        const reps1 = cluster1.representatives;
        const reps2 = cluster2.representatives;
        const metric = this._parameters.metric;

        let min_dist = Infinity;
        for (const r1 of reps1) {
            for (const r2 of reps2) {
                const dist = metric(r1, r2);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
        }
        return min_dist;
    }

    /**
     * Find the closest pair of clusters
     * @private
     * @returns {[number, number, number]} [index1, index2, distance]
     */
    _find_closest_clusters() {
        let min_dist = Infinity;
        let min_i = 0;
        let min_j = 1;

        for (let i = 0; i < this._clusters.length; ++i) {
            for (let j = i + 1; j < this._clusters.length; ++j) {
                const dist = this._cluster_distance(this._clusters[i], this._clusters[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        return [min_i, min_j, min_dist];
    }

    /**
     * Merge two clusters
     * @private
     * @param {CURECluster} cluster1
     * @param {CURECluster} cluster2
     * @returns {CURECluster}
     */
    _merge_clusters(cluster1, cluster2) {
        // Merge indices
        const merged_indices = [...cluster1.indices, ...cluster2.indices];

        // Calculate new centroid
        const size1 = cluster1.indices.length;
        const size2 = cluster2.indices.length;
        const total_size = size1 + size2;
        const D = this._D;
        const new_centroid = new Float64Array(D);

        for (let d = 0; d < D; ++d) {
            new_centroid[d] = (size1 * cluster1.centroid[d] + size2 * cluster2.centroid[d]) / total_size;
        }

        // Collect all points from both clusters
        /** @type {{index: number, point: Float64Array}[]} */
        const all_points = [];
        for (const idx of cluster1.indices) {
            all_points.push({ index: idx, point: this._matrix.row(idx) });
        }
        for (const idx of cluster2.indices) {
            all_points.push({ index: idx, point: this._matrix.row(idx) });
        }

        // Select representative points - pick points farthest from centroid
        const num_reps = Math.min(this._num_representatives, all_points.length);
        const metric = this._parameters.metric;

        // Calculate distances from centroid for all points
        const distances = all_points.map(({ point }) => metric(point, new_centroid));

        // Select num_reps points with maximum distance (farthest from centroid)
        const selected_indices = [];
        const used = new Set();

        for (let r = 0; r < num_reps; ++r) {
            let max_dist = -1;
            let max_idx = -1;

            for (let i = 0; i < distances.length; ++i) {
                if (!used.has(i) && distances[i] > max_dist) {
                    max_dist = distances[i];
                    max_idx = i;
                }
            }

            if (max_idx >= 0) {
                used.add(max_idx);
                selected_indices.push(max_idx);
            }
        }

        // Shrink representative points toward centroid
        const new_representatives = selected_indices.map((idx) => {
            const point = all_points[idx].point;
            const shrunk = new Float64Array(D);
            const alpha = this._shrink_factor;

            for (let d = 0; d < D; ++d) {
                shrunk[d] = point[d] + alpha * (new_centroid[d] - point[d]);
            }

            return shrunk;
        });

        return new CURECluster(merged_indices, new_centroid, new_representatives);
    }

    /**
     * Run CURE clustering algorithm
     * @private
     */
    _cure() {
        // Merge clusters until we have K clusters
        while (this._clusters.length > this._K) {
            const [i, j] = this._find_closest_clusters();

            // Merge clusters i and j
            const merged = this._merge_clusters(this._clusters[i], this._clusters[j]);

            // Remove the old clusters and add the merged one
            // Remove larger index first to maintain correct indices
            // min_i < min_j is always true from _find_closest_clusters
            this._clusters.splice(j, 1);
            this._clusters.splice(i, 1);

            this._clusters.push(merged);
        }

        // Build cluster list for get_cluster_list
        this._build_cluster_ids();
    }

    /**
     * Build the cluster list (point -> cluster assignment)
     * @private
     */
    _build_cluster_ids() {
        const N = this._N;
        this._cluster_ids = new Array(N).fill(-1);

        for (let c = 0; c < this._clusters.length; ++c) {
            for (const idx of this._clusters[c].indices) {
                this._cluster_ids[idx] = c;
            }
        }
    }

    /**
     * @returns {number[][]}
     */
    get_clusters() {
        return this._clusters.map((cluster) => cluster.indices);
    }

    /**
     * @returns {number[]}
     */
    get_cluster_list() {
        return this._cluster_ids;
    }
}

/**
 * @private
 * Represents a cluster in CURE algorithm
 */
class CURECluster {
    /**
     * @param {number[]} indices - Indices of points in the cluster
     * @param {Float64Array} centroid - Centroid of the cluster
     * @param {Float64Array[]} representatives - Representative points (shrunk toward centroid)
     */
    constructor(indices, centroid, representatives) {
        /** @type {number[]} */
        this.indices = indices;
        /** @type {Float64Array} */
        this.centroid = centroid;
        /** @type {Float64Array[]} */
        this.representatives = representatives;
    }
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
class HierarchicalClustering extends Clustering {
    /** @type {Cluster | null} */
    root = null;

    /**
     * @param {InputType} points - Data or distance matrix if metric is 'precomputed'
     * @param {Partial<ParametersHierarchicalClustering>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersHierarchicalClustering} */ (
                Object.assign({ linkage: "complete", metric: euclidean }, parameters)
            ),
        );
        this._id = 0;
        if (this._parameters.metric === "precomputed" && this._matrix.shape[0] !== this._matrix.shape[1]) {
            throw new Error("If metric is 'precomputed', then matrix has to be square!");
        }

        const metric = this._parameters.metric;
        const A = this._matrix;
        const N = this._N;
        this._d_min = new Float64Array(N);
        const d_min = this._d_min;
        let distance_matrix;
        if (metric !== "precomputed") {
            distance_matrix = new Matrix(N, N, Infinity);
            for (let i = 0; i < N; ++i) {
                distance_matrix.set_entry(i, i, 0);
                d_min[i] = i; // temporary
                const Ai = A.row(i);
                for (let j = i + 1; j < N; ++j) {
                    const dist = metric(Ai, A.row(j));
                    distance_matrix.set_entry(i, j, dist);
                    distance_matrix.set_entry(j, i, dist);
                }
            }
            for (let i = 0; i < N; i++) {
                let min_j = 0;
                let min_d = Infinity;
                for (let j = 0; j < N; j++) {
                    if (i === j) continue;
                    const d = distance_matrix.entry(i, j);
                    if (d < min_d) {
                        min_d = d;
                        min_j = j;
                    }
                }
                d_min[i] = min_j;
            }
        } else {
            distance_matrix = this._matrix.clone();
            for (let i = 0; i < N; ++i) {
                distance_matrix.set_entry(i, i, 0);
                d_min[i] = i === 0 ? 1 : 0;
                for (let j = 0; j < N; ++j) {
                    if (i === j) continue;
                    if (distance_matrix.entry(i, d_min[i]) > distance_matrix.entry(i, j)) {
                        d_min[i] = j;
                    }
                }
            }
        }
        this._distance_matrix = distance_matrix;
        this._clusters = new Array(N);
        const clusters = this._clusters;
        this._c_size = new Uint16Array(N);
        const c_size = this._c_size;
        for (let i = 0; i < N; ++i) {
            clusters[i] = [];
            clusters[i][0] = new Cluster(this._id++, null, null, 0, A.row(i), i, 1, 0);
            c_size[i] = 1;
        }
        const D = this._distance_matrix;
        const linkage = this._parameters.linkage;
        const p_max = N - 1;
        for (let p = 0; p < p_max; ++p) {
            let c1 = -1;
            let min_dist = Infinity;
            for (let i = 0; i < N; ++i) {
                if (D.entry(i, i) === Infinity) continue;
                const dist = D.entry(i, d_min[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    c1 = i;
                }
            }
            if (c1 === -1) break;

            const c2 = d_min[c1];
            const c1_cluster = clusters[c1][0];
            const c2_cluster = clusters[c2][0];
            const c1_cluster_indices = c1_cluster.isLeaf ? [c1_cluster.index] : c1_cluster.index;
            const c2_cluster_indices = c2_cluster.isLeaf ? [c2_cluster.index] : c2_cluster.index;
            const indices = c1_cluster_indices.concat(c2_cluster_indices);
            const new_cluster = new Cluster(this._id++, c1_cluster, c2_cluster, D.entry(c1, c2), null, indices);
            c1_cluster.parent = new_cluster;
            c2_cluster.parent = new_cluster;
            clusters[c1].unshift(new_cluster);

            const size1 = c_size[c1];
            const size2 = c_size[c2];
            c_size[c1] += size2;

            for (let j = 0; j < N; ++j) {
                if (j === c1 || j === c2 || D.entry(j, j) === Infinity) continue;
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
                        value = (size1 * D_c1_j + size2 * D_c2_j) / (size1 + size2);
                        break;
                }
                D.set_entry(j, c1, value);
                D.set_entry(c1, j, value);
            }

            D.set_entry(c2, c2, Infinity);
            for (let i = 0; i < N; ++i) {
                D.set_entry(i, c2, Infinity);
                D.set_entry(c2, i, Infinity);
            }

            // Update d_min for all rows
            for (let i = 0; i < N; i++) {
                if (D.entry(i, i) === Infinity) continue;
                if (d_min[i] === c1 || d_min[i] === c2 || i === c1) {
                    let min_j = 0;
                    let min_d = Infinity;
                    for (let j = 0; j < N; j++) {
                        if (i === j || D.entry(j, j) === Infinity) continue;
                        const d = D.entry(i, j);
                        if (d < min_d) {
                            min_d = d;
                            min_j = j;
                        }
                    }
                    d_min[i] = min_j;
                } else {
                    if (D.entry(i, c1) < D.entry(i, d_min[i])) {
                        d_min[i] = c1;
                    }
                }
            }

            this.root = new_cluster;
        }
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {Cluster[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters_raw(value, type = "distance") {
        /** @type {Cluster[][]} */
        const clusters = [];
        /** @type {(d: {dist: number, depth: number}) => number} */
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
        this._traverse(/** @type {Cluster} */ (this.root), accessor, value, clusters);
        return clusters;
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[][]} - Array of clusters with the indices of the rows in given points.
     */
    get_clusters(value, type = "distance") {
        /** @type {Cluster[][]} */
        const clusters = [];
        /** @type {(d: {dist: number, depth: number}) => number} */
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
        if (this.root) this._traverse(this.root, accessor, value, clusters);
        return clusters.map((cluster) => cluster.map((d) => d.index));
    }

    /**
     * @param {number} value - Value where to cut the tree.
     * @param {"distance" | "depth"} [type="distance"] - Type of value. Default is `"distance"`
     * @returns {number[]} - Array of clusters with the indices of the rows in given points.
     */
    get_cluster_list(value, type = "distance") {
        const clusters = this.get_clusters(value, type);
        /** @type {number[]} */
        const list = new Array(this._N).fill(0);
        for (let i = 0; i < clusters.length; ++i) {
            const cluster = clusters[i];
            for (let j = 0; j < cluster.length; ++j) {
                const index = cluster[j];
                list[index] = i;
            }
        }
        return list;
    }

    /**
     * @private
     * @param {Cluster} node
     * @param {(d: {dist: number, depth: number}) => number} f
     * @param {number} value
     * @param {Cluster[][]} result
     */
    _traverse(node, f, value, result) {
        if (f(node) <= value) {
            result.push(node.leaves());
        } else {
            if (node.left) this._traverse(node.left, f, value, result);
            if (node.right) this._traverse(node.right, f, value, result);
        }
    }
}

/** @private */
class Cluster {
    /**@type {number} */
    size;
    /**@type {number} */
    depth;
    /**@type {Cluster | null} */
    parent;

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
    constructor(id, left, right, dist, centroid, index, size, depth) {
        this.id = id;
        this.left = left;
        this.right = right;
        this.dist = dist;
        this.index = index;
        if (size) {
            this.size = size;
        } else {
            if (!left || !right) throw new Error("If size is not given, left & right cannot be null!");
            this.size = left.size + right.size;
        }

        if (depth !== undefined) {
            this.depth = depth;
        } else {
            if (!left || !right) throw new Error("If depth is not given, left & right cannot be null!");
            this.depth = Math.max(left.depth, right.depth) + 1;
        }

        if (centroid !== undefined && centroid !== null) {
            this.centroid = centroid;
        } else {
            if (!left || !right) throw new Error("If centroid is not given, left & right cannot be null!");

            this.centroid = this._calculate_centroid(left, right);
        }

        this.parent = null;
    }

    /**
     *
     * @param {Cluster} left
     * @param {Cluster} right
     * @returns {Float64Array}
     */
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

    /**
     *
     * @returns {Cluster[]}
     */
    leaves() {
        if (this.isLeaf) return [this];
        const left = this.left;
        const right = this.right;
        return (left ? (left.isLeaf ? [left] : left.leaves()) : []).concat(
            right ? (right.isLeaf ? [right] : right.leaves()) : [],
        );
    }

    /**
     *
     * @returns {Cluster[]}
     */
    descendants() {
        if (this.isLeaf) return [this];
        const left_descendants = this.left ? this.left.descendants() : [];
        const right_descendants = this.right ? this.right.descendants() : [];
        return left_descendants.concat(right_descendants).concat([this]);
    }
}

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
class DisjointSet {
    /**
     * @param {T[]?} elements
     */
    constructor(elements = null) {
        /**
         * @private
         * @type {Map<T, DisjointSetPayload<T>>}
         */
        this._list = new Map();
        if (elements) {
            for (const e of elements) {
                this.make_set(e);
            }
        }
    }

    /**
     * @private
     * @param {T} x
     * @returns {DisjointSet<T>}
     */
    make_set(x) {
        const list = this._list;
        if (!list.has(x)) {
            list.set(x, { parent: x, children: new Set([x]), size: 1 });
        }
        return this;
    }

    /**
     * @param {T} x
     * @returns
     */
    find(x) {
        const list = this._list;
        const disjoint_set = list.get(x);
        if (disjoint_set) {
            if (disjoint_set.parent !== x) {
                disjoint_set.children.add(x);
                const new_parent = this.find(disjoint_set.parent);
                if (!new_parent) throw new Error("should not happen!");
                disjoint_set.parent = new_parent;
                return disjoint_set.parent;
            } else {
                return x;
            }
        } else {
            return null;
        }
    }

    /**
     * @param {T} x
     * @param {T} y
     * @returns
     */
    union(x, y) {
        let node_x = this.find(x);
        let node_y = this.find(y);

        if (!node_x || !node_y) throw new Error("x or y not found!");

        let disjoint_set_x = this._list.get(node_x);
        let disjoint_set_y = this._list.get(node_y);

        if (!disjoint_set_x || !disjoint_set_y) throw new Error("should not happen!");

        if (node_x === node_y) return this;
        if (disjoint_set_x.size < disjoint_set_y.size) {
            [node_x, node_y] = [node_y, node_x];
            [disjoint_set_x, disjoint_set_y] = [disjoint_set_y, disjoint_set_x];
        }

        disjoint_set_y.parent = node_x;
        // keep track of children
        disjoint_set_y.children.forEach(disjoint_set_x.children.add, disjoint_set_x.children);
        disjoint_set_x.size += disjoint_set_y.size;

        return this;
    }

    /** @param {T} x */
    get_children(x) {
        const node = this._list.get(x);
        if (node) {
            return node.children;
        } else {
            return null;
        }
    }
}

/** @import { Comparator } from "./index.js" */

/**
 * @template T
 * @class
 * @category Data Structures
 */
class Heap {
    /** @type {{ element: T; value: number }[]} */
    _container;

    /** @type {Comparator} */
    _comparator;

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
    constructor(elements = null, accessor, comparator = "min") {
        /** @type {(d: T) => number} */
        this._accessor = accessor;
        this._container = [];
        if (comparator === "min") {
            this._comparator = (a, b) => a < b;
        } else if (comparator === "max") {
            this._comparator = (a, b) => a > b;
        } else {
            this._comparator = comparator;
        }
        if (elements) {
            this._container = [];
            for (const e of elements) {
                this._container.push({
                    element: e,
                    value: accessor(e),
                });
            }
            for (let i = Math.floor(elements.length / 2 - 1); i >= 0; --i) {
                this._heapify_down(i);
            }
        }
    }

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
    static heapify(elements, accessor, comparator = "min") {
        const heap = new Heap(null, accessor, comparator);
        const container = heap._container;
        for (const e of elements) {
            container.push({
                element: e,
                value: accessor(e),
            });
        }
        for (let i = Math.floor(elements.length / 2 - 1); i >= 0; --i) {
            heap._heapify_down(i);
        }
        return heap;
    }

    /**
     * Swaps elements of container array.
     *
     * @private
     * @param {number} index_a
     * @param {number} index_b
     */
    _swap(index_a, index_b) {
        const container = this._container;
        [container[index_b], container[index_a]] = [container[index_a], container[index_b]];
        return;
    }

    /** @private */
    _heapify_up() {
        const container = this._container;
        let index = container.length - 1;
        while (index > 0) {
            const parentIndex = Math.floor((index - 1) / 2);
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
     *
     * @param {T} element
     * @returns {Heap<T>}
     */
    push(element) {
        const value = this._accessor(element);
        //const node = new Node(element, value);
        const node = { element: element, value: value };
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    /**
     * @private
     * @param {Number} [start_index=0] Default is `0`
     */
    _heapify_down(start_index = 0) {
        const container = this._container;
        const comparator = this._comparator;
        const length = container.length;
        const left = 2 * start_index + 1;
        const right = 2 * start_index + 2;
        let index = start_index;
        if (index >= length) throw "index higher than length";
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
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`}).
     */
    pop() {
        const container = this._container;
        if (container.length === 0) {
            return null;
        } else if (container.length === 1) {
            const item = container.pop();
            if (!item) throw new Error("Cannot happen!");
            return item;
        }
        this._swap(0, container.length - 1);
        const item = container.pop();
        this._heapify_down();
        return item ?? null;
    }

    /**
     * Returns the top entry of the heap without removing it.
     *
     * @returns {{ element: T; value: number } | null} Object consists of the element and its value (computed by
     *   `accessor`).
     */
    get first() {
        return this._container.length > 0 ? this._container[0] : null;
    }

    /**
     * Yields the raw data
     *
     * @yields {T} Object consists of the element and its value (computed by `accessor`}).
     */
    *iterate() {
        for (let i = 0, n = this._container.length; i < n; ++i) {
            yield this._container[i].element;
        }
    }

    /**
     * Returns the heap as ordered array.
     *
     * @returns {T[]} Array consisting the elements ordered by `comparator`.
     */
    toArray() {
        return this._container.sort((a, b) => (this._comparator(a.value, b.value) ? -1 : 1)).map((d) => d.element);
    }

    /**
     * Returns elements of container array.
     *
     * @returns {T[]} Array consisting the elements.
     */
    data() {
        return this._container.map((d) => d.element);
    }

    /**
     * Returns the container array.
     *
     * @returns {{ element: T; value: number }[]} The container array.
     */
    raw_data() {
        return this._container;
    }

    /**
     * The size of the heap.
     *
     * @returns {number}
     */
    get length() {
        return this._container.length;
    }

    /**
     * Returns false if the the heap has entries, true if the heap has no entries.
     *
     * @returns {boolean}
     */
    get empty() {
        return this.length === 0;
    }
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
class KMeans extends Clustering {
    /**
     * @param {InputType} points
     * @param {Partial<ParametersKMeans>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersKMeans} */ (Object.assign({ K: 4, metric: euclidean, seed: 1212 }, parameters)),
        );

        const K = this._parameters.K;
        const seed = parameters.seed;

        // Convert points to Matrix if needed
        if (points instanceof Matrix) {
            this._matrix = points;
        } else {
            this._matrix = Matrix.from(points);
        }

        const [N, D] = this._matrix.shape;
        this._N = N;
        this._D = D;

        this._K = K > N ? N : K;
        this._randomizer = new Randomizer(seed);

        /** @type {number[]} */
        this._clusters = new Array(N).fill(0);

        this._cluster_centroids = parameters.initial_centroids
            ? parameters.initial_centroids.map((c) => new Float64Array(c))
            : this._get_random_centroids(this._K);
        let cluster_centroids = this._cluster_centroids;
        let iterations = 0;
        const max_iterations = 300;
        let clusters_changed = true;

        while (clusters_changed && iterations < max_iterations) {
            const iteration_result = this._iteration(cluster_centroids);
            cluster_centroids = iteration_result.cluster_centroids;
            clusters_changed = iteration_result.clusters_changed;
            iterations++;
        }

        this._cluster_centroids = cluster_centroids;
    }

    /** @returns {number} The number of clusters */
    get k() {
        return this._K;
    }

    /** @returns {Float64Array[]} The cluster centroids */
    get centroids() {
        return this._cluster_centroids;
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        return this._clusters;
    }

    /** @returns {number[][]} An Array of clusters with the indices of the points. */
    get_clusters() {
        const K = this._K;
        const clusters = this._clusters;
        /** @type {number[][]} */
        const result = new Array(K).fill(0).map(() => []);
        clusters.forEach((c, i) => {
            if (c >= 0 && c < K) {
                result[c].push(i);
            }
        });
        return result;
    }

    /**
     * @private
     * @param {number[]} point_indices
     * @param {number[]} candidates
     * @returns {number}
     */
    _furthest_point(point_indices, candidates) {
        const A = this._matrix;
        const metric = this._parameters.metric;

        if (point_indices.length === 0 || candidates.length === 0) {
            return candidates[0] ?? 0;
        }

        const H = Heap.heapify(
            candidates,
            (d) => {
                const Ad = A.row(d);
                let sum = 0;
                for (let j = 0; j < point_indices.length; ++j) {
                    sum += metric(Ad, A.row(point_indices[j]));
                }
                return sum;
            },
            "max",
        );

        const furthest = H.pop();
        if (!furthest) throw new Error("Should not happen!");

        return furthest.element;
    }

    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    _get_random_centroids(K) {
        const N = this._N;
        const randomizer = this._randomizer;
        const A = this._matrix;
        /** @type {Float64Array[]} */
        const cluster_centroids = new Array(K);
        const indices = linspace(0, N - 1);

        // First centroid: random selection
        const random_point = randomizer.random_int % N;
        cluster_centroids[0] = A.row(random_point);
        const init_points = [random_point];

        const sample_size = Math.max(1, Math.floor((N - K) / K));

        for (let i = 1; i < K; ++i) {
            const remaining = indices.filter((d) => !init_points.includes(d));
            if (remaining.length === 0) break;

            const sample = randomizer.choice(remaining, Math.min(sample_size, remaining.length));
            const furthest_point = this._furthest_point(init_points, sample);

            init_points.push(furthest_point);
            cluster_centroids[i] = A.row(furthest_point);
        }

        return cluster_centroids;
    }

    /**
     * @private
     * @param {Float64Array[]} cluster_centroids
     * @returns {{ clusters_changed: boolean; cluster_centroids: Float64Array[] }}
     */
    _iteration(cluster_centroids) {
        const K = cluster_centroids.length;
        const N = this._N;
        const metric = this._parameters.metric;
        const A = this._matrix;
        const clusters = this._clusters;
        let clusters_changed = false;

        // Find nearest cluster centroid for each point
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            let min_dist = Infinity;
            let min_cluster = 0;

            for (let j = 0; j < K; ++j) {
                const d = metric(cluster_centroids[j], Ai);
                if (d < min_dist) {
                    min_dist = d;
                    min_cluster = j;
                }
            }

            if (clusters[i] !== min_cluster) {
                clusters_changed = true;
                clusters[i] = min_cluster;
            }
        }

        // Update cluster centroids
        const new_centroids = this._compute_centroid(K);

        return {
            clusters_changed: clusters_changed,
            cluster_centroids: new_centroids,
        };
    }

    /**
     * @private
     * @param {number} K
     * @returns {Float64Array[]}
     */
    _compute_centroid(K) {
        const N = this._N;
        const D = this._D;
        const A = this._matrix;
        const clusters = this._clusters;

        // Initialize new centroids and counters
        /** @type {Float64Array[]} */
        const new_centroids = new Array(K);
        const cluster_counter = new Array(K).fill(0);

        for (let i = 0; i < K; ++i) {
            new_centroids[i] = new Float64Array(D);
        }

        // Sum up all points in each cluster
        for (let i = 0; i < N; ++i) {
            const Ai = A.row(i);
            const ci = clusters[i];
            if (ci >= 0 && ci < K) {
                cluster_counter[ci]++;
                const centroid = new_centroids[ci];
                for (let j = 0; j < D; ++j) {
                    centroid[j] += Ai[j];
                }
            }
        }

        // Divide by count to get mean
        for (let i = 0; i < K; ++i) {
            const n = cluster_counter[i];
            if (n > 0) {
                const centroid = new_centroids[i];
                for (let j = 0; j < D; ++j) {
                    centroid[j] /= n;
                }
            }
        }

        return new_centroids;
    }
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
class KMedoids extends Clustering {
    /**
     * @param {InputType} points - Data matrix
     * @param {Partial<ParametersKMedoids>} parameters
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(points, parameters = {}) {
        super(points, Object.assign({ K: 4, max_iter: null, metric: euclidean, seed: 1212 }, parameters));
        this._A = this._matrix.to2dArray();
        let K = this._parameters.K;
        const N = this._N;
        this._max_iter = this._parameters.max_iter ?? 10 * Math.log10(N);
        this._distance_matrix = new Matrix(N, N, "zeros");

        if (K > N) {
            this._parameters.K = K = N;
        }
        this._randomizer = new Randomizer(this._parameters.seed);
        this._clusters = new Array(N).fill(-1);
        this._cluster_medoids = this._get_random_medoids(K);
        this._is_initialized = false;
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        if (!this._is_initialized) {
            this.get_clusters();
        }
        return this._clusters;
    }

    /** @returns {number[][]} - Array of clusters with the indices of the rows in given points. */
    get_clusters() {
        const K = this._parameters.K;
        const A = this._A;
        const N = this._N;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        /** @type {number[][]} */
        const result = new Array(K).fill(0).map(() => []);
        for (let j = 0; j < N; j++) {
            const nearest = this._nearest_medoid(A[j], j);
            const cluster_idx = nearest.index_nearest;
            result[cluster_idx].push(j);
            this._clusters[j] = cluster_idx;
        }
        return result;
    }

    /** @returns {number} */
    get k() {
        return this._parameters.K;
    }

    /** @returns {number[]} */
    get medoids() {
        return this.get_medoids();
    }

    /** @returns {number[]} */
    get_medoids() {
        const K = this._parameters.K;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        return this._cluster_medoids;
    }

    async *generator() {
        const max_iter = this._max_iter;
        if (!this._is_initialized) {
            this.get_clusters();
        }
        yield this.get_clusters();
        let i = 0;
        while (i < max_iter) {
            const finish = this._iteration();
            this._update_clusters();
            yield this.get_clusters();
            if (finish) break;
            i++;
        }
    }

    /** Algorithm 1. FastPAM1: Improved SWAP algorithm */
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

    /**
     * FastPAM1: One best swap per iteration
     * @private
     * @returns {boolean}
     */
    _iteration() {
        const A = this._A;
        const K = this._parameters.K;
        const medoids = this._cluster_medoids;
        const N = this._N;

        // Precompute nearest and second nearest medoid for all points
        const cache = new Array(N);
        for (let i = 0; i < N; i++) {
            cache[i] = this._nearest_medoid(A[i], i);
        }

        let best_delta = 0;
        let best_swap = null; // { m_idx: index in medoids, x_idx: index in A }

        // For each non-medoid point j, evaluate swapping it with each medoid i
        const medoid_set = new Set(medoids);
        for (let j = 0; j < N; j++) {
            if (medoid_set.has(j)) continue;

            const x_j = A[j];
            const d_j = cache[j].distance_nearest;

            // deltaTD[i] will store the change in total distance if we swap medoid[i] with j
            const deltaTD = new Array(K).fill(-d_j);

            for (let o = 0; o < N; o++) {
                if (o === j) continue;
                const dist_o_j = this._get_distance(o, j, A[o], x_j);
                const { index_nearest: n, distance_nearest: d_n, distance_second: d_s } = cache[o];

                // If o is assigned to the current medoid being swapped out (n)
                deltaTD[n] += Math.min(dist_o_j, d_s) - d_n;

                // For all other medoids i != n, if j is closer to o than its current medoid
                if (dist_o_j < d_n) {
                    for (let i = 0; i < K; i++) {
                        if (i !== n) {
                            deltaTD[i] += dist_o_j - d_n;
                        }
                    }
                }
            }

            // Find best medoid to swap with j
            for (let i = 0; i < K; i++) {
                if (deltaTD[i] < best_delta) {
                    best_delta = deltaTD[i];
                    best_swap = { m_idx: i, x_idx: j };
                }
            }
        }

        if (best_swap && best_delta < 0) {
            medoids[best_swap.m_idx] = best_swap.x_idx;
            this._cluster_medoids = medoids;
            return false; // not finished
        }

        return true; // finished
    }

    /**
     * @private
     * Get distance between two points
     * @param {number} i
     * @param {number} j
     * @param {Float64Array?} x_i
     * @param {Float64Array?} x_j
     * @returns {number}
     */
    _get_distance(i, j, x_i = null, x_j = null) {
        if (i === j) return 0;
        const D = this._distance_matrix;
        const A = this._A;
        const metric = this._parameters.metric;
        let d_ij = D.entry(i, j);
        if (d_ij === 0) {
            d_ij = metric(x_i || A[i], x_j || A[j]);
            D.set_entry(i, j, d_ij);
            D.set_entry(j, i, d_ij);
        }
        return d_ij;
    }

    /**
     * @private
     * @param {Float64Array} x_j
     * @param {number} j
     * @returns
     */
    _nearest_medoid(x_j, j) {
        const medoids = this._cluster_medoids;
        const A = this._A;
        if (medoids.length === 0) {
            throw new Error("No medoids available. Initialization failed.");
        }

        let d_n = Infinity;
        let n = -1;
        let d_s = Infinity;
        let s = -1;

        for (let i = 0; i < medoids.length; i++) {
            const m = medoids[i];
            const d = this._get_distance(j, m, x_j, A[m]);
            if (d < d_n) {
                d_s = d_n;
                s = n;
                d_n = d;
                n = i;
            } else if (d < d_s) {
                d_s = d;
                s = i;
            }
        }

        if (s === -1) s = n;

        return {
            distance_nearest: d_n,
            index_nearest: n,
            distance_second: d_s,
            index_second: s,
        };
    }

    /**
     * @private
     */
    _update_clusters() {
        const N = this._N;
        const A = this._A;
        for (let j = 0; j < N; j++) {
            const nearest = this._nearest_medoid(A[j], j);
            this._clusters[j] = nearest.index_nearest;
        }
    }

    /**
     * Computes `K` clusters out of the `matrix`.
     * @param {number} K - Number of clusters.
     * @param {number[]} cluster_medoids
     */
    init(K, cluster_medoids) {
        if (!K) K = this._parameters.K;
        if (!cluster_medoids) cluster_medoids = this._get_random_medoids(K);
        this._cluster_medoids = cluster_medoids;
        const max_iter = this._max_iter;
        let finish = false;
        let i = 0;
        do {
            finish = this._iteration();
        } while (!finish && ++i < max_iter);
        this._update_clusters();
        this._is_initialized = true;
        return this;
    }

    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @private
     * @param {number} K - Number of clusters
     * @returns {number[]}
     */
    _get_random_medoids(K) {
        const N = this._N;
        const A = this._A;
        const indices = linspace(0, N - 1);
        const randomizer = this._randomizer;
        const n = Math.min(N, 10 + Math.ceil(Math.sqrt(N)));

        // Handle case where K >= N
        if (K >= N) {
            return indices.slice(0, N);
        }

        /** @type {number[]} */
        const medoids = [];

        // first medoid: select from a random sample of size n
        let best_j = -1;
        let min_td = Infinity;
        let S = randomizer.choice(indices, n);
        for (let j = 0; j < S.length; ++j) {
            let td = 0;
            const S_j = S[j];
            const x_j = A[S_j];
            for (let o = 0; o < S.length; ++o) {
                if (o === j) continue;
                td += this._get_distance(S_j, S[o], x_j, A[S[o]]);
            }
            if (td < min_td) {
                min_td = td;
                best_j = S_j;
            }
        }
        medoids.push(best_j);

        // other medoids: greedy additive selection (Algorithm LAB)
        for (let i = 1; i < K; ++i) {
            let best_idx = -1;
            let best_delta = Infinity;

            const remainingIndices = indices.filter((idx) => !medoids.includes(idx));
            if (remainingIndices.length === 0) break;

            S = randomizer.choice(remainingIndices, Math.min(n, remainingIndices.length));
            for (let j = 0; j < S.length; ++j) {
                let deltaTD = 0;
                const S_j = S[j];
                const x_j = A[S_j];

                // Estimate TD reduction on the sample S
                for (let o = 0; o < S.length; ++o) {
                    if (o === j) continue;
                    const S_o = S[o];
                    const x_o = A[S_o];

                    // Closest distance to current medoids
                    let min_d_existing = Infinity;
                    for (let m = 0; m < medoids.length; m++) {
                        const d = this._get_distance(S_o, medoids[m], x_o, A[medoids[m]]);
                        if (d < min_d_existing) min_d_existing = d;
                    }

                    const delta = this._get_distance(S_j, S_o, x_j, x_o) - min_d_existing;
                    if (delta < 0) {
                        deltaTD += delta;
                    }
                }

                if (deltaTD < best_delta) {
                    best_delta = deltaTD;
                    best_idx = S_j;
                }
            }
            if (best_idx !== -1) {
                medoids.push(best_idx);
            }
        }
        return medoids;
    }
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
class MeanShift extends Clustering {
    /**
     * @private
     * @type {number}
     */
    _bandwidth;
    /**
     * @private
     * @type {number}
     */
    _max_iter;
    /**
     * @private
     * @type {number}
     */
    _tolerance;
    /**
     * @private
     * @type {(dist: number) => number}
     */
    _kernel;
    /**
     * @type {Matrix}
     */
    _points;
    /**
     * @private
     * @type {number[] | undefined}
     */
    _clusters;
    /**
     * @private
     * @type {number[][] | undefined}
     */
    _cluster_list;

    /**
     *
     * @param {InputType} points
     * @param {Partial<ParametersMeanShift>} parameters
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersMeanShift} */ (
                Object.assign({ seed: 1212, metric: euclidean, bandwidth: 0, kernel: "gaussian" }, parameters)
            ),
        );

        // Ensure bandwidth is positive
        this._bandwidth = parameters.bandwidth ?? this._compute_bandwidth(this._matrix);
        this._max_iter = parameters.max_iter ?? Math.max(10, Math.floor(10 * Math.log10(this._N)));
        this._tolerance = parameters.tolerance ?? 1e-3;
        const kernel_param = parameters.kernel ?? "gaussian";
        // If kernel is a string, map to function
        if (typeof kernel_param === "string") {
            if (kernel_param === "flat") {
                this._kernel = (dist) => (dist <= this._bandwidth ? 1 : 0);
            } else {
                // gaussian (default)
                this._kernel = (dist) => Math.exp(-(dist * dist) / (2 * this._bandwidth * this._bandwidth));
            }
        } else {
            // custom function
            this._kernel = kernel_param;
        }

        // Copy points to a mutable matrix
        this._points = this._matrix.clone();

        this._mean_shift();
        this._assign_clusters();
    }

    /**
     * Helper to compute bandwidth if not provided
     * @private
     * @param {Matrix} matrix
     * @returns {number}
     */
    _compute_bandwidth(matrix) {
        const N = matrix.shape[0];
        //const D = matrix.shape[1];
        // Compute average pairwise distance
        let totalDist = 0;
        for (let i = 0; i < N; ++i) {
            const row_i = matrix.row(i);
            for (let j = i + 1; j < N; ++j) {
                const row_j = matrix.row(j);
                const dist = this._parameters.metric(row_i, row_j);
                totalDist += dist;
            }
        }
        const avgDist = totalDist / ((N * (N - 1)) / 2);
        // Use a fraction of avgDist as bandwidth
        return avgDist / 2;
    }

    /**
     * Compute kernel weight
     * @private
     * @param {number} dist
     * @returns {number}
     */
    _kernel_weight(dist) {
        return this._kernel(dist);
    }

    /**
     * Perform mean shift iterations
     * @private
     */
    _mean_shift() {
        const N = this._N;
        const D = this._D;
        const points = this._points;
        const metric = this._parameters.metric;
        //const bandwidth = this._bandwidth;
        const kernel = this._kernel_weight.bind(this);
        const tolerance = this._tolerance;

        for (let iter = 0; iter < this._max_iter; ++iter) {
            let max_shift = 0;
            // For each point compute shift
            for (let i = 0; i < N; ++i) {
                const row_i = points.row(i);
                let sum_weights = 0;
                const weighted_sum = new Float64Array(D);
                for (let j = 0; j < N; ++j) {
                    const row_j = points.row(j);
                    const dist = metric(row_i, row_j);
                    const weight = kernel(dist);
                    sum_weights += weight;
                    for (let d = 0; d < D; ++d) {
                        weighted_sum[d] += weight * row_j[d];
                    }
                }
                if (sum_weights === 0) {
                    // No neighbors within kernel, shift is zero
                    //const shift = new Float64Array(D);
                    // Compute shift magnitude
                    const shift_norm = Math.sqrt(weighted_sum.reduce((acc, v) => acc + v * v, 0));
                    max_shift = Math.max(max_shift, shift_norm);
                } else {
                    const shift = new Float64Array(D);
                    for (let d = 0; d < D; ++d) {
                        shift[d] = weighted_sum[d] / sum_weights - row_i[d];
                    }
                    const shift_norm = Math.sqrt(shift.reduce((acc, v) => acc + v * v, 0));
                    max_shift = Math.max(max_shift, shift_norm);
                    // Update point
                    for (let d = 0; d < D; ++d) {
                        row_i[d] += shift[d];
                    }
                }
            }
            if (max_shift < tolerance) {
                // Converged
                break;
            }
        }
    }

    /**
     * After convergence, assign clusters based on nearest mode
     * @private
     */
    _assign_clusters() {
        const N = this._N;
        const metric = this._parameters.metric;
        const bandwidth = this._bandwidth;

        // Group points that converged to the same mode
        // Two points are in the same mode if they're within bandwidth/2 of each other
        const mode_threshold = bandwidth * 0.5;
        /** @type {number[][]} */
        const modes = []; // Each mode contains indices of points in that mode
        const point_to_mode = new Array(N).fill(-1);

        for (let i = 0; i < N; ++i) {
            if (point_to_mode[i] !== -1) continue; // Already assigned to a mode

            const row_i = this._points.row(i);
            const mode = [i];
            point_to_mode[i] = modes.length;

            // Find all points close to this mode
            for (let j = i + 1; j < N; ++j) {
                if (point_to_mode[j] !== -1) continue;

                const row_j = this._points.row(j);
                const dist = metric(row_i, row_j);

                if (dist < mode_threshold) {
                    mode.push(j);
                    point_to_mode[j] = modes.length;
                }
            }

            modes.push(mode);
        }

        // Build final clusters - each mode becomes a cluster
        /** @type {number[][]} */
        const clusters = [];
        const cluster_ids = new Array(N).fill(-1);

        for (let mode_idx = 0; mode_idx < modes.length; ++mode_idx) {
            const mode = modes[mode_idx];
            clusters.push([...mode]);
            for (const point_idx of mode) {
                cluster_ids[point_idx] = mode_idx;
            }
        }

        this._clusters = cluster_ids;
        this._cluster_list = clusters;
    }

    /**
     * @returns {number[][]}
     */
    get_clusters() {
        // Ensure algorithm has been run
        if (!this._cluster_list) {
            this._mean_shift();
            this._assign_clusters();
        }
        return /** @type {number[][]} */ (this._cluster_list);
    }

    /**
     *
     * @returns {number[]}
     */
    get_cluster_list() {
        if (!this._clusters) {
            this._mean_shift();
            this._assign_clusters();
        }
        return /** @type {number[]} */ (this._clusters);
    }
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
class OPTICS extends Clustering {
    /**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     *
     * @param {InputType} points - The data.
     * @param {Partial<ParametersOptics>} [parameters={}]
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
    constructor(points, parameters = {}) {
        super(
            points,
            /** @type {ParametersOptics} */ (
                Object.assign({ epsilon: 1, min_points: 4, metric: euclidean }, parameters)
            ),
        );
        const matrix = this._matrix;
        /**
         * @private
         * @type {DBEntry[]}
         */
        this._ordered_list = [];
        const ordered_list = this._ordered_list;
        /** @type {number[][]} */
        this._clusters = [];
        const clusters = this._clusters;

        const N = this._N;

        /**
         * @private
         * @type {DBEntry[]}
         */
        this._DB = new Array(N).fill(0).map((_, i) => {
            return {
                element: matrix.row(i),
                index: i,
                reachability_distance: undefined,
                processed: false,
            };
        });
        const DB = this._DB;

        this._cluster_index = 0;
        let cluster_index = this._cluster_index;

        for (const p of DB) {
            if (p.processed) continue;
            p.neighbors = this._get_neighbors(p);
            p.processed = true;
            clusters.push([p.index]);
            cluster_index = clusters.length - 1;
            ordered_list.push(p);
            if (this._core_distance(p) !== undefined) {
                const seeds = new Heap(null, (d) => d.reachability_distance, "min");
                this._update(p, seeds);
                this._expand_cluster(seeds, clusters[cluster_index]);
            }
        }
    }

    /**
     * @private
     * @param {DBEntry} p - A point of the data.
     * @returns {DBEntry[]} An array consisting of the `epsilon`-neighborhood of `p`.
     */
    _get_neighbors(p) {
        if (p?.neighbors) return p.neighbors;
        const DB = this._DB;
        const metric = this._parameters.metric;
        const epsilon = this._parameters.epsilon;
        const neighbors = [];
        for (const q of DB) {
            if (q.index === p.index) continue;
            if (metric(p.element, q.element) <= epsilon) {
                neighbors.push(q);
            }
        }
        return neighbors;
    }

    /**
     * @private
     * @param {DBEntry} p - A point of `matrix`.
     * @returns {number|undefined} The distance to the `min_points`-th nearest point of `p`, or undefined if the
     *   `epsilon`-neighborhood has fewer elements than `min_points`.
     */
    _core_distance(p) {
        const min_points = this._parameters.min_points;
        const metric = this._parameters.metric;
        // Need min_points - 1 other points plus the point itself
        if (!p.neighbors || p.neighbors.length < min_points - 1) {
            return undefined;
        }
        // Sort neighbors by distance to find the MinPts-th closest
        const sortedNeighbors = p.neighbors.toSorted(
            (a, b) => metric(p.element, a.element) - metric(p.element, b.element),
        );
        // MinPts-th closest is at index min_points - 2 (0-indexed, excluding p itself)
        return metric(p.element, sortedNeighbors[min_points - 2].element);
    }

    /**
     * Updates the reachability distance of the points.
     *
     * @private
     * @param {DBEntry} p
     * @param {Heap<DBEntry>} seeds
     */
    _update(p, seeds) {
        const metric = this._parameters.metric;
        const core_distance = this._core_distance(p);
        // If p is not a core point, don't update seeds
        if (core_distance === undefined) {
            return;
        }
        const neighbors = this._get_neighbors(p); //p.neighbors;
        for (const q of neighbors) {
            if (q.processed) continue;
            const new_reachability_distance = Math.max(core_distance, metric(p.element, q.element));
            //if (q.reachability_distance == undefined) { // q is not in seeds
            if (seeds.raw_data().findIndex((d) => d.element === q) < 0) {
                q.reachability_distance = new_reachability_distance;
                seeds.push(q);
            } else {
                // q is in seeds
                if (new_reachability_distance < (q.reachability_distance ?? Infinity)) {
                    q.reachability_distance = new_reachability_distance;
                    seeds = Heap.heapify(seeds.data(), (d) => d.reachability_distance ?? Infinity, "min"); // seeds change key =/
                }
            }
        }
    }

    /**
     * Expands the `cluster` with points in `seeds`.
     *
     * @private
     * @param {Heap<DBEntry>} seeds
     * @param {number[]} cluster
     */
    _expand_cluster(seeds, cluster) {
        const ordered_list = this._ordered_list;
        while (!seeds.empty) {
            const q = /** @type {{ element: DBEntry, value: number}} */ (seeds.pop()).element;
            q.neighbors = this._get_neighbors(q);
            q.processed = true;
            cluster.push(q.index);
            ordered_list.push(q);
            if (this._core_distance(q) !== undefined) {
                this._update(q, seeds);
                // Recursive call removed - while loop handles iteration correctly
            }
        }
    }

    /**
     * Returns an array of clusters.
     *
     * @returns {number[][]} Array of clusters with the indices of the rows in given `matrix`.
     */
    get_clusters() {
        const clusters = [];
        const outliers = [];
        const min_points = this._parameters.min_points;
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
     * @returns {number[]} Returns an array, where the ith entry defines the cluster affirmation of the ith point of
     *   given data. (-1 stands for outlier)
     */
    get_cluster_list() {
        const N = this._matrix.shape[0];
        /** @type {number[]} */
        const result = new Array(N).fill(0);
        const clusters = this.get_clusters();
        for (let i = 0, n = clusters.length; i < n; ++i) {
            const cluster = clusters[i];
            for (const index of cluster) {
                result[index] = i < n - 1 ? i : -1;
            }
        }
        return result;
    }
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
class XMeans extends Clustering {
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
    constructor(points, parameters = {}) {
        const defaults = {
            K_max: 10,
            K_min: 2,
            metric: euclidean,
            seed: 1212,
            min_cluster_size: 35,
            tolerance: 0.001,
        };
        super(points, /** @type {ParametersXMeans} */ (Object.assign(defaults, parameters)));
        this._randomizer = new Randomizer(this._parameters.seed);

        /** @type {KMeans | null} */
        this._best_kmeans = null;

        // Run XMeans algorithm
        this._run();
    }

    /**
     * Run the XMeans algorithm
     *
     * @private
     */
    _run() {
        /** @type {Map<number, CandidateResult>} */
        const candidates = new Map();
        const A = this._matrix;

        // Initialize with K_min clusters
        let current_kmeans = new KMeans(this._points, {
            K: this._parameters.K_min,
            metric: this._parameters.metric,
            seed: this._parameters.seed,
        });

        let K = this._parameters.K_min;

        candidates.set(K, {
            kmeans: current_kmeans,
            score: -Infinity,
        });

        // Iteratively improve clustering
        while (K < this._parameters.K_max) {
            const clusters = current_kmeans.get_clusters();
            const centroids = current_kmeans.centroids;

            // Try splitting each cluster
            /** @type {SplitResult[]} */
            const split_results = [];

            for (let j = 0; j < clusters.length; ++j) {
                const cluster = clusters[j];

                // Skip small clusters - need enough points for reliable BIC
                if (cluster.length < this._parameters.min_cluster_size) {
                    continue;
                }

                // Get subset data for this cluster
                /** @type {number[][]} */
                const subset_points = cluster.map((idx) => {
                    const row = A.row(idx);
                    return Array.from(row);
                });

                // Calculate BIC for parent (single cluster)
                const parent_bic = this._bic([cluster], [centroids[j]]);

                // Run KMeans with K=2 on subset
                const subset_kmeans = new KMeans(subset_points, {
                    K: 2,
                    metric: this._parameters.metric,
                    seed: this._randomizer.seed,
                });

                const child_clusters_local = subset_kmeans.get_clusters();
                const child_centroids = subset_kmeans.centroids;

                // Map local indices back to global indices
                /** @type {number[][]} */
                const child_clusters_global = child_clusters_local.map((local_cluster) =>
                    local_cluster.map((local_idx) => cluster[local_idx]),
                );

                // Calculate BIC for children (split into 2 clusters)
                const children_bic = this._bic(child_clusters_global, child_centroids);

                split_results.push({
                    index: j,
                    bic_parent: parent_bic,
                    bic_children: children_bic,
                    child_clusters: child_clusters_global,
                    child_centroids: child_centroids,
                });
            }

            // Keep all splits that improve BIC (BIC_children > BIC_parent)
            /** @type {SplitResult[]} */
            const accepted_splits = split_results.filter((result) => result.bic_children > result.bic_parent);

            // If no splits improve BIC, we're done
            if (accepted_splits.length === 0) {
                break;
            }

            // Build new centroids array: keep non-split centroids + add split centroids
            /** @type {Float64Array[]} */
            const new_centroids = [];
            const split_indices = new Set();

            // Sort accepted splits by improvement (descending)
            accepted_splits.sort((a, b) => b.bic_children - b.bic_parent - (a.bic_children - a.bic_parent));

            for (const split of accepted_splits) {
                if (centroids.length + split_indices.size + 1 <= this._parameters.K_max) {
                    split_indices.add(split.index);
                } else {
                    break;
                }
            }

            for (let i = 0; i < centroids.length; ++i) {
                if (split_indices.has(i)) {
                    // This cluster was split - add both child centroids
                    const split_result = accepted_splits.find((s) => s.index === i);
                    if (split_result) {
                        new_centroids.push(...split_result.child_centroids);
                    }
                } else {
                    // This cluster wasn't split - keep its centroid
                    new_centroids.push(centroids[i]);
                }
            }

            // Run KMeans on full dataset with new centroids as initialization
            // This is crucial - we need to reassign all points to all clusters
            const newK = new_centroids.length;

            // Create a new KMeans instance with K set to new number of clusters
            current_kmeans = new KMeans(this._matrix, {
                K: newK,
                metric: this._parameters.metric,
                seed: this._randomizer.seed,
                initial_centroids: new_centroids,
            });

            // Store the candidate with the BIC of the FULL dataset
            candidates.set(newK, {
                kmeans: current_kmeans,
                score: this._bic(current_kmeans.get_clusters(), current_kmeans.centroids),
            });

            K = newK;
        }

        // Select best candidate based on BIC score
        this._best_kmeans = this._select_best_candidate(candidates);
    }

    /**
     * Select the best candidate based on BIC score
     *
     * @private
     * @param {Map<number, CandidateResult>} candidates
     * @returns {KMeans}
     */
    _select_best_candidate(candidates) {
        if (candidates.size === 0) {
            throw new Error("No candidates found");
        }

        const first_candidate = candidates.get(this._parameters.K_min);
        if (!first_candidate) {
            throw new Error("Missing initial candidate");
        }

        let best_score = first_candidate.score;
        /** @type {KMeans} */
        let best_kmeans = first_candidate.kmeans;

        for (const candidate of candidates.values()) {
            if (candidate.score > best_score) {
                best_score = candidate.score;
                best_kmeans = candidate.kmeans;
            }
        }

        return best_kmeans;
    }

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
    _bic(clusters, centroids) {
        const A = this._matrix;
        const D = this._D;
        const K = centroids.length;

        let total_variance = 0;
        let N = 0;

        // Calculate total variance (sum of squared distances)
        for (let i = 0; i < K; ++i) {
            const cluster = clusters[i];
            const centroid = centroids[i];
            N += cluster.length;

            for (let j = 0; j < cluster.length; ++j) {
                const point_idx = cluster[j];
                const point = A.row(point_idx);
                // Sum of squared distances (variance term)
                total_variance += euclidean_squared(centroid, point);
            }
        }

        // Not enough points for meaningful BIC
        if (N <= K) {
            return -Infinity;
        }

        // Estimate variance (ML estimate)
        const variance = total_variance / (N - K);

        // Handle case of zero variance (all points identical)
        if (variance <= 0) {
            return -Infinity;
        }

        // Number of free parameters: (K-1) cluster weights + K*D centroid coordinates + 1 variance
        const p = K - 1 + D * K + 1;

        // Calculate log-likelihood
        let log_likelihood = 0;
        const log_2pi = Math.log(2 * Math.PI);

        for (let i = 0; i < K; ++i) {
            const n = clusters[i].length;
            if (n <= 1) continue;

            // Log-likelihood for cluster i
            const cluster_log_likelihood =
                n * Math.log(n / N) - 0.5 * n * log_2pi - 0.5 * n * D * Math.log(variance) - 0.5 * (n - 1);

            log_likelihood += cluster_log_likelihood;
        }

        // BIC = log_likelihood - 0.5 * p * ln(N)
        return log_likelihood - 0.5 * p * Math.log(N);
    }

    /**
     * Get the computed clusters
     *
     * @returns {number[][]} Array of clusters, each containing indices of points
     */
    get_clusters() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.get_clusters();
    }

    /** @returns {number[]} The cluster list */
    get_cluster_list() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.get_cluster_list();
    }

    /**
     * Get the final centroids
     *
     * @returns {Float64Array[]} Array of centroids
     */
    get centroids() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.centroids;
    }

    /**
     * Get the optimal number of clusters found
     *
     * @returns {number} The number of clusters
     */
    get k() {
        if (!this._best_kmeans) {
            throw new Error("XMeans has not been run");
        }
        return this._best_kmeans.k;
    }
}

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
class DR {
    /** @type {number} */
    _D;
    /** @type {number} */
    _N;
    /** @type {Randomizer} */
    _randomizer;
    /** @type {boolean} */
    _is_initialized;

    /**
     * Takes the default parameters and seals them, remembers the type of input `X`, and initializes the random number
     * generator.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Para} default_parameters - Object containing default parameterization of the DR method.
     * @param {Partial<Para>} parameters - Object containing parameterization of the DR method to override defaults.
     */
    constructor(X, default_parameters, parameters = {}) {
        /** @type {T} */
        this.__input = X;

        /** @type {Para} */
        this._parameters = /** @type {Para} */ Object.seal({
            ...default_parameters,
            ...parameters,
        });
        /** @type {"array" | "matrix" | "typed"} */
        this._type;
        /** @type {Matrix} */
        this.X;
        /** @type {Matrix} */
        this.Y;

        if (Array.isArray(X)) {
            if (X[0] instanceof Float64Array) {
                this._type = "typed";
            } else {
                this._type = "array";
            }
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("No valid type for X!");
        }
        const [N, D] = this.X.shape;
        this._N = N;
        this._D = D;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._is_initialized = false;
    }

    /**
     * Get all Parameters.
     * @overload
     * @returns {Para}
     */
    /**
     * Get value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @returns {Para[K]}
     */
    /**
     * Set value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @param {Para[K]} value - Value of the parameter to set.
     * @returns {this}
     */
    /**
     * @param {keyof Para} [name] - Name of the parameter. If null, returns all parameters as an Object.
     * @param {Para[keyof Para]} [value] - Value of the parameter to set. If name is set and value is not given, returns the
     *   current value.
     * @returns {Para | Para[keyof Para] | this} On setting a parameter, returns the DR object. If name is set and value is not
     *   given, returns the parameter value. If name is null, returns all parameters. On setting a parameter, this
     *   function returns the DR object. If `name` is set and `value == null` then return actual parameter value. If
     *   `name` is not given, then returns all parameters as an Object.
     * @example
     * ```js
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for `d = 3`.
     * DR.parameter("d"); // returns 3
     * DR.parameter("d", 2); // sets parameter `d` to 2 and returns `DR`.
     * ```
     *
     */
    parameter(name, value) {
        if (name === undefined && value === undefined) {
            return Object.assign({}, this._parameters);
        }
        if (name && !Object.hasOwn(this._parameters, name)) {
            throw new Error(`${String(name)} is not a valid parameter!`);
        }
        if (name && value !== undefined) {
            this._parameters[name] = value;
            this._is_initialized = false;
            return this;
        } else if (name) {
            return this._parameters[name];
        }
        throw new Error("Should not happen!");
    }

    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {T} The projection.
     */
    transform(...args) {
        this.check_init();
        return this.projection;
    }

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
    static transform(X, parameters, ...args) {
        const dr = new DR(X, parameters, parameters);
        return /** @type {T} */ (dr.transform());
    }

    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {Generator<T, T, void>} The intermediate steps of the projection.
     */
    *generator(...args) {
        const R = this.transform(...args);
        yield R;
        return R;
    }

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
    static *generator(X, parameters, ...args) {
        const dr = new DR(X, parameters, parameters);
        const generator = dr.generator(...args);
        let result;
        do {
            result = generator.next();
            yield result.value;
        } while (!result.done);

        return result.value;
    }

    /**
     * @abstract
     * @param {...unknown} args
     */
    init(...args) {
    }

    /**
     * If the respective DR method has an `init` function, call it before `transform`.
     *
     * @returns {DR<T, Para>}
     */
    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
        return this;
    }

    /** @returns {T} The projection in the type of input `X`. */
    get projection() {
        if (Object.hasOwn(this, "Y")) {
            this.check_init();
            //return this._type === "matrix" ? this.Y : this.Y.to2dArray();
            if (this._type === "matrix") {
                return /** @type {T} */ (/** @type {any} */ (this.Y));
            } else if (this._type === "typed") {
                return /** @type {T} */ (/** @type {any} */ (this.Y.to2dArray()));
            } else {
                return /** @type {T} */ (/** @type {any} */ (this.Y.asArray()));
            }
        } else {
            throw new Error("The dataset is not transformed yet!");
        }
    }

    /**
     * Computes the projection.
     *
     * @param {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<T>} The dimensionality reduced dataset.
     */
    async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     *
     * @template {{ seed?: number }} Para
     * @param {InputType} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise<X>} A promise yielding the dimensionality reduced dataset.
     */
    static async transform_async(X, parameters, ...args) {
        return DR.transform(X, parameters, ...args);
    }
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
class FASTMAP extends DR {
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets.
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersFASTMAP>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212 }, parameters);
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     *
     * @private
     * @param {(a: number, b: number) => number} dist
     * @returns {[number, number, number]} An array consisting of first index, second index, and distance between the
     *   two points.
     */
    _choose_distant_objects(dist) {
        const X = this.X;
        const N = X.shape[0];
        let a_index = this._randomizer.random_int % N;
        /** @type {number | null} */
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        if (b_index === null) throw new Error("should not happen!");
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
     *
     * @returns {T} The `d`-dimensional projection of the data matrix `X`.
     */
    transform() {
        const X = this.X;
        const N = X.shape[0];
        const d = /** @type {number} */ (this._parameters.d);
        const metric = /** @type {typeof euclidean} */ (this._parameters.metric);
        const Y = new Matrix(N, d, 0);
        /** @type {(a: number, b: number) => number} */
        let dist = (a, b) => metric(X.row(a), X.row(b));

        for (let _col = 0; _col < d; ++_col) {
            const old_dist = dist;
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

    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        return dr.transform_async();
    }
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
class KNN {
    /** @type {T[]} */
    _elements;
    /** @type {Para} */
    _parameters;
    /** @type {"typed" | "array"} */
    _type;

    /**
     * @param {T[]} elements
     * @param {Para} parameters
     */
    constructor(elements, parameters) {
        if (elements.length === 0) throw new Error("Elements needs to contain at least one element!");
        if (elements[0] instanceof Float64Array) {
            this._type = "typed";
        } else {
            this._type = "array";
        }
        this._parameters = parameters;
        this._elements = elements;
    }

    /**
     * @abstract
     * @param {T} t
     * @param {number} k
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search(t, k) {
        throw new Error("The function search must be implemented!");
    }

    /**
     * @abstract
     * @param {number} i
     * @param {number} k
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_by_index(i, k) {
        throw new Error("The function search_by_index must be implemented!");
    }
}

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
class Annoy extends KNN {
    /**
     * Creates a new Annoy-style index with random projection trees.
     *
     * @param {T[]} elements - Elements to index
     * @param {ParametersAnnoy} [parameters={}] - Configuration parameters
     */
    constructor(
        elements,
        parameters = {
            metric: euclidean,
            numTrees: 10,
            maxPointsPerLeaf: 10,
            seed: 1212,
        },
    ) {
        // Handle empty initialization - use dummy element
        const hasElements = elements && elements.length > 0;
        const firstElement = /** @type {T} */ (hasElements ? elements[0] : new Float64Array([0]));

        super([firstElement], parameters);

        this._metric = this._parameters.metric ?? euclidean;
        this._numTrees = this._parameters.numTrees ?? 10;
        this._maxPointsPerLeaf = this._parameters.maxPointsPerLeaf ?? 10;
        this._seed = this._parameters.seed ?? 1212;
        this._randomizer = new Randomizer(this._seed);

        /**
         * @private
         * @type {AnnoyNode<T>[]}
         */
        this._trees = [];

        // Build trees
        if (hasElements) {
            // Reset elements and rebuild properly
            /** @type {T[]} */
            this._elements = [];
            this._trees = [];
            this.add(elements);
        }
    }

    /**
     * Get the number of trees in the index.
     * @returns {number}
     */
    get num_trees() {
        return this._trees.length;
    }

    /**
     * Get the total number of nodes in all trees.
     * @returns {number}
     */
    get num_nodes() {
        let total = 0;
        for (const tree of this._trees) {
            total += this._countNodes(tree);
        }
        return total;
    }

    /**
     * @private
     * @param {any} node
     * @returns {number}
     */
    _countNodes(node) {
        if (!node) return 0;
        return 1 + this._countNodes(node.left) + this._countNodes(node.right);
    }

    /**
     * Add elements to the Annoy index.
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements) {
        // Extend elements array
        this._elements = this._elements.concat(elements);

        // Rebuild all trees with new elements
        this._trees = [];
        this._buildTrees();

        return this;
    }

    /**
     * Build all random projection trees.
     * @private
     */
    _buildTrees() {
        const elements = this._elements;
        const n = elements.length;

        for (let t = 0; t < this._numTrees; t++) {
            // Create index array for this tree
            const indices = Array.from({ length: n }, (_, i) => i);
            const tree = this._buildTreeRecursive(indices);
            this._trees.push(tree);
        }
    }

    /**
     * Recursively build a random projection tree.
     * @private
     * @param {number[]} indices - Indices of elements to include
     * @returns {AnnoyNode<T>}
     */
    _buildTreeRecursive(indices) {
        const elements = this._elements;

        // Base case: small enough to be a leaf
        if (indices.length <= this._maxPointsPerLeaf) {
            return {
                isLeaf: true,
                indices: indices,
                normal: [],
                offset: 0,
                left: null,
                right: null,
            };
        }

        // Select two random points to define the splitting hyperplane
        const idx1 = indices[Math.floor(this._randomizer.random * indices.length)];
        const idx2 = indices[Math.floor(this._randomizer.random * indices.length)];

        const point1 = elements[idx1];
        const point2 = elements[idx2];

        // Compute normal vector (point2 - point1)
        const dim = point1.length;
        /** @type {number[]} */
        const normal = new Array(dim);
        for (let i = 0; i < dim; i++) {
            normal[i] = point2[i] - point1[i];
        }

        // Normalize
        let norm = 0;
        for (let i = 0; i < dim; i++) {
            norm += normal[i] * normal[i];
        }
        norm = Math.sqrt(norm);

        if (norm > 1e-10) {
            for (let i = 0; i < dim; i++) {
                normal[i] /= norm;
            }
        }

        // Compute midpoint and offset
        /** @type {number[]} */
        const midpoint = new Array(dim);
        for (let i = 0; i < dim; i++) {
            midpoint[i] = (point1[i] + point2[i]) / 2;
        }

        // Compute offset: dot(normal, midpoint)
        let offset = 0;
        for (let i = 0; i < dim; i++) {
            offset += normal[i] * midpoint[i];
        }

        // Split points based on which side of hyperplane they fall
        const leftIndices = [];
        const rightIndices = [];

        for (const idx of indices) {
            const point = elements[idx];
            let dot = 0;
            for (let i = 0; i < dim; i++) {
                dot += normal[i] * point[i];
            }

            if (dot < offset) {
                leftIndices.push(idx);
            } else {
                rightIndices.push(idx);
            }
        }

        // Handle edge case where all points fall on one side
        if (leftIndices.length === 0 || rightIndices.length === 0) {
            return {
                isLeaf: true,
                indices: indices,
                normal: [],
                offset: 0,
                left: null,
                right: null,
            };
        }

        // Recursively build subtrees
        const left = this._buildTreeRecursive(leftIndices);
        const right = this._buildTreeRecursive(rightIndices);

        return {
            isLeaf: false,
            indices: [],
            normal: normal,
            offset: offset,
            left: left,
            right: right,
        };
    }

    /**
     * Compute distance from point to hyperplane.
     * @private
     * @param {T} point
     * @param {number[]} normal
     * @param {number} offset
     * @returns {number} Signed distance (positive = right side, negative = left side)
     */
    _distanceToHyperplane(point, normal, offset) {
        let dot = 0;
        for (let i = 0; i < point.length; i++) {
            dot += normal[i] * point[i];
        }
        return dot - offset;
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

        // Collect candidates from all trees using priority queue
        const candidates = new Set();

        // Collect more candidates for better recall
        // Search at least k * numTrees * 2 candidates
        const minCandidates = Math.min(k * this._numTrees * 3, elements.length);

        for (const tree of this._trees) {
            this._searchTreePriority(tree, query, candidates, minCandidates);
        }

        // Compute exact distances for all candidates
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

        // If we still don't have enough candidates, do a linear scan fallback
        if (best.length < k) {
            for (let i = 0; i < elements.length && best.length < k; i++) {
                if (candidates.has(i)) continue;

                const element = elements[i];
                if (!element || element.length !== query.length) continue;

                const dist = metric(query, element);
                best.push({ index: i, distance: dist });
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
     * Search tree using priority queue for better recall.
     * Explores nodes in order of distance to hyperplane.
     * @private
     * @param {AnnoyNode<T>} node
     * @param {T} query
     * @param {Set<number>} candidates
     * @param {number} maxCandidates
     */
    _searchTreePriority(node, query, candidates, maxCandidates) {
        if (!node) return;

        // Priority queue entry: { node, distance }
        /** @type {Heap<{ node: AnnoyNode<T>; dist: number }>} */
        const pq = new Heap(null, (d) => d.dist, "min");
        pq.push({ node: node, dist: 0 });

        while (!pq.empty && candidates.size < maxCandidates) {
            const entry = pq.pop();
            if (!entry) continue;

            const currentNode = entry.element.node;

            // Leaf node: add all points
            if (currentNode.isLeaf) {
                for (const idx of currentNode.indices) {
                    candidates.add(idx);
                    if (candidates.size >= maxCandidates) return;
                }
                continue;
            }

            // Internal node: compute distance to hyperplane
            const dist = this._distanceToHyperplane(query, currentNode.normal, currentNode.offset);

            // Determine which side is closer
            const closerSide = dist < 0 ? currentNode.left : currentNode.right;
            const fartherSide = dist < 0 ? currentNode.right : currentNode.left;

            // Add closer side with priority 0 (explore first)
            if (closerSide) {
                pq.push({ node: closerSide, dist: 0 });
            }

            // Add farther side with priority = |dist| (explore later if needed)
            if (fartherSide && candidates.size < maxCandidates) {
                pq.push({ node: fartherSide, dist: Math.abs(dist) });
            }
        }
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

    /**
     * Alias for search_by_index for backward compatibility.
     *
     * @param {number} i - Index of the query element
     * @param {number} [k=5] - Number of nearest neighbors to return
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_index(i, k = 5) {
        return this.search_by_index(i, k);
    }
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
class BallTree extends KNN {
    /**
     * Generates a BallTree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the BallTree
     * @param {ParametersBallTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     */
    constructor(elements, parameters = { metric: euclidean, seed: 1212 }) {
        super(elements, Object.assign({ seed: 1212 }, parameters));
        /**
         * @private
         * @type {BallTreeNode<T> | BallTreeLeaf<T>}
         */
        this._root = this._construct(elements.map((element, index) => ({ index, element })));
    }

    /** @returns {Metric} */
    get _metric() {
        return this._parameters.metric;
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @returns {BallTreeNode<T> | BallTreeLeaf<T>} Root of balltree.
     */
    _construct(elements) {
        if (elements.length === 1) {
            return new BallTreeLeaf(elements);
        } else {
            const c = this._greatest_spread(elements);
            const sorted_elements = elements.sort((a, b) => a.element[c] - b.element[c]);
            const n = sorted_elements.length;
            const p_index = Math.floor(n / 2);
            const p = sorted_elements[p_index];
            const L = sorted_elements.slice(0, p_index);
            const R = sorted_elements.slice(p_index, n);
            const radius = Math.max(...elements.map((d) => this._metric(p.element, d.element)));
            let B;
            if (L.length > 0 && R.length > 0) {
                B = new BallTreeNode(p, this._construct(L), this._construct(R), radius);
            } else {
                B = new BallTreeLeaf(elements);
            }
            return B;
        }
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} B
     * @returns {number}
     */
    _greatest_spread(B) {
        const d = B[0].element.length;
        const start = new Array(d);

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
        spread = spread.map((d) => d[1] - d[0]);

        let c = 0;
        for (let i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c;
    }

    /**
     * @param {number} i
     * @param {number} k
     */
    search_by_index(i, k = 5) {
        return this.search(this._elements[i], k);
    }

    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        /** @type {Heap<ElementWithIndex<T>>} */
        const heap = new Heap(null, (d) => this._metric(d.element, t), "max");
        this._search(t, k, heap, this._root);

        // Convert heap to result array
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (heap.length > 0) {
            const item = /** @type {{ element: ElementWithIndex<T>; value: number }} */ (heap.pop());
            result.push({
                element: item.element.element,
                index: item.element.index,
                distance: item.value,
            });
        }
        return result.reverse(); // Reverse to get closest first
    }

    /**
     * @private
     * @param {T} t - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {Heap<ElementWithIndex<T>>} Q - Heap consists of the currently found `k` nearest neighbors.
     * @param {BallTreeNode<T> | BallTreeLeaf<T>} B
     */
    _search(t, k, Q, B) {
        if (!B) return;

        if (B instanceof BallTreeNode) {
            const dist_to_pivot = this._metric(t, B.pivot.element);
            if (Q.length >= k && dist_to_pivot - B.radius >= (Q.first?.value ?? -Infinity)) {
                return;
            }

            const c1 = B.child1;
            const c2 = B.child2;

            let d1 = Infinity;
            let d2 = Infinity;

            if (c1 instanceof BallTreeNode) d1 = this._metric(t, c1.pivot.element);
            else if (c1 instanceof BallTreeLeaf) d1 = this._metric(t, c1.points[0].element);

            if (c2 instanceof BallTreeNode) d2 = this._metric(t, c2.pivot.element);
            else if (c2 instanceof BallTreeLeaf) d2 = this._metric(t, c2.points[0].element);

            if (d1 < d2) {
                if (c1) this._search(t, k, Q, c1);
                if (c2) this._search(t, k, Q, c2);
            } else {
                if (c2) this._search(t, k, Q, c2);
                if (c1) this._search(t, k, Q, c1);
            }
        } else if (B instanceof BallTreeLeaf) {
            for (let i = 0, n = B.points.length; i < n; ++i) {
                const p = B.points[i];
                const dist = this._metric(p.element, t);
                if (Q.length < k) {
                    Q.push(p);
                } else if (dist < (Q.first?.value ?? Infinity)) {
                    Q.pop();
                    Q.push(p);
                }
            }
        }
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeNode {
    /**
     * @param {ElementWithIndex<T>} pivot
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child1
     * @param {BallTreeNode<T> | BallTreeLeaf<T> | null} child2
     * @param {number} radius
     */
    constructor(pivot, child1 = null, child2 = null, radius = 0) {
        this.pivot = pivot;
        this.child1 = child1;
        this.child2 = child2;
        this.radius = radius;
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class BallTreeLeaf {
    /** @param {ElementWithIndex<T>[]} points */
    constructor(points) {
        this.points = points;
    }
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
class HNSW extends KNN {
    /**
     * Creates a new HNSW index.
     *
     * @param {T[]} points - Initial points to add to the index
     * @param {ParametersHNSW} [parameters={}] - Configuration parameters
     */
    constructor(
        points,
        parameters = {
            metric: euclidean,
            heuristic: true,
            m: 16,
            ef_construction: 200,
            m0: null,
            mL: null,
            seed: 1212,
            ef: 50,
        },
    ) {
        // Handle empty initialization - use dummy element
        const hasElements = points && points.length > 0;
        let firstElement = /** @type {T} */ (hasElements ? points[0] : new Float64Array([0]));

        // Validate all points have consistent dimensions
        if (hasElements) {
            const expected_dim = firstElement.length;
            for (let i = 1; i < points.length; i++) {
                if (!points[i] || points[i].length !== expected_dim) {
                    console.warn(
                        `HNSW: Point ${i} has inconsistent dimensions (expected ${expected_dim}, got ${points[i]?.length})`,
                    );
                    // Remove invalid points
                    points = points.filter((_, idx) => idx === 0 || points[idx]?.length === expected_dim);
                    firstElement = points[0];
                }
            }
        }

        super([firstElement], parameters);

        // Store reference to elements before clearing
        const elementsToAdd = hasElements ? [...points] : [];
        /** @type {T[]} */
        this._elements = [];

        /** @type {Metric} */
        this._metric = this._parameters.metric || euclidean;

        /** @type {Function} */
        this._select = this._parameters.heuristic ? this._select_heuristic.bind(this) : this._select_simple.bind(this);

        /**
         * @private
         * @type {Map<number, Layer>}
         */
        this._graph = new Map();

        /** @type {number} */
        this._next_index = 0;

        // Validate and set parameters
        const m_param = this._parameters.m ?? 16;
        if (m_param <= 0 || !Number.isInteger(m_param)) {
            throw new Error("HNSW: parameter 'm' must be a positive integer");
        }
        /** @type {number} */
        this._m = Math.max(2, m_param);

        const ef_construction_param = this._parameters.ef_construction ?? 200;
        if (ef_construction_param <= 0 || !Number.isInteger(ef_construction_param)) {
            throw new Error("HNSW: parameter 'ef_construction' must be a positive integer");
        }
        /** @type {number} */
        this._ef_construction = ef_construction_param;

        const ef_param = this._parameters.ef ?? 50;
        if (ef_param <= 0 || !Number.isInteger(ef_param)) {
            throw new Error("HNSW: parameter 'ef' must be a positive integer");
        }
        /** @type {number} */
        this._ef = ef_param;

        const m0_param = this._parameters.m0 ?? 2 * this._m;
        if (m0_param <= 0 || !Number.isInteger(m0_param)) {
            throw new Error("HNSW: parameter 'm0' must be a positive integer");
        }
        /** @type {number} */
        this._m0 = m0_param;

        /** @type {number} */
        this._mL = this._parameters.mL ?? 1 / Math.log(this._m);

        /** @type {Randomizer} */
        this._randomizer = new Randomizer(this._parameters.seed);

        /** @type {number} - Current maximum layer in the graph */
        this._L = -1;

        /** @type {number[] | null} - Entry point indices for search */
        this._ep = null;

        // Add initial points
        if (elementsToAdd && elementsToAdd.length > 0) {
            this.add(elementsToAdd);
        }
    }

    /**
     * Add a single element to the index.
     *
     * @param {T} element - Element to add
     * @returns {HNSW<T>} This instance for chaining
     */
    addOne(element) {
        return this.add([element]);
    }

    /**
     * Add multiple elements to the index.
     *
     * @param {T[]} new_elements - Elements to add
     * @returns {HNSW<T>} This instance for chaining
     */
    add(new_elements) {
        // Handle empty array
        if (!new_elements || new_elements.length === 0) {
            return this;
        }

        const m = this._m;
        const ef_construction = this._ef_construction;
        const m0 = this._m0;
        const mL = this._mL;
        const randomizer = this._randomizer;
        const graph = this._graph;

        // Ensure _elements is a proper array that supports push
        if (!Array.isArray(this._elements)) {
            this._elements = Array.from(this._elements);
        }
        const elements = this._elements;

        // Get expected dimension from first existing element or first new element
        const expected_dim = elements.length > 0 ? elements[0].length : new_elements[0]?.length;

        for (const element of new_elements) {
            // Validate element
            if (!element || (!Array.isArray(element) && !(element instanceof Float64Array))) {
                console.warn("HNSW: Skipping invalid element (null, undefined, or not an array)");
                continue;
            }

            // Validate dimensions
            if (element.length !== expected_dim) {
                console.warn(
                    `HNSW: Skipping element with wrong dimensions (expected ${expected_dim}, got ${element.length})`,
                );
                continue;
            }

            elements.push(element);
            const global_index = elements.length - 1;

            // Assign random level to the element
            // Level is drawn from exponential distribution: l = floor(-ln(uniform(0,1)) * mL)
            const rand = Math.max(randomizer.random, 1e-10); // Avoid log(0)
            const l = Math.min(31, Math.floor(-Math.log(rand) * mL));

            let ep_indices = this._ep ? [...this._ep] : null;
            const L = this._L;

            if (L >= 0) {
                // Search from top layer down to min(L, l) + 1
                // These are the layers where element will NOT be inserted
                for (let l_c = L; l_c > l; --l_c) {
                    const search_result = this._search_layer(element, ep_indices, 1, l_c);
                    if (search_result.length > 0) {
                        ep_indices = [search_result[0].index];
                    }
                }

                // Insert element into layers l down to 0
                for (let l_c = Math.min(L, l); l_c >= 0; --l_c) {
                    const layer = graph.get(l_c);
                    if (!layer) continue;

                    layer.point_indices.push(global_index);

                    // Search for ef_construction nearest neighbors
                    let W = this._search_layer(element, ep_indices, ef_construction, l_c);

                    // If graph search returns no results (e.g., graph is empty or disconnected),
                    // fall back to linear search over all existing elements
                    if (W.length === 0 && elements.length > 1) {
                        const fallbackCandidates = [];
                        for (let i = 0; i < elements.length - 1; i++) {
                            const elem = elements[i];
                            if (elem && elem.length === element.length) {
                                fallbackCandidates.push({
                                    element: elem,
                                    index: i,
                                    distance: this._metric(element, elem),
                                });
                            }
                        }
                        fallbackCandidates.sort((a, b) => a.distance - b.distance);
                        W = fallbackCandidates.slice(0, ef_construction);
                        // Update ep_indices for next layer based on fallback results
                        if (l_c === Math.min(L, l)) {
                            ep_indices = W.map((c) => c.index);
                        }
                    }

                    // Select neighbors using heuristic or simple approach (respect heuristic setting on all layers)
                    const neighbor_indices = this._select(element, W, l_c === 0 ? m0 : m, l_c);

                    // Add bidirectional connections
                    for (const neighbor_idx of neighbor_indices) {
                        if (neighbor_idx === global_index) continue;

                        // Add connection from element to neighbor
                        if (!layer.edges.has(global_index)) {
                            layer.edges.set(global_index, []);
                        }
                        layer.edges.get(global_index)?.push(neighbor_idx);

                        // Add connection from neighbor to element
                        if (!layer.edges.has(neighbor_idx)) {
                            layer.edges.set(neighbor_idx, []);
                        }
                        const neighbor_edge_list = layer.edges.get(neighbor_idx);
                        if (neighbor_edge_list && !neighbor_edge_list.includes(global_index)) {
                            neighbor_edge_list.push(global_index);
                        }

                        // Prune connections if too many
                        const max_conn = l_c === 0 ? m0 : m;
                        const neighbor_edges = layer.edges.get(neighbor_idx);
                        if (neighbor_edges && neighbor_edges.length > max_conn) {
                            const neighbor_element = elements[neighbor_idx];
                            // Filter out self-connections before pruning
                            const valid_neighbor_edges = neighbor_edges.filter((idx) => idx !== neighbor_idx);
                            const neighbor_candidates = valid_neighbor_edges.map((idx) => ({
                                element: elements[idx],
                                index: idx,
                                distance: this._metric(neighbor_element, elements[idx]),
                            }));
                            const pruned =
                                l_c === 0
                                    ? this._select_simple(neighbor_element, neighbor_candidates, max_conn)
                                    : this._select(neighbor_element, neighbor_candidates, max_conn, l_c);
                            layer.edges.set(neighbor_idx, pruned);
                        }
                    }

                    // Use closest neighbor as entry point for next layer (following HNSW paper)
                    if (W.length > 0) {
                        ep_indices = [W[0].index];
                    }
                }
            }

            // If element's level is higher than current max, create new layers
            if (l > L) {
                for (let i = L + 1; i <= l; ++i) {
                    graph.set(i, {
                        l_c: i,
                        point_indices: [global_index],
                        edges: new Map(),
                    });
                }
                // Element becomes the new entry point
                this._ep = [global_index];
                this._L = l;
            }

            // Special case: if this is the first element (L was -1),
            // we need to ensure layer 0 has proper structure for future insertions
            if (L === -1) {
                if (!graph.has(0)) {
                    graph.set(0, {
                        l_c: 0,
                        point_indices: [global_index],
                        edges: new Map(),
                    });
                }
                const layer0 = graph.get(0);
                if (layer0 && !layer0.edges.has(global_index)) {
                    layer0.edges.set(global_index, []);
                }
            }
        }

        return this;
    }

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
    _select_heuristic(q, candidates, M, l_c, extend_candidates = true, keep_pruned_connections = true) {
        if (l_c > this._L) {
            return candidates.map((c) => c.index);
        }

        const metric = this._metric;
        const layer = this._graph.get(l_c);
        const elements = this._elements;

        // Extend candidate set with neighbors of candidates
        const W_set = new Set(candidates.map((c) => c.index));
        if (extend_candidates) {
            for (const c of candidates) {
                const edges = layer?.edges.get(c.index);
                if (edges) {
                    for (const neighbor_idx of edges) {
                        W_set.add(neighbor_idx);
                    }
                }
            }
        }

        // Create extended candidates with distances
        const W = [...W_set]
            .map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            }))
            .sort((a, b) => a.distance - b.distance);

        const R = [];
        const W_discarded = [];

        // Select neighbors: prefer points closer to query than to already selected points
        for (const e of W) {
            if (R.length >= M) break;

            let should_add = true;

            // Check if e is closer to query than to any already selected point
            for (const r of R) {
                const dist_er = metric(e.element, r.element);
                if (dist_er < e.distance) {
                    should_add = false;
                    break;
                }
            }

            if (should_add) {
                R.push(e);
            } else {
                W_discarded.push(e);
            }
        }

        // Add discarded connections if we need more
        if (keep_pruned_connections && R.length < M) {
            for (const e of W_discarded) {
                if (R.length >= M) break;
                R.push(e);
            }
        }

        return R.map((c) => c.index);
    }

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
    _select_simple(q, C, M) {
        if (C.length <= M) return C.map((c) => c.index);

        // Candidates already have distance computed, use it directly
        return C.slice()
            .sort((a, b) => a.distance - b.distance)
            .slice(0, M)
            .map((c) => c.index);
    }

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
    _search_layer(q, ep_indices, ef, l_c) {
        const metric = this._metric;
        const layer = this._graph.get(l_c);
        const elements = this._elements;

        if (!layer || layer.edges.size === 0 || !ep_indices || ep_indices.length === 0) {
            return [];
        }

        // Filter out invalid indices
        const valid_ep_indices = ep_indices.filter((idx) => elements[idx] !== undefined);
        if (valid_ep_indices.length === 0) {
            return [];
        }

        // Visited set to avoid cycles
        const visited = new Set(valid_ep_indices);

        // Candidate set (min-heap): closest unvisited candidates to expand
        const C = new Heap(
            valid_ep_indices.map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            })),
            (item) => item.distance,
            "min",
        );

        // Result set (max-heap): ef closest found neighbors
        const W = new Heap(
            valid_ep_indices.map((idx) => ({
                element: elements[idx],
                index: idx,
                distance: metric(elements[idx], q),
            })),
            (item) => item.distance,
            "max",
        );

        // Algorithm 2 stops when the distance from query to the next candidate is greater
        // than the distance to the furthest element in the result set W.
        while (!C.empty) {
            const c = C.pop();
            if (!c) break;
            const furthest_dist = W.first?.value ?? Infinity;

            // Stop if current candidate is farther than furthest result
            if (c.value > furthest_dist) {
                break;
            }

            const edges = layer.edges.get(c.element.index);
            if (!edges) continue;

            for (const neighbor_idx of edges) {
                if (!visited.has(neighbor_idx)) {
                    const neighbor_element = elements[neighbor_idx];
                    // Skip invalid elements or elements with different dimensions
                    if (!neighbor_element || neighbor_element.length !== q.length) continue;

                    // Skip self-connections
                    if (neighbor_idx === c.element.index) continue;

                    visited.add(neighbor_idx);
                    const dist_e = metric(neighbor_element, q);

                    const current_furthest = W.first?.value ?? Infinity;
                    if (dist_e < current_furthest || W.length < ef) {
                        C.push({
                            element: neighbor_element,
                            index: neighbor_idx,
                            distance: dist_e,
                        });
                        W.push({
                            element: neighbor_element,
                            index: neighbor_idx,
                            distance: dist_e,
                        });

                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }

        // Return sorted results for consistent entry point selection
        return W.data().sort((a, b) => a.distance - b.distance);
    }

    /**
     * Searches for the K nearest neighbors to a query element in the HNSW graph.
     *
     * Performs a multi-layer search starting from the entry point and traversing
     * each layer as entry points for the next.
     *
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @returns {Candidate<T>[]} K nearest neighbors with their distances
     */
    search(q, K) {
        // Validate K
        if (!Number.isInteger(K) || K <= 0) {
            throw new Error("HNSW: parameter 'K' must be a positive integer");
        }

        // Validate query dimensions
        if (!q || (!Array.isArray(q) && !(q instanceof Float64Array))) {
            throw new Error("HNSW: query must be an array");
        }

        const search_ef = this._ef;

        // Fallback to linear search if graph is not properly initialized
        if (this._L < 0 || !this._ep || this._elements.length === 0) {
            return this._linear_search(q, K);
        }

        let ep_indices = [...this._ep];

        // Search from top layer down to layer 1
        for (let l_c = this._L; l_c > 0; --l_c) {
            const result = this._search_layer(q, ep_indices, 1, l_c);
            if (result.length > 0) {
                ep_indices = [result[0].index];
            }
        }

        // Search layer 0 with ef candidates
        const result = this._search_layer(q, ep_indices, Math.max(search_ef, K), 0);

        // If graph search returns no results, fallback to linear search
        if (result.length === 0) {
            return this._linear_search(q, K);
        }

        // Return K closest
        return result.slice(0, K);
    }

    /**
     * Fallback linear search when graph search fails
     * @private
     * @param {T} q - Query element
     * @param {number} K - Number of nearest neighbors to return
     * @returns {Candidate<T>[]}
     */
    _linear_search(q, K) {
        const metric = this._metric;
        const elements = this._elements;
        const N = elements.length;

        if (N === 0) return [];

        /** @type {Candidate<T>[]} */
        const candidates = [];
        for (let i = 0; i < N; i++) {
            const element = elements[i];
            // Skip elements with different dimensions (can happen with inconsistent data)
            if (!element || element.length !== q.length) continue;

            candidates.push({
                element: element,
                index: i,
                distance: metric(q, element),
            });
        }

        candidates.sort((a, b) => a.distance - b.distance);
        return candidates.slice(0, K);
    }

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
    *search_iter(q, K, ef = null) {
        const search_ef = ef ?? this._ef;

        if (this._L < 0 || !this._ep) {
            return;
        }

        let ep_indices = [...this._ep];

        // Yield entry points at top layer instead of query itself
        const top_layer = this._graph.get(this._L);
        if (top_layer && this._ep && this._ep.length > 0) {
            const entry_candidates = this._ep
                .filter((idx) => this._elements[idx] !== undefined)
                .map((idx) => ({
                    element: this._elements[idx],
                    index: idx,
                    distance: this._metric(this._elements[idx], q),
                }));
            yield {
                layer: this._L,
                candidates: entry_candidates,
            };
        }

        for (let l_c = this._L; l_c > 0; --l_c) {
            const result = this._search_layer(q, ep_indices, 1, l_c);
            yield { layer: l_c, candidates: result };
            // Use closest candidate as entry point for next layer (following HNSW paper)
            ep_indices = result.length > 0 ? [result[0].index] : ep_indices;
        }

        const result = this._search_layer(q, ep_indices, Math.max(search_ef, K), 0);
        yield { layer: 0, candidates: result };
    }

    /**
     * Get the number of elements in the index.
     *
     * @returns {number} Number of elements
     */
    get size() {
        return this._elements?.length ?? 0;
    }

    /**
     * Get the number of layers in the graph.
     *
     * @returns {number} Number of layers
     */
    get num_layers() {
        return this._L + 1;
    }

    /**
     * Get an element by its index.
     *
     * @param {number} index - Element index
     * @returns {T} The element at the given index
     */
    get_element(index) {
        return this._elements[index];
    }

    /**
     * Search for nearest neighbors using an element index as the query.
     *
     * @param {number} i - Index of the query element
     * @param {number} [K=5] - Number of nearest neighbors to return
     * @returns {Candidate<T>[]} K nearest neighbors
     */
    search_by_index(i, K = 5) {
        const elements = this._elements;
        if (i < 0 || i >= elements.length) return [];

        const element = elements[i];
        if (!element) return [];

        return this.search(element, K);
    }
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
class KDTree extends KNN {
    /**
     * Generates a KD-Tree with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KD-Tree
     * @param {ParametersKDTree} [parameters={metric: euclidean}] Default is `{metric: euclidean}`
     */
    constructor(elements, parameters = { metric: euclidean, seed: 1212 }) {
        super(elements, Object.assign({ seed: 1212 }, parameters));
        /**
         * @private
         * @type {KDTreeNode<T> | KDTreeLeaf<T> | null}
         */
        this._root = this._construct(
            elements.map((element, index) => ({ index, element })),
            0,
        );
    }

    /** @returns {Metric} */
    get _metric() {
        return this._parameters.metric;
    }

    /**
     * @private
     * @param {ElementWithIndex<T>[]} elements
     * @param {number} depth - Current depth in the tree (determines splitting axis)
     * @returns {KDTreeNode<T> | KDTreeLeaf<T> | null} Root of KD-Tree.
     */
    _construct(elements, depth) {
        if (elements.length === 0) {
            return null;
        }

        if (elements.length === 1) {
            return new KDTreeLeaf(elements[0]);
        }

        const k = elements[0].element.length;
        const axis = depth % k;

        // Sort by the splitting axis and find median
        elements.sort((a, b) => a.element[axis] - b.element[axis]);
        const medianIndex = Math.floor(elements.length / 2);
        const medianPoint = elements[medianIndex];

        // Recursively build left and right subtrees
        const leftElements = elements.slice(0, medianIndex);
        const rightElements = elements.slice(medianIndex + 1);

        const left = this._construct(leftElements, depth + 1);
        const right = this._construct(rightElements, depth + 1);

        return new KDTreeNode(medianPoint, axis, left, right);
    }

    /**
     * @param {number} i
     * @param {number} k
     */
    search_by_index(i, k = 5) {
        return this.search(this._elements[i], k);
    }

    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        /** @type {Heap<{ point: ElementWithIndex<T>; distance: number }>} */
        const best = new Heap(null, (d) => d.distance, "max");

        this._search_recursive(t, k, this._root, best);

        // Convert heap to result array (closest first)
        /** @type {{ element: T; index: number; distance: number }[]} */
        const result = [];
        while (best.length > 0) {
            const item = /** @type {{ element: { point: ElementWithIndex<T>; distance: number }; value: number }} */ (
                best.pop()
            );
            result.push({
                element: item.element.point.element,
                index: item.element.point.index,
                distance: item.value,
            });
        }
        return result.reverse();
    }

    /**
     * @private
     * @param {T} target - Query element.
     * @param {number} k - Number of nearest neighbors to return.
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} node - Current node.
     * @param {Heap<{ point: ElementWithIndex<T>; distance: number }>} best - Heap of k best found so far.
     */
    _search_recursive(target, k, node, best) {
        if (node === null) return;

        if (node instanceof KDTreeLeaf) {
            const dist = this._metric(target, node.point.element);
            if (best.length < k) {
                best.push({ point: node.point, distance: dist });
            } else if (dist < (best.first?.value ?? Infinity)) {
                best.pop();
                best.push({ point: node.point, distance: dist });
            }
            return;
        }

        // Node is an internal node
        const axis = node.axis;
        const point = node.point;
        const pointValue = point.element[axis];
        const targetValue = target[axis];

        // Determine which subtree to search first
        const firstSubtree = targetValue < pointValue ? node.left : node.right;
        const secondSubtree = targetValue < pointValue ? node.right : node.left;

        // Search the nearer subtree
        this._search_recursive(target, k, firstSubtree, best);

        // Check if we need to search the other subtree
        // The hyperplane could contain closer points
        const distToHyperplane = Math.abs(targetValue - pointValue);
        const currentMaxDist = best.first?.value ?? Infinity;

        // Calculate distance to current point
        const distToPoint = this._metric(target, point.element);
        if (best.length < k) {
            best.push({ point: point, distance: distToPoint });
        } else if (distToPoint < currentMaxDist) {
            best.pop();
            best.push({ point: point, distance: distToPoint });
        }

        // Check if we need to explore the other side of the hyperplane
        if (best.length < k || distToHyperplane < (best.first?.value ?? Infinity)) {
            this._search_recursive(target, k, secondSubtree, best);
        }
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class KDTreeNode {
    /**
     * @param {ElementWithIndex<T>} point
     * @param {number} axis - The splitting axis
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} left
     * @param {KDTreeNode<T> | KDTreeLeaf<T> | null} right
     */
    constructor(point, axis, left = null, right = null) {
        this.point = point;
        this.axis = axis;
        this.left = left;
        this.right = right;
    }
}

/**
 * @private
 * @template {number[] | Float64Array} T
 */
class KDTreeLeaf {
    /**
     * @param {ElementWithIndex<T>} point
     */
    constructor(point) {
        this.point = point;
    }
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
class LSH extends KNN {
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
class NaiveKNN extends KNN {
    /**
     * Generates a KNN list with given `elements`.
     *
     * @param {T[]} elements - Elements which should be added to the KNN list
     * @param {ParametersNaiveKNN} parameters
     */
    constructor(elements, parameters = {}) {
        const params = Object.assign({ metric: euclidean, seed: 1212 }, parameters);
        super(elements, params);
        const N =
            this._elements instanceof Matrix ? /** @type {any} */ (this._elements).shape[0] : this._elements.length;
        if (this._parameters.metric === "precomputed") {
            this._D = Matrix.from(/** @type {number[][] | Float64Array[]} */ (/** @type {any} */ (this._elements)));
        } else {
            this._D = distance_matrix(
                /** @type {number[][] | Float64Array[]} */ (this._elements),
                this._parameters.metric,
            );
        }

        /** @type {Heap<{ value: number; index: number }>[]} */
        this.KNN = [];
        for (let row = 0; row < N; ++row) {
            const distances = this._D.row(row);
            /** @type {Heap<{ value: number; index: number }>} */
            const H = new Heap(null, (d) => d.value, "min");
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
     * @param {number} i
     * @param {number} k
     */
    search_by_index(i, k = 5) {
        if (this._parameters.metric === "precomputed") {
            const H = this.KNN[i];
            /** @type {{ element: T; index: number; distance: number }[]} */
            const result = [];
            const data = H.toArray(); // Get array representation
            const temp_heap = new Heap(data, (d) => d.value, "min");
            const N =
                this._elements instanceof Matrix ? /** @type {any} */ (this._elements).shape[0] : this._elements.length;
            for (let j = 0; j < Math.min(k, N); ++j) {
                const node = temp_heap.pop();
                if (!node) break;
                result.push({
                    element: /** @type {T} */ (
                        this._elements instanceof Matrix
                            ? /** @type {any} */ (this._elements).row(node.element.index)
                            : this._elements[node.element.index]
                    ),
                    index: /** @type {number} */ (node.element.index),
                    distance: /** @type {number} */ (node.value),
                });
            }
            return result;
        }
        return this.search(
            /** @type {T} */ (
                this._elements instanceof Matrix ? /** @type {any} */ (this._elements).row(i) : this._elements[i]
            ),
            k,
        );
    }

    /**
     * @param {T} t - Query element.
     * @param {number} [k=5] - Number of nearest neighbors to return. Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]} - List consists of the `k` nearest neighbors.
     */
    search(t, k = 5) {
        if (this._parameters.metric === "precomputed") {
            throw new Error("Search by query element is only possible when not using a precomputed distance matrix!");
        }
        /** @type {import("../metrics/index.js").Metric} */
        const metric = /** @type {any} */ (this._parameters.metric);

        const isMatrix = this._elements instanceof Matrix;
        const elementsAny = /** @type {any} */ (this._elements);
        const N = isMatrix ? elementsAny.shape[0] : this._elements.length;

        // Compute distances from query to ALL points
        const distances = [];
        for (let i = 0; i < N; i++) {
            const element = /** @type {T} */ (isMatrix ? elementsAny.row(i) : this._elements[i]);
            distances.push({
                element: element,
                index: i,
                distance: metric(t, element),
            });
        }

        // Sort by distance and return k nearest
        distances.sort((a, b) => a.distance - b.distance);
        return distances.slice(0, k);
    }
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
class NNDescent extends KNN {
    /**
     * @private
     * @type {KNNHeap<T>[]}
     */
    _B = [];
    /**
     * @private
     * @type {NNDescentNeighbor<T>[][]}
     */
    nn = [];

    /**
     * @param {T[]} elements - Called V in paper.
     * @param {Partial<ParametersNNDescent>} parameters
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     */
    constructor(elements, parameters = {}) {
        super(
            elements,
            /** @type {ParametersNNDescent} */ (
                Object.assign({ metric: euclidean, K: 10, rho: 1, delta: 1e-3, seed: 1212 }, parameters)
            ),
        );
        this._N = elements.length;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._sample_size = this._parameters.samples * this._parameters.rho;

        this._nndescent_elements = elements.map((e, i) => {
            return {
                value: e,
                index: i,
                flag: true,
            };
        });

        if (elements) {
            this.add(elements);
        }
    }

    /**
     * Samples Array A with sample size.
     *
     * @private
     * @template U
     * @param {U[]} A
     * @returns {U[]}
     */
    _sample(A) {
        const n = A.length;
        const sample_size = this._sample_size;
        if (sample_size > n) {
            return A;
        } else {
            const randomizer = this._randomizer;
            return randomizer.choice(A, sample_size);
        }
    }

    /**
     * @private
     * @param {KNNHeap<T>} B
     * @param {NNDescentNeighbor<T>} u
     * @returns {number}
     */
    _update(B, u) {
        if (B.set.has(u.index)) return 0;

        const worst = B.first;
        if (worst && B.length >= this._parameters.samples) {
            const dist = B._accessor(u);
            const worst_dist = B._accessor(worst.element);
            if (dist >= worst_dist) {
                return 0; // u is worse than the worst neighbor
            }
        }

        B.push(u);
        u.flag = true;
        if (B.length > this._parameters.samples) {
            B.pop();
        }
        return 1;
    }

    /**
     * @private
     * @param {(KNNHeap<T> | null)[]} B
     * @returns {NNDescentNeighbor<T>[][]}
     */
    _reverse(B) {
        const N = this._N;
        const R = new Array(N);
        for (let i = 0; i < N; i++) {
            R[i] = [];
        }
        for (let j = 0; j < N; j++) {
            const Bi = B[j];
            if (Bi) {
                const Bjdata = Bi.data();
                for (const neighbor of Bjdata) {
                    const v = neighbor.index;
                    R[v].push(neighbor);
                }
            }
        }
        return R;
    }

    /**
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements) {
        const randomizer = this._randomizer;
        const metric = this._parameters.metric;
        const K = this._parameters.samples;
        const delta = this._parameters.delta;
        const N = elements.length;
        this._N = N;
        /** @type {KNNHeap<T>[]} */
        const B = [];
        this._B = B;
        for (let i = 0; i < N; i++) {
            const e = elements[i];
            const sample = randomizer
                .choice(
                    elements.map((el, idx) => ({ el, idx })),
                    K,
                )
                .map((d) => {
                    return { index: d.idx, distance: metric(d.el, e), value: d.el };
                });
            const Bi = new KNNHeap(sample, (d) => d.distance, "max");
            B.push(Bi);
        }

        let c = Infinity;
        let old_c = -Infinity;
        while (c > delta * N * K && c !== old_c) {
            const old_ = new Array(N);
            const new_ = new Array(N);
            for (let i = 0; i < N; i++) {
                const Bi = B[i].data();
                const falseBs = Bi.filter((d) => !d.flag);
                const trueBs = this._sample(Bi.filter((d) => d.flag));
                for (const d of trueBs) {
                    d.flag = false;
                }
                old_[i] = new KNNHeap(falseBs, (d) => d.distance, "max");
                new_[i] = new KNNHeap(trueBs, (d) => d.distance, "max");
            }
            const old_reverse = this._reverse(old_);
            const new_reverse = this._reverse(new_);
            old_c = c;
            c = 0;
            for (let i = 0; i < N; i++) {
                for (const o of this._sample(old_reverse[i])) {
                    old_[i].push(o);
                }
                for (const n of this._sample(new_reverse[i])) {
                    new_[i].push(n);
                }

                const new_i = new_[i].data();
                const old_i = old_[i].data();
                const n1 = new_i.length;
                const n2 = old_i.length;
                for (let j = 0; j < n1; j++) {
                    const u1 = new_i[j];
                    const Bu1 = B[u1.index];
                    for (let k = 0; k < n1; k++) {
                        const u2 = new_i[k];
                        if (u1.index === u2.index) continue;
                        const Bu2 = B[u2.index];
                        c += this._update(Bu2, u1);
                        c += this._update(Bu1, u2);
                    }
                    for (let k = 0; k < n2; k++) {
                        const u2 = old_i[k];
                        if (u1.index === u2.index) continue;
                        const Bu2 = B[u2.index];
                        c += this._update(Bu2, u1);
                        c += this._update(Bu1, u2);
                    }
                }
            }
        }
        this.nn = this._B.map((heap) => heap.data());
        return this;
    }

    /**
     * @param {T} x
     * @param {number} [k=5] Default is `5`
     * @returns {{ element: T, index: number; distance: number }[]}
     */
    search(x, k = 5) {
        const metric = this._parameters.metric;
        const N = this._N;
        const elements = this._elements;

        if (N === 0) return [];
        const xLength = x.length;

        // Initialize candidate pool
        const visited = new Set();
        /** @type {{index: number, dist: number, evaluated: boolean}[]} */
        let pool = [];

        // Randomly pick initial candidates
        const randomizer = this._randomizer;
        for (let i = 0; i < Math.min(N, Math.max(k * 10, 50)); i++) {
            let rnd;
            do {
                rnd = randomizer.random_int % N;
            } while (visited.has(rnd));
            visited.add(rnd);

            const element = elements[rnd];
            if (!element || element.length !== xLength) continue;

            pool.push({
                index: rnd,
                dist: metric(x, element),
                evaluated: false,
            });
        }

        let searching = true;
        while (searching) {
            pool.sort((a, b) => a.dist - b.dist);
            // keep the top subset for exploration
            pool = pool.slice(0, Math.max(k * 5, 50));

            searching = false;
            for (let i = 0; i < pool.length; i++) {
                const candidate = pool[i];
                if (candidate.evaluated) continue;

                candidate.evaluated = true;
                searching = true;

                // get neighbors of this candidate from graph
                const neighbors = this.nn[candidate.index];
                if (!neighbors) continue;

                for (const neighbor of neighbors) {
                    const n_idx = neighbor.index;
                    if (!visited.has(n_idx)) {
                        visited.add(n_idx);
                        const element = elements[n_idx];
                        if (element && element.length === xLength) {
                            pool.push({
                                index: n_idx,
                                dist: metric(x, element),
                                evaluated: false,
                            });
                        }
                    }
                }
                // Don't break here! Look at more candidates per iteration for better convergence
                // break;
            }
        }

        pool.sort((a, b) => a.dist - b.dist);

        /** @type {{ element: T, index: number; distance: number }[]} */
        const result = [];
        for (let i = 0; i < Math.min(k, pool.length); i++) {
            const item = pool[i];
            result.push({
                element: elements[item.index],
                index: item.index,
                distance: item.dist,
            });
        }
        return result;
    }

    /**
     * @param {number} i
     * @param {number} [k=5] Default is `5`
     * @returns {{ element: T; index: number; distance: number }[]}
     */
    search_by_index(i, k = 5) {
        // Use regular search with the element at index i
        const elements = this._elements;
        if (i < 0 || i >= elements.length) return [];

        const element = elements[i];
        if (!element) return [];

        return this.search(element, k);
    }
}

/**
 * @template {number[] | Float64Array} U
 * @typedef {Object} HeapEntry
 * @property {NNDescentNeighbor<U>} element
 * @property {number} value
 */

/**
 * @template {number[] | Float64Array} U
 * @extends {Heap<NNDescentNeighbor<U>>}
 */
class KNNHeap extends Heap {
    /** @type {Set<number>} */
    set;

    /**
     * @param {NNDescentNeighbor<U>[]} elements
     * @param {(d: NNDescentNeighbor<U>) => number} accessor
     * @param {"max" | "min"} comparator
     */
    constructor(elements, accessor, comparator) {
        super(null, accessor, comparator);
        this.set = new Set();
        if (elements) {
            for (const element of elements) {
                this.push(element);
            }
        }
    }

    /**
     * @param {NNDescentNeighbor<U>} element
     * @returns {KNNHeap<U>}
     */
    push(element) {
        const set = this.set;
        if (set.has(element.index)) {
            return this;
        } else {
            set.add(element.index);
            super.push(element);
            return this;
        }
    }

    /** @returns {{ element: NNDescentNeighbor<U>; value: number } | null} */
    pop() {
        const result = super.pop();
        if (result?.element) {
            this.set.delete(result.element.index);
            return result;
        }
        return null;
    }

    /** @returns {NNDescentNeighbor<U>[]} */
    data() {
        return this._container.map((d) => d.element);
    }
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
class PCA extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersPCA>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        super(X, { d: 2, seed: 1212, eig_args: {} }, parameters);
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T} - The projected data.
     */
    transform() {
        const V = this.principal_components();
        const X = this.X;
        this.Y = X.dot(V);
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.transform();
    }

    /**
     * Computes the `d` principal components of Matrix `X`.
     *
     * @returns {Matrix}
     */
    principal_components() {
        if (this.V) {
            return this.V;
        }
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const X = this.X;
        const X_cent = X.sub(X.meanCols());
        const C = X_cent.transDot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration(C, d, eig_args);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {Matrix}
     */
    static principal_components(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.principal_components();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new PCA(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.transform_async();
    }
}

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersPaCMAP} from "./index.js" */

/**
 * Pairwise Controlled Manifold Approximation Projection (PaCMAP)
 *
 * A dimensionality reduction technique that uses three types of point pairs —
 * nearest neighbor (NN), mid-near (MN), and further (FP) pairs — with a
 * dynamic three-phase weight schedule and Adam optimization to preserve both
 * local and global structure.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersPaCMAP>
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
class PaCMAP extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersPaCMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            {
                n_neighbors: 10,
                MN_ratio: 0.5,
                FP_ratio: 2.0,
                d: 2,
                metric: euclidean,
                lr: 1.0,
                num_iters: [100, 100, 250],
                seed: 1212,
            },
            parameters,
        );
        [this._N, this._D] = this.X.shape;
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        if (n_neighbors >= this._N) {
            throw new Error(
                `Parameter n_neighbors (=${n_neighbors}) needs to be smaller than dataset size (N=${this._N})!`,
            );
        }
        this._iter = 0;
    }

    /**
     * Samples mid-near pairs for each point.
     * For each point i, repeats n_MN times: samples 6 random non-neighbor
     * candidates, picks the 2nd closest by high-dim distance.
     *
     * @protected
     * @param {Set<number>[]} nn_sets - Array of neighbor index sets per point
     * @param {number} n_MN - Number of mid-near pairs per point
     * @returns {Int32Array} Flat array of [i, j] pairs
     */
    _sample_mn_pairs(nn_sets, n_MN) {
        const N = this._N;
        const X = this.X;
        const randomizer = this._randomizer;
        const pairs = new Int32Array(N * n_MN * 2);
        let idx = 0;

        for (let i = 0; i < N; ++i) {
            const nn_set = nn_sets[i];
            const x_i = X.row(i);

            for (let k = 0; k < n_MN; ++k) {
                // Sample 6 random non-neighbor candidates
                /** @type {{idx: number, dist: number}[]} */
                const candidates = [];
                let attempts = 0;
                while (candidates.length < 6 && attempts < N * 2) {
                    const j = randomizer.random_int % N;
                    attempts++;
                    if (j !== i && !nn_set.has(j)) {
                        candidates.push({ idx: j, dist: euclidean_squared(x_i, X.row(j)) });
                    }
                }
                if (candidates.length < 2) {
                    // Fallback: use any available non-self point
                    const j = (i + 1) % N;
                    pairs[idx++] = i;
                    pairs[idx++] = j;
                    continue;
                }
                candidates.sort((a, b) => a.dist - b.dist);
                // Pick the 2nd closest (index 1)
                pairs[idx++] = i;
                pairs[idx++] = candidates[1].idx;
            }
        }
        return pairs.slice(0, idx);
    }

    /**
     * Samples further pairs for each point (random non-neighbors).
     *
     * @protected
     * @param {Set<number>[]} nn_sets - Array of neighbor index sets per point
     * @param {number} n_FP - Number of further pairs per point
     * @returns {Int32Array} Flat array of [i, j] pairs
     */
    _sample_fp_pairs(nn_sets, n_FP) {
        const N = this._N;
        const randomizer = this._randomizer;
        const pairs = new Int32Array(N * n_FP * 2);
        let idx = 0;

        for (let i = 0; i < N; ++i) {
            const nn_set = nn_sets[i];
            let count = 0;
            let attempts = 0;
            while (count < n_FP && attempts < N * 3) {
                const j = randomizer.random_int % N;
                attempts++;
                if (j !== i && !nn_set.has(j)) {
                    pairs[idx++] = i;
                    pairs[idx++] = j;
                    count++;
                }
            }
        }
        return pairs.slice(0, idx);
    }

    /**
     * Computes gradient coefficients and updates the gradient matrix for one pair type.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N×d gradient accumulator (modified in place)
     * @param {Int32Array} pairs - Flat [i, j, i, j, ...] pair array
     * @param {number} w - Weight for this pair type
     * @param {number} attr_num - Numerator constant for attractive (10 for NN, 10000 for MN); 0 for repulsive
     * @param {boolean} repulsive - Whether this is a repulsive pair type
     */
    _accumulate_gradients(grad_flat, pairs, w, attr_num, repulsive) {
        if (w === 0) return;
        const Y = this.Y;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;

        for (let p = 0; p < n_pairs; ++p) {
            const i = pairs[p * 2];
            const j = pairs[p * 2 + 1];
            const y_i = Y.row(i);
            const y_j = Y.row(j);

            // d_ij = 1 + ||y_i - y_j||²
            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                sq_dist += diff * diff;
            }
            const d_ij = 1 + sq_dist;

            /** @type {number} */
            let coeff;
            if (repulsive) {
                // FP loss: 1/(1+d_ij), gradient: -2/(1+d_ij)²
                coeff = (-w * 2) / (d_ij * d_ij);
            } else {
                // NN loss: d_ij/(attr_num+d_ij), gradient: 2*attr_num/(attr_num+d_ij)²
                const denom = attr_num + d_ij;
                coeff = (w * 2 * attr_num) / (denom * denom);
            }

            const base_i = i * d;
            const base_j = j * d;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                const g = coeff * diff;
                grad_flat[base_i + k] += g;
                grad_flat[base_j + k] -= g;
            }
        }
    }

    /**
     * Returns the weight schedule for the current iteration.
     *
     * @protected
     * @param {number} iter - Current iteration (0-indexed)
     * @returns {{ w_nn: number; w_mn: number; w_fp: number }}
     */
    _get_weights(iter) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const [p1, p2] = num_iters;
        if (iter < p1) {
            // Phase 1: MN weight linearly decays from 1000 to 3
            const t = iter / p1;
            return { w_nn: 2.0, w_mn: 1000.0 * (1 - t) + 3.0 * t, w_fp: 1.0 };
        } else if (iter < p1 + p2) {
            // Phase 2: fixed weights
            return { w_nn: 3.0, w_mn: 3.0, w_fp: 1.0 };
        } else {
            // Phase 3: MN disabled
            return { w_nn: 1.0, w_mn: 0.0, w_fp: 1.0 };
        }
    }

    /**
     * Applies Adam optimizer update to Y using accumulated gradients.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N×d gradient
     */
    _adam_update(grad_flat) {
        const lr = /** @type {number} */ (this.parameter("lr"));
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const beta1 = 0.9;
        const beta2 = 0.999;
        const eps = 1e-7;
        this._adam_t = (this._adam_t ?? 0) + 1;
        const t = /** @type {number} */ (this._adam_t);
        const bc1 = 1 - beta1 ** t;
        const bc2 = 1 - beta2 ** t;
        const Y = this.Y;
        const m = /** @type {Float64Array} */ (this._adam_m);
        const v = /** @type {Float64Array} */ (this._adam_v);

        for (let i = 0; i < N; ++i) {
            const base = i * d;
            const y_i = Y.row(i);
            for (let k = 0; k < d; ++k) {
                const g = grad_flat[base + k];
                m[base + k] = beta1 * m[base + k] + (1 - beta1) * g;
                v[base + k] = beta2 * v[base + k] + (1 - beta2) * g * g;
                const m_hat = m[base + k] / bc1;
                const v_hat = v[base + k] / bc2;
                y_i[k] -= lr * (m_hat / (Math.sqrt(v_hat) + eps));
            }
        }
    }

    /**
     * Initializes PaCMAP: PCA embedding, KNN pairs, MN pairs, FP pairs, Adam state.
     *
     * @returns {PaCMAP<T>}
     */
    init() {
        const X = this.X;
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const seed = /** @type {number} */ (this.parameter("seed"));
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const MN_ratio = /** @type {number} */ (this.parameter("MN_ratio"));
        const FP_ratio = /** @type {number} */ (this.parameter("FP_ratio"));

        // 1. PCA initialization scaled by 0.01 (X is always Matrix here)
        const pca_init = /** @type {Matrix} */ (PCA.transform(X, { d, seed }));
        this.Y = new Matrix(N, d, (i, j) => pca_init.entry(i, j) * 0.01);

        // 2. Build KNN graph for NN pairs
        const knn = new BallTree(X.to2dArray(), { metric, seed });
        const n_MN = Math.max(1, Math.round(n_neighbors * MN_ratio));
        const n_FP = Math.max(1, Math.round(n_neighbors * FP_ratio));
        /** @type {Set<number>[]} */
        const nn_sets = [];
        // NN pairs: flat [i, j] pairs
        const nn_pairs = new Int32Array(N * n_neighbors * 2);
        let nn_idx = 0;

        for (let i = 0; i < N; ++i) {
            const neighbors = knn.search(X.row(i), n_neighbors + 1);
            /** @type {number[]} */
            const idxs = [];
            for (const nb of neighbors) {
                if (nb.index !== i) idxs.push(nb.index);
                if (idxs.length >= n_neighbors) break;
            }
            nn_sets[i] = new Set(idxs);
            for (const j of idxs) {
                nn_pairs[nn_idx++] = i;
                nn_pairs[nn_idx++] = j;
            }
        }
        this._nn_pairs = nn_pairs.slice(0, nn_idx);

        // 3. MN pairs (mid-near sampling)
        this._mn_pairs = this._sample_mn_pairs(nn_sets, n_MN);

        // 4. FP pairs (random non-neighbors)
        this._fp_pairs = this._sample_fp_pairs(nn_sets, n_FP);

        // 5. Adam optimizer state
        this._adam_m = new Float64Array(N * d);
        this._adam_v = new Float64Array(N * d);
        this._adam_t = 0;

        this._iter = 0;
        return this;
    }

    /**
     * Performs one optimization step.
     *
     * @returns {Matrix}
     */
    next() {
        if (!this._nn_pairs) throw new Error("Call init() first!");
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const { w_nn, w_mn, w_fp } = this._get_weights(this._iter);

        const grad_flat = new Float64Array(N * d);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, 10, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._fp_pairs), w_fp, 0, true);
        this._adam_update(grad_flat);

        this._iter++;
        return this.Y;
    }

    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {T}
     */
    transform(iterations) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const total = iterations ?? num_iters.reduce((a, b) => a + b, 0);
        this.check_init();
        for (let i = 0; i < total; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {Generator<T, T, void>}
     */
    *generator(iterations) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const total = iterations ?? num_iters.reduce((a, b) => a + b, 0);
        this.check_init();
        for (let i = 0; i < total; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        return dr.transform_async();
    }
}

/** @import {InputType} from "../index.js" */
/** @import {ParametersLocalMAP} from "./index.js" */

/**
 * LocalMAP
 *
 * A variant of PaCMAP that improves local cluster separation by dynamically
 * resampling further pairs (FP) in phase 3 using nearby points in the current
 * low-dimensional embedding space, rather than randomly sampled non-neighbors.
 *
 * @class
 * @template {InputType} T
 * @extends PaCMAP<T>
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
class LocalMAP extends PaCMAP {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLocalMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        // Merge low_dist_thres into parameters before the seal in DR constructor.
        // DR.constructor does Object.seal({ ...defaults, ...parameters }), so
        // passing low_dist_thres here ensures it lands in the sealed object.
        super(X, { low_dist_thres: 10, ...parameters });
    }

    /**
     * Performs one optimization step.
     * In phase 3, resamples FP pairs from the current embedding space
     * and applies distance-based weight scaling.
     *
     * @returns {import("../matrix/index.js").Matrix}
     */
    next() {
        if (!this._nn_pairs) throw new Error("Call init() first!");
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const phase3_start = num_iters[0] + num_iters[1];

        if (this._iter < phase3_start) {
            // Phases 1 and 2: identical to PaCMAP
            return super.next();
        }

        // Phase 3: resample FP pairs from embedding neighbors
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const low_dist_thres = /** @type {number} */ (this._low_dist_thres ?? 10);
        const low_dist_thres_sq = low_dist_thres * low_dist_thres;
        const { w_nn, w_mn, w_fp } = this._get_weights(this._iter);
        const Y = this.Y;

        // Build a KNN structure on the current embedding to find nearby pairs
        const y_array = Y.to2dArray();
        const seed = /** @type {number} */ (this.parameter("seed"));
        const emb_knn = new BallTree(y_array, { metric: euclidean_squared, seed });
        const n_FP = /** @type {number} */ (this.parameter("FP_ratio")) *
            /** @type {number} */ (this.parameter("n_neighbors"));
        const n_FP_int = Math.max(1, Math.round(n_FP));

        // Build local FP pairs by finding nearby embedding neighbors
        // (these are not NN neighbors in high-dim space but nearby in low-dim)
        const nn_sets = this._nn_sets_cache;
        /** @type {Int32Array} */
        let local_fp_pairs;

        if (nn_sets) {
            const pair_buf = new Int32Array(N * n_FP_int * 2);
            let idx = 0;
            for (let i = 0; i < N; ++i) {
                const nn_set = nn_sets[i];
                // Search for nearby points in embedding space
                const neighbors = emb_knn.search(y_array[i], Math.min(n_FP_int * 3 + 1, N));
                let count = 0;
                for (const nb of neighbors) {
                    if (nb.index === i || (nn_set && nn_set.has(nb.index))) continue;
                    pair_buf[idx++] = i;
                    pair_buf[idx++] = nb.index;
                    count++;
                    if (count >= n_FP_int) break;
                }
            }
            local_fp_pairs = pair_buf.slice(0, idx);
        } else {
            local_fp_pairs = /** @type {Int32Array} */ (this._fp_pairs);
        }

        // Accumulate gradients with local weight scaling for FP pairs
        const grad_flat = new Float64Array(N * d);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, 10, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);

        // FP pairs with local distance-based weight scaling
        this._accumulate_gradients_local_fp(
            grad_flat,
            local_fp_pairs,
            w_fp,
            low_dist_thres,
            low_dist_thres_sq,
        );

        this._adam_update(grad_flat);
        this._iter++;
        return this.Y;
    }

    /**
     * Accumulates FP gradients with LocalMAP's distance-based weight scaling.
     * For pairs within low_dist_thres, scales w_fp by low_dist_thres / (2 * sqrt(d_ij)).
     *
     * @private
     * @param {Float64Array} grad_flat - Flat N×d gradient accumulator (modified in place)
     * @param {Int32Array} pairs - Flat [i, j, i, j, ...] pair array
     * @param {number} w_fp - Base FP weight
     * @param {number} low_dist_thres - Distance threshold
     * @param {number} low_dist_thres_sq - Squared distance threshold
     */
    _accumulate_gradients_local_fp(grad_flat, pairs, w_fp, low_dist_thres, low_dist_thres_sq) {
        if (w_fp === 0) return;
        const Y = this.Y;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;

        for (let p = 0; p < n_pairs; ++p) {
            const i = pairs[p * 2];
            const j = pairs[p * 2 + 1];
            const y_i = Y.row(i);
            const y_j = Y.row(j);

            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                sq_dist += diff * diff;
            }
            const d_ij = 1 + sq_dist;

            // Apply local weight scaling when pair is within distance threshold
            let w = w_fp;
            if (sq_dist < low_dist_thres_sq) {
                w *= low_dist_thres / (2 * Math.sqrt(d_ij));
            }

            // FP loss: 1/(1+d_ij), gradient: -2/(1+d_ij)²
            const coeff = (-w * 2) / (d_ij * d_ij);

            const base_i = i * d;
            const base_j = j * d;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                const g = coeff * diff;
                grad_flat[base_i + k] += g;
                grad_flat[base_j + k] -= g;
            }
        }
    }

    /**
     * Initializes LocalMAP (same as PaCMAP, but caches nn_sets for phase 3 resampling).
     *
     * @returns {LocalMAP<T>}
     */
    init() {
        super.init();
        // Cache low_dist_thres from sealed parameters (avoids type indexing issues)
        this._low_dist_thres = /** @type {number} */ (/** @type {any} */ (this._parameters)["low_dist_thres"] ?? 10);
        // Cache nn_sets for use in phase 3 FP resampling
        // We rebuild them from _nn_pairs
        const N = this._N;
        const nn_pairs = this._nn_pairs;
        if (!nn_pairs) return this;
        /** @type {Set<number>[]} */
        const nn_sets = Array.from({ length: N }, () => new Set());
        for (let p = 0; p < nn_pairs.length; p += 2) {
            nn_sets[nn_pairs[p]].add(nn_pairs[p + 1]);
        }
        this._nn_sets_cache = nn_sets;
        return this;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        return dr.transform_async();
    }
}

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
class SMACOF extends DR {
    /**
     * SMACOF for MDS.
     *
     * @param {T} X - The high-dimensional data or precomputed distance matrix.
     * @param {Partial<ParametersSMACOF>} [parameters] - Object containing parameterization.
     */
    constructor(X, parameters = {}) {
        super(X, { d: 2, metric: euclidean, seed: 1212, iterations: 300, epsilon: 1e-4 }, parameters);
    }

    /**
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        this.check_init();
        const X = this.X;
        const rows = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const iterations = /** @type {number} */ (this.parameter("iterations"));
        const epsilon = /** @type {number} */ (this.parameter("epsilon"));

        const target_distances = metric === "precomputed" ? X : distance_matrix(X, metric);

        let Z = new Matrix(rows, d, () => (this._randomizer.random - 0.5) * 2);

        // Center Z
        for (let j = 0; j < d; ++j) {
            const col = Z.col(j);
            const mean = col.reduce((a, b) => a + b, 0) / rows;
            for (let i = 0; i < rows; ++i) {
                Z.sub_entry(i, j, mean);
            }
        }

        this.Y = /** @type {Matrix} */ (Z); // Initial state

        let prev_stress = Infinity;

        if (!(iterations > 0)) {
            yield this.projection;
            return this.projection;
        }

        for (let iter = 0; iter < iterations; ++iter) {
            const B = new Matrix(rows, rows, 0);

            for (let i = 0; i < rows; ++i) {
                let bii = 0;
                const z_i = Z.row(i);
                for (let j = 0; j < rows; ++j) {
                    if (i === j) continue;
                    const z_j = Z.row(j);
                    const dist_Z = euclidean(z_i, z_j);
                    const dist_target = target_distances.entry(i, j);

                    let bij = 0;
                    if (dist_Z > 1e-12) {
                        bij = -dist_target / dist_Z;
                    }
                    B.set_entry(i, j, bij);
                    bii -= bij;
                }
                B.set_entry(i, i, bii);
            }

            // Z_new = 1/N * B(Z) * Z
            const Z_new = B.dot(Z)._apply(rows, (val, n) => val / n);

            this.Y = /** @type {Matrix} */ (Z_new);
            Z = /** @type {Matrix} */ (Z_new);

            // Calculate stress
            let stress_num = 0;
            let stress_den = 0;
            for (let i = 0; i < rows; ++i) {
                const z_i = Z.row(i);
                for (let j = i + 1; j < rows; ++j) {
                    const z_j = Z.row(j);
                    const dist_Y = euclidean(z_i, z_j);
                    const diff = target_distances.entry(i, j) - dist_Y;
                    stress_num += diff * diff;
                    stress_den += target_distances.entry(i, j) ** 2;
                }
            }
            const current_stress = Math.sqrt(stress_num / Math.max(stress_den, 1e-12));

            yield this.projection;

            if (Math.abs(prev_stress - current_stress) < epsilon) {
                break;
            }
            prev_stress = current_stress;
        }
        return this.projection;
    }

    /**
     * @returns {T}
     */
    transform() {
        const gen = this.generator();
        let res = /** @type {T} */ (this.X);
        for (const step of gen) {
            res = step;
        }
        return res;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new SMACOF(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new SMACOF(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSMACOF>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new SMACOF(X, parameters);
        return dr.transform_async();
    }
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
class ISOMAP extends DR {
    /**
     * Isometric feature mapping (ISOMAP).
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersISOMAP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X, parameters = {}) {
        /** @type {ParametersISOMAP} */
        const defaults = {
            neighbors: -Infinity,
            d: 2,
            metric: euclidean,
            seed: 1212,
            project: "MDS",
            eig_args: {},
        };
        super(X, defaults, parameters);

        this.defaults = defaults;

        if (this._parameters.neighbors === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N - 1));
        }

        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * @returns {T}
     */
    transform() {
        this.check_init();
        const X = this.X;
        const rows = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix(rows, rows, 0);
        D.shape = [rows, rows, (i, j) => (i <= j ? metric(X.row(i), X.row(j)) : D.entry(j, i))];

        /** @type {{ index: number; distance: number }[][]} */
        const kNearestNeighbors = [];
        const tree = new BallTree(X.to2dArray(), {
            metric,
            seed: /** @type {number} */ (this.parameter("seed")),
        });
        for (let i = 0; i < rows; ++i) {
            // BallTree search returns elements including the queried point itself (at distance 0).
            // Request neighbors + 1 and slice off the first one (which should be the query point).
            const neighborsList = tree.search_by_index(i, neighbors + 1);
            kNearestNeighbors.push(
                neighborsList.slice(1).map((n) => ({
                    index: n.index,
                    distance: n.distance,
                })),
            );
        }

        // ISOMAP requires an undirected/symmetric nearest neighbor graph.
        // If i is a nearest neighbor of j, then j should be connected to i as well.
        for (let i = 0; i < rows; ++i) {
            for (const neighbor of kNearestNeighbors[i]) {
                const j = neighbor.index;
                const d = neighbor.distance;
                const reciprocal_edge = kNearestNeighbors[j].find((n) => n.index === i);
                if (!reciprocal_edge) {
                    kNearestNeighbors[j].push({ index: i, distance: d });
                }
            }
        }

        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths using Dijkstra's algorithm
        // TODO: make extern
        const G = new Matrix(rows, rows, Infinity);

        for (let i = 0; i < rows; ++i) {
            G.set_entry(i, i, 0);
            const H = new Heap([{ index: i, distance: 0 }], (d) => d.distance, "min");

            while (!H.empty) {
                const item = H.pop();
                if (!item) break;

                const u = item.element.index;
                const dist_u = item.element.distance;

                if (dist_u > G.entry(i, u)) continue;

                for (const neighbor of kNearestNeighbors[u]) {
                    const v = neighbor.index;
                    const alt = dist_u + neighbor.distance;
                    if (alt < G.entry(i, v)) {
                        G.set_entry(i, v, alt);
                        H.push({ index: v, distance: alt });
                    }
                }
            }
        }

        let max_val = 0;
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < rows; j++) {
                const val = G.entry(i, j);
                if (val !== Infinity && val > max_val) max_val = val;
            }
        }
        const big_val = max_val * 10;

        const project = /** @type {"MDS" | "SMACOF"} */ (this.parameter("project"));

        if (project === "SMACOF") {
            // Apply SMACOF metric MDS to the distance matrix directly
            const D_matrix = new Matrix(rows, rows, (i, j) => {
                const val = G.entry(i, j);
                return val === Infinity ? big_val : val;
            });
            const smacof = new SMACOF(D_matrix, {
                metric: "precomputed",
                d,
                seed: this.parameter("seed"),
            });
            smacof.transform();
            this.Y = smacof.Y;
        } else {
            // "MDS" (Classical MDS) via Eigendecomposition of double-centered squared distance matrix
            const D_sq = new Matrix(rows, rows, (i, j) => {
                let val = G.entry(i, j);
                if (val === Infinity) val = big_val;
                return val * val;
            });

            const ai_ = D_sq.meanCols();
            const a_j = D_sq.meanRows();
            const a__ = D_sq.mean();
            const B = new Matrix(rows, rows, (i, j) => -0.5 * (D_sq.entry(i, j) - ai_[i] - a_j[j] + a__));

            // compute d eigenvectors
            const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
            this.Y = Matrix.from(V).transpose();
        }
        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersISOMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new ISOMAP(X, parameters);
        return dr.transform_async();
    }
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
class LDA extends DR {
    /**
     * Linear Discriminant Analysis.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLDA> & { labels: any[] | Float64Array }} parameters - Object containing parameterization of the DR method.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
    constructor(X, parameters) {
        super(X, { labels: parameters.labels, d: 2, seed: 1212, eig_args: {} }, parameters);
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T} - The projected data.
     */
    transform() {
        const X = this.X;
        const [rows, cols] = X.shape;
        const { d, labels, eig_args } = this._parameters;
        if (labels === null || labels.length !== rows) {
            throw new Error("LDA needs parameter label to every datapoint to work!");
        }

        /** @type {Record<string | number, { id: number; count: number; rows: Float64Array[] }>} */
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
        const X_mean = X.meanCols();
        const V_mean = new Matrix(label_id, cols);
        for (const label in unique_labels) {
            const V = Matrix.from(unique_labels[label].rows);
            const v_mean = V.meanCols();
            for (let j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }
        }
        // scatter_between
        let S_b = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const m = Matrix.from([v]).sub(Matrix.from([X_mean]));
            const N = unique_labels[label].count;
            S_b = S_b.add(m.transDot(m).mult(N));
        }

        // scatter_within
        let S_w = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const R = unique_labels[label].rows;
            for (let i = 0, n = unique_labels[label].count; i < n; ++i) {
                const row_v = Matrix.from([R[i]]).sub(Matrix.from([v]));
                S_w = S_w.add(row_v.transDot(row_v));
            }
        }

        const { eigenvectors: EV } = simultaneous_poweriteration(
            S_w.inverse().dot(S_b),
            d || Math.min(cols, label_id - 1),
            eig_args,
        );
        const V = Matrix.from(EV).transpose();
        this.Y = X.dot(V);

        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        // @ts-expect-error: LDA requires labels, but DR static transform doesn't
        const dr = new LDA(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        // @ts-expect-error: LDA requires labels, but DR static generator doesn't
        const dr = new LDA(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        // @ts-expect-error: LDA requires labels, but DR static transform doesn't
        const dr = new LDA(X, parameters);
        return dr.transform_async();
    }
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
class LLE extends DR {
    /**
     * Locally Linear Embedding.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLLE>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
                eig_args: {},
            },
            parameters,
        );
        if (this._parameters.neighbors === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }

        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T}
     */
    transform() {
        const X = this.X;
        const rows = this._N;
        const cols = this._D;
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        const nN = k_nearest_neighbors(X, neighbors, metric);
        const O = new Matrix(neighbors, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const nN_row = nN[row];
            const Z = new Matrix(neighbors, cols, (i, j) => X.entry(nN_row[i].j, j) - X.entry(row, j));
            const C = Z.dotTrans(Z);
            if (neighbors > cols) {
                const C_trace = neumair_sum(C.diag()) / 1000;
                for (let j = 0; j < neighbors; ++j) {
                    C.add_entry(j, j, C_trace);
                }
            }
            // reconstruct;
            let w = Matrix.solve_CG(C, O, this._randomizer);
            w = w.divide(w.sum());
            for (let j = 0; j < neighbors; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        const I = new Matrix(rows, rows, "identity");
        const IW = I.sub(W);
        const M = IW.transDot(IW);

        // M is symmetric positive semi-definite. Smallest eigenvalue is 0 (ones vector).
        // To find smallest eigenvalues of M, we can find largest of (C*I - M)
        // Upper bound for max eigenvalue: Frobenius norm or sum of absolute values
        const C = M.mean() * rows * 2; // Safe upper bound for a sparse-ish M in LLE
        const CI_M = new Matrix(rows, rows, (i, j) => (i === j ? C : 0) - M.entry(i, j));

        const { eigenvectors: V } = simultaneous_poweriteration(CI_M, d + 1, eig_args);
        // Skip the first eigenvector (the ones vector corresponding to eigenvalue C)
        this.Y = Matrix.from(V.slice(1, 1 + d)).T;

        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LLE(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LLE(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LLE(X, parameters);
        return dr.transform_async();
    }
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
class MDS extends DR {
    /**
     * Classical MDS.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersMDS>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        super(X, { d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T}
     */
    transform() {
        const X = this.X;
        const rows = X.shape[0];
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const A = metric === "precomputed" ? X : distance_matrix(X, metric);

        const D_sq = new Matrix(rows, rows, (i, j) => {
            const val = A.entry(i, j);
            return val * val;
        });

        const ai_ = D_sq.meanCols();
        const a_j = D_sq.meanRows();
        const a__ = D_sq.mean();

        this._d_X = A;
        const B = new Matrix(rows, rows, (i, j) => -0.5 * (D_sq.entry(i, j) - ai_[i] - a_j[j] + a__));

        const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
        this.Y = Matrix.from(V).transpose();

        return this.projection;
    }

    /** @returns {number} - The stress of the projection. */
    stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X;
        if (!d_X) throw new Error("First transform!");

        const d_Y = new Matrix(N, N, 0);
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
                top_sum += (d_X.entry(i, j) - d_Y.entry(i, j)) ** 2;
                bottom_sum += d_X.entry(i, j) ** 2;
            }
        }
        return Math.sqrt(top_sum / bottom_sum);
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new MDS(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new MDS(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new MDS(X, parameters);
        return dr.transform_async();
    }
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
class LSP extends DR {
    /**
     * Least Squares Projection.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLSP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                control_points: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
            },
            parameters,
        );
        if (this.parameter("neighbors") === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }
        if (this.parameter("control_points") === -Infinity) {
            this.parameter("control_points", Math.min(Math.ceil(Math.sqrt(this._N)), this._N - 1));
        }
        this._is_initialized = false;
    }

    /**
     * @returns {LSP<T>}
     */
    //	init(DR = MDS, DR_parameters = {}, KNN = BallTree) {
    init() {
        const DR = MDS;
        let DR_parameters = {};
        const KNN = BallTree;
        if (this._is_initialized) return this;
        const X = this.X;
        const N = this._N;
        const K = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const seed = /** @type {number} */ (this.parameter("seed"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        DR_parameters = Object.assign({ d, metric, seed }, DR_parameters);
        const nc = /** @type {number} */ (this.parameter("control_points"));
        const control_points = new KMedoids(X, { K: nc, metric }).get_medoids();
        const C = new Matrix(nc, N, "zeros");
        control_points.forEach((c_i, i) => {
            C.set_entry(i, c_i, 1);
        });

        const control_points_matrix = Matrix.from(control_points.map((c_i) => X.row(c_i)));
        const Y_C = new DR(control_points_matrix, DR_parameters).transform();

        const XA = X.to2dArray();
        const knn = new KNN(XA, { metric, seed });
        const L = new Matrix(N, N, "I");
        const alpha = -1 / K;
        XA.forEach((x_i, i) => {
            for (const { index: j } of knn.search(x_i, K)) {
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
     *
     * @returns {T} Returns the projection.
     */
    transform() {
        this.check_init();
        const A = this._A;
        const b = this._b;

        if (!A || !b) throw new Error("Call init() first!");
        const ATA = A.transDot(A);
        const ATb = A.transDot(b);
        this.Y = Matrix.solve_CG(ATA, ATb, this._randomizer);
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LSP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LSP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LSP(X, parameters);
        return dr.transform_async();
    }
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
class LTSA extends DR {
    /**
     * Local Tangent Space Alignment
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLTSA>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
                eig_args: {},
            },
            parameters,
        );
        if (this.parameter("neighbors") === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }

        const d = /** @type {number} */ (this.parameter("d"));
        if (this._D <= d) {
            throw new Error(
                `Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${d})!`,
            );
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimenionality `d`.
     *
     * @returns {T}
     */
    transform() {
        const X = this.X;
        const [rows, D] = X.shape;
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
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
            const C = X_i.dotTrans(X_i);
            const { eigenvectors: g } = simultaneous_poweriteration(C, d, eig_args);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            const G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            const W_i = G_i_t.transDot(G_i_t).add(1 / Math.sqrt(neighbors + 1));
            for (let i = 0; i < neighbors + 1; ++i) {
                for (let j = 0; j < neighbors + 1; ++j) {
                    B.add_entry(I_i[i], I_i[j], W_i.entry(i, j) - (i === j ? 1 : 0));
                }
            }
        }

        // 3. Aligning global coordinates
        const { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1, eig_args);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LTSA(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LTSA(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LTSA(X, parameters);
        return dr.transform_async();
    }
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
class SAMMON extends DR {
    /** @type {Matrix | undefined} */
    distance_matrix;

    /**
     * SAMMON's Mapping
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters] - Object containing parameterization of the DR
     *   method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                magic: 0.1,
                d: 2,
                metric: euclidean,
                seed: 1212,
                init_DR: "random",
                init_parameters: {},
            },
            parameters,
        );
    }

    /**
     * Initializes the projection.
     *
     * @param {Matrix | undefined} D
     * @returns {asserts D is Matrix}
     */
    init(D) {
        const N = this.X.shape[0];
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const init_DR = /** @type {AvailableInit} */ (this.parameter("init_DR"));
        const DR_parameters = this.parameter("init_parameters");
        if (init_DR === "random") {
            const randomizer = this._randomizer;
            this.Y = new Matrix(N, d, () => randomizer.random);
        } else if (init_DR === "PCA") {
            this.Y = Matrix.from(PCA.transform(this.X, /** @type {ParametersPCA} */ (DR_parameters)));
        } else if (init_DR === "MDS") {
            this.Y = Matrix.from(MDS.transform(this.X, /** @type {ParametersMDS} */ (DR_parameters)));
        } else {
            throw new Error('init_DR needs to be either "random" or a DR method!');
        }
        D = metric === "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        this.distance_matrix = D;
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {T} The projection of `X`.
     */
    transform(max_iter = 200) {
        this.check_init();
        if (!this.distance_matrix) this.init(this.distance_matrix);
        for (let j = 0; j < max_iter; ++j) {
            this._step();
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimenionality 2.
     *
     * @param {number} [max_iter=200] - Maximum number of iteration steps. Default is `200`
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection of
     *   `X`.
     */
    *generator(max_iter = 200) {
        this.check_init();
        if (!this.distance_matrix) this.init(this.distance_matrix);

        for (let j = 0; j < max_iter; ++j) {
            this._step();
            yield this.projection;
        }

        return this.projection;
    }

    _step() {
        if (!this.distance_matrix) this.init(this.distance_matrix);
        const MAGIC = /** @type {number} */ (this.parameter("magic"));
        const D = /** @type {Matrix} */ (this.distance_matrix);
        const N = this.X.shape[0];
        const d = /** @type {number} */ (this.parameter("d"));
        const Y = this.Y;

        const G = new Matrix(N, d, 0);

        const sum = new Float64Array(d);
        for (let i = 0; i < N; ++i) {
            const e1 = new Float64Array(d);
            const e2 = new Float64Array(d);
            const Yi = Y.row(i);
            for (let j = 0; j < N; ++j) {
                if (i === j) continue;
                const dX = D.entry(i, j);
                if (dX === 0) continue; // Skip identical points in high-dim

                const Yj = Y.row(j);
                const delta = new Float64Array(d);
                for (let k = 0; k < d; ++k) {
                    delta[k] = Yi[k] - Yj[k];
                }
                const dY = Math.max(euclidean(Yi, Yj), 1e-6);
                const dq = dX - dY;
                const dr = dX * dY;
                for (let k = 0; k < d; ++k) {
                    e1[k] += (delta[k] * dq) / dr;
                    e2[k] += (dq - (delta[k] ** 2 * (1 + dq / dY)) / dY) / dr;
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

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new SAMMON(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new SAMMON(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSAMMON<AvailableInit>>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new SAMMON(X, parameters);
        return dr.transform_async();
    }
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
class SQDMDS extends DR {
    /**
     * SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE
     * and UMAP.
     *
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
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
            },
            parameters,
        );

        this.init();
        if (this.parameter("metric") === "precomputed" && this.X.shape[0] !== this.X.shape[1]) {
            throw new Error("SQDMDS input data must be a square Matrix");
        }
    }

    init() {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));

        // initialize helpers.
        this._add = this.__add(d);
        this._sub_div = this.__sub_div(d);
        this._minus = this.__minus(d);
        this._mult = this.__mult(d);
        this._LR_init = Math.max(2, 0.005 * N);
        this._LR = this._LR_init;
        const decay_cte = /** @type {number} */ (this.parameter("decay_cte"));
        this._offset = -Math.exp(-1 / decay_cte);
        this._momentums = new Matrix(N, d, 0);
        this._grads = new Matrix(N, d, 0);
        this._indices = linspace(0, N - 1);
        // initialize projection.
        const R = this._randomizer;
        this.Y = new Matrix(N, d, () => R.random - 0.5);

        // preparing metric for optimization.
        const this_metric = /** @type {Metric | "precomputed"} */ (this.parameter("metric"));
        if (this_metric === "precomputed") {
            /** @type {(i: number, j: number, X: Matrix) => number} */
            this._HD_metric = (i, j, X) => X.entry(i, j);
            /** @type {(i: number, j: number, X: Matrix) => number} */
            this._HD_metric_exaggeration = (i, j, X) => X.entry(i, j) ** 2;
        } else {
            this._HD_metric = (i, j, X) => this_metric(X.row(i), X.row(j));
            if (this_metric === euclidean) {
                this._HD_metric_exaggeration = (i, j, X) => euclidean_squared(X.row(i), X.row(j));
            } else {
                this._HD_metric_exaggeration = (i, j, X) => this_metric(X.row(i), X.row(j)) ** 2;
            }
        }
        return;
    }

    /**
     * Computes the projection.
     *
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {T} The projection.
     */
    transform(iterations = 500) {
        this.check_init();
        const decay_start = /** @type {number} */ (this.parameter("decay_start"));
        this._decay_start = Math.round(decay_start * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
        }
        return this.projection;
    }

    /**
     * Computes the projection.
     *
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {Generator<T, T, void>} The intermediate steps of the projection.
     */
    *generator(iterations = 500) {
        this.check_init();
        const decay_start = /** @type {number} */ (this.parameter("decay_start"));
        this._decay_start = Math.round(decay_start * iterations);
        for (let i = 0; i < iterations; ++i) {
            this._step(i, iterations);
            yield this.projection;
        }

        return this.projection;
    }

    /**
     * Performs an optimization step.
     *
     * @private
     * @param {number} i - Acutal iteration.
     * @param {number} iterations - Number of iterations.
     */
    _step(i, iterations) {
        if (this._LR_init === undefined || this._offset === undefined) throw new Error("Call init() first!");

        const decay_start = /** @type {number} */ (this.parameter("decay_start"));
        if (i > decay_start) {
            const decay_cte = /** @type {number} */ (this.parameter("decay_cte"));
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
     *
     * @private
     * @returns {Uint32Array[]}
     */
    __quartets() {
        if (!this._indices) throw new Error("Call init() first!");
        if (this._offset === undefined) throw new Error("Call init() first!");
        const N = this._N;
        const max_N = N - (N % 4);
        const R = this._randomizer;
        const shuffled_indices = R.choice(this._indices, max_N);
        const result = [];
        for (let i = 0; i < max_N; i += 4) {
            result.push(
                Uint32Array.of(
                    shuffled_indices[i],
                    shuffled_indices[i + 1],
                    shuffled_indices[i + 2],
                    shuffled_indices[i + 3],
                ),
            );
        }
        return result;
    }

    /**
     * Computes and applies gradients, and updates momentum.
     *
     * @private
     * @param {boolean} distance_exaggeration
     */
    _nestrov_iteration(distance_exaggeration) {
        if (!this._momentums || !this._grads || this._LR === undefined) throw new Error("Call init() first!");
        const momentums = this._momentums.mult(0.99, { inline: true });
        const LR = this._LR;
        const grads = this._fill_MDS_grads(this.Y.add(momentums), this._grads, distance_exaggeration);
        const [n, d] = momentums.shape;
        for (let i = 0; i < n; ++i) {
            const g_i = grads.row(i);
            const g_i_norm = norm(g_i);
            if (g_i_norm === 0) continue;
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
     *
     * @param {Matrix} Y - The Projection.
     * @param {Matrix} grads - The gradients.
     * @param {boolean} [exaggeration=false] - Whether or not to use early exaggeration. Default is `false`
     * @param {boolean} [zero_grad=true] - Whether or not to reset the gradient in the beginning. Default is `true`
     * @returns {Matrix} The gradients.
     */
    _fill_MDS_grads(Y, grads, exaggeration = false, zero_grad = true) {
        if (!this._HD_metric || !this._HD_metric_exaggeration || !this._add) throw new Error("Call init() first!");
        if (zero_grad) {
            // compute new gradients
            grads.values.fill(0);
        }
        const add = this._add;
        const X = this.X;
        let HD_metric;
        if (exaggeration === true) {
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
     *
     * @private
     * @param {Matrix} Y - The acutal projection.
     * @param {number[]} quartet - The indices of the quartet.
     * @param {Float64Array} D_hd - The high-dimensional distances of the quartet.
     * @returns {Float64Array[]} The gradients for the quartet.
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
        const [gA1, gB1, gC1, gD1] = this._ABCD_grads(
            a,
            b,
            c,
            d,
            d_ab,
            d_ac,
            d_ad,
            d_bc,
            d_bd,
            d_cd,
            p_ab,
            sum_LD_dist,
        );
        const [gA2, gC2, gB2, gD2] = this._ABCD_grads(
            a,
            c,
            b,
            d,
            d_ac,
            d_ab,
            d_ad,
            d_bc,
            d_cd,
            d_bd,
            p_ac,
            sum_LD_dist,
        );
        const [gA3, gD3, gC3, gB3] = this._ABCD_grads(
            a,
            d,
            c,
            b,
            d_ad,
            d_ac,
            d_ab,
            d_cd,
            d_bd,
            d_bc,
            p_ad,
            sum_LD_dist,
        );
        const [gB4, gC4, gA4, gD4] = this._ABCD_grads(
            b,
            c,
            a,
            d,
            d_bc,
            d_ab,
            d_bd,
            d_ac,
            d_cd,
            d_ad,
            p_bc,
            sum_LD_dist,
        );
        const [gB5, gD5, gA5, gC5] = this._ABCD_grads(
            b,
            d,
            a,
            c,
            d_bd,
            d_ab,
            d_bc,
            d_ad,
            d_cd,
            d_ac,
            p_bd,
            sum_LD_dist,
        );
        const [gC6, gD6, gA6, gB6] = this._ABCD_grads(
            c,
            d,
            a,
            b,
            d_cd,
            d_ac,
            d_bc,
            d_ad,
            d_bd,
            d_ab,
            p_cd,
            sum_LD_dist,
        );

        if (!this._add) throw new Error("Call init() first!");
        const add = this._add;
        const gA = add(gA1, gA2, gA3, gA4, gA5, gA6);
        const gB = add(gB1, gB2, gB3, gB4, gB5, gB6);
        const gC = add(gC1, gC2, gC3, gC4, gC5, gC6);
        const gD = add(gD1, gD2, gD3, gD4, gD5, gD6);

        return [gA, gB, gC, gD];
    }

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
    _ABCD_grads(a, b, c, d, d_ab, d_ac, d_ad, d_bc, d_bd, d_cd, p_ab, sum_LD_dist) {
        if (!this._minus || !this._add || !this._mult || !this._sub_div) throw new Error("Call init() first!");
        const ratio = d_ab / sum_LD_dist;
        const twice_ratio = 2 * ((p_ab - ratio) / sum_LD_dist);
        const minus = this._minus;
        const add = this._add;
        const mult = this._mult;
        const sub_div = this._sub_div;
        // no side effects because sub_div creates new arrays, and the inline functions work on this new created arrays.
        const gA = mult(
            minus(mult(add(sub_div(a, b, d_ab), sub_div(a, c, d_ac), sub_div(a, d, d_ad)), ratio), sub_div(a, b, d_ab)),
            twice_ratio,
        );
        const gB = mult(
            minus(mult(add(sub_div(b, a, d_ab), sub_div(b, c, d_bc), sub_div(b, d, d_bd)), ratio), sub_div(b, a, d_ab)),
            twice_ratio,
        );
        const gC = mult(add(sub_div(c, a, d_ac), sub_div(c, b, d_bc), sub_div(c, d, d_cd)), ratio * twice_ratio);
        const gD = mult(add(sub_div(d, a, d_ad), sub_div(d, b, d_bd), sub_div(d, c, d_cd)), ratio * twice_ratio);
        return [gA, gB, gC, gD];
    }

    /**
     * Inline!
     *
     * @param {number} d
     */
    __minus(d) {
        return /** @type {(a: Float64Array, b: Float64Array) => Float64Array} */ (a, b) => {
            for (let i = 0; i < d; ++i) {
                a[i] -= b[i];
            }
            return a;
        };
    }

    /**
     * Inline!
     *
     * @param {number} d
     */
    __add(d) {
        return /** @type {(...summands: Float64Array[]) => Float64Array} */ (...summands) => {
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
     *
     * @param {number} d
     */
    __mult(d) {
        return /** @type {(a: Float64Array, v: number) => Float64Array} */ (a, v) => {
            for (let i = 0; i < d; ++i) {
                a[i] *= v;
            }
            return a;
        };
    }

    /**
     * Creates a new array `(x - y) / div`.
     *
     * @param {number} d
     */
    __sub_div(d) {
        return /** @type {(x: Float64Array, y: Float64Array, div: number) => Float64Array} */ (x, y, div) => {
            return Float64Array.from({ length: d }, (_, i) => (x[i] - y[i]) / div);
        };
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new SQDMDS(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new SQDMDS(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersSQDMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new SQDMDS(X, parameters);
        return dr.transform_async();
    }
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
class TopoMap extends DR {
    /**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTopoMap>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { metric: euclidean, seed: 1212 }, parameters);
        [this._N, this._D] = this.X.shape;
        this._distance_matrix = new Matrix(this._N, this._N, -1);
    }

    /**
     * @private
     * @param {number} i
     * @param {number} j
     * @param {import("../metrics/index.js").Metric} metric
     * @returns {number}
     */
    __lazy_distance_matrix(i, j, metric) {
        const D = this._distance_matrix;
        const X = this.X;
        const D_ij = D.entry(i, j);
        if (D_ij === -1 && i !== j) {
            const dist = metric(X.row(i), X.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
            return dist;
        }
        return i === j ? 0 : D_ij;
    }

    /**
     * Computes the minimum spanning tree, using a given metric
     *
     * @private
     * @param {import("../metrics/index.js").Metric} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    _make_minimum_spanning_tree(metric = euclidean) {
        const N = this._N;
        const X = [...this.X];

        this._disjoint_set = new DisjointSet(X);
        const disjoint_set = this._disjoint_set;
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
            if (!set_u || !set_v) throw new Error("Should not happen!");
            if (set_u !== set_v) {
                F.push([u, v, w]);
                disjoint_set.union(set_u, set_v);
            }
        }

        return F.sort((a, b) => a[2] - b[2]);
    }

    /** Initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree. */
    init() {
        const { metric } = this._parameters;
        this.Y = new Matrix(this._N, 2, 0);
        this._Emst = this._make_minimum_spanning_tree(metric);
        this._is_initialized = true;
        return this;
    }

    /**
     * Returns true if Point C is left of line AB.
     *
     * @private
     * @param {Float64Array} PointA - Point A of line AB
     * @param {Float64Array} PointB - Point B of line AB
     * @param {Float64Array} PointC - Point C
     * @returns {boolean}
     */
    __hull_cross([ax, ay], [bx, by], [sx, sy]) {
        return (bx - ax) * (sy - ay) - (by - ay) * (sx - ax) <= 0;
    }

    /**
     * Computes the convex hull of the set of Points S
     *
     * @private
     * @param {Float64Array[]} S - Set of Points.
     * @returns {Float64Array[]} Convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     */
    __hull(S) {
        const points = S.sort(([x1, y1], [x2, y2]) => y1 - y2 || x1 - x2);
        const N = points.length;
        if (N <= 2) return points;

        const lower = [];
        for (let i = 0; i < N; ++i) {
            while (
                lower.length >= 2 &&
                this.__hull_cross(lower[lower.length - 2], lower[lower.length - 1], points[i])
            ) {
                lower.pop();
            }
            lower.push(points[i]);
        }
        const upper = [];
        for (let i = N - 1; i >= 0; --i) {
            while (
                upper.length >= 2 &&
                this.__hull_cross(upper[upper.length - 2], upper[upper.length - 1], points[i])
            ) {
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
     *
     * @private
     * @param {Float64Array} PointA
     * @param {Float64Array} PointB
     * @returns {{ sin: number; cos: number }} Object containing the sinus- and cosinus-values for a rotation.
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
     * @param {Float64Array[]} hull
     * @param {Float64Array} p
     * @param {boolean} topEdge
     * @returns {{ sin: number; cos: number; tx: number; ty: number }}
     */
    __align_hull(hull, p, topEdge) {
        let v = -1;
        /** @type {number} */
        let d2 = -Infinity;
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

        const v1 = hull[v];
        let v2;
        if (topEdge) {
            v2 = hull[(v + 1) % hull.length];
        } else {
            v2 = hull[(v - 1 + hull.length) % hull.length];
        }

        /** @type {{ sin?: number; cos?: number; tx: number; ty: number }} */
        const transformation = {
            tx: -v1[0],
            ty: -v1[1],
        };

        if (hull.length >= 2) {
            const { sin, cos } = this.__findAngle(v1, v2);
            transformation.sin = sin;
            transformation.cos = cos;
        } else {
            transformation.sin = 0;
            transformation.cos = 1;
        }

        return /** @type {{ sin: number; cos: number; tx: number; ty: number }} */ (transformation);
    }

    /**
     * @private
     * @param {Float64Array} Point - The point which should get transformed.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} Transformation - Contains the values for
     *   translation and rotation.
     */
    __transform([px, py], { tx, ty, sin, cos }) {
        const x = px + tx;
        const y = py + ty;
        const xx = x * cos - y * sin;
        const yy = x * sin + y * cos;
        return [xx, yy];
    }

    /**
     * Calls `__transform` for each point in Set C
     *
     * @private
     * @param {Float64Array[]} C - Set of points.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} t - Transform object.
     * @param {number} yOffset - Value to offset set C.
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
     * @param {Float64Array} root_u - Root of component u
     * @param {Float64Array} root_v - Root of component v
     * @param {Float64Array} p_u - Point u
     * @param {Float64Array} p_v - Point v
     * @param {number} w - Edge weight w
     * @param {DisjointSet<Float64Array>} components - The disjoint set containing the components
     */
    __align_components(root_u, root_v, p_u, p_v, w, components) {
        if (!components) throw new Error("components not provided!");
        const u_children = components.get_children(root_u);
        const v_children = components.get_children(root_v);
        if (!u_children || !v_children) throw new Error("should not happen!");

        const points_u = [...u_children];
        const points_v = [...v_children];

        const hull_u = this.__hull(points_u);
        const hull_v = this.__hull(points_v);

        const t_u = this.__align_hull(hull_u, p_u, false);
        const t_v = this.__align_hull(hull_v, p_v, true);

        this.__transform_component(points_u, t_u, 0);
        this.__transform_component(points_v, t_v, w);
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {T}
     */
    transform() {
        if (!this._is_initialized) this.init();
        if (!this._Emst) throw new Error("Call init() first!");
        const Emst = this._Emst;
        const Y = this.Y.to2dArray();
        /** @type {DisjointSet<Float64Array>} */
        const components = new DisjointSet(
            Y,
            // Y.map((y, i) => {
            //     y.i = i;
            //     return y;
            // }),
        );

        for (const [u, v, w] of Emst) {
            const p_u = Y[u];
            const p_v = Y[v];
            const component_u = components.find(p_u);
            const component_v = components.find(p_v);
            if (!component_u || !component_v) throw new Error("Should not happen!");
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, p_u, p_v, w, components);
            components.union(component_u, component_v);
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {Generator<T, T, void>}
     */
    *generator() {
        if (!this._is_initialized) this.init();
        if (!this._Emst) throw new Error("call init() first!");
        const Emst = this._Emst;
        const Y = this.Y.to2dArray();
        const components = new DisjointSet(
            Y,
            // Y.map((y, i) => {
            //     y.i = i;
            //     return y;
            // }),
        );

        for (const [u, v, w] of Emst) {
            const p_u = Y[u];
            const p_v = Y[v];
            const component_u = components.find(p_u);
            const component_v = components.find(p_v);
            if (!component_u || !component_v) throw new Error("should not happen!");
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, p_u, p_v, w, components);
            components.union(component_u, component_v);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new TopoMap(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new TopoMap(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new TopoMap(X, parameters);
        return dr.transform_async();
    }
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
class TriMap extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTriMap>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                weight_adj: 500,
                n_inliers: 10,
                n_outliers: 5,
                n_random: 5,
                d: 2,
                metric: euclidean,
                tol: 1e-8,
                seed: 1212,
            },
            parameters,
        );
    }

    /**
     * @param {Matrix | null} [pca=null] - Initial Embedding (if null then PCA gets used). Default is `null`
     * @param {import("../knn/KNN.js").KNN<number[] | Float64Array, any> | null} [knn=null] - KNN Object (if null then BallTree gets used). Default is `null`
     */
    init(pca = null, knn = null) {
        const X = this.X;
        const N = X.shape[0];
        //const c = /** @type {number} */ (this._parameters.c);
        const d = /** @type {number} */ (this._parameters.d);
        const metric = /** @type {Metric} */ (this._parameters.metric);
        const seed = /** @type {number} */ (this._parameters.seed);
        this.n_inliers = /** @type {number} */ (this._parameters.n_inliers);
        this.n_outliers = /** @type {number} */ (this._parameters.n_outliers);
        this.n_random = /** @type {number} */ (this._parameters.n_random);
        this.Y = pca ?? PCA.transform(X, { d, seed });
        this.knn = knn ?? new BallTree(X.to2dArray(), { metric, seed });
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
     *
     * @param {number} n_inliers
     * @param {number} n_outliers
     * @param {number} n_random
     */
    _generate_triplets(n_inliers, n_outliers, n_random) {
        const metric = /** @type {Metric} */ (this._parameters.metric);
        const weight_adj = /** @type {number} */ (this._parameters.weight_adj);
        const X = this.X;
        const N = X.shape[0];
        const knn = this.knn;
        if (!knn) throw new Error("Call init() first!");
        const n_extra = Math.min(n_inliers + 20, N);
        const nbrs = new Matrix(N, n_extra);
        const knn_distances = new Matrix(N, n_extra);
        for (let i = 0; i < N; ++i) {
            const results = knn
                .search(X.row(i), n_extra + 1)
                .filter((d) => d.distance !== 0)
                .sort((a, b) => a.distance - b.distance);

            results.forEach((d, j) => {
                if (j < n_extra) {
                    nbrs.set_entry(i, j, d.index);
                    knn_distances.set_entry(i, j, d.distance);
                }
            });
        }
        // scale parameter
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            sig[i] = Math.max(
                (knn_distances.entry(i, 3) + knn_distances.entry(i, 4) + knn_distances.entry(i, 5)) / 3,
                1e-10,
            );
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
            if (Number.isNaN(weights[i])) {
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
     *
     * @private
     * @param {Matrix} knn_distances - Matrix of pairwise knn distances
     * @param {Float64Array} sig - Scaling factor for the distances
     * @param {Matrix} nbrs - Nearest neighbors
     * @returns {Matrix} Pairwise similarity matrix
     */
    _find_p(knn_distances, sig, nbrs) {
        const [N, n_neighbors] = knn_distances.shape;
        return new Matrix(N, n_neighbors, (i, j) => {
            return Math.exp(-(knn_distances.entry(i, j) ** 2 / sig[i] / sig[nbrs.entry(i, j)]));
        });
    }

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
    _sample_knn_triplets(P, nbrs, n_inliers, n_outliers) {
        const N = nbrs.shape[0];
        const triplets_list = [];
        for (let i = 0; i < N; ++i) {
            const sort_indices = this.__argsort(P.row(i));
            for (let j = 0; j < n_inliers; ++j) {
                const sim = nbrs.entry(i, sort_indices[sort_indices[j] === i ? j + 1 : j]);
                const rejects = [i, ...Array.from(sort_indices.slice(0, j + 2)).map((idx) => nbrs.entry(i, idx))];
                const samples = this._rejection_sample(n_outliers, N, rejects);
                for (let k = 0; k < samples.length; ++k) {
                    const out = samples[k];
                    triplets_list.push([i, sim, out]);
                }
            }
        }
        const triplets = new Matrix(triplets_list.length, 3);
        for (let t = 0; t < triplets_list.length; ++t) {
            triplets.set_entry(t, 0, triplets_list[t][0]);
            triplets.set_entry(t, 1, triplets_list[t][1]);
            triplets.set_entry(t, 2, triplets_list[t][2]);
        }
        return triplets;
    }

    /**
     * Should do the same as np.argsort()
     *
     * @private
     * @param {Float64Array | number[]} A
     */
    __argsort(A) {
        return linspace(0, A.length - 1).sort((i, j) => A[j] - A[i]);
    }

    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are
     * in the {@link rejects}.
     *
     * @private
     * @param {number} n_samples
     * @param {number} max_int
     * @param {number[]} rejects
     */
    _rejection_sample(n_samples, max_int, rejects) {
        const randomizer = this._randomizer;
        const interval = linspace(0, max_int - 1).filter((d) => rejects.indexOf(d) < 0);
        return randomizer.choice(interval, Math.min(n_samples, interval.length));
    }

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
     *
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    _sample_random_triplets(X, n_random, sig) {
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const randomizer = this._randomizer;
        const N = X.shape[0];
        const random_triplets = new Matrix(N * n_random, 3);
        const random_weights = new Float64Array(N * n_random);
        for (let i = 0; i < N; ++i) {
            const n_i = i * n_random;
            const indices = Array.from({ length: N }, (_, idx) => idx).filter((idx) => idx !== i);
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
                random_weights[index] = 0.1 * (p_sim / p_out);
            }
        }
        return {
            random_triplets: random_triplets,
            random_weights: random_weights,
        };
    }

    /**
     * Computes the gradient for updating the embedding.
     *
     * @param {Matrix} Y - The embedding
     */
    _grad(Y) {
        const n_inliers = this.n_inliers;
        const n_outliers = this.n_outliers;
        const triplets = this.triplets;
        const weights = this.weights;
        if (!triplets || n_inliers === undefined || n_outliers === undefined || !weights)
            throw new Error("Call init() first!");
        const [N, dim] = Y.shape;
        const n_triplets = triplets.shape[0];
        const grad = new Matrix(N, dim, 0);
        const y_ij = new Float64Array(dim);
        const y_ik = new Float64Array(dim);
        let d_ij = 1;
        let d_ik = 1;
        let n_viol = 0;
        let loss = 0;
        const n_knn_triplets = N * n_inliers * n_outliers;

        for (let t = 0; t < n_triplets; ++t) {
            const [i, j, k] = triplets.row(t);
            // update y_ij, y_ik, d_ij, d_ik
            if (t % n_outliers === 0 || t >= n_knn_triplets) {
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
            const w = weights[t] / (d_ij + d_ik) ** 2;
            for (let d = 0; d < dim; ++d) {
                const gs = y_ij[d] * d_ik * w;
                const go = y_ik[d] * d_ij * w;
                grad.add_entry(i, d, gs - go);
                grad.sub_entry(j, d, gs);
                grad.add_entry(k, d, go);
            }
        }
        return { grad, loss, n_viol };
    }

    /**
     * @param {number} max_iteration
     * @returns {T}
     */
    transform(max_iteration = 800) {
        this.check_init();
        for (let iter = 0; iter < max_iteration; ++iter) {
            this._next(iter);
        }
        return this.projection;
    }

    /**
     * @param {number} max_iteration
     * @returns {Generator<T, T, void>}
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
     *
     * @private
     * @param {number} iter
     */
    _next(iter) {
        const gamma = iter > 250 ? 0.5 : 0.3;
        const old_C = this.C;
        const vel = this.vel;
        if (!vel || old_C === undefined || this.lr === undefined) throw new Error("Call init() first!");
        const Y = this.Y.add(vel.mult(gamma));
        const { grad, loss } = this._grad(Y);
        this.C = loss;
        this.Y = this._update_embedding(Y, iter, grad);
        const tol = /** @type {number} */ (this.parameter("tol"));
        this.lr *= old_C > loss + tol ? 1.01 : 0.9;
        return this.Y;
    }

    /**
     * Updates the embedding.
     *
     * @private
     * @param {Matrix} Y
     * @param {number} iter
     * @param {Matrix} grad
     */
    _update_embedding(Y, iter, grad) {
        const [N, dim] = Y.shape;
        const gamma = iter > 250 ? 0.8 : 0.5; // moment parameter
        const min_gain = 0.01;
        const gain = this.gain;
        const vel = this.vel;
        const lr = this.lr;
        if (!vel || !gain || lr === undefined) throw new Error("Call init() first!");
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const new_gain =
                    Math.sign(vel.entry(i, d)) !== Math.sign(grad.entry(i, d))
                        ? gain.entry(i, d) + 0.2
                        : Math.max(gain.entry(i, d) * 0.8, min_gain);
                gain.set_entry(i, d, new_gain);
                vel.set_entry(i, d, gamma * vel.entry(i, d) - lr * gain.entry(i, d) * grad.entry(i, d));
                Y.set_entry(i, d, Y.entry(i, d) + vel.entry(i, d));
            }
        }
        return Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new TriMap(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new TriMap(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new TriMap(X, parameters);
        return dr.transform_async();
    }
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
class TSNE extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTSNE>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            {
                perplexity: 50,
                epsilon: 10,
                d: 2,
                metric: euclidean_squared,
                seed: 1212,
            },
            parameters,
        );
        [this._N, this._D] = this.X.shape;
        this._iter = 0;
        const d = /** @type {number} */ (this.parameter("d"));
        this.Y = new Matrix(this._N, d, () => this._randomizer.gauss_random() * 1e-4);
    }

    init() {
        // init
        const perplexity = /** @type {number} */ (this.parameter("perplexity"));
        const Htarget = Math.log(perplexity);
        const N = this._N;
        const D = this._D;
        const metric = /** @type {Metric | "precomputed"} */ (this._parameters.metric);
        const X = this.X;
        let Delta;
        if (metric === "precomputed") {
            Delta = Matrix.from(X);
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

        const P = new Matrix(N, N, 0);

        this._ystep = new Matrix(N, D, 0);
        this._gains = new Matrix(N, D, 1);

        // search for fitting sigma
        const tol = 1e-4;
        const maxtries = 50;
        for (let i = 0; i < N; ++i) {
            const dist_i = Delta.row(i);
            const prow = P.row(i);
            let betamin = -Infinity;
            let betamax = Infinity;
            let beta = 1;
            let cnt = maxtries;
            let done = false;
            let psum = 0;

            while (!done && cnt--) {
                // compute entropy and kernel row with beta precision
                psum = 0;
                let dp_sum = 0;
                for (let j = 0; j < N; ++j) {
                    const dist = dist_i[j];
                    const pj = i !== j ? Math.exp(-dist * beta) : 0;
                    dp_sum += dist * pj;
                    prow[j] = pj;
                    psum += pj;
                }
                // compute entropy
                const H = psum > 0 ? Math.log(psum) + (beta * dp_sum) / psum : 0;
                if (H > Htarget) {
                    betamin = beta;
                    beta = betamax === Infinity ? beta * 2 : (beta + betamax) / 2;
                } else {
                    betamax = beta;
                    beta = betamin === -Infinity ? beta / 2 : (beta + betamin) / 2;
                }
                done = Math.abs(H - Htarget) < tol;
            }
            // normalize p
            for (let j = 0; j < N; ++j) {
                prow[j] /= psum;
            }
        }

        // compute probabilities
        const N2 = N * 2;
        for (let i = 0; i < N; ++i) {
            for (let j = i; j < N; ++j) {
                const p = Math.max((P.entry(i, j) + P.entry(j, i)) / N2, 1e-100);
                P.set_entry(i, j, p);
                P.set_entry(j, i, p);
            }
        }
        this._P = P;
        return this;
    }

    /**
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {T} The projection.
     */
    transform(iterations = 500) {
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     * @param {number} [iterations=500] - Number of iterations. Default is `500`
     * @returns {Generator<T, T, void>} - The projection.
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
     * Performs a optimization step
     *
     * @private
     * @returns {Matrix}
     */
    next() {
        const iter = ++this._iter;
        if (!this._P || !this._ystep || !this._gains) throw new Error("Call init() first!");
        const P = this._P;
        const ystep = this._ystep;
        const gains = this._gains;
        const N = this._N;
        const dim = /** @type {number} */ (this._parameters.d);
        const epsilon = /** @type {number} */ (this._parameters.epsilon);
        const Y = this.Y;

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
                    grad.add_entry(i, d, premult * (Y.entry(i, d) - Y.entry(j, d)));
                }
            }
        }

        // perform gradient step
        const ymean = new Float64Array(dim);
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

                Y.add_entry(i, d, newsid);
                ymean[d] += Y.entry(i, d);
            }
        }

        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                Y.sub_entry(i, d, ymean[d] / N);
            }
        }

        return this.Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new TSNE(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new TSNE(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTSNE>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new TSNE(X, parameters);
        return dr.transform_async();
    }
}

/**
 * @template {Float64Array | number[]} T
 * @category Optimization
 * @param {(d: T) => number} f
 * @param {T} x0
 * @param {number} [max_iter=300] Default is `300`
 * @returns {T}
 * @see http://optimization-js.github.io/optimization-js/optimization.js.html#line438
 */
function powell(f, x0, max_iter = 300) {
    const epsilon = 1e-2;
    const n = x0.length;
    let alpha = 1e-3;
    let pfx = 10000;
    const x = /** @type {T} */ (x0.slice());
    let fx = f(x);
    let convergence = false;

    while (max_iter-- >= 0 && !convergence) {
        convergence = true;
        for (let i = 0; i < n; ++i) {
            x[i] += 1e-6;
            const fxi = f(x);
            x[i] -= 1e-6;
            const dx = (fxi - fx) / 1e-6;
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
class UMAP extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersUMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            {
                n_neighbors: 15,
                local_connectivity: 1,
                min_dist: 1,
                d: 2,
                metric: euclidean,
                seed: 1212,
                _spread: 1,
                _set_op_mix_ratio: 1,
                _repulsion_strength: 1,
                _negative_sample_rate: 5,
                _n_epochs: 350,
                _initial_alpha: 1,
            },
            parameters,
        );
        [this._N, this._D] = this.X.shape;
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const local_connectivity = /** @type {number} */ (this.parameter("local_connectivity"));
        const d = /** @type {number} */ (this.parameter("d"));
        /* let n_neighbors = Math.min(this._N - 1, parameters.n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(this.parameter("local_connectivity"), n_neighbors - 1)); */
        if (n_neighbors > this._N) {
            throw new Error(
                `Parameter n_neighbors (=${n_neighbors}) needs to be smaller than dataset size (N=${this._N})!`,
            );
        }
        if (local_connectivity > n_neighbors) {
            throw new Error(
                `Parameter local_connectivity (=${local_connectivity}) needs to be smaller than parameter n_neighbors (=${n_neighbors})`,
            );
        }
        this._iter = 0;
        const randomizer = this._randomizer;
        this.Y = new Matrix(this._N, d, () => randomizer.random);
    }

    /**
     * @private
     * @param {number} spread
     * @param {number} min_dist
     * @returns {number[]}
     */
    _find_ab_params(spread, min_dist) {
        /** @type {(x: number, a: number, b: number) => number} */
        const curve = (x, a, b) => 1 / (1 + a * x ** (2 * b));
        const xv = linspace(0, spread * 3, 300);
        const yv = linspace(0, spread * 3, 300);

        for (let i = 0, n = xv.length; i < n; ++i) {
            const xv_i = xv[i];
            yv[i] = xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread);
        }

        /** @type {(p: [number, number]) => number} */
        const err = (p) => {
            const error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map((e) => e * e)));
        };

        return powell(err, [1, 1]);
    }

    /**
     * @private
     * @param {{ element: Float64Array; index: number; distance: number }[][]} distances
     * @param {number[]} sigmas
     * @param {number[]} rhos
     * @returns {{ element: Float64Array; index: number; distance: number }[][]}
     */
    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            const rho = rhos[i];
            const curr_dist = distances[i];
            for (let j = 0, m = curr_dist.length; j < m; ++j) {
                const v = curr_dist[j].distance - rho;
                curr_dist[j].distance = v > 0 ? Math.exp(-v / sigmas[i]) : 1.0;
            }
        }
        return distances;
    }

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
    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = /** @type {number} */ (this._parameters.local_connectivity);
        const metric = /** @type {Metric | "precomputed"} */ (this._parameters.metric);
        const target = Math.log2(k);
        const rhos = [];
        const sigmas = [];
        const X = this.X;
        const N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());

        /** @type {{ element: Float64Array; index: number; distance: number }[][]} */
        const distances = [];
        if (metric === "precomputed" || knn instanceof NaiveKNN) {
            for (let i = 0; i < N; ++i) {
                distances.push(knn.search_by_index(i, k).reverse());
            }
        } else {
            for (const x_i of X) {
                distances.push(knn.search(x_i, k).reverse());
            }
        }

        const index = Math.floor(local_connectivity);
        const interpolation = local_connectivity - index;
        for (let i = 0; i < N; ++i) {
            let lo = 0;
            let hi = Infinity;
            let mid = 1;
            let rho = 0;

            const search_result = distances[i];
            const non_zero_dist = search_result.filter((d) => d.distance > 0);
            const non_zero_dist_length = non_zero_dist.length;
            if (non_zero_dist_length >= local_connectivity) {
                if (index > 0) {
                    rho = non_zero_dist[index - 1].distance;
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rho += interpolation * (non_zero_dist[index].distance - non_zero_dist[index - 1].distance);
                    }
                } else {
                    rho = interpolation * non_zero_dist[0].distance;
                }
            } else if (non_zero_dist_length > 0) {
                rho = non_zero_dist[non_zero_dist_length - 1].distance;
            }
            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    const d = search_result[j].distance - rho;
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

            //let mean_d = null;
            if (rho > 0) {
                const mean_ithd = search_result.reduce((a, b) => a + b.distance, 0) / search_result.length;
                if (mid < MIN_K_DIST_SCALE * mean_ithd) {
                    mid = MIN_K_DIST_SCALE * mean_ithd;
                }
            } else {
                const mean_d = distances.reduce(
                    (acc, res) => acc + res.reduce((a, b) => a + b.distance, 0) / res.length,
                    0,
                );
                if (mid < MIN_K_DIST_SCALE * mean_d) {
                    mid = MIN_K_DIST_SCALE * mean_d;
                }
            }
            rhos[i] = rho;
            sigmas[i] = mid;
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
     * @param {number} n_neighbors
     * @returns {Matrix}
     */
    _fuzzy_simplicial_set(X, n_neighbors) {
        const N = X.shape[0];
        const metric = /** @type {Metric | "precomputed"} */ (this._parameters.metric);
        const _set_op_mix_ratio = /** @type {number} */ (this._parameters._set_op_mix_ratio);

        const knn =
            metric === "precomputed"
                ? new NaiveKNN(X.to2dArray(), {
                      metric: "precomputed",
                      seed: /** @type {number} */ (this._parameters.seed),
                  })
                : new BallTree(X.to2dArray(), {
                      metric,
                      seed: /** @type {number} */ (this._parameters.seed),
                  });
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        const result = new Matrix(N, N, "zeros");
        for (let i = 0; i < N; ++i) {
            const distances_i = distances[i];
            for (let j = 0; j < distances_i.length; ++j) {
                result.set_entry(i, distances_i[j].index, distances_i[j].distance);
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
     * @param {number} n_epochs
     * @returns {Float32Array}
     */
    _make_epochs_per_sample(n_epochs) {
        if (!this._weights) throw new Error("Call init() first!");
        const weights = this._weights;
        const result = new Float32Array(weights.length).fill(-1);
        const weight_scl = n_epochs / max(weights);
        weights.forEach((w, i) => {
            const sample = w * weight_scl;
            if (sample > 0) result[i] = Math.round(n_epochs / sample);
        });
        return result;
    }

    /**
     * @private
     * @param {Matrix} graph
     * @returns {{ rows: number[]; cols: number[]; data: number[] }}
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
     *
     * @returns {UMAP<T>}
     */
    init() {
        const _spread = /** @type {number} */ (this._parameters._spread);
        const min_dist = /** @type {number} */ (this._parameters.min_dist);
        const n_neighbors = /** @type {number} */ (this._parameters.n_neighbors);
        const _n_epochs = /** @type {number} */ (this._parameters._n_epochs);
        const _negative_sample_rate = /** @type {number} */ (this._parameters._negative_sample_rate);
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
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {T}
     */
    transform(iterations = 350) {
        if (this.parameter("_n_epochs") !== iterations) {
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
     * @param {number} [iterations=350] - Number of iterations. Default is `350`
     * @returns {Generator<T, T, void>}
     */
    *generator(iterations = 350) {
        if (this.parameter("_n_epochs") !== iterations) {
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
     * @param {number} x
     * @returns {number}
     */
    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

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
    _optimize_layout(head_embedding, tail_embedding, head, tail) {
        const randomizer = this._randomizer;
        const _repulsion_strength = /** @type {number} */ (this.parameter("_repulsion_strength"));
        const dim = /** @type {number} */ (this.parameter("d"));
        const {
            _alpha: alpha,
            _a: a,
            _b: b,
            _epochs_per_sample: epochs_per_sample,
            _epochs_per_negative_sample: epochs_per_negative_sample,
            _epoch_of_next_negative_sample: epoch_of_next_negative_sample,
            _epoch_of_next_sample: epoch_of_next_sample,
            _clip: clip,
        } = this;
        if (
            alpha === undefined ||
            a === undefined ||
            b === undefined ||
            epochs_per_sample === undefined ||
            epochs_per_negative_sample === undefined ||
            epoch_of_next_negative_sample === undefined ||
            epoch_of_next_sample === undefined ||
            clip === undefined
        ) {
            throw new Error("call init() first!");
        }
        const tail_length = tail.length;

        for (let i = 0, n = epochs_per_sample.length; i < n; ++i) {
            if (epoch_of_next_sample[i] <= this._iter) {
                const j = head[i];
                const k = tail[i];
                const current = head_embedding.row(j);
                const other = tail_embedding.row(k);
                const dist = euclidean_squared(current, other);
                if (dist > 0) {
                    const grad_coeff = (-2 * a * b * dist ** (b - 1)) / (a * dist ** b + 1);
                    for (let d = 0; d < dim; ++d) {
                        const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                        current[d] += grad_d;
                        other[d] -= grad_d;
                    }
                }
                epoch_of_next_sample[i] += epochs_per_sample[i];
                const n_neg_samples = (this._iter - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i];
                for (let p = 0; p < n_neg_samples; ++p) {
                    const k = randomizer.random_int % tail_length;
                    const other = tail_embedding.row(tail[k]);
                    const dist = euclidean_squared(current, other);
                    if (dist > 0) {
                        const grad_coeff = (2 * _repulsion_strength * b) / ((0.01 + dist) * (a * dist ** b + 1));
                        for (let d = 0; d < dim; ++d) {
                            const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                            current[d] += grad_d;
                            other[d] -= grad_d;
                        }
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
        if (!this._head || !this._tail) throw new Error("Call init() first!");
        const iter = ++this._iter;
        const Y = this.Y;
        const _initial_alpha = /** @type {number} */ (this._parameters._initial_alpha);
        const _n_epochs = /** @type {number} */ (this._parameters._n_epochs);
        this._alpha = _initial_alpha * (1 - iter / _n_epochs);
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new UMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new UMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersUMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new UMAP(X, parameters);
        return dr.transform_async();
    }
}

const version = pkg.version;

export { Annoy, BallTree, CURE, DisjointSet, FASTMAP, HNSW, Heap, HierarchicalClustering, ISOMAP, KDTree, KMeans, KMedoids, LDA, LLE, LSH, LSP, LTSA, LocalMAP, MDS, Matrix, MeanShift, NNDescent, NaiveKNN, OPTICS, PCA, PaCMAP, Randomizer, SAMMON, SMACOF, SQDMDS, TSNE, TopoMap, TriMap, UMAP, XMeans, bray_curtis, canberra, chebyshev, cosine, distance_matrix, euclidean, euclidean_squared, goodman_kruskal, hamming, haversine, inner_product, jaccard, k_nearest_neighbors, kahan_sum, linspace, manhattan, max, min, neumair_sum, norm, normalize, powell, qr, qr_householder, simultaneous_poweriteration, sokal_michener, version, wasserstein, yule };
//# sourceMappingURL=druid.js.map
