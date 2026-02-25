import { linspace } from "../matrix/index.js";

/**
 * @category Utils
 * @class
 */
export class Randomizer {
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
