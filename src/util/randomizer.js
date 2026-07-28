/**
 * Expands a single 32-bit seed into a well-mixed sequence, used to fill the generator state.
 *
 * A generator seeded with a small integer (1, 2, 42 …) must not start in a state correlated with
 * that integer, or nearby seeds would produce correlated streams. splitmix32 decorrelates them.
 *
 * @private
 * @param {number} seed
 * @returns {() => number} Successive 32-bit values.
 */
function splitmix32(seed) {
    let a = seed | 0;
    return () => {
        a = (a + 0x9e3779b9) | 0;
        let t = a ^ (a >>> 16);
        t = Math.imul(t, 0x21f0aaad);
        t = t ^ (t >>> 15);
        t = Math.imul(t, 0x735a2d97);
        t = t ^ (t >>> 15);
        return t >>> 0;
    };
}

/**
 * Seeded pseudo-random number generator.
 *
 * Implements **sfc32** (Small Fast Counting, Doty-Humphrey), a 128-bit-state counter-based
 * generator that passes TestU01 BigCrush. It is preferred here over the more familiar Mersenne
 * Twister, which needs a 2.5 KB state, is markedly slower to seed, and whose GF(2)-linear output
 * fails BigCrush's matrix-rank and linear-complexity tests.
 *
 * Every operation is 32-bit integer arithmetic (`^ << >>> +` and `Math.imul`), all of which
 * ECMAScript specifies exactly. The stream is therefore identical on every engine for a given seed —
 * unlike floating point transcendentals, which are only "implementation-approximated". See
 * {@link DR} for what that does and does not guarantee about the output of an algorithm.
 *
 * @category Utils
 * @class
 * @see {@link https://pracrand.sourceforge.net/|PractRand, where sfc32 originates}
 * @example
 * const R = new Randomizer(1212);
 * R.random;      // float in [0, 1)
 * R.random_int;  // uint32
 */
export class Randomizer {
    /** @type {number} */
    _a = 0;
    /** @type {number} */
    _b = 0;
    /** @type {number} */
    _c = 0;
    /** @type {number} */
    _d = 0;
    /** @type {number} */
    _seed;
    /** @type {number | null} */
    _val = null;

    /**
     * @param {number} [_seed=new Date().getTime()] - The seed for the random number generator. If `_seed == null` then
     *   the actual time gets used as seed. Default is `new Date().getTime()`
     */
    constructor(_seed) {
        this._seed = _seed ?? Date.now();
        this.seed = this._seed;
    }

    /** @type {number} seed */
    set seed(_seed) {
        this._seed = _seed;
        const mix = splitmix32(_seed);
        this._a = mix();
        this._b = mix();
        this._c = mix();
        this._d = mix();
        this._val = null;
        // sfc32 needs a short run-in before its output is well distributed.
        for (let i = 0; i < 12; ++i) this.random_int;
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
        const a = this._a | 0;
        const b = this._b | 0;
        const c = this._c | 0;
        const d = this._d | 0;

        const t = (((a + b) | 0) + d) | 0;
        this._d = (d + 1) | 0;
        this._a = b ^ (b >>> 9);
        this._b = (c + (c << 3)) | 0;
        this._c = ((((c << 21) | (c >>> 11)) + t) | 0) >>> 0;

        return t >>> 0;
    }

    /**
     * Returns a normally distributed number with mean 0 and standard deviation 1.
     *
     * Uses the Marsaglia polar method, which yields two values per iteration; the spare is cached
     * and returned by the following call.
     *
     * @returns {number} A standard normal variate.
     */
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
     * Returns `n` samples drawn from `A` without replacement.
     *
     * Uses a partial Fisher-Yates shuffle over a scratch index array, which runs in O(n): `A` itself
     * is never touched, and exactly `n` random values are consumed. Removing each picked index with
     * `Array.prototype.splice` instead would shift the tail on every draw and make the call O(n²) —
     * measurably worse from a few hundred elements upward, and 27× slower at n = 4000.
     *
     * @template T Returns samples from an input Matrix or Array.
     * @param {T[]} A - The input Matrix or Array.
     * @param {number} n - The number of samples.
     * @returns {T[]} A random selection form `A` of `n` samples.
     */
    choice(A, n) {
        if (!Array.isArray(A)) throw new Error("A must be an Array!");
        const rows = A.length;
        // Guard explicitly: `new Array(n)` below turns a non-numeric `n` into a one-element array
        // of `undefined` rather than failing, which would hand the caller silent garbage.
        if (!Number.isInteger(n) || n < 0) {
            throw new Error("n must be a non-negative integer!");
        }
        if (n > rows) {
            throw new Error("n bigger than A!");
        }
        const index_list = new Array(rows);
        for (let i = 0; i < rows; ++i) index_list[i] = i;

        const sample = new Array(n);
        for (let i = 0; i < n; ++i) {
            const j = i + (this.random_int % (rows - i));
            const swap = index_list[i];
            index_list[i] = index_list[j];
            index_list[j] = swap;
            sample[i] = A[index_list[i]];
        }
        return sample;
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
