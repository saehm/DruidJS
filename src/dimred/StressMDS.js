import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { WASM_MIN_ROWS } from "../wasm/thresholds.js";
import { DR } from "./DR.js";
import { PCA } from "./PCA.js";
import { wasmStressMDSEnergy, wasmStressMDSStep } from "./StressMDS.wasm.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersStressMDS, WeightSpec} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/** @typedef {"MDS" | "PCA" | "random"} AvailableInit */

/** Raw stress: every pair counts equally. Recovers the objective {@link SMACOF} minimises. */
export const WEIGHTS_UNIFORM = 0;
/** Sammon stress. Recovers the objective {@link SAMMON} minimises, up to a constant factor. */
export const WEIGHTS_SAMMON = -1;
/** Elastic scaling, also known as the Kamada-Kawai energy. See {@link KKMDS}. */
export const WEIGHTS_ELASTIC = -2;

/**
 * Weighted metric MDS (stress majorization family)
 *
 * Minimises `σ(Y) = Σ_{i<j} w_ij ⋅ (‖y_i - y_j‖ - d_ij)²` for a weighting you choose. An exponent
 * `q` gives `w_ij = d_ij^q`, recovering the classical family: `0` is raw stress ({@link SMACOF}'s
 * objective), `-1` Sammon stress ({@link SAMMON}'s), `-2` elastic scaling ({@link KKMDS}). A matrix
 * or a function works too, where a zero weight drops the pair from the objective — the usual way to
 * express missing dissimilarities.
 *
 * It does not replace {@link SMACOF} or {@link SAMMON}: those minimise the same objectives with
 * their own historical algorithms and reach different local minima.
 *
 * Optimised by Jacobi-preconditioned gradient descent with a backtracking line search, warm-started
 * from classical MDS. O(N²·d) per iteration.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersStressMDS>
 * @category Dimensionality Reduction
 * @see {@link KKMDS} for the `weights: -2` preset.
 * @example
 * // Sammon stress, converged further than SAMMON's own optimizer takes it
 * const Y = new StressMDS(X, { weights: -1 }).transform();
 * @example
 * // Ignore pairs whose dissimilarity was never measured
 * const W = new Matrix(N, N, (i, j) => (observed(i, j) ? 1 : 0));
 * const Y = new StressMDS(D, { metric: "precomputed", weights: W }).transform();
 */
export class StressMDS extends DR {
    /**
     * The target distances. Named apart from the base class's `_D`, which is the *input*
     * dimensionality of `X`, not a matrix.
     *
     * @type {Matrix | undefined}
     */
    _target_distances;

    /** @type {Matrix | undefined} */
    _weights;

    /** @type {number} */
    _energy = Infinity;

    /** @protected */
    get _wasm_session_keys() {
        return ["stress_mds"];
    }

    /**
     * Weighted metric MDS.
     *
     * @param {T} X - The high-dimensional data, or a precomputed distance matrix.
     * @param {Partial<ParametersStressMDS>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        super(
            X,
            {
                d: 2,
                metric: euclidean,
                weights: WEIGHTS_ELASTIC,
                iterations: 300,
                epsilon: 1e-6,
                learning_rate: 0.1,
                init_DR: "MDS",
                seed: 1212,
                eig_args: {},
            },
            parameters,
        );
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }
    }

    /**
     * The weighted stress `σ(Y)` of the current embedding.
     *
     * @returns {number}
     */
    get energy() {
        this.check_init();
        return /** @type {number} */ (this._energy);
    }

    /**
     * The weight matrix actually in use, whatever form `weights` was given in.
     *
     * @returns {Matrix}
     */
    get weights() {
        this.check_init();
        return /** @type {Matrix} */ (this._weights);
    }

    /**
     * Materialises `weights` into an N ⨯ N matrix.
     *
     * Non-positive and non-finite weights are normalised to exactly 0, which is the kernel's "skip
     * this pair" signal. That is what makes a zero distance harmless under a negative exponent: `0^-2`
     * is `Infinity`, and coincident rows are ordinary in real data, so the pair is dropped rather
     * than given infinite stiffness.
     *
     * @private
     * @param {Matrix} D
     * @returns {Matrix}
     */
    _build_weights(D) {
        const N = this._N;
        const spec = /** @type {WeightSpec} */ (this.parameter("weights"));

        /** @type {(i: number, j: number) => number} */
        let weight_of;
        if (typeof spec === "number") {
            const q = spec;
            weight_of = q === 0 ? () => 1 : (i, j) => D.entry(i, j) ** q;
        } else if (typeof spec === "function") {
            weight_of = (i, j) => spec(D.entry(i, j), i, j);
        } else if (spec instanceof Matrix || Array.isArray(spec)) {
            const W = spec instanceof Matrix ? spec : Matrix.from(spec);
            const [rows, cols] = W.shape;
            if (rows !== N || cols !== N) {
                throw new Error(`weights matrix is ${rows}⨯${cols}, but X has ${N} rows!`);
            }
            weight_of = (i, j) => W.entry(i, j);
        } else {
            throw new Error("weights must be a number, a matrix, or a function!");
        }

        return new Matrix(N, N, (i, j) => {
            if (i === j) return 0;
            const w = weight_of(i, j);
            return Number.isFinite(w) && w > 0 ? w : 0;
        });
    }

    /**
     * Classical MDS on the target distances, scaled by the square roots of the eigenvalues.
     *
     * The scaling matters: the descent starts from here, and an unscaled eigenvector basis would put
     * the configuration at a completely different magnitude from `d_ij`, so the first steps would be
     * spent on scale rather than structure.
     *
     * @private
     * @param {Matrix} D
     * @param {number} d
     * @returns {Matrix}
     */
    _classical_mds(D, d) {
        const N = this._N;
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const D_sq = new Matrix(N, N, (i, j) => D.entry(i, j) ** 2);
        const ai_ = D_sq.meanCols();
        const a_j = D_sq.meanRows();
        const a__ = D_sq.mean();
        const B = new Matrix(N, N, (i, j) => -0.5 * (D_sq.entry(i, j) - ai_[i] - a_j[j] + a__));

        const { eigenvalues, eigenvectors } = simultaneous_poweriteration(B, d, eig_args);
        const scales = Array.from({ length: d }, (_, k) => Math.sqrt(Math.max(eigenvalues[k], 0)));
        return new Matrix(N, d, (i, k) => eigenvectors[k][i] * scales[k]);
    }

    /** Computes the target distances, the weights, and the starting configuration. */
    init() {
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const init_DR = /** @type {AvailableInit} */ (this.parameter("init_DR"));

        this._target_distances = metric === "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);
        const D = this._target_distances;
        this._weights = this._build_weights(D);

        if (init_DR === "random") {
            const randomizer = this._randomizer;
            // Scaled to the target distances so the descent does not start by fixing magnitude.
            let scale = 0;
            for (let i = 0; i < this._N; ++i) {
                for (let j = i + 1; j < this._N; ++j) {
                    const value = D.entry(i, j);
                    if (Number.isFinite(value) && value > scale) scale = value;
                }
            }
            if (scale === 0) scale = 1;
            this.Y = new Matrix(this._N, d, () => (randomizer.random - 0.5) * scale);
        } else if (init_DR === "PCA") {
            if (metric === "precomputed") {
                throw new Error('init_DR "PCA" needs the original data, not a precomputed distance matrix!');
            }
            this.Y = Matrix.from(PCA.transform(this.X, { d, seed: this.parameter("seed") }));
        } else if (init_DR === "MDS") {
            this.Y = this._classical_mds(D, d);
        } else {
            throw new Error('init_DR needs to be "MDS", "PCA" or "random"!');
        }

        this._energy = this._compute_energy(this.Y, D, this._weights);
        this._is_initialized = true;
        return this;
    }

    /**
     * @private
     * @param {Matrix} Y
     * @param {Matrix} D
     * @param {Matrix} W
     * @returns {number}
     */
    _compute_energy(Y, D, W) {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));

        if (N >= WASM_MIN_ROWS) {
            const energy = wasmStressMDSEnergy(Y.values, D.values, W.values, N, d);
            if (energy !== null) return energy;
        }

        let energy = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                const w_ij = W.entry(i, j);
                if (!(w_ij > 0)) continue;
                let sum_sq = 0;
                for (let k = 0; k < d; ++k) {
                    const diff = Y.entry(i, k) - Y.entry(j, k);
                    sum_sq += diff * diff;
                }
                const residual = Math.sqrt(sum_sq) - D.entry(i, j);
                energy += w_ij * residual * residual;
            }
        }
        return energy;
    }

    /**
     * One preconditioned gradient step with a backtracking line search, in JS.
     *
     * Mirrors `stress_mds_step_f64` in `StressMDS.as.ts` loop for loop, so the two paths agree to
     * floating point tolerance. See that file for why the gradient is divided by the weighted degree.
     *
     * @private
     * @param {Matrix} Y - Updated in place on an accepted step.
     * @param {Matrix} D
     * @param {Matrix} W
     * @param {Float64Array} G - Gradient scratch, `N * d`.
     * @param {number} step
     * @param {number} energy_current
     * @returns {{ energy: number; step: number; accepted: boolean }}
     */
    _step_js(Y, D, W, G, step, energy_current) {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        G.fill(0);

        for (let i = 0; i < N; ++i) {
            let degree = 0;
            for (let j = 0; j < N; ++j) {
                if (i === j) continue;
                const w_ij = W.entry(i, j);
                if (!(w_ij > 0)) continue;
                degree += w_ij;
                let sum_sq = 0;
                for (let k = 0; k < d; ++k) {
                    const diff = Y.entry(i, k) - Y.entry(j, k);
                    sum_sq += diff * diff;
                }
                const dist = Math.sqrt(sum_sq);
                if (dist < 1e-12) continue;
                const coeff = (2 * w_ij * (dist - D.entry(i, j))) / dist;
                for (let k = 0; k < d; ++k) {
                    G[i * d + k] += coeff * (Y.entry(i, k) - Y.entry(j, k));
                }
            }
            if (degree > 0) {
                const inv_degree = 1 / degree;
                for (let k = 0; k < d; ++k) G[i * d + k] *= inv_degree;
            }
        }

        let g_norm = 0;
        for (let idx = 0; idx < G.length; ++idx) g_norm += G[idx] * G[idx];
        g_norm = Math.sqrt(g_norm);
        if (g_norm < 1e-12) return { energy: energy_current, step, accepted: false };

        const trial = new Matrix(N, d);
        for (let attempt = 0; attempt < 12; ++attempt) {
            for (let idx = 0; idx < G.length; ++idx) {
                trial.values[idx] = Y.values[idx] - step * G[idx];
            }
            const energy = this._compute_energy(trial, D, W);
            if (energy < energy_current) {
                Y.values.set(trial.values);
                return { energy, step: step * 1.2, accepted: true };
            }
            step *= 0.5;
        }
        return { energy: energy_current, step, accepted: false };
    }

    /**
     * Computes the projection.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        this.check_init();
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const iterations = /** @type {number} */ (this.parameter("iterations"));
        const epsilon = /** @type {number} */ (this.parameter("epsilon"));
        const learning_rate = /** @type {number} */ (this.parameter("learning_rate"));
        const D = /** @type {Matrix} */ (this._target_distances);
        const W = /** @type {Matrix} */ (this._weights);
        const Y = this.Y;

        // The warm start is a meaningful configuration in its own right, and yielding it here is
        // what guarantees the generator produces at least one value: the line search can reject the
        // very first step — on an input classical MDS already solves, it usually does — and a
        // consumer that saw nothing at all would have no projection to fall back on.
        yield this.projection;

        if (N < 2 || !(iterations > 0)) {
            return this.projection;
        }

        // Dimensionless: the preconditioned gradient is already in configuration units, so
        // `learning_rate` means the same thing whether the distances run 0-1 or 0-1000, and across
        // weightings. A value of 1 takes the full preconditioned step.
        let step = learning_rate;
        let energy = /** @type {number} */ (this._energy);
        const G = new Float64Array(N * d);

        try {
            for (let iter = 0; iter < iterations; ++iter) {
                /** @type {{ energy: number; step: number; accepted: boolean } | null} */
                let result = null;
                if (N >= WASM_MIN_ROWS) {
                    result = wasmStressMDSStep(Y.values, D.values, W.values, N, d, step, energy);
                }
                if (result === null) {
                    result = this._step_js(Y, D, W, G, step, energy);
                }

                step = result.step;
                if (!result.accepted) break;

                const improvement = (energy - result.energy) / Math.max(energy, 1e-12);
                energy = result.energy;
                this._energy = energy;
                yield this.projection;

                if (improvement < epsilon) break;
            }
            return this.projection;
        } finally {
            this._release_wasm();
        }
    }

    /**
     * Computes the projection.
     *
     * @returns {T}
     */
    transform() {
        // Drained for its effect on `Y`, and the projection read back afterwards rather than taken
        // from the last yielded value — the loop body must never decide what `transform` returns.
        for (const _ of this.generator()) {
            // no-op
        }
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new StressMDS(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new StressMDS(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersStressMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new StressMDS(X, parameters);
        return dr.transform_async();
    }
}
