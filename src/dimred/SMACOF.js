import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { WASM_MIN_ROWS } from "../wasm/thresholds.js";
import { DR } from "./DR.js";
import { wasmSmacofStep } from "./SMACOF.wasm.js";

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
export class SMACOF extends DR {
    /**
     * @private
     * @type {Matrix | undefined}
     */
    _target_distances;

    /** @protected */
    get _wasm_session_keys() {
        return ["smacof"];
    }

    /**
     * SMACOF for MDS.
     *
     * @param {T} X - The high-dimensional data or precomputed distance matrix.
     * @param {Partial<ParametersSMACOF>} [parameters] - Object containing parameterization.
     */
    constructor(X, parameters = {}) {
        super(
            X,
            { d: 2, metric: euclidean, seed: 1212, iterations: 300, epsilon: 1e-4, init_Y: null },
            parameters,
        );
    }

    /**
     * Computes the target distances and the starting configuration: `init_Y` if given, otherwise
     * random and centered. Runs automatically before the first step; {@link SMACOF#generator} and
     * {@link SMACOF#transform} continue from the current `Y` afterwards.
     *
     * @returns {this}
     */
    init() {
        const rows = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const init_Y = /** @type {Matrix | number[][] | Float64Array[] | null} */ (this.parameter("init_Y"));

        this._target_distances = metric === "precomputed" ? Matrix.from(this.X) : distance_matrix(this.X, metric);

        if (init_Y) {
            const Y = Matrix.from(init_Y);
            if (Y.shape[0] !== rows || Y.shape[1] !== d) {
                throw new Error(`init_Y needs shape ${rows} × ${d}, but has shape ${Y.shape[0]} × ${Y.shape[1]}!`);
            }
            this.Y = Y.clone();
        } else {
            const Z = new Matrix(rows, d, () => (this._randomizer.random - 0.5) * 2);
            // Center Z
            for (let j = 0; j < d; ++j) {
                const col = Z.col(j);
                const mean = col.reduce((a, b) => a + b, 0) / rows;
                for (let i = 0; i < rows; ++i) {
                    Z.sub_entry(i, j, mean);
                }
            }
            this.Y = Z;
        }
        this._is_initialized = true;
        return this;
    }

    /**
     * Runs SMACOF from the current configuration (the start from {@link SMACOF#init} on the first
     * call), so calling it again continues where the last run stopped.
     *
     * @param {number} [iterations] - Maximum number of steps. Defaults to the `iterations` parameter.
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator(iterations = /** @type {number} */ (this.parameter("iterations"))) {
        this.check_init();
        const rows = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const epsilon = /** @type {number} */ (this.parameter("epsilon"));
        const target_distances = /** @type {Matrix} */ (this._target_distances);

        let Z = this.Y;
        let prev_stress = Infinity;

        if (!(iterations > 0)) {
            yield this.projection;
            return this.projection;
        }

        try {
            for (let iter = 0; iter < iterations; ++iter) {
                const Z_new = new Matrix(rows, d);
                let current_stress = null;

                if (rows >= WASM_MIN_ROWS) {
                    current_stress = wasmSmacofStep(Z.values, target_distances.values, Z_new.values, rows, d);
                }

                if (current_stress === null) {
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
                    const Z_new_js = B.dot(Z)._apply(rows, (val, n) => val / n);
                    Z_new.values.set(Z_new_js.values);

                    // Calculate stress
                    let stress_num = 0;
                    let stress_den = 0;
                    for (let i = 0; i < rows; ++i) {
                        const z_i = Z_new.row(i);
                        for (let j = i + 1; j < rows; ++j) {
                            const z_j = Z_new.row(j);
                            const dist_Y = euclidean(z_i, z_j);
                            const diff = target_distances.entry(i, j) - dist_Y;
                            stress_num += diff * diff;
                            stress_den += target_distances.entry(i, j) ** 2;
                        }
                    }
                    current_stress = Math.sqrt(stress_num / Math.max(stress_den, 1e-12));
                }

                this.Y = Z_new;
                Z = Z_new;

                yield this.projection;

                if (Math.abs(prev_stress - current_stress) < epsilon) {
                    break;
                }
                prev_stress = current_stress;
            }
            return this.projection;
        } finally {
            this._release_wasm();
        }
    }

    /**
     * @param {number} [iterations] - Maximum number of steps. Defaults to the `iterations` parameter.
     * @returns {T}
     */
    transform(iterations) {
        const gen = this.generator(iterations);
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
