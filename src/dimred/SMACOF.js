import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { wasmSmacofStep } from "../wasm/index.js";
import { WASM_MIN_ROWS } from "../wasm/thresholds.js";
import { DR } from "./DR.js";

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
