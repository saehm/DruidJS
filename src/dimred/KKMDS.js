import { StressMDS, WEIGHTS_ELASTIC } from "./StressMDS.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersKKMDS} from "./index.js" */

/**
 * Kamada-Kawai Multidimensional Scaling (KKMDS)
 *
 * {@link StressMDS} fixed at `weights: -2`, weighting each pair by `1 / d_ij²`. That is the
 * Kamada-Kawai energy from graph drawing, known in the MDS literature as *elastic scaling*.
 *
 * Its classic use is laying out a graph from its shortest-path distances — pass those as a
 * `"precomputed"` matrix, as {@link MINFOTree} does.
 *
 * @class
 * @template {InputType} T
 * @extends StressMDS<T>
 * @category Dimensionality Reduction
 * @see {@link https://doi.org/10.1016/0020-0190(89)90102-6}
 * @see {@link StressMDS} for other weightings.
 */
export class KKMDS extends StressMDS {
    /**
     * Kamada-Kawai weighted MDS.
     *
     * @param {T} X - The high-dimensional data, or a precomputed distance matrix.
     * @param {Partial<ParametersKKMDS>} [parameters] - Object containing parameterization of the DR
     *   method. `weights` is not among them; it is fixed at `-2`.
     */
    constructor(X, parameters = {}) {
        super(X, { ...parameters, weights: WEIGHTS_ELASTIC });
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new KKMDS(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new KKMDS(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersKKMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new KKMDS(X, parameters);
        return dr.transform_async();
    }
}
