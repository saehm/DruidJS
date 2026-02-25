import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { DR } from "./DR.js";

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
export class LDA extends DR {
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
