/**
 * @class
 * @alias TriMap
 * @extends DR
 */
export class TriMap extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.weight_adj = 500] - scaling factor.
     * @param {Number} [parameters.c = 5] - number of triplets multiplier.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.tol = 1e-8] -
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X: Matrix, parameters: {
        weight_adj?: number;
        c?: number;
        d?: number;
        tol?: number;
        metric?: Function;
        seed?: number;
    });
    /**
     *
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used).
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used).
     */
    init(pca?: Matrix, knn?: KNN): TriMap;
    n_inliers: number;
    n_outliers: number;
    n_random: number;
    Y: any;
    knn: any;
    triplets: Matrix;
    weights: Float64Array;
    lr: number;
    C: number;
    vel: Matrix;
    gain: Matrix;
    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {Number} n_inliers
     * @param {Number} n_outliers
     * @param {Number} n_random
     */
    _generate_triplets(n_inliers: number, n_outliers: number, n_random: number): {
        triplets: Matrix;
        weights: Float64Array;
    };
    /**
     * Calculates the similarity matrix P
     * @private
     * @param {Matrix} knn_distances - matrix of pairwise knn distances
     * @param {Float64Array} sig - scaling factor for the distances
     * @param {Matrix} nbrs - nearest neighbors
     * @returns {Matrix} pairwise similarity matrix
     */
    private _find_p;
    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix {@link P}. Row i corresponds to the i-th point.
     * @param {Number} n_inliers - Number of inlier points.
     * @param {Number} n_outliers - Number of outlier points.
     *
     */
    private _sample_knn_triplets;
    /**
     * Should do the same as np.argsort()
     * @private
     * @param {Array} A
     */
    private __argsort;
    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are in the {@link rejects}.
     * @private
     * @param {*} n_samples
     * @param {*} max_int
     * @param {*} rejects
     */
    private _rejection_sample;
    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - scaling factor for the distances.
     */
    private _find_weights;
    /**
     * Sample uniformly ranom triplets
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {Number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    private _sample_random_triplets;
    /**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */
    _grad(Y: Matrix): {
        grad: Matrix;
        loss: number;
        n_viol: number;
    };
    /**
     * Does the iteration step.
     * @private
     * @param {Number} iter
     */
    private _next;
    /**
     * Updates the embedding.
     * @private
     * @param {Matrix} Y
     * @param {Number} iter
     * @param {Matrix} grad
     */
    private _update_embedding;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=TriMap.d.ts.map