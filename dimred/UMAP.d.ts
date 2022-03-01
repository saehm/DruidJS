/**
 * @class
 * @alias UMAP
 * @extends DR
 */
export class UMAP extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias UMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [parameters.local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [parameters.min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.
     * @param {Number} [parameters._spread = 1] - The effective scale of embedded points. (In combination with {@link parameters.min_dist})
     * @param {Number} [parameters._set_op_mix_ratio = 1] - Interpolate between union and intersection.
     * @param {Number} [parameters._repulsion_strength = 1]  - Weighting applied to negative samples.
     * @param {Number} [parameters._negative_sample_rate = 5] - The number of negative samples per positive sample.
     * @param {Number} [parameters._n_epochs = 350] - The number of training epochs.
     * @param {Number} [parameter._initial_alpha = 1] - The initial learning rate for the optimization.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {UMAP}
     */
    constructor(X: Matrix, parameters: {
        n_neighbors?: number;
        local_connectivity?: number;
        min_dist?: number;
        d?: number;
        metric?: Function;
        _spread?: number;
        _set_op_mix_ratio?: number;
        _repulsion_strength?: number;
        _negative_sample_rate?: number;
        _n_epochs?: number;
    });
    _iter: number;
    Y: Matrix;
    /**
     * @private
     * @param {Number} spread
     * @param {Number} min_dist
     * @returns {Array}
     */
    private _find_ab_params;
    /**
     * @private
     * @param {Array<Array>} distances
     * @param {Array<Number>} sigmas
     * @param {Array<Number>} rhos
     * @returns {Array}
     */
    private _compute_membership_strengths;
    /**
     * @private
     * @param {KNN|BallTree} knn
     * @param {Number} k
     * @returns {Object}
     */
    private _smooth_knn_dist;
    /**
     * @private
     * @param {Matrix} X
     * @param {Number} n_neighbors
     * @returns {Matrix}
     */
    private _fuzzy_simplicial_set;
    /**
     * @private
     * @param {Number} n_epochs
     * @returns {Array}
     */
    private _make_epochs_per_sample;
    /**
     * @private
     * @param {Matrix} graph
     * @returns {Object}
     */
    private _tocoo;
    /**
     * Computes all necessary
     * @returns {UMAP}
     */
    init(): UMAP;
    _a: any;
    _b: any;
    _graph: Matrix;
    _head: any;
    _tail: any;
    _weights: any;
    _epochs_per_sample: any[];
    _epochs_per_negative_sample: number[];
    _epoch_of_next_sample: any[];
    _epoch_of_next_negative_sample: number[];
    graph(): {
        cols: any;
        rows: any;
        weights: any;
    };
    /**
     * @private
     * @param {Number} x
     * @returns {Number}
     */
    private _clip;
    /**
     * performs the optimization step.
     * @private
     * @param {Matrix} head_embedding
     * @param {Matrix} tail_embedding
     * @param {Matrix} head
     * @param {Matrix} tail
     * @returns {Matrix}
     */
    private _optimize_layout;
    /**
     * @private
     * @returns {Matrix}
     */
    private next;
    _alpha: number;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=UMAP.d.ts.map