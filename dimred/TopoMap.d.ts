/**
 * @class
 * @alias TopoMap
 * @memberof module:dimensionality_reduction
 * @extends DR
 */
export class TopoMap extends DR {
    /**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TopoMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X: Matrix, parameters: {
        metric?: Function;
        seed?: number;
    });
    _distance_matrix: Matrix;
    /**
     * @private
     */
    private __lazy_distance_matrix;
    /**
     * Computes the minimum spanning tree, using a given metric
     * @private
     * @param {Function} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    private _make_minimum_spanning_tree;
    /**
     * initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init(): TopoMap;
    Y: Matrix;
    _Emst: any[][];
    /**
     * Returns true if Point C is left of line AB.
     * @private
     * @param {Array} PointA - Point A of line AB
     * @param {Array} PointB - Point B of line AB
     * @param {Array} PointC - Point C
     * @returns {Boolean}
     */
    private __hull_cross;
    /**
     * Computes the convex hull of the set of Points S
     * @private
     * @param {Array} S - Set of Points.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     * @returns {Array} convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     */
    private __hull;
    /**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     * @private
     * @param {Array} PointA
     * @param {Array} PointB
     * @return {Object} Object containing the sinus- and cosinus-values for a rotation.
     */
    private __findAngle;
    /**
     * @private
     * @param {Array} hull
     * @param {Array} p
     * @param {Bool} topEdge
     */
    private __align_hull;
    /**
     * @private
     * @param {Array} Point - The point which should get transformed.
     * @param {Object} Transformation - contains the values for translation and rotation.
     */
    private __transform;
    /**
     * Calls {@link __transform} for each point in Set C
     * @private
     * @param {Array} C - Set of points.
     * @param {Object} t - Transform object.
     * @param {Number} yOffset - value to offset set C.
     */
    private __transform_component;
    /**
     * @private
     * @param {Array} u - point u
     * @param {Array} v - point v
     * @param {Number} w - edge weight w
     */
    private __align_components;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=TopoMap.d.ts.map