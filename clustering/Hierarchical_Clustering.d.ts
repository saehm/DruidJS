/**
 * @class
 * @alias Hierarchical_Clustering
 */
export class Hierarchical_Clustering {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias Hierarchical_Clustering
     * @todo needs restructuring.
     * @param {Matrix} - Data or distance matrix if metric is 'precomputed'
     * @param {("single"|"complete"|"average")} [linkage = "complete"]
     * @param {Function|"precomputed"} [metric = euclidean]
     * @returns {Hierarchical_Clustering}
     */
    constructor(matrix: any, linkage?: ("single" | "complete" | "average"), metric?: Function | "precomputed");
    _id: number;
    _matrix: Matrix;
    _metric: Function | "precomputed";
    _linkage: "complete" | "single" | "average";
    root: Cluster;
    /**
     *
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(value: number, type?: ("distance" | "depth")): Array<any[]>;
    /**
     * @private
     * @param {} node
     * @param {*} f
     * @param {*} value
     * @param {*} result
     */
    private _traverse;
    /**
     * computes the tree.
     */
    init(): Hierarchical_Clustering;
    _n: any;
    _d_min: Float64Array;
    _distance_matrix: Matrix;
    _clusters: any[];
    _c_size: Uint16Array;
    /**
     * computes the tree.
     */
    do(): Cluster;
}
import { Matrix } from "../matrix/index.js";
declare class Cluster {
    constructor(id: any, left: any, right: any, dist: any, centroid: any, index: any, size: any, depth: any);
    id: any;
    left: any;
    right: any;
    dist: any;
    index: any;
    size: any;
    depth: any;
    centroid: any;
    parent: any;
    _calculate_centroid(left: any, right: any): Float64Array;
    get isLeaf(): boolean;
    leaves(): any;
    descendants(): any;
}
export {};
//# sourceMappingURL=Hierarchical_Clustering.d.ts.map