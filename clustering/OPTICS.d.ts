/**
 * @class
 * @alias OPTICS
 */
export class OPTICS {
    /**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     * @constructor
     * @memberof module:clustering
     * @alias OPTICS
     * @todo needs restructuring.
     * @param {Matrix} matrix - the data.
     * @param {Number} epsilon - the minimum distance which defines whether a point is a neighbor or not.
     * @param {Number} min_points - the minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
     * @param {Function} [metric = euclidean] - the distance metric which defines the distance between two points of the {@link matrix}.
     * @returns {OPTICS}
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
    constructor(matrix: Matrix, epsilon: number, min_points: number, metric?: Function);
    _matrix: Matrix;
    _epsilon: number;
    _min_points: number;
    _metric: Function;
    _ordered_list: any[];
    _clusters: any[];
    _DB: any[];
    /**
     * Computes the clustering.
     */
    init(): OPTICS;
    _cluster_index: number;
    /**
     *
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Array} An array consisting of the {@link epsilon}-neighborhood of {@link p}.
     */
    private _get_neighbors;
    /**
     *
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Number} The distance to the {@link min_points}-th nearest point of {@link p}, or undefined if the {@link epsilon}-neighborhood has fewer elements than {@link min_points}.
     */
    private _core_distance;
    /**
     * Updates the reachability distance of the points.
     * @private
     * @param {Object} p
     * @param {Heap} seeds
     */
    private _update;
    /**
     * Expands the {@link cluster} with points in {@link seeds}.
     * @private
     * @param {Heap} seeds
     * @param {Array} cluster
     */
    private _expand_cluster;
    /**
     * Returns an array of clusters.
     * @returns {Array<Array>} Array of clusters with the indices of the rows in given {@link matrix}.
     */
    get_clusters(): Array<any[]>;
    /**
     * @returns {Array} Returns an array, where the ith entry defines the cluster affirmation of the ith point of {@link matrix}. (-1 stands for outlier)
     */
    get_cluster_affirmation(): any[];
}
//# sourceMappingURL=OPTICS.d.ts.map