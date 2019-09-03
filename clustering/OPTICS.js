import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/index";
import { Heap } from "../datastructure/index";
import { linspace } from "../matrix/index";
import { ENETUNREACH } from "constants";

/**
 * @class
 * @alias OPTICS
 */
export class OPTICS {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias OPTICS
     * @todo needs restructuring. 
     * @param {Matrix} matrix 
     * @param {Number} epsilon
     * @param {Number} min_points 
     * @param {Function} [metric = euclidean]
     * @returns {OPTICS}
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     */
    constructor(matrix, epsilon, min_points, metric = euclidean) {
        this._matrix = matrix;
        this._epsilon = epsilon;
        this._min_points = min_points;
        this._metric = metric;
        this._ordered_list = [];
        this._clusters = [];
        this._DB = new Array(N).fill();
        this.init();
        return this;
    }

    /**
     * 
     */
    init() {
        const ordered_list = this._ordered_list;
        const matrix = this._matrix;
        const N = matrix.shape[0];
        const DB = this._DB;
        for (let i = 0; i < N; ++i) {
            const e = matrix.row(i);
            DB[i] = {
                "element": e,
                "index": i,
                "reachability_distance": undefined,
                "processed": false,
            }
        }
        
        const clusters = this._clusters;
        let cluster_index = this._cluster_index = 0;

        for (const p of DB) {
            if (p.processed) continue;
            p.neighbors = this._get_neighbors(p);
            p.processed = true;
            clusters.push([p.index])
            cluster_index = clusters.length - 1;
            ordered_list.push(p);
            if (this._core_distance(p) != undefined) {
                const seeds = new Heap(null, d => d.reachability_distance, "min")
                this._update(p, seeds);
                this._expand_cluster(seeds, clusters[cluster_index]);
            }
        }
    }

    _get_neighbors(p) {
        //if ("neighbors" in p) return p.neighbors;
        const DB = this._DB;
        const metric = this._metric;
        const epsilon = this._epsilon;
        const neighbors = [];
        /*const p_neighbors = Heap.heapify(DB, d => metric(d.element, p.element), "min");
        let searching = true;
        while (searching && !p_neighbors.empty) {
            const { "element": neighbor, "value": distance } = p_neighbors.pop();
            if (distance < epsilon) {
                if (p != neighbor) {
                    neighbors.push(neighbor)
                }    
            } else {
                searching = false;
            }
        }
        return neighbors.reverse();*/

        for (const q of DB) {
            if (q.index == p.index) continue;
            if (metric(p.element, q.element) < epsilon) {
                neighbors.push(q);
            }
        }
        return neighbors;
    }

    _core_distance(p) {
        const min_points = this._min_points;
        const metric = this._metric;
        if (p.neighbors && p.neighbors.length <= min_points) {
            return undefined;
        }
        return metric(p.element, p.neighbors[min_points].element);
    }

    _update(p, seeds) {
        const metric = this._metric;
        const core_distance = this._core_distance(p);
        const neighbors = this._get_neighbors(p);//p.neighbors;
        for (const o of neighbors) {
            if (o.processed) continue;
            const new_reachability_distance = Math.max(core_distance, metric(p.element, o.element));
            if (o.reachability_distance == undefined) { // o is not in seeds
                o.reachability_distance = new_reachability_distance;
                seeds.push(o);
            } else { // o is in seeds
                if (new_reachability_distance < o.reachability_distance) {
                    o.reachability_distance = new_reachability_distance;
                    seeds = Heap.heapify(seeds.data(), d => d.reachability_distance, "min"); // seeds change key =/
                }
            }
            
        }
    }

    _expand_cluster(seeds, cluster) {
        const ordered_list = this._ordered_list;
        while (!seeds.empty) {
            const q = seeds.pop().element;
            q.neighbors = this._get_neighbors(q);
            q.processed = true;
            cluster.push(q.index);
            ordered_list.push(q);
            if (this._core_distance(q) != undefined) {
                this._update(q, seeds);
                this._expand_cluster(seeds, cluster);
            }
        }
    }

    get_clusters(threshold) {
        /*const ordered_list = this._ordered_list;
        console.log(ordered_list.map(e => e.reachability_distance))
        const N = ordered_list.length;
        const clusters = [];
        const outliers = [];
        let switched = true;
        let cluster_index = -1;
        for (const e of ordered_list) {
            if (e.reachability_distance < threshold) {
                if (switched) {
                    switched = false;
                    cluster_index++;
                    clusters.push([]);
                } else {

                }
                clusters[cluster_index].push(e.index);
            } else {
                if (switched) {
                    outliers.push(e.index);
                } else {
                    outliers.push(e.index);
                    switched = true;
                }
            }
        }
        //console.log(clusters, outliers)
        clusters.push(outliers)
        return clusters;*/

        const clusters = [];
        const outliers = [];
        const min_points = this._min_points;
        for (const cluster of this._clusters) {
            if (cluster.length < min_points) {
                outliers.push(...cluster);
            } else {
                clusters.push(cluster);
            }
        }
        clusters.push(outliers);
        return clusters;

        //return this._clusters
    }

    
}
