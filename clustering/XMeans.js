import { euclidean, euclidean_squared } from "../metrics/index.js";
import { Matrix, linspace } from "../matrix/index.js";
import { KMeans } from "../clustering/index.js";
/**
 * @class
 * @alias XMeans
 */
export class XMeans{
    /**
     * @constructor
     * @memberof module:clustering
     * @alias XMeans
     * @todo needs restructuring and repairing!!
     * @param {Matrix} matrix 
     * @param {Numbers} K_max
     * @param {Numbers} K_min
     * @param {Function} [metric = euclidean] 
     * @param {Number} [seed = 1987]
     * @returns {XMeans}
     * @see {@link https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf}
     * @see {@link https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py}
     * @see {@link https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java}
     */
    constructor(matrix, K_max = 10, K_min = 2, metric = euclidean, seed=1987) {
        //const first = super(matrix, K_min, metric, seed, false);
        const first = new KMeans(matrix, K_min, metric, seed, false);
        this._K_max = K_max;
        this._K_min = K_min;
        this._metric = metric;
        first.init(K_min, first._get_random_centroids(K_min));
        const randomizer = this._randomizer = first._randomizer;
        const centroids = first._cluster_centroids;
        const candidates = this._candidates = {};
        let K = K_min
        candidates[K] = {
            "kmeans": first,
            "cluster_centroids": centroids,
            "score": null,
        };
        const A = this._matrix = matrix;
        const N = A.shape[0];
        // foreach K in [K_min, K_max];
        do {
            console.log(K, candidates)
            const candidate = candidates[K];
            const clusters = candidate.kmeans.get_clusters();
            const parent_bic = this._bic(clusters, centroids, linspace(0, N - 1))
            candidate.score = parent_bic;
            const child_bic = [];
            const child_kmeans = [];
            const child_indices = [];
            // foreach cluster
            for (let j = 0; j < K; ++j) {
                const cluster = clusters[j];
                console.log(cluster.length)
                if (cluster.length < K_max) continue;
                const subset = Matrix.from(cluster.map(d => A.row(d)));
                const subset_kmeans = new KMeans(subset, 2, metric, 1987, false);
                subset_kmeans._randomizer = randomizer;
                subset_kmeans.init()
                const subset_cluster = subset_kmeans.get_clusters();
                const subset_centroids = subset_kmeans._cluster_centroids;
                const bic = this._bic(subset_cluster, subset_centroids, cluster)
                child_bic.push(bic);
                child_kmeans.push(subset_kmeans);
                child_indices.push(j);
            }
            //if (child_bic.length < (K )) break;
            //choose best
            let best_split = child_indices[0];
            let best_bic = child_bic[0];
            for (let i = 0; i < child_bic.length; ++i) {
                if (best_bic > child_bic[i]) {
                    best_split = child_indices[i];
                    best_bic = child_bic[i];
                }
            }
            const best_cluster_centroids = candidate.cluster_centroids.splice(best_split, 1, ...child_kmeans[best_split]._cluster_centroids);
            console.log(candidate.cluster_centroids, child_kmeans[best_split]._cluster_centroids)
            
            const parent_clusters = candidate.kmeans._clusters;
            const best_candidate_clusters = child_kmeans[best_split]._clusters;
            const best_candidate = new KMeans(A, K + 1, metric, 1987, false)
            best_candidate._randomizer = randomizer;
            // set clusters and centroids
            let counter = 0;
            //let cluster_number = best_cluster_centroids.length;
            best_candidate._clusters = parent_clusters.map(c => {
                if (c == best_split) {
                    return c + best_candidate_clusters[counter++];
                } else if (c > best_split) {
                    return c + 1;
                }
                return c;
            })
            //best_candidate._K = K + 1;
            console.log(best_candidate.get_clusters())
            //best_candidate.init(K + 1, best_cluster_centroids);
            console.log(best_candidate)
            //save best candidate.
            candidates[K + 1] = {
                "kmeans": best_candidate,
                "cluster_centroids": best_cluster_centroids,
                "score": child_bic[best_split],
            }
            

        } while (++K < K_max)


        // return best candidate.
        return this;
    }

    get_clusters() {
        let K_min = this._K_min;
        let K_max = this._K_max;
        const candidates = this._candidates;
        let best_score = candidates[K_min].score;
        let best_candidate = candidates[K_min].kmeans;
        for (let i = K_min + 1; i < K_max; ++i) {
            if (!(i in candidates)) break;
            const candidate = candidates[i];
            const score = candidate.score
            if (best_score < score) {
                best_score = score;
                best_candidate = candidate.kmeans;
            }
        }
        return best_candidate.get_clusters();
    }
    
    _bic(clusters, centroids, indices) {
        const A = this._matrix;
        const D = this._matrix.shape[1];
        const K = centroids.length;
        //const result = new Array(K).fill();
        let result = 0;

        let variance = 0;
        for (let i = 0; i < K; ++i) {
            const cluster = clusters[i];
            const centroid = centroids[i];
            const n = cluster.length;
            for (let j = 0; j < n; ++j) {
                variance += euclidean_squared(centroid, A.row(indices[cluster[j]])) ** 2;
            }
        }
        const N = clusters.reduce((a, b) => a + b.length, 0);
        const p = (K - 1) + D * K + 1;
        variance /= (N - K);

        for (let i = 0; i < K; ++i) {
            const n = clusters[i].length;
            const log_likelihood = 
                (n * Math.log(2 * Math.PI) -
                n * D * Math.log(variance) -
                (n - K)) * .5 + 
                n * Math.log(n) - 
                n * Math.log(N);
            result += log_likelihood - p * .5 * Math.log(N);
        }
        return result;
    }
    
}
