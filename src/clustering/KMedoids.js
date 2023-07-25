import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { linspace, Matrix } from "../matrix/index.js";
import { min } from "../util/index.js";
/**
 * @class
 * @alias KMedoids
 */
export class KMedoids {
    /**
     * @constructor
     * @memberof module:clustering
     * @alias KMedoids
     * @todo needs restructuring. 
     * @param {Matrix} matrix - data matrix
     * @param {Numbers} K - number of clusters
     * @param {number} [max_iter=null] - maximum number of iterations. Default is 10 * Math.log10(N)
     * @param {Function} [metric = euclidean] - metric defining the dissimilarity 
     * @param {Number} [seed = 1212] - seed value for random number generator
     * @returns {KMedoids}
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
    constructor(matrix, K, max_iter=null, metric = euclidean, seed=1212) {
        this._metric = metric;
        this._matrix = matrix;
        this._A = this._matrix.to2dArray;
        this._K = K;
        const [N, D] = matrix.shape;
        this._N = N;
        this._D = D;
        this._max_iter = max_iter || 10 * Math.log10(N) 
        this._distance_matrix = new Matrix(N, N, "zeros");
        /* for (let i = 1; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dist = metric(this._A[i], this._A[j]);
                this._distance_matrix.set_entry(i, j, dist);
                this._distance_matrix.set_entry(j, i, dist)
            }
        } */
        if (K > N) K = N;
        this._randomizer = new Randomizer(seed);
        this._clusters = new Array(N).fill(undefined);
        this._cluster_medoids = this._get_random_medoids(K);
        //if (init) this.init(K, this._cluster_medoids);
        this._is_initialized = false;
        return this;
    }

    /**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */
    get_clusters() {
        const K = this._K;
        const A = this._A;
        if (!this._is_initialized) {
            this.init(K, this._cluster_medoids);
        }
        const result = new Array(K).fill().map(() => new Array());
        A.forEach((x_j, j) => {
            result[this._nearest_medoid(x_j, j).index_nearest].push(j);
        })
        result.medoids = this._cluster_medoids;
        return result;
    }

    async* generator() {
        const max_iter = this._max_iter;
        yield this.get_clusters()
        let finish = false;
        let i = 0
        do {
            finish = this._iteration();
            yield this.get_clusters();
        } while (!finish && ++i < max_iter)
    }

    /**
     * Algorithm 1. FastPAM1: Improved SWAP algorithm
     */
    /* _iteration_1() {
        const A = this._A;
        const N = this._N;
        const K = this._K;
        const medoids = this._cluster_medoids;
        let DeltaTD = 0;
        let m0 = null;
        let x0 = null;
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const nearest_medoid = this._nearest_medoid(x_j, j);
                const d_j = nearest_medoid.distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    // disance to new medoid
                    const d_oj = this._get_distance(o, j, x_o, x_j);
                    const {
                        "index_nearest": n,
                        "distance_nearest": d_n,
                        "distance_second": d_s,
                    } = this._nearest_medoid(x_o, o); 
                    this._clusters[o] = n; // cached values
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change
                    if (d_oj < d_n) { // reassignment check
                        deltaTD.forEach((d_i, i) => {
                            if (n !== i) {
                                deltaTD[i] = d_i + d_oj - d_n; // update loss change
                            }
                        });
                    }
                });
                // choose best medoid i;
                const i = deltaTD
                    .map((d, i) => [d, i])
                    .sort((d1, d2) => d1[0] - d2[0])[0][1];
                const deltaTD_i = deltaTD[i];
                // store
                if (deltaTD_i < DeltaTD) {
                    DeltaTD = deltaTD_i;
                    m0 = i;
                    x0 = j;
                }
            }
        });

        if (DeltaTD >= 0) {
            return true // break loop if DeltaTD >= 0
        }
        // swap roles of medoid m and non-medoid x;
        medoids[m0] = x0;
        this._cluster_medoids = medoids;
        return false
    } */

    /** Algorithm 2. FastPAM2: SWAP with multiple candidates
     * 
     */
    _iteration() {
        const A = this._A;
        const K = this._K;
        const medoids = this._cluster_medoids;
        const cache = A.map((x_o, o) => this._nearest_medoid(x_o, o));
        // empty best candidates array
        const DeltaTD = new Array(K).fill(0);
        const xs = new Array(K).fill(null);
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const d_j = cache[j].distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    if (j === o) return;
                    const d_oj = this._get_distance(o, j, x_o, x_j); // distance to new medoid
                    const {"index_nearest": n, "distance_nearest": d_n, "distance_second": d_s} = cache[o]; // cached
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change for x_o
                    // Reassignment check
                    if (d_oj < d_n) { 
                        // update loss change
                        for (let i = 0; i < K; ++i) {
                            if (i !== n) deltaTD[i] += d_oj - d_n;
                        }
                    }
                });
                // remember best swap for i;
                deltaTD
                    .map((d, i) => [d, i])
                    .filter(([d, i]) => d < DeltaTD[i])
                    .forEach(([d, i]) => {
                        if (d < DeltaTD[i]) {
                            DeltaTD[i] = d;
                            xs[i] = j;
                        }
                    })
            }
        })
        // stop if no improvements were found
        if (min(DeltaTD) >= 0) return true; 

        // execute all improvements
        while (min(DeltaTD) < 0) {
            // swap roles of medoid m_i and non_medoid xs_i
            const i = DeltaTD
                .map((d, i) => [d, i])
                .sort(([a], [b]) => a - b)[0][1];
            if (medoids.filter(m => m == xs[i]).length == 0) {
                medoids[i] = xs[i];
            }
            // disable the swap just performed
            DeltaTD[i] = 0; 
            // recompute TD for remaining swap candidates
            DeltaTD
                .map((d_j, j) => [d_j, j])
                .filter(([d_j]) => d_j < 0)
                .forEach(([_, j]) => {
                    const x_j = A[j];
                    let sum = 0;
                    A.forEach((x_o, o) => {
                        if (medoids.findIndex(m => m != j && m == o) >= 0) return;
                        if (i == j) return;
                        if (cache[o].index_nearest === medoids[j])
                            sum += (Math.min(this._get_distance(o, j, x_o, x_j), cache[o].distance_second) - cache[o].distance_nearest); 
                        else {
                            sum += (Math.min(this._get_distance(o, j, x_o, x_j) - cache[o].distance_nearest, 0));
                        }
                    });
                    DeltaTD[j] = sum;
                })
        }
        this._cluster_medoids = medoids;
        return false;
    }

    _get_distance(i, j, x_i=null, x_j=null) {
        if (i === j) return 0;
        const D = this._distance_matrix;
        const A = this._A;
        const metric = this._metric;
        let d_ij = D.entry(i, j);
        if (d_ij === 0) {
            d_ij = metric(x_i || A[i], x_j || A[j]);
            D.set_entry(i, j, d_ij);
            D.set_entry(j, i, d_ij);
        }
        return d_ij;
    }

    _nearest_medoid(x_j, j) {
        const medoids = this._cluster_medoids;
        const A = this._A;
        const [nearest, second] = medoids
            .map((m, i) => {
                const x_m = A[m]; 
                return [this._get_distance(j, m, x_j, x_m), i];
            })
            .sort((m1, m2) => m1[0] - m2[0]);
        
        return { 
            "distance_nearest": nearest[0], 
            "index_nearest": nearest[1],
            "distance_second": second[0],
            "index_second": second[1],
        };
    }

    /**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */
    init(K, cluster_medoids) {
        if (!K) K = this._K;
        if (!cluster_medoids) cluster_medoids = this._get_random_medoids(K);
        const max_iter = this._max_iter;
        let finish = false;
        let i = 0
        do {
            finish = this._iteration();
        } while (!finish && ++i < max_iter)
        return this;
    }

    /**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @param {number} K - number of clusters
     * 
     */
    _get_random_medoids(K) {
        const N = this._N;
        const A = this._A;
        const indices = linspace(0, N - 1);
        const randomizer = this._randomizer;
        const n = Math.min(N, 10 + Math.ceil(Math.sqrt(N)));
        const TD = new Array(n).fill(Infinity);
        const medoids = [];
        // first medoid
        let TD0 = Infinity;
        let S = randomizer.choice(indices, n);
        for (let j = 0; j < n; ++j) {
            const S_j = S[j];
            const x_j = A[S_j];
            for (let o = 0; o < n; ++o) {
                if (o === j) continue;
                const x_o = A[S[o]];
                TD[j] += this._get_distance(j, o, x_j, x_o);
            }
            if (TD[j] < TD0) {
                TD0 = TD[j]; // smallest distance sum
                medoids.push(S_j);
            }
        }
        // other medoids
        for (let i = 1; i < K; ++i) {
            let DeltaTD = Infinity;
            S = randomizer.choice(indices.filter(index => medoids.findIndex(d => d === index) < 0), n);
            for (let j = 0; j < n; ++j) {
                let deltaTD = 0;
                const S_j = S[j];
                const x_j = A[S_j];
                for (let o = 0; o < n; ++o) {
                    if (o === j) continue;
                    const S_o = S[o];
                    const x_o = A[S_o];
                    let delta = this._get_distance(S_j, S_o, x_j, x_o) - min(medoids.map(m => this._get_distance(S_o, m, x_o)));
                    if (delta < 0) {
                        deltaTD = deltaTD + delta;
                    }
                }
                // best reduction
                if (deltaTD < DeltaTD) {
                    DeltaTD = deltaTD;
                    medoids.push(S_j);
                }
            }
            TD0 += DeltaTD;
        }
        return medoids.slice(0, K);
    }
    
}
