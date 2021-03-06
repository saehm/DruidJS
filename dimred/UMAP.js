import { Matrix } from "../matrix/index";
import { euclidean, euclidean_squared } from "../metrics/index";
import { BallTree } from "../knn/index";
import { neumair_sum } from "../numerical/index";
import { linspace } from "../matrix/index";
import { powell } from "../optimization/index";
import { DR } from "./DR.js";
import { max } from "../util/index";
import { KNN } from "../knn/KNN";

/**
 * @class
 * @alias UMAP
 */
export class UMAP extends DR {

    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias UMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {UMAP}
     */
    constructor(X, n_neighbors=15, local_connectivity=1, min_dist=1, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed)
        super.parameter_list = ["n_neighbors", "local_connectivity", "min_dist"];
        [ this._N, this._D ] = this.X.shape;
        n_neighbors = Math.min(this._N - 1, n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(local_connectivity, n_neighbors - 1));
        this.parameter("min_dist", min_dist);
        this._iter = 0;
        this._spread = 1;
        this._set_op_mix_ratio = 1;
        this._repulsion_strength = 1;
        this._negative_sample_rate = 5;
        this._n_epochs = 350;
        this._initial_alpha = 1;
        this.Y = new Matrix(this._N, this._d, () => this._randomizer.random);
        return this;
    }

    /**
     * @private
     * @param {Number} spread 
     * @param {Number} min_dist 
     * @returns {Array}
     */
    _find_ab_params(spread, min_dist) {
        const curve = (x, a, b) => 1 / (1 + a * Math.pow(x, 2 * b));
        const xv = linspace(0, spread * 3, 300);
        const yv = linspace(0, spread * 3, 300);
        
        for (let i = 0, n = xv.length; i < n; ++i) {
            const xv_i = xv[i];
            yv[i] = (xv_i < min_dist ? 1 : Math.exp(-(xv_i - min_dist) / spread));
        }
      
        const err = (p) => {
            const error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map(e => e * e)));
        }
      
        return powell(err, [1, 1]);
    }

    /**
     * @private
     * @param {Array<Array>} distances 
     * @param {Array<Number>} sigmas 
     * @param {Array<Number>} rhos 
     * @returns {Array}
     */
    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            for (let j = 0, m = distances[i].length; j < m; ++j) {
                const v = distances[i][j].value - rhos[i];
                distances[i][j].value = v > 0 ? Math.exp(-v / sigmas[i]) : 1;
            }
        }
        return distances;
    }

    /**
     * @private
     * @param {KNN|BallTree} knn 
     * @param {Number} k 
     * @returns {Object}
     */
    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = this._local_connectivity;
        const target = Math.log2(k);
        const rhos = [];
        const sigmas = [];
        const X = this.X;
        const N = X.shape[0];
        //const distances = [...X].map(x_i => knn.search(x_i, k).raw_data().reverse());

        const distances = [];
        if (this._metric === "precomputed") {
            for (let i = 0; i < N; ++i) {
                distances.push(knn.search(i, k).reverse())
            }
        } else {
           for (const x_i of X) {
                distances.push(knn.search(x_i, k).raw_data().reverse())
            }
        }

        for (let i = 0; i < N; ++i) {
            let lo = 0;
            let hi = Infinity;
            let mid = 1;

            const search_result = distances[i]
            const non_zero_dist = search_result.filter(d => d.value > 0);
            const non_zero_dist_length = non_zero_dist.length;
            if (non_zero_dist_length >= local_connectivity) {
                const index = Math.floor(local_connectivity);
                const interpolation = local_connectivity - index;
                if (index > 0) {
                    rhos.push(non_zero_dist[index - 1]);
                    if (interpolation > SMOOTH_K_TOLERANCE) {
                        rhos[i].value += interpolation * (non_zero_dist[index].value - non_zero_dist[index - 1]);
                    }
                } else {
                    rhos[i].value = interpolation * non_zero_dist[0].value;
                }
            } else if (non_zero_dist_length > 0) {
                rhos[i] = non_zero_dist[non_zero_dist_length - 1].value;
            }
            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    const d = search_result[j].value - rhos[i];
                    psum += (d > 0 ? Math.exp(-(d / mid)) : 1);
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    [hi, mid] = [mid, (lo + hi) / 2];
                } else {
                    if (hi === Infinity) {
                        [lo, mid] = [mid, mid * 2];
                    } else {
                        [lo, mid] = [mid, (lo + hi) / 2];
                    }
                }
            }
            sigmas[i] = mid;

            const mean_ithd = search_result.reduce((a, b) => a + b.value, 0) / search_result.length;
            //let mean_d = null;
            if (rhos[i] > 0) {
                if (sigmas[i] < MIN_K_DIST_SCALE * mean_ithd) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_ithd;
                }
            } else {
                const mean_d = distances.reduce((acc, res) => acc + res.reduce((a, b) => a + b.value, 0) / res.length);
                if (sigmas[i] > MIN_K_DIST_SCALE * mean_d) {
                    sigmas[i] = MIN_K_DIST_SCALE * mean_d;
                }
                
            }
        }
        return {
            "distances": distances, 
            "sigmas": sigmas, 
            "rhos": rhos
        }
    }

    /**
     * @private
     * @param {Matrix} X 
     * @param {Number} n_neighbors 
     * @returns {Matrix}
     */
    _fuzzy_simplicial_set(X, n_neighbors) {
        const N = X.shape[0];
        const metric = this._metric;
        const knn = metric === "precomputed" ? new KNN(X, "precomputed") : new BallTree(X.to2dArray, metric);
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        const result = new Matrix(N, N, "zeros");
        for (let i = 0; i < N; ++i) {
            const distances_i = distances[i];
            for (let j = 0; j < distances_i.length; ++j) {
                result.set_entry(i, distances_i[j].element.index, distances_i[j].value);
            }
        }
        const transposed_result = result.T;
        const prod_matrix = result.mult(transposed_result);
        return result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(this._set_op_mix_ratio)
            .add(prod_matrix.mult(1 - this._set_op_mix_ratio));
    }

    /**
     * @private
     * @param {Number} n_epochs 
     * @returns {Array}
     */
    _make_epochs_per_sample(n_epochs) {
        const weights = this._weights;
        const result = new Float32Array(weights.length).fill(-1);
        const weights_max = max(weights);
        const n_samples = weights.map(w => n_epochs * (w / weights_max));
        for (let i = 0; i < result.length; ++i) 
          if (n_samples[i] > 0) result[i] = Math.round(n_epochs / n_samples[i]);
        return result;
    }

    /**
     * @private
     * @param {Matrix} graph 
     * @returns {Object}
     */
    _tocoo(graph) {
        const rows = [];
        const cols = [];
        const data = [];
        const [ rows_n, cols_n ] = graph.shape;
        for (let row = 0; row < rows_n; ++row) {
            for (let col = 0; col < cols_n; ++col) {
                const entry = graph.entry(row, col);
                if (entry !== 0) {
                    rows.push(row);
                    cols.push(col);
                    data.push(entry);
                }
            }
        }
        return {
            "rows": rows, 
            "cols": cols, 
            "data": data
        };
    }

    /**
     * Computes all necessary 
     * @returns {UMAP}
     */
    init() {
        const [ a, b ] = this._find_ab_params(this._spread, this._min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, this._n_neighbors);
        const { rows, cols, data: weights } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        this._weights = weights;
        this._epochs_per_sample = this._make_epochs_per_sample(this._n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map(d => d * this._negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        return this;
    }

    set local_connectivity(value) {
        this._local_connectivity = value;
    }

    get local_connectivity() {
        return this._local_connectivity;
    }

    set min_dist(value) {
        this._min_dist = value;
    }

    get min_dist() {
        return this._min_dist;
    }

    graph() {
        this.check_init();
        return { cols: this._head, rows: this._tail, weights: this._weights };
    }

    /**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    transform(iterations=350) {
        if (this._n_epochs != iterations) {
            this._n_epochs = iterations;
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }


    /**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */
    * generator(iterations=350) {
        if (this._n_epochs != iterations) {
            this._n_epochs = iterations;
            this.init();
        }
        this.check_init();
        for (let i = 0; i < iterations; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @private
     * @param {Number} x 
     * @returns {Number}
     */
    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

    /**
     * performs the optimization step.
     * @private
     * @param {Matrix} head_embedding 
     * @param {Matrix} tail_embedding 
     * @param {Matrix} head 
     * @param {Matrix} tail 
     * @returns {Matrix}
     */
    _optimize_layout(head_embedding, tail_embedding, head, tail) {
        const { 
            _d: dim, 
            _alpha: alpha, 
            _repulsion_strength: repulsion_strength, 
            _a: a, 
            _b: b,
            _epochs_per_sample: epochs_per_sample,
            _epochs_per_negative_sample: epochs_per_negative_sample,
            _epoch_of_next_negative_sample: epoch_of_next_negative_sample,
            _epoch_of_next_sample: epoch_of_next_sample,
            _clip: clip
        } = this;
        const tail_length = tail.length;

        for (let i = 0, n = epochs_per_sample.length; i < n; ++i) {
            if (epoch_of_next_sample[i] <= this._iter) {
                const j = head[i];
                const k = tail[i];
                const current = head_embedding.row(j);
                const other = tail_embedding.row(k);
                const dist = euclidean_squared(current, other);
                let grad_coeff = 0;
                if (dist > 0) {
                    grad_coeff = (-2 * a * b * Math.pow(dist, b - 1)) / (a * Math.pow(dist, b) + 1);
                }
                for (let d = 0; d < dim; ++d) {
                    const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                    const c = current[d] + grad_d;
                    const o = other[d] - grad_d;
                    current[d] = c;
                    other[d] = o;
                    head_embedding.set_entry(j, d, c);
                    tail_embedding.set_entry(k, d, o);
                }
                epoch_of_next_sample[i] += epochs_per_sample[i];
                const n_neg_samples = (this._iter - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i];
                for (let p = 0; p < n_neg_samples; ++p) {
                    const k = Math.floor(this._randomizer.random * tail_length);
                    const other = tail_embedding.row(tail[k]);
                    const dist = euclidean_squared(current, other);
                    let grad_coeff = 0;
                    if (dist > 0) {
                        grad_coeff = (2 * repulsion_strength * b) / ((.01 + dist) * (a * Math.pow(dist, b) + 1));
                    } else if (j === k) {
                        continue;
                    }
                    for (let d = 0; d < dim; ++d) {
                        const grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
                        const c = current[d] + grad_d;
                        const o = other[d] - grad_d;
                        current[d] = c;
                        other[d] = o;
                        head_embedding.set_entry(j, d, c);
                        tail_embedding.set_entry(tail[k], d, o);
                    }
                }
                epoch_of_next_negative_sample[i] += (n_neg_samples * epochs_per_negative_sample[i]);
            }
        }
        return head_embedding;
    }

    /**
     * @private
     * @returns {Matrix}
     */
    next() {
        let iter = ++this._iter;
        let Y = this.Y;

        this._alpha = (this._initial_alpha * (1 - iter / this._n_epochs));
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }
} 