import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { BallTree } from "../knn/index";
import { neumair_sum } from "../numerical/index";
import { linspace } from "../matrix/index";
import { powell } from "../optimization/index";
import { DR } from "./DR.js";

export class UMAP extends DR {
    static parameter_list = ["local_connectivity", "min_dist"];

    constructor(X, local_connectivity=1, min_dist=1, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed)
        super.parameter_list = UMAP.parameter_list;
        [ this._N, this._D ] = X.shape;
        this.parameter("local_connectivity", local_connectivity);
        this.parameter("min_dist", min_dist);
        this._iter = 0;
        this._n_neighbors = 11;
        this._spread = 1;
        this._set_op_mix_ratio = 1;
        this._repulsion_strength = 1;
        this._negative_sample_rate = 5;
        this._n_epochs = 350;
        this._initial_alpha = 1;
        this.Y = new Matrix(this._N, this._d, () => this._randomizer.random);
    }

    _find_ab_params(spread, min_dist) {
        function curve(x, a, b) {
            return 1 / (1 + a * Math.pow(x, 2 * b));
        }
      
        var xv = linspace(0, spread * 3, 300);
        var yv = linspace(0, spread * 3, 300);
        
        for ( var i = 0, n = xv.length; i < n; ++i ) {
            if (xv[i] < min_dist) {
                yv[i] = 1
            } else {
                yv[i] = Math.exp(-(xv[i] - min_dist) / spread)
            }
        }
      
        function err(p) {
            var error = linspace(1, 300).map((_, i) => yv[i] - curve(xv[i], p[0], p[1]));
            return Math.sqrt(neumair_sum(error.map(e => e * e)));
        }
      
        var [ a, b ] = powell(err, [1, 1])
        return [ a, b ]
    }

    _compute_membership_strengths(distances, sigmas, rhos) {
        for (let i = 0, n = distances.length; i < n; ++i) {
            for (let j = 0, m = distances[i].length; j < m; ++j) {
                let v = distances[i][j].value - rhos[i];
                let value = 1;
                if (v > 0) {
                    value = Math.exp(-v / sigmas[i]);
                }
                distances[i][j].value = value;
            }
        }
        return distances;
    }

    _smooth_knn_dist(knn, k) {
        const SMOOTH_K_TOLERANCE = 1e-5;
        const MIN_K_DIST_SCALE = 1e-3;
        const n_iter = 64;
        const local_connectivity = this._local_connectivity
        const bandwidth = 1
        const target = Math.log2(k) * bandwidth
        const rhos = []
        const sigmas = []
        const X = this.X

        let distances = [];
        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            let x_i = X.row(i);
            distances.push(knn.search(x_i, Math.max(local_connectivity, k)).raw_data().reverse())
        }

        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            let search_result = distances[i]
            rhos.push(search_result[0].value);

            let lo = 0;
            let hi = Infinity;
            let mid = 1;

            for (let x = 0; x < n_iter; ++x) {
                let psum = 0;
                for (let j = 0; j < k; ++j) {
                    let d = search_result[j].value - rhos[i];
                    psum += (d > 0 ? Math.exp(-(d / mid)) : 1);
                }
                if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) {
                    break;
                }
                if (psum > target) {
                    //[hi, mid] = [mid, (lo + hi) / 2];
                    hi = mid;
                    mid = (lo + hi) / 2; // PROBLEM mit hi?????
                } else {
                    lo = mid;
                    if (hi === Infinity) {
                        mid *= 2;
                    } else {
                        mid = (lo + hi) / 2;
                    }
                }
            }
            sigmas[i] = mid

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

    _fuzzy_simplicial_set(X, n_neighbors) {
        const knn = new BallTree(X.to2dArray, euclidean);
        let { distances, sigmas, rhos } = this._smooth_knn_dist(knn, n_neighbors);
        distances = this._compute_membership_strengths(distances, sigmas, rhos);
        let result = new Matrix(X.shape[0], X.shape[0], "zeros");
        for (let i = 0, n = X.shape[0]; i < n; ++i) {
            for (let j = 0; j < n_neighbors; ++j) {
                result.set_entry(i, distances[i][j].element.index, distances[i][j].value);
            }
        }
        const transposed_result = result.T;
        const prod_matrix = result.mult(transposed_result);
        result = result
            .add(transposed_result)
            .sub(prod_matrix)
            .mult(this._set_op_mix_ratio)
            .add(prod_matrix.mult(1 - this._set_op_mix_ratio));
        return result;
    }

    _make_epochs_per_sample(graph, n_epochs) {
        const { data: weights } = this._tocoo(graph);
        let result = new Array(weights.length).fill(-1);
        const weights_max = Math.max(...weights);
        const n_samples = weights.map(w => n_epochs * (w / weights_max));
        result = result.map((d, i) => (n_samples[i] > 0) ? Math.round(n_epochs / n_samples[i]) : d);
        return result;
    }

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
        return {rows: rows, cols: cols, data: data};
    }

    init() {
        const [ a, b ] = this._find_ab_params(this._spread, this._min_dist);
        this._a = a;
        this._b = b;
        this._graph = this._fuzzy_simplicial_set(this.X, this._n_neighbors);
        this._epochs_per_sample = this._make_epochs_per_sample(this._graph, this._n_epochs);
        this._epochs_per_negative_sample = this._epochs_per_sample.map(d => d * this._negative_sample_rate);
        this._epoch_of_next_sample = this._epochs_per_sample.slice();
        this._epoch_of_next_negative_sample = this._epochs_per_negative_sample.slice();
        const { rows, cols } = this._tocoo(this._graph);
        this._head = rows;
        this._tail = cols;
        return this
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

    transform(iterations) {
        this.check_init();
        iterations = iterations || this._n_epochs;
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this.projection;
    }

    * generator() {
        this.check_init();
        this._iter = 0
        while (this._iter < this._n_epochs) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    _clip(x) {
        if (x > 4) return 4;
        if (x < -4) return -4;
        return x;
    }

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
                const dist = euclidean(current, other)//this._metric(current, other);
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
                    const dist = euclidean(current, other)//this._metric(current, other);
                    let grad_coeff = 0;
                    if (dist > 0) {
                        grad_coeff = (2 * repulsion_strength * b) / ((.01 + dist) * (a * Math.pow(dist, b) + 1));
                    } else if (j == k) {
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

    next() {
        let iter = ++this._iter;
        let Y = this.Y;

        this._alpha = (this._initial_alpha * (1 - iter / this._n_epochs));
        this.Y = this._optimize_layout(Y, Y, this._head, this._tail);

        return this.Y;
    }
} 