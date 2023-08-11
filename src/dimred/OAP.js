import { euclidean, chebyshev } from "../metrics/index.js";
import { MDS } from "../dimred/MDS.js";
import { Randomizer } from "../util/index.js";
import { BallTree } from "../knn/BallTree.js";
import { Matrix } from "../matrix/index.js";
import { neumair_sum } from "../numerical/index.js";

/**
 *
 */
export class OAP {
    constructor(X, depth_field_lag, step_size, depth_weight, d = 2, metric = euclidean, seed = 1212) {
        this._X = X;
        this._d = d;
        [this._N, this._D] = X.shape;
        this._depth_field_lag = depth_field_lag;
        this._step_size = step_size;
        this._depth_weight = depth_weight;
        this._J = 3;
        this._max_iter = 1;
        this._metric = metric;
        this._seed = seed;
        this._randomizer = new Randomizer(seed);
    }

    _data_depth(technique = "chebyshev") {
        const X = this._X;
        const N = this._N;
        const h = new Float32Array(N);
        let deepest_point = 0;
        if (technique === "mdb") {
            h.fill(1);

            /*
            // Modified Band Depth 
            // https://www.tandfonline.com/doi/pdf/10.1198/jasa.2009.0108?casa_token=E1Uzntgs-5AAAAAA:Eo8mUpJDhpLQ5RHBkCB3Mdz0tbGM3Q0v78bwyCIAv7-peLGwfG3TcXLqShIaYuJLEqKc7GvaKlgvUg 
            const randomizer = this._randomizer;
            const h = new Float32Array(this._N);
            const J = this._J;
            const N = this._N;
            const D = this._D;
            const X = this._X;

            const one_div_by_n_choose_j = 1;
            for (let row = 0; row < N; ++row) {
                const x = X.row(row);
                const B_min = new Float32Array(D).fill(Infinity);
                const B_max = new Float32Array(D).fill(-Infinity);
                let r = Math.floor(randomizer.random * N);
                for (let i = 0; i < J; ++i) {
                    const x_j = X.row(r);
                    for (let d = 0; d < D; ++d) {
                        const x_jd = x_j[d]
                        B_min[d] = Math.min(B_min[d], x_jd);
                        B_max[d] = Math.max(B_max[d], x_jd);
                    }
                    r += Math.floor(randomizer.random * (N - 1));
                    r = r % N;
                }
                for (let d = 0; d < D; ++d) {
                    const x_d = x[d];
                    if (x_d >= B_min[d] && x_d <= B_max[d]) {
                        ++h[row]
                    }
                }
            }
            this._h = h;*/
        } else if (technique === "chebyshev") {
            // L∞ Depth
            // https://arxiv.org/pdf/1506.01332.pdf
            for (let i = 0; i < N; ++i) {
                let x = X.row(i);
                let sum = 0;
                for (let j = 0; j < N; ++j) {
                    if (i !== j) {
                        sum += chebyshev(x, X.row(j));
                    }
                }
                h[i] = 1 / (1 + sum / N);
                if (h[deepest_point] < h[i]) {
                    deepest_point = i;
                }
            }
        }
        this._h = h;
        this._deepest_point = deepest_point;
    }

    init() {
        this._iter = 0;
        // init with MDS
        const init_MDS = new MDS(this._X, this._d, this._metric);
        //console.log(init_MDS)
        this._Y = init_MDS.transform();

        // try häääh?
        this._X_distances = init_MDS._d_X;
        /*let max = -Infinity
        init_MDS._d_X._data.forEach(dx => max = Math.max(dx, max));
        this._X_distances = init_MDS._d_X.divide(max);*/
        // end try hääääh?

        // compute order statistics
        this._data_depth();
        this._M = this._monotonic_field(this._Y);
        //
        return this;
    }

    set depth_field_lag(value) {
        this._depth_field_lag = value;
    }

    get depth_field_lag() {
        return this._depth_field_lag;
    }

    set step_size(value) {
        this._step_size = value;
    }

    get step_size() {
        return this._step_size;
    }

    set depth_weight(value) {
        this._depth_weight = value;
    }

    get depth_weight() {
        return this._depth_weight;
    }

    transform(iterations = this._max_iter) {
        for (let i = 0; i < iterations; ++i) {
            this.next();
        }
        return this._Y;
    }

    *transform_iter() {
        while (true) {
            this.next();
            yield this._Y;
        }
    }

    _monotonic_field(Y) {
        const h = this._h;
        const Y_ = this._Y_;
        Y, h, Y_;
        const nn = new BallTree();
        nn.add(Y.to2dArray);

        const N = 5;
        let M = (x) => {
            let neighbors = nn.search(x, N).toArray();
            let d_sum = 0; //neighbors.reduce((a, b) => a + b.value, 0);
            let m = 0;
            for (let i = 0; i < N; ++i) {
                d_sum += neighbors[i].value;
                m += h[neighbors[i].element.index] * neighbors[i].value;
            }
            //console.log(m, d_sum)
            m /= d_sum;
            return m;
        };
        return M;
    }

    next() {
        const iter = ++this._iter;
        const l = this._depth_field_lag;
        const step_size = this._step_size;
        const w = this._depth_weight;
        const N = this._N;
        const dim = this._d;
        const d_X = this._X_distances;
        const h = this._h;
        let Y = this._Y;

        if (iter % l === 1) {
            // compute monotonic field
            this._Y_ = this._Y.clone();
            this._M = this._monotonic_field(Y);
        }
        const M = this._M;
        // perform gradient step

        // MDS stress step
        /*for (let i = 0; i < N; ++i) {
            const d_x = d_X.row(i);
            const y_i = Y.row(i)
            const delta_mds_stress = new Float32Array(dim);
            for (let j = 0; j < N; ++j) {
                if (i !== j) {
                    const y_j = Y.row(j)
                    const d_y = metric(y_i, y_j);
                    const d_x_j = d_x[j] === 0 ? 1e-2 : d_x[j]
                    const mult = 1 - (d_x_j / d_y)
                    for (let d = 0; d < dim; ++d) {
                        delta_mds_stress[d] += (mult * (y_i[d] - y_j[d]));
                    }
                }
            }
            for (let d = 0; d < dim; ++d) {
                Y.set_entry(i, d, Y.entry(i, d) - step_size * delta_mds_stress[d] / N)
            }
        }*/

        // MDS stress step
        const d_Y = new Matrix();
        d_Y.shape = [
            N,
            N,
            (i, j) => {
                return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
            },
        ];
        const ratio = new Matrix(); //d_X.divide(d_Y).mult(-1);
        ratio.shape = [
            N,
            N,
            (i, j) => {
                if (i === j) return 1e-8;
                return i < j ? -d_X.entry(i, j) / d_Y.entry(i, j) : ratio.entry(j, i);
            },
        ];
        for (let i = 0; i < N; ++i) {
            ratio.sub_entry(i, i, neumair_sum(ratio.row(i)));
        }
        const mds_Y = ratio.dot(Y).divide(N);

        // Data depth step
        const diff_Y = new Matrix(N, dim, (i, j) => mds_Y.entry(i, j) - Y.entry(i, j));

        for (let i = 0; i < N; ++i) {
            const m = M(Y.row(i));
            const dm = M(mds_Y.row(i));
            const h_i = h[i];
            for (let d = 0; d < dim; ++d) {
                Y.add_entry(i, d, step_size * (diff_Y.entry(i, d) + w * 2 * (m - h_i) * dm));
            }
        }

        this._Y = Y;

        return this._Y;
    }

    get projection() {
        return this._Y;
    }
}
