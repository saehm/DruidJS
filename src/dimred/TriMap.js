import { Matrix, linspace } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { PCA } from "./PCA.js";
import { BallTree } from "../knn/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias TriMap
 * @extends DR
 */
export class TriMap extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.weight_adj = 500] - scaling factor.
     * @param {number} [parameters.c = 5] - number of triplets multiplier.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {number} [parameters.tol = 1e-8] -
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X, parameters) {
        super(X, { weight_adj: 500, c: 5, d: 2, metric: euclidean, tol: 1e-8, seed: 1212 }, parameters);
        return this;
    }

    /**
     *
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used).
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used).
     */
    init(pca = null, knn = null) {
        const X = this.X;
        const N = X.shape[0];
        const { c, d, metric, seed } = this._parameters;
        this.n_inliers = 2 * c;
        this.n_outliers = 1 * c;
        this.n_random = 1 * c;
        this.Y = pca || new PCA(X, { d, seed }).transform();
        this.knn = knn || new BallTree(X.to2dArray, metric);
        const { triplets, weights } = this._generate_triplets(this.n_inliers, this.n_outliers, this.n_random);
        this.triplets = triplets;
        this.weights = weights;
        this.lr = (1000 * N) / triplets.shape[0];
        this.C = Infinity;
        this.vel = new Matrix(N, d, 0);
        this.gain = new Matrix(N, d, 1);
        return this;
    }

    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {number} n_inliers
     * @param {number} n_outliers
     * @param {number} n_random
     */
    _generate_triplets(n_inliers, n_outliers, n_random) {
        const { metric, weight_adj } = this._parameters;
        const X = this.X;
        const N = X.shape[0];
        const knn = this.knn;
        const n_extra = Math.min(n_inliers + 20, N);
        const nbrs = new Matrix(N, n_extra);
        const knn_distances = new Matrix(N, n_extra);
        for (let i = 0; i < N; ++i) {
            knn.search(X.row(i), n_extra + 1)
                .raw_data()
                .filter((d) => d.value != 0)
                .sort((a, b) => a.value - b.value)
                .forEach((d, j) => {
                    nbrs.set_entry(i, j, d.element.index);
                    knn_distances.set_entry(i, j, d.value);
                });
        }
        // scale parameter
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            sig[i] = Math.max((knn_distances.entry(i, 3) + knn_distances.entry(i, 4) + knn_distances.entry(i, 5) + knn_distances.entry(i, 6)) / 4, 1e-10);
        }

        const P = this._find_p(knn_distances, sig, nbrs);

        let triplets = this._sample_knn_triplets(P, nbrs, n_inliers, n_outliers);
        let n_triplets = triplets.shape[0];
        const outlier_distances = new Float64Array(n_triplets);
        for (let i = 0; i < n_triplets; ++i) {
            const j = triplets.entry(i, 0);
            const k = triplets.entry(i, 2);
            outlier_distances[i] = metric(X.row(j), X.row(k));
        }
        let weights = this._find_weights(triplets, P, nbrs, outlier_distances, sig);

        if (n_random > 0) {
            const { random_triplets, random_weights } = this._sample_random_triplets(X, n_random, sig);
            triplets = triplets.concat(random_triplets, "vertical");
            weights = Float64Array.from([...weights, ...random_weights]);
        }
        n_triplets = triplets.shape[0];
        let max_weight = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            if (isNaN(weights[i])) {
                weights[i] = 0;
            }
            if (max_weight < weights[i]) max_weight = weights[i];
        }
        let max_weight_2 = -Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight;
            weights[i] += 0.0001;
            weights[i] = Math.log(1 + weight_adj * weights[i]);
            if (max_weight_2 < weights[i]) max_weight_2 = weights[i];
        }
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] /= max_weight_2;
        }
        return {
            triplets: triplets,
            weights: weights,
        };
    }

    /**
     * Calculates the similarity matrix P
     * @private
     * @param {Matrix} knn_distances - matrix of pairwise knn distances
     * @param {Float64Array} sig - scaling factor for the distances
     * @param {Matrix} nbrs - nearest neighbors
     * @returns {Matrix} pairwise similarity matrix
     */
    _find_p(knn_distances, sig, nbrs) {
        const [N, n_neighbors] = knn_distances.shape;
        return new Matrix(N, n_neighbors, (i, j) => {
            return Math.exp(-(knn_distances.entry(i, j) ** 2 / sig[i] / sig[nbrs.entry(i, j)]));
        });
    }

    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix {@link P}. Row i corresponds to the i-th point.
     * @param {number} n_inliers - number of inlier points.
     * @param {number} n_outliers - number of outlier points.
     *
     */
    _sample_knn_triplets(P, nbrs, n_inliers, n_outliers) {
        const N = nbrs.shape[0];
        const triplets = new Matrix(N * n_inliers * n_outliers, 3);
        for (let i = 0; i < N; ++i) {
            let n_i = i * n_inliers * n_outliers;
            const sort_indices = this.__argsort(P.row(i));
            for (let j = 0; j < n_inliers; ++j) {
                let n_j = j * n_outliers;
                const sim = nbrs.entry(i, sort_indices[j]);
                const samples = this._rejection_sample(n_outliers, N, sort_indices.slice(0, j + 1));
                for (let k = 0; k < n_outliers; ++k) {
                    const index = n_i + n_j + k;
                    const out = samples[k];
                    triplets.set_entry(index, 0, i);
                    triplets.set_entry(index, 1, sim);
                    triplets.set_entry(index, 2, out);
                }
            }
        }
        return triplets;
    }

    /**
     * Should do the same as np.argsort()
     * @private
     * @param {number[]} A
     */
    __argsort(A) {
        return linspace(0, A.length - 1).sort((i, j) => A[j] - A[i]);
    }

    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are in the {@link rejects}.
     * @private
     * @param {*} n_samples
     * @param {*} max_int
     * @param {*} rejects
     */
    _rejection_sample(n_samples, max_int, rejects) {
        const randomizer = this._randomizer;
        const interval = linspace(0, max_int - 1).filter((d) => rejects.indexOf(d) < 0);
        return randomizer.choice(interval, Math.min(n_samples, interval.length - 2));
    }

    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - scaling factor for the distances.
     */
    _find_weights(triplets, P, nbrs, outlier_distances, sig) {
        const n_triplets = triplets.shape[0];
        const weights = new Float64Array(n_triplets);
        for (let t = 0; t < n_triplets; ++t) {
            const i = triplets.entry(t, 0);
            const sim = nbrs.row(i).indexOf(triplets.entry(t, 1));
            const p_sim = P.entry(i, sim);
            let p_out = Math.exp(-(outlier_distances[t] ** 2 / (sig[i] * sig[triplets.entry(t, 2)])));
            if (p_out < 1e-20) p_out = 1e-20;
            weights[t] = p_sim / p_out;
        }
        return weights;
    }

    /**
     * Sample uniformly ranom triplets
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {number} n_random - number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    _sample_random_triplets(X, n_random, sig) {
        const metric = this.parameter("metric");
        const randomizer = this._randomizer;
        const N = X.shape[0];
        const random_triplets = new Matrix(N * n_random, 3);
        const random_weights = new Float64Array(N * n_random);
        for (let i = 0; i < N; ++i) {
            const n_i = i * n_random;
            const indices = [...linspace(0, i - 1), ...linspace(i + 1, N - 1)];
            for (let j = 0; j < n_random; ++j) {
                let [sim, out] = randomizer.choice(indices, 2);
                let p_sim = Math.exp(-(metric(X.row(i), X.row(sim)) ** 2 / (sig[i] * sig[sim])));
                if (p_sim < 1e-20) p_sim = 1e-20;
                let p_out = Math.exp(-(metric(X.row(i), X.row(out)) ** 2 / (sig[i] * sig[out])));
                if (p_out < 1e-20) p_out = 1e-20;

                if (p_sim < p_out) {
                    [sim, out] = [out, sim];
                    [p_sim, p_out] = [p_out, p_sim];
                }
                const index = n_i + j;
                random_triplets.set_entry(index, 0, i);
                random_triplets.set_entry(index, 1, sim);
                random_triplets.set_entry(index, 2, out);
                random_weights[index] = p_sim / p_out;
            }
        }
        return {
            random_triplets: random_triplets,
            random_weights: random_weights,
        };
    }

    /**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */
    _grad(Y) {
        const n_inliers = this.n_inliers;
        const n_outliers = this.n_outliers;
        const triplets = this.triplets;
        const weights = this.weights;
        const [N, dim] = Y.shape;
        const n_triplets = triplets.shape[0];
        const grad = new Matrix(N, dim, 0);
        let y_ij = new Float64Array(dim);
        let y_ik = new Float64Array(dim);
        let d_ij = 1;
        let d_ik = 1;
        let n_viol = 0;
        let loss = 0;
        const n_knn_triplets = N * n_inliers * n_outliers;

        for (let t = 0; t < n_triplets; ++t) {
            const [i, j, k] = triplets.row(t);
            // update y_ij, y_ik, d_ij, d_ik
            if (t % n_outliers == 0 || t >= n_knn_triplets) {
                d_ij = 1;
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_jd = Y.entry(j, d);
                    const Y_kd = Y.entry(k, d);
                    y_ij[d] = Y_id - Y_jd;
                    y_ik[d] = Y_id - Y_kd;
                    d_ij += y_ij[d] ** 2;
                    d_ik += y_ik[d] ** 2;
                }
                // update y_ik and d_ik only
            } else {
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_kd = Y.entry(k, d);
                    y_ik[d] = Y_id - Y_kd;
                    d_ik += y_ik[d] ** 2;
                }
            }

            if (d_ij > d_ik) ++n_viol;
            loss += weights[t] / (1 + d_ik / d_ij);
            const w = (weights[t] / (d_ij + d_ik)) ** 2;
            for (let d = 0; d < dim; ++d) {
                const gs = y_ij[d] * d_ik * w;
                const go = y_ik[d] * d_ij * w;
                grad.add_entry(i, d, gs - go);
                grad.sub_entry(j, d, gs);
                grad.add_entry(k, d, go);
            }
        }
        return { grad, loss, n_viol };
    }

    /**
     *
     * @param {number} max_iteration
     * @returns {Matrix|number[][]}
     */
    transform(max_iteration = 400) {
        this.check_init();
        for (let iter = 0; iter < max_iteration; ++iter) {
            this._next(iter);
        }
        return this.projection;
    }

    /**
     * @param {number} max_iteration
     * @yields {Matrix|number[][]}
     * @returns {Matrix|number[][]}
     */
    *generator(max_iteration = 800) {
        this.check_init();
        for (let iter = 0; iter < max_iteration; ++iter) {
            this._next(iter);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * Does the iteration step.
     * @private
     * @param {number} iter
     */
    _next(iter) {
        const gamma = iter > 150 ? 0.5 : 0.3;
        const old_C = this.C;
        const vel = this.vel;
        const Y = this.Y.add(vel.mult(gamma));
        const { grad, loss, n_viol } = this._grad(Y);
        this.C = loss;
        this.Y = this._update_embedding(Y, iter, grad);
        this.lr *= old_C > loss + this._parameters.tol ? 1.01 : 0.9;
        return this.Y;
    }

    /**
     * Updates the embedding.
     * @private
     * @param {Matrix} Y
     * @param {number} iter
     * @param {Matrix} grad
     */
    _update_embedding(Y, iter, grad) {
        const [N, dim] = Y.shape;
        const gamma = iter > 150 ? 0.9 : 0.5; // moment parameter
        const min_gain = 0.01;
        const gain = this.gain;
        const vel = this.vel;
        const lr = this.lr;
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const new_gain = Math.sign(vel.entry(i, d)) != Math.sign(grad.entry(i, d)) ? gain.entry(i, d) + 0.2 : Math.max(gain.entry(i, d) * 0.8, min_gain);
                gain.set_entry(i, d, new_gain);
                vel.set_entry(i, d, gamma * vel.entry(i, d) - lr * gain.entry(i, d) * grad.entry(i, d));
                Y.set_entry(i, d, Y.entry(i, d) + vel.entry(i, d));
            }
        }
        return Y;
    }
}
