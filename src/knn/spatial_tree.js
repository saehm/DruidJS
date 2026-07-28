import { chebyshev, euclidean, manhattan } from "../metrics/index.js";
import { BallTree } from "./BallTree.js";
import { KDTree } from "./KDTree.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * Metrics for which a KD-Tree may prune soundly.
 *
 * `KDTree` discards a subtree when the distance to the splitting hyperplane, `|a[axis] - b[axis]|`,
 * already exceeds the current k-th best distance. That is only valid for metrics which are bounded
 * below by the per-coordinate difference — true for the L_p family, but false for e.g. `cosine`,
 * `canberra`, `haversine`, or `euclidean_squared`, where it silently drops true neighbors.
 *
 * @type {Set<Metric>}
 */
const KDTREE_SAFE_METRICS = new Set([euclidean, manhattan, chebyshev]);

/**
 * Picks the fastest spatial index that is still correct for the given metric.
 *
 * Returns a {@link KDTree} for metrics it can prune soundly, and a {@link BallTree} — which only
 * relies on the triangle inequality and therefore works for any metric — otherwise.
 *
 * @template {number[] | Float64Array} T
 * @param {T[]} elements - Elements to index.
 * @param {{ metric: Metric; seed: number }} parameters - Metric and seed for the tree.
 * @returns {KDTree<T> | BallTree<T>} The constructed index.
 */
export function spatial_tree(elements, parameters) {
    if (KDTREE_SAFE_METRICS.has(parameters.metric)) {
        return new KDTree(elements, parameters);
    }
    return new BallTree(elements, parameters);
}
