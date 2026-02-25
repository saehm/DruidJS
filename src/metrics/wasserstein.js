/**
 * Computes the 1D Wasserstein distance (Earth Mover's Distance) between two distributions.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First distribution (histogram or probability mass)
 * @param {number[] | Float64Array} b - Second distribution (histogram or probability mass)
 * @returns {number} The Wasserstein/EMD distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Wasserstein_metric}
 */
export function wasserstein(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sumA = 0;
    let sumB = 0;
    for (let i = 0; i < n; i++) {
        sumA += a[i];
        sumB += b[i];
    }

    // Fallback if sums are 0
    if (sumA === 0 && sumB === 0) return 0;
    if (sumA === 0 || sumB === 0) return Infinity;

    let distance = 0;
    let cumA = 0;
    let cumB = 0;
    for (let i = 0; i < n; i++) {
        cumA += a[i] / sumA;
        cumB += b[i] / sumB;
        distance += Math.abs(cumA - cumB);
    }
    return distance;
}
