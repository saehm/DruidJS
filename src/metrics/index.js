/** @category Metrics */
/**
 * @callback Metric
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number}
 */

/** @module metric */

export { bray_curtis } from "./bray_curtis.js";
export { canberra } from "./canberra.js";
export { chebyshev } from "./chebyshev.js";
export { cosine } from "./cosine.js";
export { euclidean } from "./euclidean.js";
export { euclidean_squared } from "./euclidean_squared.js";
export { goodman_kruskal } from "./goodman_kruskal.js";
export { hamming } from "./hamming.js";
export { haversine } from "./haversine.js";
export { jaccard } from "./jaccard.js";
export { manhattan } from "./manhattan.js";
export { sokal_michener } from "./sokal_michener.js";
export { wasserstein } from "./wasserstein.js";
export { yule } from "./yule.js";
