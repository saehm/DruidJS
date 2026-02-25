/**
 * Computes the Haversine distance between two points on a sphere of unit length 1. Multiply the result with the radius of the sphere. (For instance Earth's radius is 6371km)
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - Point [lat1, lon1] in radians
 * @param {number[] | Float64Array} b - Point [lat2, lon2] in radians
 * @returns {number} The Haversine distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Haversine_formula}
 */
export function haversine(a, b) {
    if (a.length !== 2 || b.length !== 2)
        throw new Error("Haversine distance requires exactly 2 coordinates [lat, lon] for each point!");
    const lat1 = a[0];
    const lon1 = a[1];
    const lat2 = b[0];
    const lon2 = b[1];

    const dlat = lat2 - lat1;
    const dlon = lon2 - lon1;

    const sin_dlat2 = Math.sin(dlat / 2);
    const sin_dlon2 = Math.sin(dlon / 2);

    const x = sin_dlat2 * sin_dlat2 + Math.cos(lat1) * Math.cos(lat2) * sin_dlon2 * sin_dlon2;
    const c = 2 * Math.atan2(Math.sqrt(x), Math.sqrt(1 - x));

    return c;
}
