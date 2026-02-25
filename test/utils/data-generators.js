/**
 * Helper function to generate simple test data (small for fast tests)
 * @param {number} n - Number of points
 * @param {number} d - Number of dimensions
 * @returns {number[][]}
 */
export function generateTestData(n = 10, d = 4) {
    const data = [];
    for (let i = 0; i < n; i++) {
        const row = [];
        for (let j = 0; j < d; j++) {
            row.push(Math.sin(i * 0.5 + j) + Math.cos(i * 0.3 - j));
        }
        data.push(row);
    }
    return data;
}

/**
 * Helper function to generate clustered test data
 * @param {number} seed - Seed for reproducibility
 * @returns {number[][]}
 */
export function generateClusteredData(seed = 42) {
    const data = [];
    // Cluster 1 around (0, 0, 0, 0)
    for (let i = 0; i < 5; i++) {
        data.push([
            Math.sin(seed + i) * 0.5,
            Math.cos(seed + i) * 0.5,
            Math.sin(seed + i * 2) * 0.5,
            Math.cos(seed + i * 2) * 0.5,
        ]);
    }
    // Cluster 2 around (10, 10, 10, 10)
    for (let i = 0; i < 5; i++) {
        data.push([
            10 + Math.sin(seed + i) * 0.5,
            10 + Math.cos(seed + i) * 0.5,
            10 + Math.sin(seed + i * 2) * 0.5,
            10 + Math.cos(seed + i * 2) * 0.5,
        ]);
    }
    return data;
}
