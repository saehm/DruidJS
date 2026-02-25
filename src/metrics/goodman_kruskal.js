/**
 * Computes the Goodman-Kruskal gamma coefficient for ordinal association.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a - First categorical/ordinal variable
 * @param {number[] | Float64Array} b - Second categorical/ordinal variable
 * @returns {number} The Goodman-Kruskal gamma coefficient between `a` and `b` (-1 to 1).
 * @see {@link https://en.wikipedia.org/wiki/Goodman_and_Kruskal%27s_gamma}
 */
export function goodman_kruskal(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    if (n < 2) return 0;

    let concordant = 0;
    let discordant = 0;
    let tie_a = 0;
    let tie_b = 0;

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const a_diff = a[i] - a[j];
            const b_diff = b[i] - b[j];
            const a_tied = a_diff === 0;
            const b_tied = b_diff === 0;

            if (a_tied && b_tied) {
            } else if (a_tied) {
                tie_a++;
            } else if (b_tied) {
                tie_b++;
            } else if (a_diff * b_diff > 0) {
                concordant++;
            } else {
                discordant++;
            }
        }
    }

    const denominator = concordant + discordant + tie_a + tie_b;
    if (denominator === 0) return 0;

    const numerator = concordant + discordant;
    if (numerator === 0) return 0;

    return (concordant - discordant) / numerator;
}
