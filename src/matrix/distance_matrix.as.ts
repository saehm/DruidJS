// src/matrix/distance_matrix.as.ts

import { sqdist_simd_f64 } from "../wasm/shared.as";

/**
 * Pairwise Euclidean distance matrix computation for X (n x d).
 * Computes symmetric n x n distance matrix D.
 *
 * Only the upper triangle is computed; each distance is mirrored into the lower triangle, which
 * halves the work compared to visiting every ordered pair.
 */
export function euclidean_distance_matrix_f64(
    x_ptr: usize,
    out_ptr: usize,
    n: i32,
    d: i32
): void {
    for (let i = 0; i < n; ++i) {
        const row_i = x_ptr + ((i * d) << 3);
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = i + 1; j < n; ++j) {
            const dist = Math.sqrt(sqdist_simd_f64(row_i, x_ptr + ((j * d) << 3), d));
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
            store<f64>(out_ptr + ((j * n + i) << 3), dist);
        }
    }
}

/**
 * Swaps entries `a` and `b` of the parallel distance and index arrays.
 */
function swap_pair(dist_base: usize, idx_base: usize, a: i32, b: i32): void {
    const a_off = a << 3;
    const b_off = b << 3;

    const dist_a = load<f64>(dist_base + a_off);
    store<f64>(dist_base + a_off, load<f64>(dist_base + b_off));
    store<f64>(dist_base + b_off, dist_a);

    const idx_a = load<f64>(idx_base + a_off);
    store<f64>(idx_base + a_off, load<f64>(idx_base + b_off));
    store<f64>(idx_base + b_off, idx_a);
}

/**
 * Restores the max-heap property from `root` downwards, over the first `len` entries.
 */
function sift_down(dist_base: usize, idx_base: usize, root: i32, len: i32): void {
    let parent = root;
    while (true) {
        const left = 2 * parent + 1;
        const right = left + 1;
        let largest = parent;

        if (left < len && load<f64>(dist_base + (left << 3)) > load<f64>(dist_base + (largest << 3))) {
            largest = left;
        }
        if (right < len && load<f64>(dist_base + (right << 3)) > load<f64>(dist_base + (largest << 3))) {
            largest = right;
        }
        if (largest == parent) break;

        swap_pair(dist_base, idx_base, parent, largest);
        parent = largest;
    }
}

/**
 * Exact k-nearest-neighbour search for rows [start_row, end_row), selection included.
 *
 * The whole search loop lives here rather than just the metric, so the JS/WASM boundary is paid
 * once per block instead of once per pair -- the only shape in which these kernels beat a JS KNN
 * structure. Nothing proportional to n*n crosses back: a bounded max-heap keeps the k best
 * candidates per row and only those are written out.
 *
 * The heap is built directly in the output slot, so the kernel allocates nothing. Distances and
 * indices share one f64 buffer -- indices are exact in f64 well past any workable n -- because the
 * worker pool splits a single output array. Rows are written at their absolute offset so a worker
 * copying back only its own range lands in the right place.
 *
 * Output row `i` is `2 * k` wide: k ascending distances, then their k indices. `i` itself is never
 * a candidate. Rows with fewer than k other points are padded with an index of -1.
 */
export function euclidean_knn_block_f64(
    x_ptr: usize,
    out_ptr: usize,
    n: i32,
    d: i32,
    k: i32,
    start_row: i32,
    end_row: i32
): void {
    const row_len = k << 1;

    for (let i = start_row; i < end_row; ++i) {
        const row_i = x_ptr + ((i * d) << 3);
        const dist_base = out_ptr + ((i * row_len) << 3);
        const idx_base = dist_base + (k << 3);

        let size = 0;

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;

            const dist = Math.sqrt(sqdist_simd_f64(row_i, x_ptr + ((j * d) << 3), d));

            if (size < k) {
                // Grow the heap, then sift the new leaf up.
                store<f64>(dist_base + (size << 3), dist);
                store<f64>(idx_base + (size << 3), <f64>j);
                let child = size++;
                while (child > 0) {
                    const parent = (child - 1) >> 1;
                    if (load<f64>(dist_base + (parent << 3)) >= load<f64>(dist_base + (child << 3))) break;
                    swap_pair(dist_base, idx_base, parent, child);
                    child = parent;
                }
            } else if (dist < load<f64>(dist_base)) {
                // Replace the worst and sift down.
                store<f64>(dist_base, dist);
                store<f64>(idx_base, <f64>j);
                sift_down(dist_base, idx_base, 0, k);
            }
        }

        for (let p = size; p < k; ++p) {
            store<f64>(dist_base + (p << 3), Infinity);
            store<f64>(idx_base + (p << 3), -1.0);
        }

        // Heapsort in place: repeatedly move the max to the back, leaving ascending order.
        for (let end = size - 1; end > 0; --end) {
            swap_pair(dist_base, idx_base, 0, end);
            sift_down(dist_base, idx_base, 0, end);
        }
    }
}

/**
 * Multi-threaded Range-based Pairwise Euclidean Distance Matrix Kernel.
 * Allows worker threads to compute distance matrix rows [start_row, end_row) in parallel.
 */
export function euclidean_distance_matrix_range_f64(
    x_ptr: usize,
    out_ptr: usize,
    n: i32,
    d: i32,
    start_row: i32,
    end_row: i32
): void {
    for (let i = start_row; i < end_row; ++i) {
        const row_i = x_ptr + ((i * d) << 3);
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const dist = Math.sqrt(sqdist_simd_f64(row_i, x_ptr + ((j * d) << 3), d));
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
        }
    }
}
