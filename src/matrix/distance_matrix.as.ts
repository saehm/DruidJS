// src/matrix/distance_matrix.as.ts

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
    const d_simd = d - 1;

    for (let i = 0; i < n; ++i) {
        const i_d = i * d;
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = i + 1; j < n; ++j) {
            const j_d = j * d;
            let sum: f64 = 0.0;
            let k = 0;

            let sum_v = f64x2.splat(0.0);

            for (; k < d_simd; k += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + k) << 3));
                const xj_v = v128.load(x_ptr + ((j_d + k) << 3));
                const diff = f64x2.sub(xi_v, xj_v);
                const diff_sq = f64x2.mul(diff, diff);
                sum_v = f64x2.add(sum_v, diff_sq);
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < d; ++k) {
                const diff = load<f64>(x_ptr + ((i_d + k) << 3)) - load<f64>(x_ptr + ((j_d + k) << 3));
                sum += diff * diff;
            }

            const dist = Math.sqrt(sum);
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
            store<f64>(out_ptr + ((j * n + i) << 3), dist);
        }
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
    const d_simd = d - 1;
    const row_len = k << 1;

    for (let i = start_row; i < end_row; ++i) {
        const i_d = i * d;
        const dist_base = out_ptr + ((i * row_len) << 3);
        const idx_base = dist_base + (k << 3);

        let size = 0;

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;

            const j_d = j * d;
            let sum: f64 = 0.0;
            let c = 0;

            let sum_v = f64x2.splat(0.0);

            for (; c < d_simd; c += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + c) << 3));
                const xj_v = v128.load(x_ptr + ((j_d + c) << 3));
                const diff = f64x2.sub(xi_v, xj_v);
                sum_v = f64x2.add(sum_v, f64x2.mul(diff, diff));
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; c < d; ++c) {
                const diff = load<f64>(x_ptr + ((i_d + c) << 3)) - load<f64>(x_ptr + ((j_d + c) << 3));
                sum += diff * diff;
            }

            const dist = Math.sqrt(sum);

            if (size < k) {
                // Grow the heap, then sift the new leaf up.
                store<f64>(dist_base + (size << 3), dist);
                store<f64>(idx_base + (size << 3), <f64>j);
                let ch = size++;
                while (ch > 0) {
                    const pa = (ch - 1) >> 1;
                    if (load<f64>(dist_base + (pa << 3)) >= load<f64>(dist_base + (ch << 3))) break;
                    const td = load<f64>(dist_base + (pa << 3));
                    store<f64>(dist_base + (pa << 3), load<f64>(dist_base + (ch << 3)));
                    store<f64>(dist_base + (ch << 3), td);
                    const ti = load<f64>(idx_base + (pa << 3));
                    store<f64>(idx_base + (pa << 3), load<f64>(idx_base + (ch << 3)));
                    store<f64>(idx_base + (ch << 3), ti);
                    ch = pa;
                }
            } else if (dist < load<f64>(dist_base)) {
                // Replace the worst and sift down.
                store<f64>(dist_base, dist);
                store<f64>(idx_base, <f64>j);
                let pa = 0;
                while (true) {
                    const l = 2 * pa + 1;
                    const r = l + 1;
                    let m = pa;
                    if (l < k && load<f64>(dist_base + (l << 3)) > load<f64>(dist_base + (m << 3))) m = l;
                    if (r < k && load<f64>(dist_base + (r << 3)) > load<f64>(dist_base + (m << 3))) m = r;
                    if (m == pa) break;
                    const td = load<f64>(dist_base + (pa << 3));
                    store<f64>(dist_base + (pa << 3), load<f64>(dist_base + (m << 3)));
                    store<f64>(dist_base + (m << 3), td);
                    const ti = load<f64>(idx_base + (pa << 3));
                    store<f64>(idx_base + (pa << 3), load<f64>(idx_base + (m << 3)));
                    store<f64>(idx_base + (m << 3), ti);
                    pa = m;
                }
            }
        }

        for (let p = size; p < k; ++p) {
            store<f64>(dist_base + (p << 3), Infinity);
            store<f64>(idx_base + (p << 3), -1.0);
        }

        // Heapsort in place: repeatedly move the max to the back, leaving ascending order.
        for (let end = size - 1; end > 0; --end) {
            const td = load<f64>(dist_base);
            store<f64>(dist_base, load<f64>(dist_base + (end << 3)));
            store<f64>(dist_base + (end << 3), td);
            const ti = load<f64>(idx_base);
            store<f64>(idx_base, load<f64>(idx_base + (end << 3)));
            store<f64>(idx_base + (end << 3), ti);

            let pa = 0;
            while (true) {
                const l = 2 * pa + 1;
                const r = l + 1;
                let m = pa;
                if (l < end && load<f64>(dist_base + (l << 3)) > load<f64>(dist_base + (m << 3))) m = l;
                if (r < end && load<f64>(dist_base + (r << 3)) > load<f64>(dist_base + (m << 3))) m = r;
                if (m == pa) break;
                const sd = load<f64>(dist_base + (pa << 3));
                store<f64>(dist_base + (pa << 3), load<f64>(dist_base + (m << 3)));
                store<f64>(dist_base + (m << 3), sd);
                const si = load<f64>(idx_base + (pa << 3));
                store<f64>(idx_base + (pa << 3), load<f64>(idx_base + (m << 3)));
                store<f64>(idx_base + (m << 3), si);
                pa = m;
            }
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
    const d_simd = d - 1;

    for (let i = start_row; i < end_row; ++i) {
        const i_d = i * d;
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const j_d = j * d;
            let sum: f64 = 0.0;
            let k = 0;

            let sum_v = f64x2.splat(0.0);

            for (; k < d_simd; k += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + k) << 3));
                const xj_v = v128.load(x_ptr + ((j_d + k) << 3));
                const diff = f64x2.sub(xi_v, xj_v);
                const diff_sq = f64x2.mul(diff, diff);
                sum_v = f64x2.add(sum_v, diff_sq);
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < d; ++k) {
                const diff = load<f64>(x_ptr + ((i_d + k) << 3)) - load<f64>(x_ptr + ((j_d + k) << 3));
                sum += diff * diff;
            }

            const dist = Math.sqrt(sum);
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
        }
    }
}
