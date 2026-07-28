// src/dimred/ISOMAP.as.ts

/**
 * Native AssemblyScript Min-Heap for Dijkstra Shortest Path Search.
 */
class HeapItem {
    index: i32;
    dist: f64;
    constructor(idx: i32, d: f64) {
        this.index = idx;
        this.dist = d;
    }
}

class MinHeap {
    private data: Array<HeapItem> = new Array<HeapItem>();

    push(item: HeapItem): void {
        this.data.push(item);
        this.up(this.data.length - 1);
    }

    pop(): HeapItem | null {
        if (this.data.length == 0) return null;
        const top = this.data[0];
        const last = this.data.pop();
        if (this.data.length > 0 && last != null) {
            this.data[0] = last;
            this.down(0);
        }
        return top;
    }

    get isEmpty(): boolean {
        return this.data.length == 0;
    }

    private up(idx: i32): void {
        while (idx > 0) {
            const parent = (idx - 1) >> 1;
            if (this.data[idx].dist < this.data[parent].dist) {
                const tmp = this.data[idx];
                this.data[idx] = this.data[parent];
                this.data[parent] = tmp;
                idx = parent;
            } else {
                break;
            }
        }
    }

    private down(idx: i32): void {
        const len = this.data.length;
        while ((idx << 1) + 1 < len) {
            let left = (idx << 1) + 1;
            let right = left + 1;
            let smallest = idx;

            if (left < len && this.data[left].dist < this.data[smallest].dist) {
                smallest = left;
            }
            if (right < len && this.data[right].dist < this.data[smallest].dist) {
                smallest = right;
            }
            if (smallest != idx) {
                const tmp = this.data[idx];
                this.data[idx] = this.data[smallest];
                this.data[smallest] = tmp;
                idx = smallest;
            } else {
                break;
            }
        }
    }
}

/**
 * Repeated Dijkstra Shortest Path WASM SIMD Kernel for ISOMAP.
 * Operates on k-NN adjacency list buffers.
 */
export function dijkstra_apsp_f64(
    neighbors_ptr: usize, // N x k (i32)
    distances_ptr: usize, // N x k (f64)
    out_d_ptr: usize,     // N x N (f64)
    n: i32,
    k: i32
): void {
    dijkstra_apsp_range_f64(neighbors_ptr, distances_ptr, out_d_ptr, n, k, 0, n);
}

/**
 * Range-based Parallel WASM Min-Heap Repeated Dijkstra Kernel for ISOMAP.
 * Allows worker threads to process a slice of source nodes [start_src, end_src) concurrently.
 */
export function dijkstra_apsp_range_f64(
    neighbors_ptr: usize,
    distances_ptr: usize,
    out_d_ptr: usize,
    n: i32,
    k: i32,
    start_src: i32,
    end_src: i32
): void {
    for (let src = start_src; src < end_src; ++src) {
        const src_n = src * n;
        for (let j = 0; j < n; ++j) {
            store<f64>(out_d_ptr + ((src_n + j) << 3), F64.POSITIVE_INFINITY);
        }
        store<f64>(out_d_ptr + ((src_n + src) << 3), 0.0);

        const heap = new MinHeap();
        heap.push(new HeapItem(src, 0.0));

        while (!heap.isEmpty) {
            const top = heap.pop();
            if (top == null) break;

            const u = top.index;
            const dist_u = top.dist;
            const current_dist = load<f64>(out_d_ptr + ((src_n + u) << 3));

            if (dist_u > current_dist) continue;

            const u_k = u * k;
            for (let neighbor_idx = 0; neighbor_idx < k; ++neighbor_idx) {
                const v = load<i32>(neighbors_ptr + ((u_k + neighbor_idx) << 2));
                if (v < 0 || v >= n) continue;
                const edge_weight = load<f64>(distances_ptr + ((u_k + neighbor_idx) << 3));
                const alt = dist_u + edge_weight;

                const v_offset = (src_n + v) << 3;
                const current_v_dist = load<f64>(out_d_ptr + v_offset);

                if (alt < current_v_dist) {
                    store<f64>(out_d_ptr + v_offset, alt);
                    heap.push(new HeapItem(v, alt));
                }
            }
        }
    }
}
