import { WASM_BASE64 } from "./wasm_bytes.js";

/**
 * Worker pool that splits a row-range kernel across CPU cores.
 *
 * Only Node is supported. Blocking a caller until the workers finish needs `Atomics.wait`, which
 * throws on a browser main thread, and the shared output needs `SharedArrayBuffer`, which browsers
 * gate behind COOP/COEP headers. {@link parallel_available} reports false there and every caller
 * falls back to its single-threaded path.
 *
 * @module wasm/worker_pool
 */

/**
 * Runs inside each worker. Kept as a string with no imports so that bundling the library does not
 * have to emit a separate worker chunk — the source travels inside the bundle and is spawned with
 * `eval: true`.
 *
 * @type {string}
 */
const WORKER_SOURCE = `
import { parentPort, workerData } from "node:worker_threads";

const { wasm_base64, control } = workerData;
const flags = new Int32Array(control);

const binary = Uint8Array.from(atob(wasm_base64), (c) => c.charCodeAt(0));
const instance = new WebAssembly.Instance(new WebAssembly.Module(binary), {
    env: { abort: () => { throw new Error("wasm abort"); } },
});
const exports = instance.exports;

parentPort.on("message", (task) => {
    const ptrs = [];
    try {
        const memory = exports.memory;
        const args = [];
        for (const input of task.inputs) {
            const ptr = exports.allocate(input.data.byteLength);
            ptrs.push(ptr);
            if (input.kind === "i32") {
                new Int32Array(memory.buffer, ptr, input.data.length).set(input.data);
            } else {
                new Float64Array(memory.buffer, ptr, input.data.length).set(input.data);
            }
            args.push(ptr);
        }

        const out_ptr = exports.allocate(task.out_bytes);
        ptrs.push(out_ptr);
        args.push(out_ptr);

        const returned = exports[task.kernel](...args, ...task.scalars, task.start, task.end);

        // Only this worker's rows: the shared output is written by every worker at once.
        const offset = task.start * task.row_length;
        const count = (task.end - task.start) * task.row_length;
        new Float64Array(task.out).set(new Float64Array(memory.buffer, out_ptr + offset * 8, count), offset);

        // Kernels that reduce (a max shift, a loss) hand their partial result back per slot.
        if (task.returns) new Float64Array(task.returns)[task.slot] = typeof returned === "number" ? returned : 0;
    } catch (error) {
        Atomics.store(flags, 1, 1);
        parentPort.postMessage(String(error));
    } finally {
        for (const ptr of ptrs) exports.free(ptr);
        Atomics.add(flags, 0, 1);
        Atomics.notify(flags, 0);
    }
});

Atomics.add(flags, 2, 1);
Atomics.notify(flags, 2);
`;

/** @type {any[] | null} */
let workers = null;
/** @type {Int32Array | null} */
let control = null;
/** @type {boolean} */
let unavailable = false;

/**
 * Number of workers the pool will start.
 *
 * @returns {number}
 */
function worker_count() {
    try {
        const os = /** @type {any} */ (globalThis).process?.getBuiltinModule?.("node:os");
        // One worker per core minus the caller's own thread, which blocks rather than computing.
        return Math.max(1, Math.min(8, (os?.availableParallelism?.() ?? 4) - 1));
    } catch {
        return 4;
    }
}

/**
 * Whether row-range kernels can be split across workers in this environment.
 *
 * @returns {boolean} True only on Node, where `Atomics.wait` may block the calling thread.
 * @example
 * import { parallel_available } from "@saehrimnir/druidjs";
 * parallel_available(); // false in a browser
 */
export function parallel_available() {
    if (unavailable) return false;
    try {
        const p = /** @type {any} */ (globalThis).process;
        if (!p?.versions?.node || typeof SharedArrayBuffer === "undefined") return false;
        // Opt-out, so the single-threaded path can be exercised without changing any code.
        const off = p.env?.DRUID_DISABLE_PARALLEL;
        if (off !== undefined && off !== "" && off !== "0" && off !== "false") return false;
        // `Atomics.wait` throws on a browser main thread; the pool has no way to block there.
        return p.getBuiltinModule("node:worker_threads").isMainThread === true;
    } catch {
        return false;
    }
}

/**
 * Starts the workers, or returns the running ones. Blocks until every worker has instantiated its
 * own copy of the module, so the first task does not pay for start-up.
 *
 * @private
 * @returns {any[] | null} The workers, or null if the pool could not start.
 */
function ensure_pool() {
    if (workers) return workers;
    if (!parallel_available()) return null;

    try {
        const { Worker } = /** @type {any} */ (globalThis).process.getBuiltinModule("node:worker_threads");
        const n = worker_count();
        control = new Int32Array(new SharedArrayBuffer(12));

        const started = [];
        for (let i = 0; i < n; ++i) {
            started.push(
                new Worker(WORKER_SOURCE, {
                    eval: true,
                    workerData: { wasm_base64: WASM_BASE64, control: control.buffer },
                }),
            );
        }

        // Wait for every worker to report ready, so start-up is not billed to the first task.
        const deadline = Date.now() + 10000;
        while (Atomics.load(control, 2) < n) {
            if (Date.now() > deadline) {
                for (const w of started) w.terminate();
                unavailable = true;
                return null;
            }
            Atomics.wait(control, 2, Atomics.load(control, 2), 50);
        }

        for (const w of started) w.unref();
        workers = started;
        return workers;
    } catch {
        unavailable = true;
        return null;
    }
}

/**
 * Splits a row-range kernel across the pool and blocks until every worker is done.
 *
 * The kernel is called as `(...input_ptrs, out_ptr, ...scalars, start, end)`, and each worker copies
 * back only the rows it produced. `out` must be backed by a `SharedArrayBuffer` so the workers write
 * into the caller's memory.
 *
 * @param {string} kernel - Name of the exported range kernel.
 * @param {{ data: Float64Array | Int32Array; kind: "f64" | "i32" }[]} inputs - Copied into each worker.
 * @param {Float64Array} out - Shared output, `total_rows * row_length` long.
 * @param {number[]} scalars - Kernel arguments between the output pointer and the row range.
 * @param {number} total_rows - Rows to divide between workers.
 * @param {number} row_length - Output entries per row.
 * @param {Float64Array} [returns] - Shared, one slot per worker, for kernels that return a value.
 * @returns {boolean} True if the pool ran the kernel; false if the caller must fall back to JS.
 */
export function run_row_range(kernel, inputs, out, scalars, total_rows, row_length, returns) {
    const pool = ensure_pool();
    if (!pool || !control) return false;
    if (!(out.buffer instanceof SharedArrayBuffer)) return false;

    const n = Math.min(pool.length, total_rows);
    if (n < 1) return false;

    Atomics.store(control, 0, 0);
    Atomics.store(control, 1, 0);

    const chunk = Math.ceil(total_rows / n);
    let dispatched = 0;
    for (let i = 0; i < n; ++i) {
        const start = i * chunk;
        const end = Math.min((i + 1) * chunk, total_rows);
        if (start >= end) break;
        pool[i].postMessage({
            kernel,
            inputs,
            out: out.buffer,
            out_bytes: out.byteLength,
            scalars,
            start,
            end,
            row_length,
            returns: returns?.buffer,
            slot: i,
        });
        ++dispatched;
    }

    const deadline = Date.now() + 120000;
    while (Atomics.load(control, 0) < dispatched) {
        if (Date.now() > deadline) return false;
        Atomics.wait(control, 0, Atomics.load(control, 0), 100);
    }

    return Atomics.load(control, 1) === 0;
}

/**
 * Stops the workers. They are unreferenced, so this is only needed to release them early.
 *
 * @returns {void}
 */
export function terminate_pool() {
    if (workers) {
        for (const w of workers) w.terminate();
        workers = null;
        control = null;
    }
}
