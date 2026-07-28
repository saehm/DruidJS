import { WASM_BASE64 } from "./wasm_bytes.js";

/** @type {WebAssembly.Instance | null} */
let wasmInstance = null;

/**
 * Set once the module traps. An AssemblyScript `abort` unwinds straight back to the JS caller and
 * leaves the instance in an undefined state, so it must never be used again — every kernel then
 * reports failure and the caller falls back to its JS implementation.
 *
 * @type {boolean}
 */
let wasmAborted = false;

/**
 * Opt-out switch, see {@link setWasmEnabled}. Reads `DRUID_DISABLE_WASM` from the environment so
 * the JS paths can be exercised without changing any code.
 *
 * @type {boolean}
 */
let wasmEnabled = !read_env_flag("DRUID_DISABLE_WASM");

/**
 * Reads a boolean flag from `process.env`, tolerating environments without a `process` global.
 *
 * @param {string} name
 * @returns {boolean}
 */
function read_env_flag(name) {
    try {
        const env = /** @type {any} */ (globalThis).process?.env;
        const value = env?.[name];
        return value !== undefined && value !== "" && value !== "0" && value !== "false";
    } catch {
        return false;
    }
}

/**
 * Enables or disables WASM acceleration at runtime.
 *
 * With WASM disabled every accelerated function takes its JS fallback path. This exists so the two
 * implementations can be compared in tests and benchmarks.
 *
 * For most of the library the two paths agree to within floating point tolerance — the kernels
 * accumulate in a different order than the scalar loops, and `Math.pow`/`Math.exp` are only
 * "implementation-approximated" in ECMAScript, so the last bit may differ. `PCA`, `MDS`, `TSNE` and
 * `TriMap` come out bit-identical in practice.
 *
 * {@link UMAP} and {@link SAMMON} are the exception: their optimisers are chaotic, so a last-bit
 * difference grows into a visibly different (though equally valid) layout. See the note on
 * reproducibility in those classes.
 *
 * @param {boolean} enabled - Whether kernels may run in WASM.
 * @returns {boolean} The new state.
 * @example
 * import { setWasmEnabled } from "@saehrimnir/druidjs";
 * setWasmEnabled(false); // force the pure JS implementations
 */
export function setWasmEnabled(enabled) {
    wasmEnabled = Boolean(enabled);
    return wasmEnabled;
}

/**
 * Reads an AssemblyScript string out of linear memory, for the abort handler.
 *
 * AssemblyScript strings are UTF-16 with their byte length in the header word immediately before
 * the payload. A null pointer is normal — `assert(cond)` without a message aborts with one.
 *
 * @private
 * @param {number} ptr
 * @returns {string | null} The string, or null if there is nothing readable there.
 */
function read_wasm_string(ptr) {
    if (!ptr) return null;
    try {
        const memory = /** @type {any} */ (wasmInstance).exports.memory;
        const length = new Uint32Array(memory.buffer, ptr - 4, 1)[0] >>> 1;
        return String.fromCharCode(...new Uint16Array(memory.buffer, ptr, length));
    } catch {
        // Reporting the abort matters more than decorating it.
        return null;
    }
}

export function initWasm() {
    if (!wasmEnabled || wasmAborted) return null;
    if (wasmInstance) return wasmInstance;
    try {
        if (typeof WebAssembly === "undefined") return null;
        let binaryString = "";
        const g = /** @type {any} */ (globalThis);
        if (typeof atob === "function") {
            binaryString = atob(WASM_BASE64);
        } else if (g?.Buffer) {
            binaryString = g.Buffer.from(WASM_BASE64, "base64").toString("binary");
        } else {
            return null;
        }
        const bytes = new Uint8Array(binaryString.length);
        for (let i = 0; i < binaryString.length; i++) {
            bytes[i] = binaryString.charCodeAt(i);
        }
        const module = new WebAssembly.Module(bytes);
        wasmInstance = new WebAssembly.Instance(module, {
            env: {
                abort: (
                    /** @type {number} */ msg,
                    /** @type {number} */ _file,
                    /** @type {number} */ line,
                    /** @type {number} */ column,
                ) => {
                    // The instance cannot be trusted after a trap; retire it and fall back to JS.
                    // The location is worth decoding: a bare `assert` passes a null message, so
                    // "~lib/rt/tlsf.ts" is the only thing distinguishing heap corruption caused by
                    // a bad buffer size from an ordinary kernel bug.
                    // Read the location out of the still-live memory before retiring the instance.
                    const where = read_wasm_string(_file) ?? "kernel";
                    const reason = read_wasm_string(msg);
                    retire_instance();
                    console.error(
                        `DruidJS: WASM aborted in ${where} at ${line}:${column}${reason ? `: ${reason}` : ""}, falling back to JS.`,
                    );
                },
            },
        });
        return wasmInstance;
    } catch (_e) {
        return null;
    }
}

/**
 * Returns whether WASM acceleration is available in current environment.
 * @returns {boolean}
 */
export function isWasmAvailable() {
    return initWasm() !== null;
}

/**
 * Returns whether SharedArrayBuffer and WASM multi-threading is supported in current environment (Node.js, Bun, Deno, or Web Browser with COOP/COEP headers).
 * @returns {boolean}
 */
export function isWasmThreadsSupported() {
    try {
        if (typeof SharedArrayBuffer === "undefined") return false;
        const mem = new WebAssembly.Memory({ initial: 1, maximum: 1, shared: true });
        return mem.buffer instanceof SharedArrayBuffer;
    } catch {
        return false;
    }
}

/**
 * Allocates a buffer in the WASM heap and records the pointer in `ptrs`.
 *
 * Every kernel allocates through this and releases through {@link free_all} in a `finally`, so that
 * a failing allocation part-way through a list still frees the ones that already succeeded — WASM
 * linear memory is never reclaimed by the garbage collector, so a leaked pointer is permanent.
 *
 * @private
 * @param {any} exports - The WASM instance exports.
 * @param {number[]} ptrs - Accumulator holding every pointer allocated so far.
 * @param {number} size - Size in bytes.
 * @returns {number} Pointer to the allocated buffer.
 */
export function alloc(exports, ptrs, size) {
    const ptr = exports.allocate(size);
    ptrs.push(ptr);
    return ptr;
}

/**
 * Frees every pointer collected by {@link alloc}.
 *
 * @private
 * @param {any} exports - The WASM instance exports.
 * @param {number[]} ptrs - Pointers to release.
 */
export function free_all(exports, ptrs) {
    for (let i = 0; i < ptrs.length; ++i) {
        exports.free(ptrs[i]);
    }
    ptrs.length = 0;
}

/**
 * Buffers held across calls to an iterative kernel.
 *
 * The optimisers call their kernel once per iteration with operands that mostly do not change.
 * t-SNE is the clearest case: `P`, the joint probability matrix, is built once by `init()` and read
 * by all 500 iterations, and the N ⨯ N scratch matrices are never read from JS at all. Allocating
 * and copying those per call moves roughly 15 GB over a run at N = 2000 to no purpose.
 *
 * A session keeps the pointers alive and re-copies an operand only when the caller hands over a
 * different array. Held buffers are not garbage: WASM linear memory is never reclaimed, so a
 * session is released when the shape changes, and {@link release_wasm_buffers} exists for callers
 * that want the memory back before then.
 *
 * One slot per key rather than one slot overall: SQDMDS calls two kernels per iteration and shares
 * `grads` between them, so a single slot would have them evicting each other every iteration and
 * end up slower than copying.
 *
 * @private
 * @type {Map<string, { inst: object, ptrs: Record<string, number>, sizes: Record<string, number>, sources: Record<string, Float64Array | null> }>}
 */
const sessions = new Map();

/**
 * Retires the current instance after a trap and drops everything that references it.
 *
 * The abort handler calls this the moment a kernel traps; it is also the state a test needs to
 * reach without corrupting the heap, so it is exported (privately) rather than inlined.
 *
 * Clearing `sessions` is the part that is easy to forget and expensive to omit. Every session entry
 * holds `inst`, so leaving them in the Map would pin the dead instance's whole linear memory — tens
 * to hundreds of MB after a large run — for the life of the process: once `wasmAborted` is set,
 * `initWasm` returns null and no later kernel call ever reaches `get_session`/`release_sessions` to
 * clear them. The pointers are deliberately not freed — the heap they came from dies with the
 * instance.
 *
 * @private
 * @returns {void}
 */
export function retire_instance() {
    wasmAborted = true;
    wasmInstance = null;
    sessions.clear();
}

/**
 * Returns the session for `key`, allocating it if the shape changed or nothing is cached yet.
 *
 * The session is tied to the instance it was allocated from. A kernel trap retires the instance, so
 * comparing identity is what stops a pointer from a dead instance ever being reused — in that case
 * the old buffers are dropped rather than freed, since the heap they came from is gone with it.
 *
 * @private
 * @param {object} inst - The WASM instance the buffers belong to.
 * @param {any} exports - The WASM instance exports.
 * @param {string} key - Identifies the kernel and its operand shape.
 * @param {Record<string, number>} sizes - Byte size per buffer name.
 * @returns {{ ptrs: Record<string, number>, sources: Record<string, Float64Array | null> }}
 */
export function get_session(inst, exports, key, sizes) {
    // Every byte size is compared, not just the `key`. Sizing a buffer from anything other than the
    // array that will be copied into it is how this went wrong once already: t-SNE's `ystep` and
    // `gains` are N ⨯ D in the *input* dimensionality, not the N ⨯ dim the kernel indexes, so a
    // buffer sized from `n * dim` took a view sized from `ystep.length` and wrote through the next
    // TLSF block header. Reusing a buffer only when its size is unchanged makes that unrepresentable.
    const names = Object.keys(sizes);
    const existing = sessions.get(key);
    if (existing && existing.inst === inst) {
        let fits = Object.keys(existing.sizes).length === names.length;
        if (fits) {
            for (const name of names) {
                if (existing.sizes[name] !== sizes[name]) {
                    fits = false;
                    break;
                }
            }
        }
        if (fits) return existing;
    }

    if (existing) {
        // Only free when the buffers came from the instance still in use.
        if (existing.inst === inst) {
            for (const name of Object.keys(existing.ptrs)) exports.free(existing.ptrs[name]);
        }
        sessions.delete(key);
    }

    /** @type {Record<string, number>} */
    const ptrs = {};
    for (const name of names) ptrs[name] = exports.allocate(sizes[name]);
    const created = { inst, ptrs, sizes: { ...sizes }, sources: {} };
    sessions.set(key, created);
    return created;
}

/**
 * Copies `source` into the session buffer `name`, unless the same array is already there.
 *
 * Identity, not contents, is the test: the iterative callers build an operand once and hand over
 * the same `Float64Array` every iteration, and a caller that rebuilds it produces a new array. An
 * operand mutated in place behind a stable identity would be missed, so this is only used for
 * operands the callers treat as constant for the run.
 *
 * The view is built from the source's own constructor, so this works for the `Int32Array` and
 * `Float32Array` operands too — UMAP's constant edge arrays are of all three kinds.
 *
 * @private
 * @param {any} exports - The WASM instance exports.
 * @param {{ ptrs: Record<string, number>, sources: Record<string, any> }} s
 * @param {string} name
 * @param {Float64Array | Float32Array | Int32Array | Uint32Array} source
 */
export function copy_once(exports, s, name, source) {
    if (s.sources[name] === source) return;
    const View = /** @type {any} */ (source.constructor);
    new View(exports.memory.buffer, s.ptrs[name], source.length).set(source);
    s.sources[name] = source;
}

/**
 * Releases the named sessions, if they are still held.
 *
 * The optimisers call this from a `finally` when a run ends, which is what keeps the retained
 * buffers — around 96 MB for t-SNE at N = 2000 — from outliving the projection that needed them.
 * Freeing a session is never a correctness matter: the next call simply reallocates.
 *
 * @private
 * @param {string[]} keys
 * @returns {void}
 */
export function release_sessions(keys) {
    const inst = /** @type {any} */ (wasmInstance);
    for (const key of keys) {
        const held = sessions.get(key);
        if (!held) continue;
        if (inst && held.inst === inst) {
            for (const name of Object.keys(held.ptrs)) inst.exports.free(held.ptrs[name]);
        }
        sessions.delete(key);
    }
}

/**
 * How many buffer sessions are currently held.
 *
 * Internal, and only used by the tests: that the optimisers release their buffers when a run ends
 * is otherwise invisible from outside, and an unreleased session is a leak rather than a failure.
 *
 * @private
 * @returns {number}
 */
export function held_session_count() {
    return sessions.size;
}

/**
 * Releases every session.
 *
 * Internal: the optimisers release their own sessions when a run ends, so nothing outside needs
 * this. It remains for tests and for the rare caller driving `next()` by hand.
 *
 * @private
 * @returns {void}
 */
export function release_wasm_buffers() {
    const inst = /** @type {any} */ (wasmInstance);
    for (const held of sessions.values()) {
        if (inst && held.inst === inst) {
            for (const name of Object.keys(held.ptrs)) inst.exports.free(held.ptrs[name]);
        }
    }
    sessions.clear();
}
