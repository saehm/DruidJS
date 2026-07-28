import { describe, expect, it } from "vitest";
import { TSNE } from "../../src/dimred/TSNE.js";
import { held_session_count, isWasmAvailable, release_wasm_buffers, retire_instance } from "../../src/wasm/index.js";

/**
 * A kernel trap retires the WASM instance for the life of the process. That poisons every later
 * WASM path in the same module, which is why this lives in its own file: Vitest isolates test files,
 * so retiring the instance here does not affect the rest of the suite.
 *
 * The regression under test: the abort handler used to null the instance but leave the persistent
 * buffer sessions in place. Each session holds a reference to the retired instance, so the Map alone
 * kept its whole linear memory alive — and nothing ever cleared it, because `initWasm` returns null
 * once aborted. `retire_instance()` is the handler's own retire path, exercised directly so the test
 * does not have to corrupt the heap to reach it.
 */

/**
 * @param {number} n
 * @param {number} d
 * @returns {number[][]}
 */
function data(n, d) {
    let x = 42;
    const random = () => {
        x ^= x << 13;
        x ^= x >>> 17;
        x ^= x << 5;
        return (x >>> 0) / 4294967296;
    };
    return Array.from({ length: n }, () => Array.from({ length: d }, () => random()));
}

describe.runIf(isWasmAvailable())("WASM abort handling", () => {
    it("drops held buffer sessions when the instance is retired", () => {
        release_wasm_buffers();

        // A hand-driven generator leaves a session held (see the Symbol.dispose case in
        // session.test.js), which is exactly the state that would leak across an abort.
        const tsne = new TSNE(data(40, 5), { d: 2, seed: 42, perplexity: 10 });
        const steps = tsne.generator(50);
        steps.next();
        steps.next();
        expect(held_session_count(), "a half-driven generator holds its buffers").toBe(1);

        retire_instance();

        expect(held_session_count(), "retiring the instance must not leave the session pinning its memory").toBe(0);
        expect(isWasmAvailable(), "a retired instance stays retired").toBe(false);
    });
});
