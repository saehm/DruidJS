//import "../dist/druid.js";
import * as druid from "./test_index.js";
import * as assert from "assert";

describe("eig", () => {
    const N = 20;
    const R = new druid.Randomizer(1212);
    let M = new druid.Matrix(N, N, () => (R.random - .5)  * 5);
    M = M.dot(M.T);

    it("qr", () => {
        assert.ok(druid.qr(M));
        assert.ok(druid.qr_householder(M));
    });
    it("simultanious poweriteration", () => {
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr}));
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr_householder}));
        let eigs = druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, N)
        assert.equal(eigs.eigenvectors[0].length, 2)

        eigs = druid.simultaneous_poweriteration(M, 6, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, N)
        assert.equal(eigs.eigenvectors[0].length, 6)
    }).timeout(10000);
});
