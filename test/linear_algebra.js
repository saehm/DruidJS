//import "../dist/druid.js";
import * as druid from "./test_index.js";
import * as assert from "assert";

describe("eig", () => {
    const R = new druid.Randomizer(1212);
    let M = new druid.Matrix(20, 20, () => (R.random - .5)  * 5);
    M = M.dot(M.T);

    it("qr", () => {
        assert.ok(druid.qr(M));
        assert.ok(druid.qr_householder(M));
        //assert.ok(druid.qr_givens(M));
    });
    it("simultanious poweriteration", () => {
        assert.ok(druid.simultaneous_poweriteration(M, 2, 100, 1212, druid.qr));
        assert.ok(druid.simultaneous_poweriteration(M, 2, 100, 1212, druid.qr_householder));
        //assert.ok(druid.simultaneous_poweriteration(M, 2, 100, 1212, druid.qr_givens));
        const eigs = druid.simultaneous_poweriteration(M, 2, 100, 1212);
        //console.log(eigs);
    });
});
