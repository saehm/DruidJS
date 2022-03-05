import * as druid from "./test_index.js";
//import "../dist/druid.js";
import * as assert from "assert";

describe("SQMDS", () => {
    const R = new druid.Randomizer(1212);
    const X = new druid.Matrix(100, 10, () => R.random - 0.5);

    it("DR", async () => {
        assert.ok(new druid.SQDMDS(X, { metric: druid.euclidean, seed: 2323 }));
        const dr = new druid.SQDMDS(X, { metric: druid.euclidean });
        assert.ok(dr.para("metric", druid.euclidean));
        let Y;
        assert.ok((Y = dr.transform(100)));
        let generator;
        assert.ok((generator = dr.generator(100)));
        for (const Y of generator) {
            for (const v of Y.values) {
                assert.ok(!isNaN(v));
            }
        }
        assert.ok(dr.transform_async(1));
        assert.equal(dr.p("metric"), druid.euclidean);

        assert.ok((Y = druid.SQDMDS.transform(X, { seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.SQDMDS.transform_async(X, { metric: druid.manhattan })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("helpers", async () => {
        const dr = new druid.SQDMDS(X);

        const sub_div = dr.__sub_div(3);
        const add = dr.__add(3);
        const mult = dr.__mult(3);
        const minus = dr.__minus(3);

        let a = [1, 2, 3];
        let b = [3, 4, 5];
        // sub div is not inline
        assert.deepEqual([...sub_div(b, a, 2)], [1, 1, 1]);
        assert.deepEqual([...sub_div(a, b, 2)], [-1, -1, -1]);

        // add is inline! adds up into the first summand!
        // oh my.. this side effects :D
        assert.deepEqual([...add(a, b)], [4, 6, 8]);
        a = [1, 2, 3];
        assert.deepEqual([...add(a, b, a.slice(), b)], [8, 12, 16]);
        a = [1, 2, 3];
        assert.deepEqual([...add(a, b)], [4, 6, 8]);
        a = [1, 2, 3];

        // mult is inline
        assert.deepEqual([...mult(a, 3)], [3, 6, 9]);
        a = [1, 2, 3];
        assert.deepEqual([...mult(a, 0)], [0, 0, 0]);
        a = [1, 2, 3];

        // minus is inline!
        assert.deepEqual([...minus(a, a)], [0, 0, 0]);
        a = [1, 2, 3];
        assert.deepEqual([...minus(a, b)], [-2, -2, -2]);
        a = [1, 2, 3];
    });

    it("metric helper", () => {
        const dr = new druid.SQDMDS(X, { metric: druid.euclidean });
        dr.init();
        let metric = dr._HD_metric;
        let metric_ex = dr._HD_metric_exaggeration;
        assert.deepEqual(metric(1, 2, X), druid.euclidean(X.row(1), X.row(2)));
        assert.deepEqual(metric_ex(1, 2, X), druid.euclidean_squared(X.row(1), X.row(2)));

        dr.parameter("metric", druid.canberra);
        dr.init();
        metric = dr._HD_metric;
        metric_ex = dr._HD_metric_exaggeration;
        assert.deepEqual(metric(1, 2, X), druid.canberra(X.row(1), X.row(2)));
        assert.deepEqual(metric_ex(1, 2, X), druid.canberra(X.row(1), X.row(2)) ** 2);
    });
});
