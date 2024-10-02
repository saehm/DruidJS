import * as druid from "./test_index.js";
//import "../dist/druid.js";
import * as assert from "assert";
/* 
describe("DR", () => {
    describe("transforms", () => {
        it("PCA: static transform", () => {
            const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
            const X2 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
            const X = X1.concat(X2, "vertical");
            const Y = druid.PCA.transform(X);
            for (const y of Y.values) {
                assert.ok(!Number.isNaN(y));
            }
        });
        it("PCA: using transform with parameter", () => {
            const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
            const X2 = new druid.Matrix(8, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
            const DR = new druid.PCA(X1);
            const Y1 = DR.transform();
            const Y2 = DR.transform(X2);
            assert.deepEqual(Y1.shape, [10, 2]);
            assert.deepEqual(Y2.shape, [8, 2]);
            for (const Y of [Y1, Y2]) {
                for (const y of Y.values) {
                    assert.ok(!Number.isNaN(y));
                }
            }
            const V = DR.principal_components();
            assert.deepEqual(V.shape, [4, 2]);
        });

        it("UMAP: static transform.", () => {
            const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
            const X2 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
            const X = X1.concat(X2, "vertical");
            const Y = druid.UMAP.transform(X);
            for (const y of Y.values) {
                assert.ok(!Number.isNaN(y));
            }
        });
    });

    const X = new druid.Matrix(50, 5, () => Math.random());
    describe("async", () => {
        it("UMAP: static transform_async.", async () => {
            const Y = await druid.UMAP.transform_async(X);
            assert.deepEqual(Y.shape, [X.shape[0], 2]);
            for (const y of Y.values) {
                assert.ok(!Number.isNaN(y));
            }
        });
        it("UMAP: transform_async", async () => {
            const DR = new druid.UMAP(X);
            const Y = await DR.transform_async();
            assert.deepEqual(Y.shape, [X.shape[0], 2]);
            for (const y of Y.values) {
                assert.ok(!Number.isNaN(y));
            }
        });
    });
}); */

const methods = ["UMAP", "TSNE", "TriMap", "PCA", "LLE", "LTSA", "ISOMAP", "FASTMAP", "MDS", "LSP", "LDA", "TopoMap", "SAMMON"];
describe("DR techniques", () => {
    const R = new druid.Randomizer(1212);
    const X = new druid.Matrix(100, 10, () => R.random - 0.5);
    const L = Float64Array.from({ length: 100 }, (_, i) => (i < 50 ? 0 : 1));
    it("FASTMAP", async () => {
        assert.ok(new druid.FASTMAP(X, { d: 1, metric: druid.manhattan }));
        const dr = new druid.FASTMAP(X, { d: 5, metric: druid.canberra });
        assert.ok(dr.parameter("d", 8));
        assert.ok(dr.para("d", 3));
        assert.ok(dr.p("d", 2));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.canberra);

        let Y;
        assert.ok((Y = druid.FASTMAP.transform(X, { d: 1 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 1], Y.shape);
        for (const v of Y.values) assert.ok(!isNaN(v));
        assert.ok((Y = druid.FASTMAP.transform_async(X, { metric: druid.manhattan })));
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("MDS", async () => {
        assert.ok(new druid.MDS(X, { d: 1, metric: druid.manhattan }));
        const dr = new druid.MDS(X, { d: 5, metric: druid.canberra });
        assert.ok(dr.parameter("d", 8));
        assert.ok(dr.para("d", 3));
        assert.ok(dr.p("d", 2));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.canberra);

        let Y;
        assert.ok((Y = druid.MDS.transform(X, { d: 1 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 1], Y.shape);
        for (const v of Y.values) assert.ok(!isNaN(v));
        assert.ok((Y = druid.MDS.transform_async(X, { metric: druid.manhattan, eig_args: { max_iterations: 200, qr: druid.qr, tol: 1e-2 } })));
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("ISOMAP", async () => {
        assert.ok(new druid.ISOMAP(X, { d: 3, metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.ISOMAP(X, { d: 2, metric: druid.manhattan });
        assert.ok(dr.para("neighbors", 8));
        assert.ok(dr.p("neighbors", 12));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.manhattan);

        let Y;
        assert.ok((Y = druid.ISOMAP.transform(X, { d: 2 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.ISOMAP.transform_async(X, { metric: druid.manhattan, neighbors: 20 })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("LLE", async () => {
        assert.ok(new druid.LLE(X, { d: 3, metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.LLE(X, { d: 2, metric: druid.manhattan });
        assert.ok(dr.para("neighbors", 8));
        assert.ok(dr.p("neighbors", 12));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.manhattan);

        let Y;
        assert.ok((Y = druid.LLE.transform(X, { d: 2 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.LLE.transform_async(X, { metric: druid.manhattan, neighbors: 20 })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("LTSA", async () => {
        assert.ok(new druid.LTSA(X, { d: 3, metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.LTSA(X, { d: 2, metric: druid.manhattan });
        assert.ok(dr.para("neighbors", 8));
        assert.ok(dr.p("neighbors", 12));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.manhattan);

        let Y;
        assert.ok((Y = druid.LTSA.transform(X, { d: 2 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.LTSA.transform_async(X, { metric: druid.manhattan, neighbors: 20 })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("LSP", async () => {
        assert.ok(new druid.LSP(X, { neighbors: 15, control_points: 10, d: 3, metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.LSP(X, { d: 2, metric: druid.manhattan });
        assert.ok(dr.para("neighbors", 8));
        assert.ok(dr.p("neighbors", 12));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.manhattan);

        let Y;
        assert.ok((Y = druid.LSP.transform(X, { d: 2, seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.LSP.transform_async(X, { metric: druid.manhattan, neighbors: 20 })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("LDA", async () => {
        assert.ok(new druid.LDA(X, { d: 4, seed: 1212, labels: L }));
        assert.ok(new druid.LDA(X, { labels: L }));
        const dr = new druid.LDA(X, { labels: L });
        assert.ok(dr.para("labels", L));
        assert.ok(dr.p("labels", L));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());
        let Y;
        assert.ok((Y = druid.LDA.transform(X, { d: 3, labels: L, eig_args: { max_iterations: 200, qr: druid.qr, tol: 1e-2 } })));
        assert.deepEqual([100, 3], Y.shape);
        assert.ok(Y instanceof druid.Matrix);
        for (const v of Y.values) assert.ok(!isNaN(v));
        assert.ok((Y = druid.LDA.transform_async(X, { labels: L })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
    });

    it("PCA", async () => {
        assert.ok(new druid.PCA(X, { d: 1, seed: 12121212 }));
        const dr = new druid.PCA(X, { d: 5 });
        assert.ok(dr.parameter("d", 8));
        assert.ok(dr.para("d", 3));
        assert.ok(dr.p("d", 2));
        assert.ok(dr.transform());
        assert.ok(dr.generator());
        assert.ok(dr.transform_async());

        let Y;
        assert.ok((Y = druid.PCA.transform(X, { d: 1 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 1], Y.shape);
        for (const v of Y.values) assert.ok(!isNaN(v));
        assert.ok((Y = druid.PCA.transform_async(X, { eig_args: { max_iterations: 200, qr: druid.qr, tol: 1e-2 } })));
        assert.deepEqual([100, 2], (await Y).shape);

        let pc;
        assert.ok((pc = druid.PCA.principal_components(X, { seed: 12 })));
        for (const v of pc.values) assert.ok(!isNaN(v));
    });



    it("SAMMON", async () => {
        assert.ok(new druid.SAMMON(X, { magic: 0.2, d: 3, metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.SAMMON(X, { d: 2, metric: druid.manhattan, init_DR: "PCA" });
        assert.ok(dr.para("magic", 0.6));
        assert.ok(dr.p("magic", 0.5));
        assert.equal(dr._is_initialized, false);
        assert.ok(dr.init());
        assert.equal(dr._is_initialized, true);
        assert.ok(dr.transform());

        let generator;
        assert.ok(generator = dr.generator(20));
        for (const Y of generator) {
                for (const v of Y.values) {
                    assert.ok(!isNaN(v));
                }
        }
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.manhattan);

        let Y;
        assert.ok((Y = druid.SAMMON.transform(X, { d: 2, seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.SAMMON.transform_async(X, { metric: druid.manhattan, magic: 0.2 })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);

    });


    it("TopoMap", async () => {
        assert.ok(new druid.TopoMap(X, { metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.TopoMap(X, { metric: druid.manhattan });
        assert.ok(dr.para("metric", druid.canberra));
        assert.ok(dr.transform());

        let generator;
        assert.ok(generator = dr.generator());
        for (const Y of generator) {
                for (const v of Y.values) {
                    assert.ok(!isNaN(v));
                }
        }
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.canberra);

        let Y;
        assert.ok((Y = druid.TopoMap.transform(X, { seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.TopoMap.transform_async(X, { metric: druid.manhattan })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);

    });

    it("TSNE", async () => {
        assert.ok(new druid.TSNE(X, { metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.TSNE(X, { metric: druid.manhattan });
        assert.ok(dr.para("metric", druid.canberra));
        assert.ok(dr.transform());

        let generator;
        assert.ok(generator = dr.generator(200));
        for (const Y of generator) {
                for (const v of Y.values) {
                    assert.ok(!isNaN(v));
                }
        }
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.canberra);

        let Y;
        assert.ok((Y = druid.TSNE.transform(X, { seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.TSNE.transform_async(X, { metric: druid.manhattan })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);

    });

    it("UMAP", async () => {
        assert.ok(new druid.UMAP(X, { metric: druid.manhattan, seed: 2323 }));
        const dr = new druid.UMAP(X, { metric: druid.manhattan });
        assert.ok(dr.para("metric", druid.canberra));
        assert.ok(dr.transform());
        let generator;
        assert.ok(generator = dr.generator(200));
        for (const Y of generator) {
                for (const v of Y.values) {
                    assert.ok(!isNaN(v));
                }
        }
        assert.ok(dr.transform_async());
        assert.equal(dr.p("metric"), druid.canberra);

        let Y;
        assert.ok((Y = druid.UMAP.transform(X, { seed: 1111 })));
        assert.ok(Y instanceof druid.Matrix);
        assert.deepEqual([100, 2], Y.shape);
        for (const v of Y.values) {
            assert.ok(!isNaN(v));
        }
        assert.ok((Y = druid.UMAP.transform_async(X, { metric: druid.manhattan })));
        assert.ok(Y instanceof Promise);
        assert.deepEqual([100, 2], (await Y).shape);
        //done();
    }).timeout(10000);
});
