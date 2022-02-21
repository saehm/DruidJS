//import * as druid from "./test_index.js";
import "../dist/druid.js"
import * as assert from "assert";

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
});
