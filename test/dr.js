import * as druid from "./test_index.js";
import * as assert from "assert";

describe("DR", () => {
    it("PCA", () => {
        const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
        const X2 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
        const X = X1.concat(X2, "vertical");
        const Y = druid.PCA.transform(X);
        for (const y of Y.values) {
            assert.ok(!Number.isNaN(y));
        }
    });
    it("PCA 2", () => {
        const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
        const X2 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
        const DR = new druid.PCA(X1);
        const Y1 = DR.transform();
        const Y2 = DR.transform(X2);
        assert.deepEqual(Y1.shape, [10, 2]);
        assert.deepEqual(Y2.shape, [10, 2]);
        for (const Y of [Y1, Y2]) {
            for (const y of Y.values) {
                assert.ok(!Number.isNaN(y));
            }
        }
        const V = DR.principal_components();
        assert.deepEqual(V.shape, [4, 2]);
    });
    it("UMAP", () => {
        const X1 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? 3 : 0));
        const X2 = new druid.Matrix(10, 4, (i, j) => Math.random() + (j == 1 ? -3 : 0));
        const X = X1.concat(X2, "vertical");
        const Y = druid.UMAP.transform(X);
        for (const y of Y.values) {
            assert.ok(!Number.isNaN(y));
        }
    });
});
