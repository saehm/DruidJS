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
});
