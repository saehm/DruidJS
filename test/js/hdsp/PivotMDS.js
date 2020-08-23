"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.PivotMDS = void 0;
var PowerMethod_1 = require("./PowerMethod");
var Utils_1 = require("./Utils");
var PivotMDS = (function () {
    function PivotMDS() {
    }
    PivotMDS.project = function (featureVectors, K, D) {
        if (featureVectors == null) {
            return null;
        }
        var N = featureVectors.length;
        if (N == 0) {
            return [];
        }
        K = Math.min(N, K);
        var distances = Utils_1.Utils.distance(featureVectors);
        var result = Utils_1.Utils.array2D(N, D);
        var C = Utils_1.Utils.array2D(K, N);
        for (var k = 0; k < K; k++) {
            for (var n = 0; n < N; n++) {
                C[k][n] = distances[k][n] * distances[k][n];
            }
        }
        var rMeans = Utils_1.Utils.array(K);
        var cMeans = Utils_1.Utils.array(N);
        var mean = 0;
        for (var k = 0; k < K; k++) {
            for (var n = 0; n < N; n++) {
                rMeans[k] += C[k][n];
                cMeans[n] += C[k][n];
                mean += C[k][n];
            }
            rMeans[k] /= N;
        }
        for (var n = 0; n < N; n++) {
            cMeans[n] /= K;
        }
        mean /= K * N;
        for (var k = 0; k < K; k++) {
            for (var n = 0; n < N; n++) {
                C[k][n] = -.5 * (C[k][n] - rMeans[k] - cMeans[n] + mean);
            }
        }
        var decomposition = PowerMethod_1.PowerMethod.singularValueDecomposition(C, D);
        var eVals = decomposition.values;
        var eVecs = decomposition.vectors;
        for (var i = 0; i < eVecs.length; i++) {
            var scale = Math.sqrt(eVals[i]);
            if (isNaN(scale)) {
                scale = 0;
            }
            for (var j = 0; j < eVecs[0].length; j++) {
                if (isNaN(eVecs[i][j])) {
                    eVecs[i][j] = 0;
                }
                eVecs[i][j] *= scale;
            }
        }
        for (var n = 0; n < N; n++) {
            for (var d = 0; d < D; d++) {
                result[n][d] = eVecs[d][n];
            }
        }
        return result;
    };
    return PivotMDS;
}());
exports.PivotMDS = PivotMDS;
//# sourceMappingURL=PivotMDS.js.map