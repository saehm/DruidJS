"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Decomposition = exports.PowerMethod = void 0;
var Utils_1 = require("./Utils");
var PowerMethod = (function () {
    function PowerMethod() {
    }
    PowerMethod.singularValueDecomposition = function (U, K) {
        var MMT = PowerMethod.selfProd(U);
        var eigen = PowerMethod.eigenDecomposition(MMT, K);
        for (var k = 0; k < K; k++) {
            eigen.values[k] = Math.sqrt(eigen.values[k]);
        }
        var eVecs = PowerMethod.product(U, eigen.vectors);
        for (var k = 0; k < K; k++) {
            PowerMethod.normalize(eVecs[k]);
        }
        return new Decomposition(eigen.values, eVecs);
    };
    PowerMethod.product = function (U, V) {
        var result = Utils_1.Utils.array2D(V.length, U[0].length);
        for (var m = 0; m < result.length; m++) {
            for (var i = 0; i < U[0].length; i++) {
                result[m][i] = 0;
                for (var j = 0; j < U.length; j++) {
                    result[m][i] += U[j][i] * V[m][j];
                }
            }
        }
        return result;
    };
    PowerMethod.eigenDecomposition = function (U, K) {
        var N = U.length;
        var eVecs = PowerMethod.matrix(K, N);
        var eVals = PowerMethod.vector(K);
        var temp = Utils_1.Utils.array2D(K, N);
        var epsilon = 0.0000000001;
        var r, sumEVals = 0, sumEValsLast, fac;
        do {
            for (var k = 0; k < K; k++) {
                for (var n = 0; n < N; n++) {
                    temp[k][n] = eVecs[k][n];
                    eVecs[k][n] = 0;
                }
            }
            for (var k = 0; k < K; k++) {
                for (var i = 0; i < N; i++) {
                    for (var j = 0; j < N; j++) {
                        eVecs[k][j] += U[i][j] * temp[k][i];
                    }
                }
            }
            for (var k = 0; k < K; k++) {
                for (var p = 0; p < k; p++) {
                    fac = PowerMethod.scalar(eVecs[k], eVecs[p]) / PowerMethod.scalar(eVecs[p], eVecs[p]);
                    for (var n = 0; n < N; n++) {
                        eVecs[k][n] -= fac * eVecs[p][n];
                    }
                }
            }
            for (var k = 0; k < K; k++) {
                eVals[k] = PowerMethod.normalize(eVecs[k]);
            }
            r = 1;
            sumEValsLast = sumEVals;
            sumEVals = 0;
            for (var k = 0; k < K; k++) {
                r = Math.min(Math.abs(PowerMethod.scalar(eVecs[k], temp[k])), r);
                sumEVals += eVals[k];
            }
        } while (r < 1 - epsilon && Math.abs(sumEVals - sumEValsLast) > epsilon);
        return new Decomposition(eVals, eVecs);
    };
    PowerMethod.selfProd = function (U) {
        var N = U.length;
        var M = U[0].length;
        var result = Utils_1.Utils.array2D(N, M);
        var sum;
        for (var i = 0; i < N; i++) {
            for (var j = 0; j <= i; j++) {
                sum = 0;
                for (var k = 0; k < M; k++) {
                    sum += U[i][k] * U[j][k];
                }
                result[i][j] = result[j][i] = sum;
            }
        }
        return result;
    };
    PowerMethod.scalar = function (v, u) {
        var s = 0;
        for (var i = 0; i < v.length; i++) {
            s += v[i] * u[i];
        }
        return s;
    };
    PowerMethod.normalize = function (v) {
        var norm = Math.sqrt(PowerMethod.scalar(v, v));
        for (var i = 0; i < v.length; i++) {
            v[i] /= norm;
        }
        return norm;
    };
    PowerMethod.vector = function (length) {
        var v = Utils_1.Utils.array(length);
        for (var i = 0; i < length; i++) {
            v[i] = Math.random();
        }
        return v;
    };
    PowerMethod.matrix = function (i, j) {
        var m = Utils_1.Utils.array2D(i, j);
        for (var k = 0; k < i; k++) {
            for (var l = 0; l < j; l++) {
                m[k][l] = Math.random();
            }
        }
        return m;
    };
    return PowerMethod;
}());
exports.PowerMethod = PowerMethod;
var Decomposition = (function () {
    function Decomposition(values, vectors) {
        this.values = values;
        this.vectors = vectors;
    }
    return Decomposition;
}());
exports.Decomposition = Decomposition;
//# sourceMappingURL=PowerMethod.js.map