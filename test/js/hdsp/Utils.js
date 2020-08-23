"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Utils = void 0;
var Utils = (function () {
    function Utils() {
    }
    Utils.array = function (N) {
        return Utils.fillArray(N, 0);
    };
    ;
    Utils.fillArray = function (N, value) {
        var result = [];
        for (var n = 0; n < N; n++) {
            result.push(value);
        }
        return result;
    };
    ;
    Utils.randomArray = function (N) {
        var result = [];
        for (var n = 0; n < N; n++) {
            result.push(Math.random());
        }
        return result;
    };
    ;
    Utils.array2D = function (M, N) {
        var result = [];
        for (var m = 0; m < M; m++) {
            result.push(Utils.fillArray(N, 0));
        }
        return result;
    };
    Utils.random2D = function (M, N) {
        var result = [];
        for (var m = 0; m < M; m++) {
            result.push(Utils.randomArray(N));
        }
        return result;
    };
    Utils.distance = function (featureVectors) {
        var N = featureVectors.length;
        var result = Utils.array2D(N, N);
        for (var n = 0; n < N; n++) {
            for (var m = 0; m <= n; m++) {
                result[n][m] = result[m][n] = Utils.euclidean(featureVectors[n], featureVectors[m]);
            }
        }
        return result;
    };
    Utils.euclidean = function (u, v) {
        var sum = 0;
        for (var i = 0; i < u.length; i++) {
            sum += Math.pow(u[i] - v[i], 2);
        }
        return Math.sqrt(sum);
    };
    Utils.sub = function (u, v) {
        var sub = [];
        for (var i = 0; i < u.length; i++) {
            sub.push(u[i] - v[i]);
        }
        return sub;
    };
    Utils.add = function (u, v) {
        var add = [];
        for (var i = 0; i < u.length; i++) {
            add.push(u[i] + v[i]);
        }
        return add;
    };
    Utils.mult = function (u, v) {
        var mult = [];
        for (var i = 0; i < u.length; i++) {
            mult.push(u[i] * v[i]);
        }
        return mult;
    };
    Utils.norm = function (u) {
        var squaredSum = 0;
        for (var i = 0; i < u.length; i++) {
            squaredSum += Math.pow(u[i], 2);
        }
        return Math.sqrt(squaredSum);
    };
    Utils.shuffle = function (u) {
        return u.sort(function () { return Math.random() - 0.5; });
    };
    return Utils;
}());
exports.Utils = Utils;
//# sourceMappingURL=Utils.js.map