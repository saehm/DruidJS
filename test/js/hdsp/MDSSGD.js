"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.MDSSGD = void 0;
var Utils_1 = require("./Utils");
var MDSSGD = (function () {
    function MDSSGD() {
    }
    MDSSGD.project = function (featureVectors, D) {
        if (featureVectors == null) {
            return null;
        }
        var N = featureVectors.length;
        if (N == 0) {
            return [];
        }
        var distances = Utils_1.Utils.distance(featureVectors);
        var result = Utils_1.Utils.random2D(N, D);
        var constraints = [];
        var w;
        var w_min = Number.MAX_VALUE;
        var w_max = Number.MIN_VALUE;
        for (var i = 0; i < distances.length; i++) {
            for (var j = 0; j < i; j++) {
                w = 1 / Math.pow(distances[i][j], 2);
                w_min = Math.min(w_min, w);
                w_max = Math.max(w_max, w);
                constraints.push([i, j, w]);
            }
        }
        var num_iter = 60;
        var epsilon = 0.1;
        var eta_max = 1 / w_min;
        var eta_min = epsilon / w_max;
        var lambda = Math.log(eta_min / eta_max) / (num_iter - 1);
        var eta = function (t) {
            return eta_max * Math.pow(lambda, t);
        };
        var schedule = [];
        for (var i = 0; i < num_iter; i++) {
            schedule.push(eta(i));
        }
        var wc;
        var pq;
        var mag;
        var r;
        var m;
        for (var _i = 0, schedule_1 = schedule; _i < schedule_1.length; _i++) {
            var c = schedule_1[_i];
            constraints = Utils_1.Utils.shuffle(constraints);
            for (var _a = 0, constraints_1 = constraints; _a < constraints_1.length; _a++) {
                var _b = constraints_1[_a], i = _b[0], j = _b[1], w_1 = _b[2];
                wc = w_1 * c;
                if (wc > 1) {
                    wc = 1;
                }
                pq = Utils_1.Utils.sub(result[i], result[j]);
                mag = Utils_1.Utils.norm(pq);
                r = (distances[i][j] - mag) / 2;
                m = Utils_1.Utils.mult(pq, Utils_1.Utils.fillArray(D, wc * r / mag));
                result[i] = Utils_1.Utils.add(result[i], m);
                result[j] = Utils_1.Utils.sub(result[j], m);
            }
        }
        return result;
    };
    return MDSSGD;
}());
exports.MDSSGD = MDSSGD;
//# sourceMappingURL=MDSSGD.js.map