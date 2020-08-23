"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Scaling = void 0;
var Utils_1 = require("./Utils");
var Scaling = (function () {
    function Scaling() {
    }
    Scaling.normalize = function (coordinates) {
        var D = coordinates[0].length;
        var min = Utils_1.Utils.fillArray(D, Number.MAX_VALUE);
        var max = Utils_1.Utils.fillArray(D, Number.MIN_VALUE);
        var x;
        for (var i = 0; i < coordinates.length; i++) {
            for (var d = 0; d < D; d++) {
                x = coordinates[i][d];
                if (x < min[d]) {
                    min[d] = x;
                }
                if (x > max[d]) {
                    max[d] = x;
                }
            }
        }
        var length, longestLength = 0;
        for (var d = 0; d < D; d++) {
            length = max[d] - min[d];
            if (length > longestLength) {
                longestLength = length;
            }
        }
        for (var d = 0; d < D; d++) {
            for (var i = 0; i < coordinates.length; i++) {
                x = coordinates[i][d];
                coordinates[i][d] = (x - min[d]) / longestLength;
            }
        }
    };
    return Scaling;
}());
exports.Scaling = Scaling;
//# sourceMappingURL=Scaling.js.map