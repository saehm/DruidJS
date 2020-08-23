"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Procrustes = void 0;
var Procrustes = (function () {
    function Procrustes() {
    }
    Procrustes.toPoint2D = function (A, mirror) {
        var points = [];
        A.forEach(function (a, index) {
            if (mirror) {
                points.push(new Point2D(a[0], -a[1]));
            }
            else {
                points.push(new Point2D(a[0], a[1]));
            }
        });
        return points;
    };
    Procrustes.toArray = function (A) {
        var points = [];
        A.forEach(function (a, index) {
            points.push([a.x, a.y]);
        });
        return points;
    };
    Procrustes.mean = function (A) {
        var mean = new Point2D(0, 0);
        A.forEach(function (a, index) {
            mean.x += a.x;
            mean.y += a.y;
        });
        mean.x /= A.length;
        mean.y /= A.length;
        return mean;
    };
    Procrustes.center = function (A) {
        var mean = Procrustes.mean(A);
        var scale = 0;
        A.forEach(function (a, index) {
            a.x -= mean.x;
            a.y -= mean.y;
            scale += a.x * a.x;
            scale += a.y * a.y;
        });
        return Math.sqrt(scale / A.length);
    };
    Procrustes.scale = function (A, scale) {
        A.forEach(function (a, index) {
            a.x /= scale;
            a.y /= scale;
        });
    };
    Procrustes.rotate = function (A, B, theta) {
        if (theta == null) {
            var n = 0, d = 0;
            for (var i = 0; i < A.length; i++) {
                n += A[i].x * B[i].y - A[i].y * B[i].x;
                d += A[i].x * B[i].x + A[i].y * B[i].y;
            }
            theta = Math.atan(n / d);
        }
        var x = 0, y = 0;
        for (var i = 0; i < A.length; i++) {
            x = Math.cos(theta) * A[i].x - Math.sin(theta) * A[i].y;
            y = Math.sin(theta) * A[i].x + Math.cos(theta) * A[i].y;
            A[i].x = x;
            A[i].y = y;
        }
        return theta;
    };
    Procrustes.transform = function (A, B) {
        var A_p = Procrustes.toPoint2D(A, false);
        var A_p_m = Procrustes.toPoint2D(A, true);
        var B_p = Procrustes.toPoint2D(B, false);
        Procrustes.scale(A_p, Procrustes.center(A_p));
        Procrustes.scale(B_p, Procrustes.center(B_p));
        Procrustes.scale(A_p_m, Procrustes.center(A_p_m));
        var A_p_r = Procrustes.copy(A_p);
        var theta = Procrustes.rotate(A_p, B_p);
        Procrustes.rotate(A_p_r, B_p, theta - Math.PI);
        var A_p_m_r = Procrustes.copy(A_p_m);
        theta = Procrustes.rotate(A_p_m, B_p);
        Procrustes.rotate(A_p_m_r, B_p, Math.PI + theta);
        var d = Procrustes.distance(A_p, B_p);
        var d_r = Procrustes.distance(A_p_r, B_p);
        var d_m = Procrustes.distance(A_p_m, B_p);
        var d_m_r = Procrustes.distance(A_p_m_r, B_p);
        var min = Math.min(d, Math.min(d_r, Math.min(d_m, d_m_r)));
        if (d == min) {
            return Procrustes.toArray(A_p);
        }
        else if (d_r == min) {
            return Procrustes.toArray(A_p_r);
        }
        else if (d_m == min) {
            return Procrustes.toArray(A_p_m);
        }
        else if (d_m_r == min) {
            return Procrustes.toArray(A_p_m_r);
        }
        else {
            return null;
        }
    };
    Procrustes.copy = function (A) {
        var copy = [];
        for (var _i = 0, A_1 = A; _i < A_1.length; _i++) {
            var a = A_1[_i];
            copy.push(new Point2D(a.x, a.y));
        }
        return copy;
    };
    Procrustes.distance = function (A, B) {
        var k = A.length;
        var e = 0;
        for (var i = 0; i < k; i++) {
            e += Math.sqrt(Math.pow(A[i].x - B[i].x, 2) + Math.pow(A[i].y - B[i].y, 2));
        }
        return e;
    };
    return Procrustes;
}());
exports.Procrustes = Procrustes;
var Point2D = (function () {
    function Point2D(x, y) {
        this.x = x;
        this.y = y;
    }
    return Point2D;
}());
//# sourceMappingURL=Procrustes.js.map