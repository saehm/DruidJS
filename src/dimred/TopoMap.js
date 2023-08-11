import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { DisjointSet } from "../datastructure/index.js";

/**
 * @class
 * @alias TopoMap
 * @memberof module:dimensionality_reduction
 * @extends DR
 */
export class TopoMap extends DR {
    /**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TopoMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { metric: euclidean, seed: 1212 }, parameters);
        [this._N, this._D] = this.X.shape;
        this._distance_matrix = new Matrix(this._N, this._N, 0);
        return this;
    }

    /**
     * @private
     */
    __lazy_distance_matrix(i, j, metric) {
        const D = this._distance_matrix;
        const X = this.X;
        const D_ij = D.entry(i, j);
        if (D_ij === 0) {
            let dist = metric(X.row(i), X.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
            return dist;
        }
        return D_ij;
    }

    /**
     * Computes the minimum spanning tree, using a given metric
     * @private
     * @param {function} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    _make_minimum_spanning_tree(metric = euclidean) {
        const N = this._N;
        const X = [...this.X];

        let disjoint_set = new DisjointSet(X);
        const F = [];
        let E = [];
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                E.push([i, j, this.__lazy_distance_matrix(i, j, metric)]);
            }
        }
        E = E.sort((a, b) => a[2] - b[2]);

        for (const [u, v, w] of E) {
            const set_u = disjoint_set.find(X[u]);
            const set_v = disjoint_set.find(X[v]);
            if (set_u !== set_v) {
                F.push([u, v, w]);
                disjoint_set.union(set_u, set_v);
            }
        }

        return F.sort((a, b) => a[2] - b[2]);
    }

    /**
     * initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.
     */
    init() {
        const { metric } = this._parameters;
        this.Y = new Matrix(this._N, 2, 0);
        this._Emst = this._make_minimum_spanning_tree(metric);
        this._is_initialized = true;
        return this;
    }

    /**
     * Returns true if Point C is left of line AB.
     * @private
     * @param {number[][]} PointA - Point A of line AB
     * @param {number[][]} PointB - Point B of line AB
     * @param {number[][]} PointC - Point C
     * @returns {boolean}
     */
    __hull_cross([ax, ay], [bx, by], [sx, sy]) {
        return (bx - ax) * (sy - ay) - (by - ay) * (sx - ax) <= 0;
    }

    /**
     * Computes the convex hull of the set of Points S
     * @private
     * @param {number[][]} S - Set of Points.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     * @returns {number[][]} convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     */
    __hull(S) {
        const points = S.sort(([x1, y1], [x2, y2]) => y1 - y2 || x1 - x2);
        const N = points.length;
        if (N <= 2) return points;

        const lower = [];
        for (let i = 0; i < N; ++i) {
            while (lower.length >= 2 && this.__hull_cross(lower[lower.length - 2], lower[lower.length - 1], points[i])) {
                lower.pop();
            }
            lower.push(points[i]);
        }
        const upper = [];
        for (let i = N - 1; i >= 0; --i) {
            while (upper.length >= 2 && this.__hull_cross(upper[upper.length - 2], upper[upper.length - 1], points[i])) {
                upper.pop();
            }
            upper.push(points[i]);
        }
        upper.pop();
        lower.pop();
        return lower.concat(upper);
    }

    /**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     * @private
     * @param {number[]} PointA
     * @param {number[]} PointB
     * @return {object} Object containing the sinus- and cosinus-values for a rotation.
     */
    __findAngle([p1x, p1y], [p2x, p2y]) {
        const n = euclidean([p1x, p1y], [p2x, p2y]);
        if (n === 0)
            return {
                sin: 0,
                cos: 1,
            };
        const vec = [(p2x - p1x) / n, (p2y - p1y) / n];
        const cos = vec[0];
        let sin = Math.sqrt(1 - cos * cos);
        sin = vec[1] >= 0 ? -sin : sin;
        return {
            sin: sin,
            cos: cos,
        };
    }

    /**
     * @private
     * @param {number[][]} hull
     * @param {number[]} p
     * @param {boolean} topEdge
     */
    __align_hull(hull, p, topEdge) {
        let v = -1;
        let d2;
        for (let i = 0; i < hull.length; ++i) {
            const d = euclidean(hull[i], p);
            if (v === -1) {
                d2 = d;
                v = i;
            } else {
                if (d2 > d) {
                    d2 = d;
                    v = i;
                }
            }
        }

        let v1;
        let v2;
        if (topEdge) {
            v1 = hull[v];
            v2 = hull[(v + 1) % hull.length];
        } else {
            if (v == 0) v = hull.length - 1;
            v1 = hull[v];
            v2 = hull[(v - 1) % hull.length];
        }

        const transformation = {
            tx: -hull[v][0],
            ty: -hull[v][1],
        };

        if (hull.length >= 2) {
            const { sin, cos } = this.__findAngle(v1, v2);
            transformation.sin = sin;
            transformation.cos = cos;
        } else {
            transformation.sin = 0;
            transformation.cos = 1;
        }

        return transformation;
    }

    /**
     * @private
     * @param {number[][]} Point - The point which should get transformed.
     * @param {object} Transformation - contains the values for translation and rotation.
     */
    __transform([px, py], { tx, ty, sin, cos }) {
        let x = px + tx;
        let y = py + ty;
        let xx = x * cos - y * sin;
        let yy = x * sin + y * cos;
        return [xx, yy];
    }

    /**
     * Calls {@link __transform} for each point in Set C
     * @private
     * @param {number[][]} C - Set of points.
     * @param {object} t - Transform object.
     * @param {number} yOffset - value to offset set C.
     */
    __transform_component(C, t, yOffset) {
        const N = C.length;
        for (let i = 0; i < N; ++i) {
            const c = C[i];
            const [cx, cy] = this.__transform(c, t);
            c[0] = cx;
            c[1] = cy + yOffset;
        }
    }

    /**
     * @private
     * @param {number[]} u - point u
     * @param {number[]} v - point v
     * @param {number} w - edge weight w
     */
    __align_components(u, v, w) {
        const points_u = [...u.__disjoint_set.children];
        const points_v = [...v.__disjoint_set.children];

        const hull_u = this.__hull(points_u);
        const hull_v = this.__hull(points_v);

        const t_u = this.__align_hull(hull_u, u, false);
        const t_v = this.__align_hull(hull_v, v, true);

        this.__transform_component(points_u, t_u, 0);
        this.__transform_component(points_v, t_v, w);
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality 2.
     */
    transform() {
        if (!this._is_initialized) this.init();
        const Emst = this._Emst;
        const Y = this.Y.to2dArray;
        const components = new DisjointSet(
            Y.map((y, i) => {
                y.i = i;
                return y;
            })
        );

        for (const [u, v, w] of Emst) {
            const component_u = components.find(Y[u]);
            const component_v = components.find(Y[v]);
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, w);
            components.union(component_u, component_v);
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality 2.
     * @yields {Matrix|number[][]}
     */
    *generator() {
        if (!this._is_initialized) this.init();
        const Emst = this._Emst;
        const Y = this.Y.to2dArray;
        const components = new DisjointSet(
            Y.map((y, i) => {
                y.i = i;
                return y;
            })
        );

        for (const [u, v, w] of Emst) {
            const component_u = components.find(Y[u]);
            const component_v = components.find(Y[v]);
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, w);
            components.union(component_u, component_v);
            yield this.projection;
        }
        return this.projection;
    }
}
