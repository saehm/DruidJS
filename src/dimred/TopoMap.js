import { DisjointSet } from "../datastructure/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersTopoMap} from "./index.js" */

/**
 * TopoMap
 *
 * A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
 * It aims to preserve the topological structure of the data by maintaining
 * the connectivity of a minimum spanning tree.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersTopoMap>
 * @category Dimensionality Reduction
 */
export class TopoMap extends DR {
    /**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTopoMap>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X, parameters) {
        super(X, { metric: euclidean, seed: 1212 }, parameters);
        [this._N, this._D] = this.X.shape;
        this._distance_matrix = new Matrix(this._N, this._N, -1);
    }

    /**
     * @private
     * @param {number} i
     * @param {number} j
     * @param {import("../metrics/index.js").Metric} metric
     * @returns {number}
     */
    __lazy_distance_matrix(i, j, metric) {
        const D = this._distance_matrix;
        const X = this.X;
        const D_ij = D.entry(i, j);
        if (D_ij === -1 && i !== j) {
            const dist = metric(X.row(i), X.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
            return dist;
        }
        return i === j ? 0 : D_ij;
    }

    /**
     * Computes the minimum spanning tree, using a given metric
     *
     * @private
     * @param {import("../metrics/index.js").Metric} metric
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */
    _make_minimum_spanning_tree(metric = euclidean) {
        const N = this._N;
        const X = [...this.X];

        this._disjoint_set = new DisjointSet(X);
        const disjoint_set = this._disjoint_set;
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
            if (!set_u || !set_v) throw new Error("Should not happen!");
            if (set_u !== set_v) {
                F.push([u, v, w]);
                disjoint_set.union(set_u, set_v);
            }
        }

        return F.sort((a, b) => a[2] - b[2]);
    }

    /** Initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree. */
    init() {
        const { metric } = this._parameters;
        this.Y = new Matrix(this._N, 2, 0);
        this._Emst = this._make_minimum_spanning_tree(metric);
        this._is_initialized = true;
        return this;
    }

    /**
     * Returns true if Point C is left of line AB.
     *
     * @private
     * @param {Float64Array} PointA - Point A of line AB
     * @param {Float64Array} PointB - Point B of line AB
     * @param {Float64Array} PointC - Point C
     * @returns {boolean}
     */
    __hull_cross([ax, ay], [bx, by], [sx, sy]) {
        return (bx - ax) * (sy - ay) - (by - ay) * (sx - ax) <= 0;
    }

    /**
     * Computes the convex hull of the set of Points S
     *
     * @private
     * @param {Float64Array[]} S - Set of Points.
     * @returns {Float64Array[]} Convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     */
    __hull(S) {
        const points = S.sort(([x1, y1], [x2, y2]) => y1 - y2 || x1 - x2);
        const N = points.length;
        if (N <= 2) return points;

        const lower = [];
        for (let i = 0; i < N; ++i) {
            while (
                lower.length >= 2 &&
                this.__hull_cross(lower[lower.length - 2], lower[lower.length - 1], points[i])
            ) {
                lower.pop();
            }
            lower.push(points[i]);
        }
        const upper = [];
        for (let i = N - 1; i >= 0; --i) {
            while (
                upper.length >= 2 &&
                this.__hull_cross(upper[upper.length - 2], upper[upper.length - 1], points[i])
            ) {
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
     *
     * @private
     * @param {Float64Array} PointA
     * @param {Float64Array} PointB
     * @returns {{ sin: number; cos: number }} Object containing the sinus- and cosinus-values for a rotation.
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
     * @param {Float64Array[]} hull
     * @param {Float64Array} p
     * @param {boolean} topEdge
     * @returns {{ sin: number; cos: number; tx: number; ty: number }}
     */
    __align_hull(hull, p, topEdge) {
        let v = -1;
        /** @type {number} */
        let d2 = -Infinity;
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

        const v1 = hull[v];
        let v2;
        if (topEdge) {
            v2 = hull[(v + 1) % hull.length];
        } else {
            v2 = hull[(v - 1 + hull.length) % hull.length];
        }

        /** @type {{ sin?: number; cos?: number; tx: number; ty: number }} */
        const transformation = {
            tx: -v1[0],
            ty: -v1[1],
        };

        if (hull.length >= 2) {
            const { sin, cos } = this.__findAngle(v1, v2);
            transformation.sin = sin;
            transformation.cos = cos;
        } else {
            transformation.sin = 0;
            transformation.cos = 1;
        }

        return /** @type {{ sin: number; cos: number; tx: number; ty: number }} */ (transformation);
    }

    /**
     * @private
     * @param {Float64Array} Point - The point which should get transformed.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} Transformation - Contains the values for
     *   translation and rotation.
     */
    __transform([px, py], { tx, ty, sin, cos }) {
        const x = px + tx;
        const y = py + ty;
        const xx = x * cos - y * sin;
        const yy = x * sin + y * cos;
        return [xx, yy];
    }

    /**
     * Calls `__transform` for each point in Set C
     *
     * @private
     * @param {Float64Array[]} C - Set of points.
     * @param {{ sin: number; cos: number; tx: number; ty: number }} t - Transform object.
     * @param {number} yOffset - Value to offset set C.
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
     * @param {Float64Array} root_u - Root of component u
     * @param {Float64Array} root_v - Root of component v
     * @param {Float64Array} p_u - Point u
     * @param {Float64Array} p_v - Point v
     * @param {number} w - Edge weight w
     * @param {DisjointSet<Float64Array>} components - The disjoint set containing the components
     */
    __align_components(root_u, root_v, p_u, p_v, w, components) {
        if (!components) throw new Error("components not provided!");
        const u_children = components.get_children(root_u);
        const v_children = components.get_children(root_v);
        if (!u_children || !v_children) throw new Error("should not happen!");

        const points_u = [...u_children];
        const points_v = [...v_children];

        const hull_u = this.__hull(points_u);
        const hull_v = this.__hull(points_v);

        const t_u = this.__align_hull(hull_u, p_u, false);
        const t_v = this.__align_hull(hull_v, p_v, true);

        this.__transform_component(points_u, t_u, 0);
        this.__transform_component(points_v, t_v, w);
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {T}
     */
    transform() {
        if (!this._is_initialized) this.init();
        if (!this._Emst) throw new Error("Call init() first!");
        const Emst = this._Emst;
        const Y = this.Y.to2dArray();
        /** @type {DisjointSet<Float64Array>} */
        const components = new DisjointSet(
            Y,
            // Y.map((y, i) => {
            //     y.i = i;
            //     return y;
            // }),
        );

        for (const [u, v, w] of Emst) {
            const p_u = Y[u];
            const p_v = Y[v];
            const component_u = components.find(p_u);
            const component_v = components.find(p_v);
            if (!component_u || !component_v) throw new Error("Should not happen!");
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, p_u, p_v, w, components);
            components.union(component_u, component_v);
        }
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimensionality 2.
     *
     * @returns {Generator<T, T, void>}
     */
    *generator() {
        if (!this._is_initialized) this.init();
        if (!this._Emst) throw new Error("call init() first!");
        const Emst = this._Emst;
        const Y = this.Y.to2dArray();
        const components = new DisjointSet(
            Y,
            // Y.map((y, i) => {
            //     y.i = i;
            //     return y;
            // }),
        );

        for (const [u, v, w] of Emst) {
            const p_u = Y[u];
            const p_v = Y[v];
            const component_u = components.find(p_u);
            const component_v = components.find(p_v);
            if (!component_u || !component_v) throw new Error("should not happen!");
            if (component_u === component_v) continue;
            this.__align_components(component_u, component_v, p_u, p_v, w, components);
            components.union(component_u, component_v);
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new TopoMap(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new TopoMap(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTopoMap>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new TopoMap(X, parameters);
        return dr.transform_async();
    }
}
