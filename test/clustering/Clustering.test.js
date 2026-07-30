import { describe, expect, it } from "vitest";
import { Clustering } from "../../src/clustering/Clustering.js";
import { KMeans } from "../../src/clustering/KMeans.js";
import { Matrix } from "../../src/matrix/Matrix.js";

/**
 * The shared base every clustering method extends. Its input handling and the abstract methods it
 * defines are only reached through subclasses, which all override them — so the base's own
 * behaviour, including the errors it raises when a subclass forgets to implement something, was
 * never executed.
 */

const points = [
    [0, 0],
    [0.1, 0.1],
    [5, 5],
    [5.1, 5.1],
];

describe("Clustering (base class)", () => {
    describe("input handling", () => {
        it("should accept a nested array and record its shape", () => {
            const c = new Clustering(points);
            expect(c._N).toBe(4);
            expect(c._D).toBe(2);
            expect(c._matrix).toBeInstanceOf(Matrix);
        });

        it("should keep a Matrix as given rather than copying it", () => {
            const M = Matrix.from(points);
            expect(new Clustering(M)._matrix).toBe(M);
        });

        it("should accept Float64Array rows", () => {
            const c = new Clustering(points.map((p) => Float64Array.from(p)));
            expect(c._N).toBe(4);
            expect(c._D).toBe(2);
        });

        it("should agree with a subclass on the shape it derived", () => {
            const km = new KMeans(points, { K: 2, seed: 1212 });
            const base = new Clustering(points);
            expect([km._N, km._D]).toEqual([base._N, base._D]);
        });
    });

    describe("abstract methods", () => {
        it("should refuse get_clusters until a subclass implements it", () => {
            expect(() => new Clustering(points).get_clusters()).toThrow(
                "The function get_clusters must be implemented!",
            );
        });

        it("should refuse get_cluster_list until a subclass implements it", () => {
            expect(() => new Clustering(points).get_cluster_list()).toThrow(
                "The function get_cluster_list must be implemented!",
            );
        });

        it("should be satisfied by a subclass that implements them", () => {
            const km = new KMeans(points, { K: 2, seed: 1212 });
            expect(() => km.get_clusters()).not.toThrow();
            expect(() => km.get_cluster_list()).not.toThrow();
        });
    });
});
