import { describe, expect, it } from "vitest";
import { DisjointSet, Heap } from "../../src/datastructure/index.js";

describe("Heap", () => {
    describe("constructor", () => {
        it("should create an empty heap", () => {
            const heap = new Heap(null, (x) => x, "min");
            expect(heap.length).toBe(0);
            expect(heap.empty).toBe(true);
        });

        it("should create a min heap from elements", () => {
            const heap = new Heap([5, 3, 8, 1, 2], (x) => x, "min");
            expect(heap.length).toBe(5);
            expect(heap.first.value).toBe(1);
        });

        it("should create a max heap from elements", () => {
            const heap = new Heap([5, 3, 8, 1, 2], (x) => x, "max");
            expect(heap.length).toBe(5);
            expect(heap.first.value).toBe(8);
        });
    });

    describe("push", () => {
        it("should maintain heap property after push", () => {
            const heap = new Heap([5, 3, 8], (x) => x, "min");
            heap.push(1);
            expect(heap.first.value).toBe(1);
        });
    });

    describe("pop", () => {
        it("should pop elements in sorted order from min heap", () => {
            const heap = new Heap([5, 3, 8, 1, 2, 7, 4, 6], (x) => x, "min");
            const result = [];
            while (!heap.empty) {
                result.push(heap.pop().value);
            }
            expect(result).toEqual([1, 2, 3, 4, 5, 6, 7, 8]);
        });
    });

    describe("toArray", () => {
        it("should return sorted array", () => {
            const heap = new Heap([5, 3, 8, 1, 2], (x) => x, "min");
            expect(heap.toArray()).toEqual([1, 2, 3, 5, 8]);
        });
    });
});

describe("DisjointSet", () => {
    it("should merge sets", () => {
        const ds = new DisjointSet([1, 2, 3]);
        ds.union(1, 2);
        expect(ds.find(1)).toBe(ds.find(2));
        expect(ds.find(1)).not.toBe(ds.find(3));
    });

    it("should merge multiple sets", () => {
        const ds = new DisjointSet([1, 2, 3, 4]);
        ds.union(1, 2);
        ds.union(3, 4);
        ds.union(1, 3);
        const root = ds.find(1);
        expect(ds.find(4)).toBe(root);
    });

    it("should handle path compression", () => {
        const ds = new DisjointSet([1, 2, 3, 4, 5]);
        ds.union(1, 2);
        ds.union(2, 3);
        ds.union(3, 4);
        ds.union(4, 5);
        const root = ds.find(5);
        expect(ds.find(1)).toBe(root);
    });
});
