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

    describe("comparator", () => {
        it("should order a max heap largest first", () => {
            const heap = new Heap([5, 3, 8, 1, 2], (x) => x, "max");
            expect(heap.first.value).toBe(8);
            expect(heap.toArray()).toEqual([8, 5, 3, 2, 1]);
        });

        it("should accept a custom comparator function", () => {
            // Ordering by distance from 4, which neither "min" nor "max" expresses.
            const heap = new Heap(
                [0, 3, 5, 9],
                (x) => Math.abs(x - 4),
                (a, b) => a < b,
            );
            expect(heap.first.element).toBe(3);
            expect(heap.toArray()).toEqual([3, 5, 0, 9]);
        });
    });

    describe("pop", () => {
        it("should return null when empty", () => {
            const heap = new Heap(null, (x) => x, "min");
            expect(heap.pop()).toBe(null);
        });

        it("should drain a single-element heap", () => {
            const heap = new Heap([42], (x) => x, "min");
            expect(heap.pop().element).toBe(42);
            expect(heap.empty).toBe(true);
            expect(heap.pop()).toBe(null);
        });

        it("should keep the heap ordered while draining", () => {
            const heap = new Heap([9, 4, 7, 1, 8, 2], (x) => x, "min");
            const drained = [];
            while (!heap.empty) drained.push(heap.pop().element);
            expect(drained).toEqual([1, 2, 4, 7, 8, 9]);
        });
    });

    describe("first", () => {
        it("should return null when empty", () => {
            expect(new Heap(null, (x) => x, "min").first).toBe(null);
        });

        it("should not remove the entry it returns", () => {
            const heap = new Heap([5, 3, 8], (x) => x, "min");
            expect(heap.first.element).toBe(3);
            expect(heap.length).toBe(3);
        });
    });

    describe("iterate", () => {
        it("should yield every element", () => {
            const heap = new Heap([5, 3, 8, 1], (x) => x, "min");
            expect([...heap.iterate()].sort((a, b) => a - b)).toEqual([1, 3, 5, 8]);
        });

        it("should yield nothing for an empty heap", () => {
            expect([...new Heap(null, (x) => x, "min").iterate()]).toEqual([]);
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

    it("should return null when finding an element it does not hold", () => {
        const ds = new DisjointSet([1, 2, 3]);
        expect(ds.find(99)).toBe(null);
    });

    it("should throw when uniting an element it does not hold", () => {
        const ds = new DisjointSet([1, 2, 3]);
        expect(() => ds.union(1, 99)).toThrow("x or y not found!");
        expect(() => ds.union(99, 1)).toThrow("x or y not found!");
    });

    it("should be a no-op when uniting an element with its own set", () => {
        const ds = new DisjointSet([1, 2, 3]);
        ds.union(1, 2);
        const root = ds.find(1);
        expect(ds.union(1, 2)).toBe(ds);
        expect(ds.find(1)).toBe(root);
        expect(ds.find(2)).toBe(root);
    });

    it("should attach the smaller set under the larger one", () => {
        // Union by size: merging a 3-set with a 1-set must keep the 3-set's root, whichever way
        // round the arguments come.
        const ds = new DisjointSet([1, 2, 3, 4]);
        ds.union(1, 2);
        ds.union(2, 3);
        const big_root = ds.find(1);

        ds.union(4, 1);
        expect(ds.find(4)).toBe(big_root);
    });

    it("should track children of a root", () => {
        const ds = new DisjointSet([1, 2, 3]);
        ds.union(1, 2);
        const root = ds.find(1);
        const children = ds.get_children(root);
        expect(children).toBeInstanceOf(Set);
        expect(children.size).toBeGreaterThan(0);
    });

    it("should return null for the children of an element it does not hold", () => {
        expect(new DisjointSet([1, 2]).get_children(99)).toBe(null);
    });

    it("should default to an empty set when constructed without elements", () => {
        const ds = new DisjointSet();
        expect(ds.find(1)).toBe(null);
    });
});
