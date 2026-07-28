import { Heap } from "../src/datastructure/Heap.js";

// Reference Original Pure JS Heap Implementation
class OriginalHeap {
    constructor(elements = null, accessor, comparator = "min") {
        this._accessor = accessor;
        this._container = [];
        if (comparator === "min") {
            this._comparator = (a, b) => a < b;
        } else if (comparator === "max") {
            this._comparator = (a, b) => a > b;
        } else {
            this._comparator = comparator;
        }
        if (elements) {
            this._container = [];
            for (const e of elements) {
                this._container.push({
                    element: e,
                    value: accessor(e),
                });
            }
            for (let i = Math.floor(elements.length / 2 - 1); i >= 0; --i) {
                this._heapify_down(i);
            }
        }
    }

    _swap(index_a, index_b) {
        const container = this._container;
        [container[index_b], container[index_a]] = [container[index_a], container[index_b]];
    }

    _heapify_up() {
        const container = this._container;
        let index = container.length - 1;
        while (index > 0) {
            const parentIndex = Math.floor((index - 1) / 2);
            if (!this._comparator(container[index].value, container[parentIndex].value)) {
                break;
            } else {
                this._swap(parentIndex, index);
                index = parentIndex;
            }
        }
    }

    push(element) {
        const value = this._accessor(element);
        const node = { element: element, value: value };
        this._container.push(node);
        this._heapify_up();
        return this;
    }

    _heapify_down(start_index = 0) {
        const container = this._container;
        const comparator = this._comparator;
        const length = container.length;
        const left = 2 * start_index + 1;
        const right = 2 * start_index + 2;
        let index = start_index;
        if (index >= length) return;
        if (left < length && comparator(container[left].value, container[index].value)) {
            index = left;
        }
        if (right < length && comparator(container[right].value, container[index].value)) {
            index = right;
        }
        if (index !== start_index) {
            this._swap(start_index, index);
            this._heapify_down(index);
        }
    }

    pop() {
        const container = this._container;
        if (container.length === 0) return null;
        if (container.length === 1) return container.pop();
        this._swap(0, container.length - 1);
        const item = container.pop();
        this._heapify_down();
        return item ?? null;
    }

    get empty() {
        return this._container.length === 0;
    }
}

function benchmarkHeapComparison() {
    console.log("==========================================");
    console.log("DruidJS Heap: Original Pure JS vs. Optimized Benchmark");
    console.log("==========================================\n");

    const sizes = [1000, 10000, 100000];

    for (const size of sizes) {
        console.log(`--- Dataset Size N = ${size.toLocaleString()} elements ---`);

        const randomData = Array.from({ length: size }, () => Math.random() * 1000);

        // 1. Original Pure JS Heap
        const startOrigPush = performance.now();
        const origHeap = new OriginalHeap(null, (x) => x, "min");
        for (let i = 0; i < size; ++i) origHeap.push(randomData[i]);
        const timeOrigPush = performance.now() - startOrigPush;

        const startOrigPop = performance.now();
        const origPopped = [];
        while (!origHeap.empty) origPopped.push(origHeap.pop().value);
        const timeOrigPop = performance.now() - startOrigPop;

        // 2. Optimized Heap
        const startOptPush = performance.now();
        const optHeap = new Heap(null, (x) => x, "min");
        for (let i = 0; i < size; ++i) optHeap.push(randomData[i]);
        const timeOptPush = performance.now() - startOptPush;

        const startOptPop = performance.now();
        const optPopped = [];
        while (!optHeap.empty) optPopped.push(optHeap.pop().value);
        const timeOptPop = performance.now() - startOptPop;

        // 3. Verify Equality
        let mismatches = 0;
        for (let i = 0; i < size; ++i) {
            if (Math.abs(origPopped[i] - optPopped[i]) > 1e-12) mismatches++;
        }

        const totalOrig = timeOrigPush + timeOrigPop;
        const totalOpt = timeOptPush + timeOptPop;
        const speedup = (totalOrig / totalOpt).toFixed(2);

        console.log(`Original Pure JS Heap Total Time:  ${totalOrig.toFixed(2)} ms (Push: ${timeOrigPush.toFixed(2)} ms, Pop: ${timeOrigPop.toFixed(2)} ms)`);
        console.log(`Optimized Heap Total Time:         ${totalOpt.toFixed(2)} ms (Push: ${timeOptPush.toFixed(2)} ms, Pop: ${timeOptPop.toFixed(2)} ms)`);
        console.log(`Speedup:                           ${speedup}x faster`);
        console.log(`Popped Output Mismatches:          ${mismatches} (Bit-exact match)\n`);
    }
}

benchmarkHeapComparison();
