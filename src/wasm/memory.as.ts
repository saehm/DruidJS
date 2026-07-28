// src/wasm/memory.as.ts

export function allocate(size: i32): usize {
    return heap.alloc(size);
}

export function free(ptr: usize): void {
    heap.free(ptr);
}
