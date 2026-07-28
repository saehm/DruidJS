/**
 * QuickSelect module providing O(N) average selection and median partitioning.
 *
 * The pivot is drawn from a caller-supplied seeded {@link Randomizer} rather than `Math.random`.
 * A random pivot keeps the expected runtime linear on any input, including the sorted and organ-pipe
 * arrangements that defeat a fixed pivot rule, while seeding keeps every result reproducible — the
 * rest of the library derives all randomness from a `Randomizer` and nothing may depend on
 * `Math.random`.
 *
 * Partitioning is three-way (Dutch national flag), which keeps the runtime linear on inputs with
 * many duplicate keys instead of degrading to O(N²).
 *
 * @module util/quickselect
 */

/** @import { Randomizer } from "./randomizer.js" */

/**
 * Default numeric compare helper.
 * @param {any} a
 * @param {any} b
 * @returns {number}
 */
function defaultCompare(a, b) {
    return Number(a) - Number(b);
}

/**
 * In-place QuickSelect algorithm to partition an array around the k-th smallest element.
 * Runs in O(N) average time complexity (compared to O(N log N) for full Array.prototype.sort).
 *
 * After the call `arr[k]` holds the k-th smallest element, every element left of `k` compares
 * less than or equal to it, and every element right of `k` compares greater than or equal to it.
 *
 * @template T
 * @param {T[]} arr - Array to partition in-place
 * @param {Randomizer} randomizer - Seeded source of randomness for pivot selection.
 * @param {number} k - Target 0-indexed rank to select
 * @param {(a: T, b: T) => number} [compareFn] - Comparison function
 * @param {number} [start_left] - Start index (inclusive)
 * @param {number} [start_right] - End index (inclusive)
 * @returns {T} The k-th smallest element in the array
 */
export function quickselect(
    arr,
    randomizer,
    k,
    compareFn = defaultCompare,
    start_left = 0,
    start_right = arr.length - 1,
) {
    if (arr.length === 0 || k < 0 || k >= arr.length) {
        return arr[k];
    }

    let left = start_left;
    let right = start_right;

    while (left < right) {
        const pivot = arr[left + Math.floor(randomizer.random * (right - left + 1))];
        const [lt, gt] = partition(arr, left, right, pivot, compareFn);
        if (k < lt) {
            right = lt - 1;
        } else if (k > gt) {
            left = gt + 1;
        } else {
            // arr[lt..gt] all compare equal to the pivot, so arr[k] is already in place.
            return arr[k];
        }
    }
    return arr[k];
}

/**
 * QuickSelect by specific dimension axis for spatial trees (KDTree, BallTree).
 * Partitions array in-place along `arr[i].element[axis]`.
 *
 * @template {{ element: number[] | Float64Array }} T
 * @param {T[]} arr - Array to partition in-place
 * @param {Randomizer} randomizer - Seeded source of randomness for pivot selection.
 * @param {number} k - Target 0-indexed rank to select
 * @param {number} axis - Dimension coordinate index
 * @param {number} [start_left] - Start index (inclusive)
 * @param {number} [start_right] - End index (inclusive)
 * @returns {T} The k-th element along specified axis
 */
export function quickselectByAxis(arr, randomizer, k, axis, start_left = 0, start_right = arr.length - 1) {
    return quickselect(arr, randomizer, k, (a, b) => a.element[axis] - b.element[axis], start_left, start_right);
}

/**
 * In-place three-way partition (Dutch national flag) around `pivot`.
 *
 * @private
 * @template T
 * @param {T[]} arr
 * @param {number} left
 * @param {number} right
 * @param {T} pivot - The pivot element to partition around.
 * @param {(a: T, b: T) => number} compareFn
 * @returns {[number, number]} `[lt, gt]` such that `arr[lt..gt]` all compare equal to `pivot`.
 */
function partition(arr, left, right, pivot, compareFn) {
    let lt = left;
    let i = left;
    let gt = right;

    while (i <= gt) {
        const cmp = compareFn(arr[i], pivot);
        if (cmp < 0) {
            swap(arr, i, lt);
            ++lt;
            ++i;
        } else if (cmp > 0) {
            swap(arr, i, gt);
            --gt;
        } else {
            ++i;
        }
    }

    return [lt, gt];
}

/**
 * Swaps two entries of an array in-place.
 *
 * @private
 * @template T
 * @param {T[]} arr
 * @param {number} a
 * @param {number} b
 */
function swap(arr, a, b) {
    const temp = arr[a];
    arr[a] = arr[b];
    arr[b] = temp;
}
