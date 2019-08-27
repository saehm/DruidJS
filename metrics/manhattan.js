export default function(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i])
    }
    return sum
}