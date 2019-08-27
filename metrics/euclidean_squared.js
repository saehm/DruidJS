import { neumair_sum } from "../numerical/index";

export default function(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        let x = a[i];
        let y = b[i]
        s[i] = ((x - y) * (x - y))
    }
    return neumair_sum(s);
}