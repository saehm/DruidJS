[@saehrimnir/druidjs](../globals.md) / norm

# Function: norm()

> **norm**(`v`, `metric?`): `number`

Defined in: [matrix/norm.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/norm.js#L14)

Computes the norm of a vector, by computing its distance to **0**.

## Parameters

### v

Vector.

`number`[] | `Float64Array`\<`ArrayBufferLike`\> | [`Matrix`](../classes/Matrix.md)

### metric?

[`Metric`](../type-aliases/Metric.md) = `euclidean`

Which metric should be used to compute the norm. Default is `euclidean`

## Returns

`number`

- The norm of `v`.
