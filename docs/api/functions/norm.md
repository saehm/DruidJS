[@saehrimnir/druidjs](../globals.md) / norm

# Function: norm()

> **norm**(`v`, `metric?`): `number`

Defined in: [matrix/norm.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/matrix/norm.js#L14)

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
