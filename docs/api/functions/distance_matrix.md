[@saehrimnir/druidjs](../globals.md) / distance\_matrix

# Function: distance\_matrix()

> **distance\_matrix**(`A`, `metric?`): [`Matrix`](../classes/Matrix.md)

Defined in: [matrix/distance\_matrix.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/distance_matrix.js#L22)

Computes the distance matrix of datamatrix `A`.

## Parameters

### A

Matrix.

[`Matrix`](../classes/Matrix.md) | `Float64Array`\<`ArrayBufferLike`\>[] | `number`[][]

### metric?

[`Metric`](../type-aliases/Metric.md) = `euclidean`

The diistance metric. Default is `euclidean`

## Returns

[`Matrix`](../classes/Matrix.md)

The distance matrix of `A`.
