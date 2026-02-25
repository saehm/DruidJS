[@saehrimnir/druidjs](../globals.md) / distance_matrix

# Function: distance_matrix()

> **distance_matrix**(`A`, `metric?`): [`Matrix`](../classes/Matrix.md)

Defined in: [matrix/distance_matrix.js:22](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/matrix/distance_matrix.js#L22)

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
