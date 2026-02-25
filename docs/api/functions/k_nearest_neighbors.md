[@saehrimnir/druidjs](../globals.md) / k_nearest_neighbors

# Function: k_nearest_neighbors()

> **k_nearest_neighbors**(`A`, `k`, `metric?`): `object`[][]

Defined in: [matrix/k_nearest_neighbors.js:17](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/matrix/k_nearest_neighbors.js#L17)

Computes the k-nearest neighbors of each row of `A`.

## Parameters

### A

[`Matrix`](../classes/Matrix.md)

Either the data matrix, or a distance matrix.

### k

`number`

The number of neighbors to compute.

### metric?

Default is `euclidean`

[`Metric`](../type-aliases/Metric.md) | `"precomputed"`

## Returns

`object`[][]

The kNN graph.
