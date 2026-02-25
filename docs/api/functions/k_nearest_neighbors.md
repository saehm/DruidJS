[@saehrimnir/druidjs](../globals.md) / k\_nearest\_neighbors

# Function: k\_nearest\_neighbors()

> **k\_nearest\_neighbors**(`A`, `k`, `metric?`): `object`[][]

Defined in: [matrix/k\_nearest\_neighbors.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/k_nearest_neighbors.js#L17)

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
