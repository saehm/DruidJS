[@saehrimnir/druidjs](../globals.md) / ParametersLLE

# Interface: ParametersLLE

Defined in: [dimred/index.js:59](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L59)

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:61](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L61)

the dimensionality of the projection.

---

### eig_args?

> `optional` **eig_args**: `Partial`\<[`EigenArgs`](EigenArgs.md)\>

Defined in: [dimred/index.js:64](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L64)

Parameters for the eigendecomposition algorithm.

---

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [dimred/index.js:62](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L62)

the metric which defines the distance between two points.

---

### neighbors?

> `optional` **neighbors**: `number`

Defined in: [dimred/index.js:60](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L60)

The number of neighbors for LLE.

---

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:63](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L63)

the seed for the random number generator.
