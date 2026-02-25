[@saehrimnir/druidjs](../globals.md) / ParametersISOMAP

# Interface: ParametersISOMAP

Defined in: [dimred/index.js:41](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L41)

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:43](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L43)

the dimensionality of the projection.

---

### eig_args?

> `optional` **eig_args**: `Partial`\<[`EigenArgs`](EigenArgs.md)\>

Defined in: [dimred/index.js:47](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L47)

Parameters for the eigendecomposition algorithm.

---

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [dimred/index.js:44](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L44)

the metric which defines the distance between two points.

---

### neighbors?

> `optional` **neighbors**: `number`

Defined in: [dimred/index.js:42](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L42)

The number of neighbors ISOMAP should use to project the data.

---

### project?

> `optional` **project**: `"MDS"` \| `"SMACOF"`

Defined in: [dimred/index.js:45](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L45)

Whether to use classical MDS or SMACOF for the final DR.

---

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:46](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L46)

the seed for the random number generator.
