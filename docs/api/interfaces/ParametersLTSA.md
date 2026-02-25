[@saehrimnir/druidjs](../globals.md) / ParametersLTSA

# Interface: ParametersLTSA

Defined in: [dimred/index.js:68](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L68)

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:70](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L70)

the dimensionality of the projection.

---

### eig_args?

> `optional` **eig_args**: `Partial`\<[`EigenArgs`](EigenArgs.md)\>

Defined in: [dimred/index.js:73](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L73)

Parameters for the eigendecomposition algorithm.

---

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [dimred/index.js:71](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L71)

the metric which defines the distance between two points.

---

### neighbors?

> `optional` **neighbors**: `number`

Defined in: [dimred/index.js:69](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L69)

The number of neighbors for LTSA.

---

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:72](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L72)

the seed for the random number generator.
