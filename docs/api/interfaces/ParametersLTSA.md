[@saehrimnir/druidjs](../globals.md) / ParametersLTSA

# Interface: ParametersLTSA

Defined in: [dimred/index.js:68](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L68)

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:70](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L70)

the dimensionality of the projection.

***

### eig\_args?

> `optional` **eig\_args**: `Partial`\<[`EigenArgs`](EigenArgs.md)\>

Defined in: [dimred/index.js:73](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L73)

Parameters for the eigendecomposition algorithm.

***

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [dimred/index.js:71](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L71)

the metric which defines the distance between two points.

***

### neighbors?

> `optional` **neighbors**: `number`

Defined in: [dimred/index.js:69](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L69)

The number of neighbors for LTSA.

***

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:72](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L72)

the seed for the random number generator.
