[@saehrimnir/druidjs](../globals.md) / ParametersMDS

# Interface: ParametersMDS

Defined in: [dimred/index.js:77](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L77)

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:78](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L78)

the dimensionality of the projection.

***

### eig\_args?

> `optional` **eig\_args**: `Partial`\<[`EigenArgs`](EigenArgs.md)\>

Defined in: [dimred/index.js:81](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L81)

Parameters for the eigendecomposition algorithm.

***

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md) \| `"precomputed"`

Defined in: [dimred/index.js:79](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L79)

the metric which defines the distance between two points.

***

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L80)

the seed for the random number generator.
