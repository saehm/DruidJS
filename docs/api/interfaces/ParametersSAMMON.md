[@saehrimnir/druidjs](../globals.md) / ParametersSAMMON

# Interface: ParametersSAMMON\<K\>

Defined in: [dimred/index.js:93](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L93)

## Type Parameters

### K

`K` *extends* keyof `ChooseDR`

## Properties

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:94](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L94)

the dimensionality of the projection.

***

### init\_DR?

> `optional` **init\_DR**: `K`

Defined in: [dimred/index.js:96](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L96)

Either "PCA" or "MDS", with which SAMMON initialiates the projection.

***

### init\_parameters?

> `optional` **init\_parameters**: `ChooseDR`\[`K`\]

Defined in: [dimred/index.js:97](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L97)

Parameters for the "init"-DR method.

***

### magic?

> `optional` **magic**: `number`

Defined in: [dimred/index.js:98](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L98)

learning rate for gradient descent.

***

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md) \| `"precomputed"`

Defined in: [dimred/index.js:95](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L95)

the metric which defines the distance between two points.

***

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:99](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/index.js#L99)

the seed for the random number generator.
