[@saehrimnir/druidjs](../globals.md) / ParametersHNSW

# Interface: ParametersHNSW

Defined in: [knn/index.js:30](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L30)

## Properties

### ef

> **ef**: `number`

Defined in: [knn/index.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L38)

Size of candidate list during search. Default is `50`

***

### ef\_construction

> **ef\_construction**: `number`

Defined in: [knn/index.js:34](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L34)

Size of candidate list during construction. Default is `200`

***

### heuristic

> **heuristic**: `boolean`

Defined in: [knn/index.js:32](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L32)

Use heuristics or naive selection. Default is `true`

***

### m

> **m**: `number`

Defined in: [knn/index.js:33](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L33)

Max number of connections per element (excluding ground layer). Default is `16`

***

### m0

> **m0**: `number` \| `null`

Defined in: [knn/index.js:35](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L35)

Max number of connections for ground layer (layer 0). Default is `2 * m`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/index.js:31](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L31)

Metric to use: (a, b) => distance. Default is `euclidean`

***

### mL

> **mL**: `number` \| `null`

Defined in: [knn/index.js:36](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L36)

Normalization factor for level generation. Default is `1 / Math.log(m)`

***

### seed

> **seed**: `number`

Defined in: [knn/index.js:37](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L37)

Seed for random number generator. Default is `1212`
