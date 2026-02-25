[@saehrimnir/druidjs](../globals.md) / ParametersMeanShift

# Interface: ParametersMeanShift

Defined in: [clustering/index.js:52](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L52)

## Properties

### bandwidth

> **bandwidth**: `number`

Defined in: [clustering/index.js:53](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L53)

bandwidth

***

### kernel

> **kernel**: `"flat"` \| `"gaussian"` \| (`dist`) => `number`

Defined in: [clustering/index.js:56](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L56)

Kernel function. Default is `gaussian`

***

### max\_iter?

> `optional` **max\_iter**: `number`

Defined in: [clustering/index.js:57](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L57)

Maximum number of iterations. Default is `Math.max(10, Math.floor(10 * Math.log10(N)))`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:54](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L54)

Metric defining the dissimilarity. Default is `euclidean`

***

### seed

> **seed**: `number`

Defined in: [clustering/index.js:55](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L55)

Seed value for random number generator. Default is `1212`

***

### tolerance?

> `optional` **tolerance**: `number`

Defined in: [clustering/index.js:58](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L58)

Convergence tolerance. Default is `1e-3`
