[@saehrimnir/druidjs](../globals.md) / ParametersKMedoids

# Interface: ParametersKMedoids

Defined in: [clustering/index.js:28](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L28)

## Properties

### K

> **K**: `number`

Defined in: [clustering/index.js:29](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L29)

Number of clusters

***

### max\_iter

> **max\_iter**: `number` \| `null`

Defined in: [clustering/index.js:30](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L30)

Maximum number of iterations. Default is 10 * Math.log10(N). Default is `null`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:31](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L31)

Metric defining the dissimilarity. Default is `euclidean`

***

### seed

> **seed**: `number`

Defined in: [clustering/index.js:32](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L32)

Seed value for random number generator. Default is `1212`
