[@saehrimnir/druidjs](../globals.md) / ParametersKMedoids

# Interface: ParametersKMedoids

Defined in: [clustering/index.js:28](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L28)

## Properties

### K

> **K**: `number`

Defined in: [clustering/index.js:29](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L29)

Number of clusters

---

### max_iter

> **max_iter**: `number` \| `null`

Defined in: [clustering/index.js:30](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L30)

Maximum number of iterations. Default is 10 \* Math.log10(N). Default is `null`

---

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:31](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L31)

Metric defining the dissimilarity. Default is `euclidean`

---

### seed

> **seed**: `number`

Defined in: [clustering/index.js:32](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L32)

Seed value for random number generator. Default is `1212`
