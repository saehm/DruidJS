[@saehrimnir/druidjs](../globals.md) / ParametersXMeans

# Interface: ParametersXMeans

Defined in: [clustering/index.js:42](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L42)

## Properties

### K_max

> **K_max**: `number`

Defined in: [clustering/index.js:44](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L44)

Maximum number of clusters. Default is `10`

---

### K_min

> **K_min**: `number`

Defined in: [clustering/index.js:43](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L43)

Minimum number of clusters. Default is `2`

---

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:45](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L45)

Distance metric function. Default is `euclidean`

---

### min_cluster_size

> **min_cluster_size**: `number`

Defined in: [clustering/index.js:47](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L47)

Minimum points required to consider splitting a cluster. Default is `25`

---

### seed

> **seed**: `number`

Defined in: [clustering/index.js:46](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L46)

Random seed. Default is `1212`

---

### tolerance

> **tolerance**: `number`

Defined in: [clustering/index.js:48](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L48)

Convergence tolerance for KMeans. Default is `0.001`
