[@saehrimnir/druidjs](../globals.md) / ParametersCURE

# Interface: ParametersCURE

Defined in: [clustering/index.js:62](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L62)

## Properties

### K

> **K**: `number`

Defined in: [clustering/index.js:63](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L63)

Target number of clusters. Default is `2`

---

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:66](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L66)

Distance metric function. Default is `euclidean`

---

### num_representatives

> **num_representatives**: `number`

Defined in: [clustering/index.js:64](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L64)

Number of representative points per cluster. Default is `5`

---

### seed

> **seed**: `number`

Defined in: [clustering/index.js:67](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L67)

Random seed. Default is `1212`

---

### shrink_factor

> **shrink_factor**: `number`

Defined in: [clustering/index.js:65](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/index.js#L65)

Factor to shrink representatives toward centroid (0-1). Default is `0.5`
