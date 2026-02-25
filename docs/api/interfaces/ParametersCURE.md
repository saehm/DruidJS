[@saehrimnir/druidjs](../globals.md) / ParametersCURE

# Interface: ParametersCURE

Defined in: [clustering/index.js:62](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L62)

## Properties

### K

> **K**: `number`

Defined in: [clustering/index.js:63](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L63)

Target number of clusters. Default is `2`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:66](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L66)

Distance metric function. Default is `euclidean`

***

### num\_representatives

> **num\_representatives**: `number`

Defined in: [clustering/index.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L64)

Number of representative points per cluster. Default is `5`

***

### seed

> **seed**: `number`

Defined in: [clustering/index.js:67](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L67)

Random seed. Default is `1212`

***

### shrink\_factor

> **shrink\_factor**: `number`

Defined in: [clustering/index.js:65](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L65)

Factor to shrink representatives toward centroid (0-1). Default is `0.5`
