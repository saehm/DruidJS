[@saehrimnir/druidjs](../globals.md) / ParametersNNDescent

# Interface: ParametersNNDescent

Defined in: [knn/index.js:61](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L61)

## Properties

### delta

> **delta**: `number`

Defined in: [knn/index.js:65](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L65)

= 0.0001 - Precision parameter. Default is `0.0001`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/index.js:62](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L62)

Called sigma in paper. Default is `euclidean`

***

### rho

> **rho**: `number`

Defined in: [knn/index.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L64)

= .8 - Sample rate. Default is `.8`

***

### samples

> **samples**: `number`

Defined in: [knn/index.js:63](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L63)

=10 - Number of samples. Default is `10`

***

### seed

> **seed**: `number`

Defined in: [knn/index.js:66](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L66)

= 1212 - Seed for the random number generator. Default is `1212`
