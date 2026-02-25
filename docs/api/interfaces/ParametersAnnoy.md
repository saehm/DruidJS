[@saehrimnir/druidjs](../globals.md) / ParametersAnnoy

# Interface: ParametersAnnoy

Defined in: [knn/index.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L17)

## Properties

### maxPointsPerLeaf

> **maxPointsPerLeaf**: `number`

Defined in: [knn/index.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L20)

Maximum points per leaf node. Default is `10`

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/index.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L18)

Metric to use: (a, b) => distance. Default is `euclidean`

***

### numTrees

> **numTrees**: `number`

Defined in: [knn/index.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L19)

Number of random projection trees to build. Default is `10`

***

### seed

> **seed**: `number`

Defined in: [knn/index.js:21](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/index.js#L21)

Seed for random number generator. Default is `1212`
