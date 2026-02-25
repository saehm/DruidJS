[@saehrimnir/druidjs](../globals.md) / ParametersLSH

# Interface: ParametersLSH

Defined in: [knn/index.js:48](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/index.js#L48)

## Properties

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/index.js:49](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/index.js#L49)

Metric to use: (a, b) => distance. Default is `euclidean`

---

### numHashFunctions

> **numHashFunctions**: `number`

Defined in: [knn/index.js:51](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/index.js#L51)

Number of hash functions per table. Default is `10`

---

### numHashTables

> **numHashTables**: `number`

Defined in: [knn/index.js:50](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/index.js#L50)

Number of hash tables. Default is `10`

---

### seed

> **seed**: `number`

Defined in: [knn/index.js:52](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/index.js#L52)

Seed for random number generator. Default is `1212`
