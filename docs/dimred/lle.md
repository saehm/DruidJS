<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LLE(mistle.IRIS.values, {neighbors: 60});
    const lle = dr.transform()
</script>

# LLE

Locally Linear Embedding (LLE) seeks a lower-dimensional projection of the data which preserves distances within local neighborhoods.

## How It Works

Locally Linear Embedding (LLE) represents each data point as a linear combination of its nearest neighbors. It then finds a low-dimensional embedding that preserves these local linear relationships.

## Why or When to Use

Use LLE when the data lies on a non-linear manifold and you want to preserve local neighborhood distances (local geometry) while unfolding the manifold.

## Example

<RoundScatterplot :data=lle :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const lle = new druid.LLE(data);

// 2. Compute the projection
const projection = lle.transform();
```

## Larger Datasets

LLE builds a nearest-neighbor graph before solving for the embedding. By default it uses an exact
index — a [KD-Tree](/api/classes/KDTree) or [BallTree](/api/classes/BallTree), whichever suits the
metric — which costs O(N²) on large inputs. Pass any KNN index as `knn` to swap that out for an
approximate one:

```javascript
const knn = new druid.HNSW(data, { metric: druid.euclidean, ef: 100 });
const lle = new druid.LLE(data, { neighbors: 10, knn });
```

The eigenproblem LLE solves tolerates a few incorrect edges, so an approximate graph is usually
enough. [HNSW](/api/classes/HNSW), [Annoy](/api/classes/Annoy) and
[NNDescent](/api/classes/NNDescent) all work here.
