# Automatic Cluster Discovery

One of the biggest challenges in unsupervised learning is choosing the number of clusters ($K$).

While standard [K-Means](/api/classes/KMeans) requires you to specify $K$ upfront, [**X-Means**](/api/classes/XMeans) (available in DruidJS) automatically discovers the optimal number of clusters using the **Bayesian Information Criterion (BIC)**.

<script setup>
import XMeansShowcase from '../components/XMeansShowcase.vue'
</script>

<XMeansShowcase />

## How it works

1. **Initial K**: X-Means starts with a minimum number of clusters (typically $K=2$).
2. **Recursive Splitting**: For each cluster, the algorithm attempts to split it into two sub-clusters using K-Means.
3. **Statistical Validation**: It compares the BIC score of the original cluster versus the split. If the split provides a significantly better fit, the new clusters are kept.
4. **Iterative Refinement**: This process repeats until it reaches a maximum $K$ or no further improvements are found.

## Benefits

- **Zero Configuration**: No need to run "elbow plots" or manual searches for $K$.
- **Objective Fit**: Uses a robust statistical measure (BIC) to prevent over-fitting.
- **Speed**: Built on top of K-Means, it remains efficient for large datasets.

## How-to (Code)

Use [`XMeans`](/api/classes/XMeans) when you want the algorithm to determine the cluster count for you. You only need to provide a search range (`K_min` and `K_max`).

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... high-dimensional data ... */
];

// Initialize XMeans with a search range for clusters
const xmeans = new druid.XMeans(data, {
  K_min: 2,
  K_max: 10,
});

// Run the discovery algorithm
// It returns an array of cluster indices for each point
const clusters = xmeans.get_cluster_list();

// You can check how many clusters were found
console.log(`Optimal number of clusters: ${xmeans.k}`);
```
