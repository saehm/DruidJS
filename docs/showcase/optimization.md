# Interactive Optimization

Many modern DR algorithms ([t-SNE](/api/classes/TSNE), [UMAP](/api/classes/UMAP), [TriMap](/api/classes/TriMap), [Sammon](/api/classes/SAMMON), [SMACOF](/api/classes/SMACOF), [SQDMDS](/api/classes/SQDMDS)) are iterative. DruidJS provides `generator` support, allowing you to observe the optimization process as it happens. Watch how clusters form and separate over time.

<script setup>
import AnimatedProjection from '../components/AnimatedProjection.vue'
</script>

<AnimatedProjection />

## How-to (Code)

To animate the optimization process, use the `generator()` method instead of `transform()`. This returns an iterator that yields intermediate projections.

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... data ... */
];
const tsne = new druid.TSNE(data);

// Use generator for many iterations
const iterations = 500;
for (const projection of tsne.generator(iterations)) {
  // Update your UI/Plot with the current 'projection'
  updatePlot(projection);
}
```
