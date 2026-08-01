# Interactive Optimization

Many modern DR algorithms ([t-SNE](/api/classes/TSNE), [UMAP](/api/classes/UMAP), [PaCMAP](/api/classes/PaCMAP), [LocalMAP](/api/classes/LocalMAP), [TriMap](/api/classes/TriMap), [Sammon](/api/classes/SAMMON), [SMACOF](/api/classes/SMACOF), [SQDMDS](/api/classes/SQDMDS)) are iterative. DruidJS provides `generator` support, allowing you to observe the optimization process as it happens. Watch how clusters form and separate over time.

<script setup>
import AnimatedProjection from '../components/AnimatedProjection.vue'
</script>

# Interactive Optimization

Many modern DR algorithms ([t-SNE](/api/classes/TSNE), [UMAP](/api/classes/UMAP), [PaCMAP](/api/classes/PaCMAP), [LocalMAP](/api/classes/LocalMAP), [TriMap](/api/classes/TriMap), [Sammon](/api/classes/SAMMON), [SMACOF](/api/classes/SMACOF), [SQDMDS](/api/classes/SQDMDS)) are iterative. DruidJS provides `generator` support, allowing you to observe the optimization process as it happens. Watch how clusters form and separate over time.

[PaCMAP](/dimred/pacmap) and [LocalMAP](/dimred/localmap) are worth watching in particular: their optimization runs in three phases, and the mid-near weight decays from 1000 to 3 across the first one. The global arrangement is settled early and barely moves afterwards, which is visibly different from t-SNE or UMAP pulling structure out gradually. The phase lengths here scale with the iteration slider so the whole schedule stays visible at any setting.

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
