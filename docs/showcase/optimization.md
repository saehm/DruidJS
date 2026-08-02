# Interactive Optimization

Many modern DR algorithms ([t-SNE](/api/classes/TSNE), [UMAP](/api/classes/UMAP), [PaCMAP](/api/classes/PaCMAP), [LocalMAP](/api/classes/LocalMAP), [TriMap](/api/classes/TriMap), [Sammon](/api/classes/SAMMON), [SMACOF](/api/classes/SMACOF), [SQDMDS](/api/classes/SQDMDS)) are iterative. DruidJS provides `generator` support, allowing you to observe the optimization process as it happens. Watch how clusters form and separate over time.

<script setup>
import AnimatedProjection from '../components/AnimatedProjection.vue'
</script>

# Interactive Optimization

Many modern DR algorithms ([t-SNE](/api/classes/TSNE), [UMAP](/api/classes/UMAP), [PaCMAP](/api/classes/PaCMAP), [LocalMAP](/api/classes/LocalMAP), [TriMap](/api/classes/TriMap), [Sammon](/api/classes/SAMMON), [SMACOF](/api/classes/SMACOF), [StressMDS](/api/classes/StressMDS), [SQDMDS](/api/classes/SQDMDS)) are iterative. DruidJS provides `generator` support, allowing you to observe the optimization process as it happens. Watch how clusters form and separate over time.

[PaCMAP](/dimred/pacmap) and [LocalMAP](/dimred/localmap) are worth watching in particular: their optimization runs in three phases, and the mid-near weight decays from 1000 to 3 across the first one. The global arrangement is settled early and barely moves afterwards, which is visibly different from t-SNE or UMAP pulling structure out gradually. The phase lengths here scale with the iteration slider so the whole schedule stays visible at any setting.

[StressMDS](/dimred/stressmds) is the opposite kind of thing to watch: it starts from a classical MDS solution rather than from noise, so the first frame is already a reasonable layout and the run only refines it. Sliding the **weight exponent** from `0` to `-2` changes which objective it is refining — `0` is the stress [SMACOF](/dimred/smacof) minimizes, `-1` is Sammon's, `-2` is Kamada-Kawai. It also stops as soon as the improvement falls below its tolerance, so a short run is convergence rather than a truncated animation.

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
