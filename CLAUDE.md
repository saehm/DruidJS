# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Build:**
```sh
npm run build        # TypeScript compilation (tsc) + Rollup bundling
npm run dev          # Watch mode for development
```

**Test:**
```sh
npm test                  # Run full test suite (Node + Browser via Playwright)
npm run test:coverage     # Run tests with v8 coverage reporting
npm run test:types        # TypeScript type checking
npm run test:types:strict # Strict TypeScript type checking
```

To run a single test file:
```sh
npx vitest run test/dimred/PCA.test.js
```

**Code Quality:**
```sh
npm run lint    # Biomejs lint on ./src
npm run format  # Biomejs format ./src and ./test
npm run check   # Biomejs check with auto-write on ./src
```

**Docs:**
```sh
npm run docs:dev    # VitePress dev server
npm run docs:build  # Build TypeDoc API docs + VitePress static site
```

## Architecture

DruidJS is a dimensionality reduction library. Source lives in `src/` and distributes as both ESM (`dist/druid.js`) and CJS (`dist/druid.cjs`). The main entry point is `src/index.js`, which re-exports all submodules.

**Module breakdown:**
- `src/dimred/` — 17 DR algorithms (PCA, MDS, TSNE, UMAP, ISOMAP, LDA, LLE, etc.) plus a base DR class
- `src/clustering/` — 8 clustering algorithms (KMeans, KMedoids, OPTICS, CURE, MeanShift, etc.)
- `src/knn/` — 7 nearest-neighbor implementations (BallTree, KDTree, HNSW, NNDescent, Annoy, LSH, NaiveKNN)
- `src/matrix/` — Core Matrix class with linear algebra operations (transpose, dot, inverse, QR, eigen)
- `src/linear_algebra/` — QR decomposition, eigendecomposition via simultaneous power iteration
- `src/metrics/` — 17 distance metrics (Euclidean, Cosine, Jaccard, Wasserstein, etc.)
- `src/datastructure/` — Heap, DisjointSet
- `src/numerical/` — Numerically stable summation (Kahan, Neumair)
- `src/optimization/` — Powell optimizer
- `src/util/` — Randomizer, min/max utilities

Each subdirectory exports its public API through a local `index.js`.

**Type system:** The project uses plain JavaScript with JSDoc annotations. TypeScript (`tsc`) reads those annotations to emit declaration files into `dist/types/`, which Rollup then bundles into `dist/druid.d.ts` (ESM) and `dist/druid.d.cts` (CJS). Type tests live in `test/types/`.

**Testing:** Vitest runs two project environments simultaneously — standard Node and multi-browser via Playwright (Chromium, Firefox, WebKit). Test files mirror the `src/` directory layout under `test/`.

**Package manager:** pnpm (v10+).
