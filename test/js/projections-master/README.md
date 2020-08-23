# High-Dimensional Space Projections
TypeScript library for sparse and dense high-dimensional space projections.

## Usage

### TypeDocs
https://ibm.github.io/projections/

### Examples

#### Multi-Dimensional Scaling (PivotMDS)
```typescript
import * as hdsp from "hdsp";
hdsp.PivotMDS.project(
    [
        [1,0,0,1,0,0,0,1],  // feature vector 1
        [0,1,1,0,1,1,1,1],  // feature vector 2
                            // ...
        [0,0,0,1,0,1,0,1]   // feature vector n
    ],
    4, // number of pivots
    2  // output dimensionality
);
```

#### Multi-dimensional Scaling by Stochastic Gradient Descent (MDSSGD)
```typescript
import * as hdsp from "hdsp";
hdsp.MDSSGD.project(
    [
        [1,0,0,1,0,0,0,1],  // feature vector 1
        [0,1,1,0,1,1,1,1],  // feature vector 2
                            // ...
        [0,0,0,1,0,1,0,1]   // feature vector n
    ],
    2  // output dimensionality
);
```

#### Procrustes Analysis (2D)
```typescript
import * as hdsp from "hdsp";
hdsp.Procrustes.transform(
    [[0,0],[0,-2],[-1,-1]],  // input
    [[0,0],[0,2],[1,1]]      // reference
);
```

## Developers
### Code repository
https://github.com/IBM/projections

### Build / Release Instructions
1. `tsc`
2. `typedoc --out docs`
3. `touch docs/.nojekyll`
4. `npm version patch`
5. `npm publish`
6. `git push`
