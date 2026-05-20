# Dimensionality Reduction

**Feature Extraction** creates entirely new features (dimensions) by combining the original ones mathematically. It doesn't just pick the "best" variables; it compresses the information from _all_ variables into a more compact, efficient form.

Think of it as **data synthesis**. Instead of throwing away data points, you are melting them down and recasting them into a stronger, leaner structure.

## 1. Linear Feature Extraction

These methods work best when the data follows simple, straight-line patterns. They project the data onto lower-dimensional planes.

- **Principal Component Analysis (PCA):** As we discussed, this focuses on finding the directions of maximum variance. It is unsupervised (it doesn't care about labels like "cat" or "dog") and is the "go-to" for general data compression.
- **Linear Discriminant Analysis (LDA):** unlike PCA, this is **supervised**. It tries to reduce dimensions in a way that maximizes the separation between different classes.
- _Analogy:_ If PCA tries to spread everyone out in a room to see them clearly, LDA tries to group all the "Team Red" people in one corner and "Team Blue" in the other.

## 2. Non-Linear Feature Extraction (Manifold Learning)

Real-world data is rarely a straight line; it's often twisted and curved (like a Swiss roll). Linear methods fail here because they squash the curves flat, destroying the structure. Non-linear methods "unroll" or "unfold" the data.

- **t-SNE (t-Distributed Stochastic Neighbor Embedding):** This is excellent for **visualization**. It focuses on keeping similar data points close together in the new lower-dimensional space. It's famous for taking complex high-dimensional data (like pixels of handwritten digits) and clustering them perfectly in 2D.
- **UMAP (Uniform Manifold Approximation and Projection):** Similar to t-SNE but generally faster and better at preserving the "global structure" (the relationship between distant clusters), not just local neighbors.
- **PaCMAP / LocalMAP:** Newer alternatives that use explicit pair types (nearest-neighbor, mid-near, and further pairs) with a three-phase optimization schedule, often producing cleaner cluster separation than UMAP with less hyperparameter sensitivity.

## Why Focus on Extraction?

- **Uncovering Latent Features:** Sometimes the "real" driver of your data isn't in the columns you measured, but in a combination of them.
- _Example:_ You measure a house's _Length_ and _Width_. Feature extraction might combine them to create a new feature representing _Area_, which is far more predictive of price than length or width alone.

- **De-noising:** By reconstructing the data in fewer dimensions, extraction methods often leave out the random noise found in the original high-dimensional space.
- **Solving Multicollinearity:** If two features are highly correlated (e.g., "salary in USD" and "salary in Euros"), models can get confused. Feature extraction merges these correlated features into a single component, removing the redundancy.
