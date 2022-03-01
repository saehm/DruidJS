# druid.DR

This module provides dimensionality reduction methods. There are several 

## API Reference

### <code>druid._DR_</code>

The class `druid.DR` extends all of the dimensionality reduction methods implemented in Druid<sub>JS</sub> with some basic functionalities every method should have. The constructor creates a <code>_DR_</code> object.

<a href="#DR.parameter" name="DR.parameter">#</a> <code>_DR_.**parameter**(name, [value])</code>

If *value* is specified, sets the value of the specified parameter and returns the <code>_DR_</code> object. If *value* is not specified, returns the current value of the parameter.

<a href="#DR.para" name="DR.para">#</a> <code>_DR_.**para**(name, [value])</code>

Is an shorthand alias for <code>_DR_.parameter</code>.

<a href="#DR.p" name="DR.p">#</a> <code>_DR_.**para**(name, [value])</code>

Is an shorthand alias for <code>_DR_.parameter</code>.

<a href="#DR.transform" name="DR.transform">#</a> <code>_DR_.**transform**(X, [parameters])</code>

Computes a projection from dataset `X` with given parameterization `parameters` and returns the result.

The dataset `X` can be either `Number[][]` or `druid.Matrix`, which also defines the type of the returned projection. If `X` is type of `Number[][]`, then the return value will be of type `Number[][]`. If `X` is type of `druid.Matrix`, then the returned projection will be of type `druid.Matrix`.

The second argument `parameters` is an object for parameterizing the used dimensionality reduction method.

<a href="#DR.transform_async" name="DR.transform_async">#</a> <code> _DR_.**transform_async**()</code>

Works like [DR.transform](#DR.transform) but returns a Promise.

<a href="#DR.generator" name="DR.generator">#</a> <code>_DR_.**generator**(X, [parameters])</code>

Returns a generator, which yields either the intermediate step of the dimensionality reduction methods if the method is iterative, or yields the final result if it is computed in one step.

```javascript
const DR = new druid.DR(X, parameters);
const generator = DR.generator();
for (const intermediate_result of generator) {
    // do something
}
```

<a href="#DR.projection" name="DR.projection">#</a> <code>get _DR_.**projection**</code>
Returns the computed projection of the <code>_DR_</code> in the type of the given high-dimensional data.

<a href="#DR.static_transform" name="DR.static_transform">#</a> <code>static druid.DR.**transform**(X, [parameters])</code>

Is the static variant of the [`transform`](#DR.transform) function. 

In example, 
```javascript
const DR = new druid.DR(X, parameters);
const projection = DR.transform();
```
creates first a <code>_DR_</code> object which computes the projection with the `transform` call in the next line.
This two lines can be compressed to:
```javascript
const projection = druid.DR.transform(X, parameters);
```
which does the previous two lines interally.

<a href="#DR.static_transform_async" name="DR.static_transform_async">#</a> <code> static druid.DR.**transform_async**(X, [parameters])</code>

Is the static variant of the [`transform_async`](#DR.transform_async) function. 

<a href="#DR.static_generator" name="DR.static_generator">#</a> <code>static druid.DR.**generator**(X, [parameters])</code>

Is the static variant of the [`generator`](#DR.generator) function.

---

### <code>druid.PCA</code>

The dimensionality reduction method "Principal Component Analysis" is a very old technique. It projects high-dimensional data, so that the projection will have the biggest variance of the original high-dimensional datapoints. The high-dimensional data gets centered, from that a covariance matrix gets derived. The $d$ first eigenvectors of that covariance matrix are the so called *principal components*. A $d$-dimensional projection is then computed with the dot-product of the original high-dimensional data and the $d$ biggest principal components.

This technqiue is linear, which results in interpretable axis of the projection.

<a href="#DR.PCA" name="DR.PCA">#</a> <code>druid.PCA(X, [parameters])</code>

Creates an <code>_PCA_</code> object for projecting the data `X` with given `parameters`, which can be
 - `parameters.d` - The dimensionality of the projection,
 - `parameters.seed` - The seed for the random number generator,
 - `parameters.eig_args` - The parameters for the eigendecomposition algorithm.

<a href="#PCA.principal_components" name="PCA.principal_components">#</a> <code>_PCA_.**principal_components**()</code>
Computes the $d$ principal components - if not already computed - and returns them.

<a href="#PCA.transform" name="PCA.transform">#</a> <code>_PCA_.**transform**([A])</code>

If `A` is not given, returns the projection of `X` (see [DR.transform](#DR.transform)). If `A` is given, returns the projection of `A` with the $d$ principal components of `X`.

---

### <code>druid.MDS</code>

There is a whole class of "Multi Dimensional Scaling" techniques. In Druid<sub>JS</sub> is *classical*-MDS implemented. 
A MDS projection tries to maintain the distance - given by a metric - between all points.
This technique works on a distance matrix, which is a quadratic matrix, consisting of the distances between all pairs of points of `X`.

<a href="#DR.MDS" name="DR.MDS">#</a> <code>druid.MDS(X, [parameters])</code>

Creates an <code>_MDS_</code> object for projecting the data `X` with given `parameters`, which can be
 - `parameters.d` - The dimensionality of the projection,
 - `parameters.metric` - The metric which defines the distance between two points of `X`,
 - `parameters.seed` - The seed for the random number generator,
 - `parameters.eig_args` - The parameters for the eigendecomposition algorithm.

If the value of `parameters.metric` is `"precomputed"`, then `X` gets used as the distance matrix.

<a href="#MDS.transform" name="MDS.transform">#</a> <code>_MDS_.**transform**()</code>

See [<code>_DR_.transform</code>](DR.transform).

<a href="#MDS.stress" name="MDS.stress">#</a> <code>get _MDS_.stress</code>

Returns the stress of the MDS projection.







