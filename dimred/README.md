
<a href="#module_dimensionality_reduction" name="module_dimensionality_reduction">#</a> <code>**dimensionality_reduction**</code>



* [dimensionality_reduction](#module_dimensionality_reduction)

    * [TopoMap](#TopoMap)

        * [new exports.TopoMap(X, parameters)](#new_TopoMap_new)

        * [.projection](#DR+projection)

        * [.init()](#TopoMap+init)

        * [.transform()](#TopoMap+transform)

        * [.parameter(name, [value])](#DR+parameter)

        * [.generator()](#DR+generator)

        * [.check_init()](#DR+check_init)

        * [.transform_async(...args)](#DR+transform_async)

    * _global_
        * [DR](#DR)

            * [new exports.DR(X, parameters)](#new_DR_new)

            * _instance_
                * [.projection](#DR+projection)

                * [.parameter(name, [value])](#DR+parameter)

                * [.transform()](#DR+transform)

                * [.generator()](#DR+generator)

                * [.check_init()](#DR+check_init)

                * [.transform_async(...args)](#DR+transform_async)

                * [.para(name, [value])](#DR+para)

                * [.p(name, [value])](#DR+p)

                * [.para(name, [value])](#DR+para)

                * [.p(name, [value])](#DR+p)

            * _static_
                * [.transform(...args)](#DR.transform)

                * [.transform_async(...args)](#DR.transform_async)

                * [.generator(...args)](#DR.generator)

        * [FASTMAP](#FASTMAP)

            * [new exports.FASTMAP(X, parameters)](#new_FASTMAP_new)

            * [.projection](#DR+projection)

            * [.transform()](#FASTMAP+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [HIPP](#HIPP)

            * [new exports.HIPP(X, [d], [metric])](#new_HIPP_new)

            * [.projection](#HIPP+projection)

            * [.transform()](#HIPP+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [ISOMAP](#ISOMAP)

            * [new exports.ISOMAP(X, parameters)](#new_ISOMAP_new)

            * [.projection](#DR+projection)

            * [.transform()](#ISOMAP+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [LDA](#LDA)

            * [new exports.LDA(X, parameters)](#new_LDA_new)

            * [.projection](#DR+projection)

            * [.transform()](#LDA+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [LLE](#LLE)

            * [new exports.LLE(X, parameters, neighbors, [d], [metric], [seed])](#new_LLE_new)

            * [.projection](#DR+projection)

            * [.transform()](#LLE+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [LSP](#LSP)

            * [new exports.LSP(X, parameters)](#new_LSP_new)

            * [.projection](#DR+projection)

            * [.init(DR, DR_parameters)](#LSP+init)

            * [.transform()](#LSP+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [LTSA](#LTSA)

            * [new exports.LTSA(X, parameters)](#new_LTSA_new)

            * [.projection](#DR+projection)

            * [.transform()](#LTSA+transform)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [MDS](#MDS)

            * [new exports.MDS(X, parameters)](#new_MDS_new)

            * [.projection](#DR+projection)

            * [.transform()](#MDS+transform)

            * [.stress()](#MDS+stress)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [PCA](#PCA)

            * [new exports.PCA(X, parameters)](#new_PCA_new)

            * [.projection](#DR+projection)

            * [.transform([A])](#PCA+transform)

            * [.principal_components()](#PCA+principal_components)

            * [.parameter(name, [value])](#DR+parameter)

            * [.generator()](#DR+generator)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [SAMMON](#SAMMON)

            * [new exports.SAMMON(X, parameters)](#new_SAMMON_new)

            * [.projection](#DR+projection)

            * [.transform([max_iter])](#SAMMON+transform)

            * [.generator([max_iter])](#SAMMON+generator)

            * [.parameter(name, [value])](#DR+parameter)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [TriMap](#TriMap)

            * [new exports.TriMap(X, parameters)](#new_TriMap_new)

            * [.projection](#DR+projection)

            * [.init([pca], [knn])](#TriMap+init)

            * [._generate_triplets(n_inliers, n_outliers, n_random)](#TriMap+_generate_triplets)

            * [._grad(Y)](#TriMap+_grad)

            * [.transform(max_iteration)](#TriMap+transform)

            * [.generator(max_iteration)](#TriMap+generator)

            * [.parameter(name, [value])](#DR+parameter)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [TSNE](#TSNE)

            * [new exports.TSNE(X, parameters)](#new_TSNE_new)

            * [.projection](#DR+projection)

            * [.init(distance_matrix)](#TSNE+init)

            * [.transform([iterations])](#TSNE+transform)

            * [.generator([iterations])](#TSNE+generator)

            * [.parameter(name, [value])](#DR+parameter)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)

        * [UMAP](#UMAP)

            * [new exports.UMAP(X, parameters)](#new_UMAP_new)

            * [.projection](#DR+projection)

            * [.init()](#UMAP+init)

            * [.transform([iterations])](#UMAP+transform)

            * [.generator([iterations])](#UMAP+generator)

            * [.parameter(name, [value])](#DR+parameter)

            * [.check_init()](#DR+check_init)

            * [.transform_async(...args)](#DR+transform_async)



<a href="#TopoMap" name="TopoMap">#</a> <code>*dimensionality_reduction***TopoMap**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://arxiv.org/pdf/2009.01512.pdf](https://arxiv.org/pdf/2009.01512.pdf)  

* [TopoMap](#TopoMap)

    * [new exports.TopoMap(X, parameters)](#new_TopoMap_new)

    * [.projection](#DR+projection)

    * [.init()](#TopoMap+init)

    * [.transform()](#TopoMap+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_TopoMap_new" name="new_TopoMap_new">#</a> new <code>**exports.TopoMap**</code>
(X, parameters)

TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*topoMap*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#TopoMap+init" name="TopoMap+init">#</a> <code>*topoMap*.**init**</code>
()

initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.


<a href="#TopoMap+transform" name="TopoMap+transform">#</a> <code>*topoMap*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimensionality 2.


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*topoMap*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*topoMap*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*topoMap*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*topoMap*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#DR" name="DR">#</a> <code>*dimensionality_reduction***DR**</code>



* [DR](#DR)

    * [new exports.DR(X, parameters)](#new_DR_new)

    * _instance_
        * [.projection](#DR+projection)

        * [.parameter(name, [value])](#DR+parameter)

        * [.transform()](#DR+transform)

        * [.generator()](#DR+generator)

        * [.check_init()](#DR+check_init)

        * [.transform_async(...args)](#DR+transform_async)

        * [.para(name, [value])](#DR+para)

        * [.p(name, [value])](#DR+p)

        * [.para(name, [value])](#DR+para)

        * [.p(name, [value])](#DR+p)

    * _static_
        * [.transform(...args)](#DR.transform)

        * [.transform_async(...args)](#DR.transform_async)

        * [.generator(...args)](#DR.generator)



<a href="#new_DR_new" name="new_DR_new">#</a> new <code>**exports.DR**</code>
(X, parameters)

Takes the default parameters and seals them, remembers the type of input [X](X), and initializes the random number generator.


- X <code>Matrix</code> | <code>Array.&lt;Array.&lt;Number&gt;&gt;</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed value for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*dR*.**projection**</code>


**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*dR*.**parameter**</code>
(name, [value])

Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+transform" name="DR+transform">#</a> <code>*dR*.**transform**</code>
()

Computes the projection.

**Returns**: <code>Matrix</code> - - Returns the projection.  

<a href="#DR+generator" name="DR+generator">#</a> <code>*dR*.**generator**</code>
()

Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*dR*.**check_init**</code>
()

If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*dR*.**transform_async**</code>
(...args)


- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#DR+para" name="DR+para">#</a> <code>*dR*.**para**</code>
(name, [value])

Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+p" name="DR+p">#</a> <code>*dR*.**p**</code>
(name, [value])

Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+para" name="DR+para">#</a> <code>*dR*.**para**</code>
(name, [value])

Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+p" name="DR+p">#</a> <code>*dR*.**p**</code>
(name, [value])

Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR.transform" name="DR.transform">#</a> <code>*DR*.**transform**</code>
(...args)


- ...args <code>any</code> - Takes the same arguments of the constructor of the respective DR method.

**Returns**: <code>Matrix</code> \| <code>Array</code> - - The dimensionality reduced dataset.  

<a href="#DR.transform_async" name="DR.transform_async">#</a> <code>*DR*.**transform_async**</code>
(...args)


- ...args <code>any</code> - Takes the same arguments of the constructor of the respective DR method.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#DR.generator" name="DR.generator">#</a> <code>*DR*.**generator**</code>
(...args)


- ...args <code>any</code> - Takes the same arguments of the constructor of the respective DR method.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#FASTMAP" name="FASTMAP">#</a> <code>*dimensionality_reduction***FASTMAP**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://doi.org/10.1145/223784.223812](https://doi.org/10.1145/223784.223812)  

* [FASTMAP](#FASTMAP)

    * [new exports.FASTMAP(X, parameters)](#new_FASTMAP_new)

    * [.projection](#DR+projection)

    * [.transform()](#FASTMAP+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_FASTMAP_new" name="new_FASTMAP_new">#</a> new <code>**exports.FASTMAP**</code>
(X, parameters)

FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the dimensionality of the projection.


<a href="#DR+projection" name="DR+projection">#</a> <code>*fastmaP*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#FASTMAP+transform" name="FASTMAP+transform">#</a> <code>*fastmaP*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Computes the projection.

**Returns**: <code>Matrix</code> - The [d](d)-dimensional projection of the data matrix [X](X).  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*fastmaP*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*fastmaP*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*fastmaP*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*fastmaP*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#HIPP" name="HIPP">#</a> <code>*dimensionality_reduction***HIPP**</code>


**Extends**: [<code>DR</code>](#DR)  

* [HIPP](#HIPP)

    * [new exports.HIPP(X, [d], [metric])](#new_HIPP_new)

    * [.projection](#HIPP+projection)

    * [.transform()](#HIPP+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_HIPP_new" name="new_HIPP_new">#</a> new <code>**exports.HIPP**</code>
(X, [d], [metric])


- X <code>Matrix</code> - the high-dimensional data.
- [d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
- [metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.


<a href="#HIPP+projection" name="HIPP+projection">#</a> <code>*hipP*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> - Returns the projection.  

<a href="#HIPP+transform" name="HIPP+transform">#</a> <code>*hipP*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Computes the projection.

**Returns**: <code>Matrix</code> - Returns the projection.  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*hipP*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*hipP*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*hipP*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*hipP*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#ISOMAP" name="ISOMAP">#</a> <code>*dimensionality_reduction***ISOMAP**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://doi.org/10.1126/science.290.5500.2319](https://doi.org/10.1126/science.290.5500.2319)  

* [ISOMAP](#ISOMAP)

    * [new exports.ISOMAP(X, parameters)](#new_ISOMAP_new)

    * [.projection](#DR+projection)

    * [.transform()](#ISOMAP+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_ISOMAP_new" name="new_ISOMAP_new">#</a> new <code>**exports.ISOMAP**</code>
(X, parameters)

Isometric feature mapping (ISOMAP).


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - .neighbors <code>Number</code> - the number of neighbors [ISOMAP](#ISOMAP) should use to project the data.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*isomaP*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#ISOMAP+transform" name="ISOMAP+transform">#</a> <code>*isomaP*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Computes the projection.

**Returns**: <code>Matrix</code> - Returns the projection.  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*isomaP*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*isomaP*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*isomaP*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*isomaP*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#LDA" name="LDA">#</a> <code>*dimensionality_reduction***LDA**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x](https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x)  

* [LDA](#LDA)

    * [new exports.LDA(X, parameters)](#new_LDA_new)

    * [.projection](#DR+projection)

    * [.transform()](#LDA+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_LDA_new" name="new_LDA_new">#</a> new <code>**exports.LDA**</code>
(X, parameters)

Linear Discriminant Analysis.


- X <code>Matrix</code> - The high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - .labels <code>Array</code> - The labels / classes for each data point.
    - [.d] <code>number</code> <code> = 2</code> - The dimensionality of the projection.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*ldA*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#LDA+transform" name="LDA+transform">#</a> <code>*ldA*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimenionality [d](d).


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*ldA*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*ldA*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*ldA*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*ldA*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#LLE" name="LLE">#</a> <code>*dimensionality_reduction***LLE**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://doi.org/10.1126/science.290.5500.2323](https://doi.org/10.1126/science.290.5500.2323)  

* [LLE](#LLE)

    * [new exports.LLE(X, parameters, neighbors, [d], [metric], [seed])](#new_LLE_new)

    * [.projection](#DR+projection)

    * [.transform()](#LLE+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_LLE_new" name="new_LLE_new">#</a> new <code>**exports.LLE**</code>
(X, parameters, neighbors, [d], [metric], [seed])

Locally Linear Embedding.


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
- neighbors <code>Number</code> - the label / class of each data point.
- [d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
- [metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
- [seed] <code>Number</code> <code> = 1212</code> - the dimensionality of the projection.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*llE*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#LLE+transform" name="LLE+transform">#</a> <code>*llE*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimenionality [d](d).


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*llE*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*llE*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*llE*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*llE*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#LSP" name="LSP">#</a> <code>*dimensionality_reduction***LSP**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://ieeexplore.ieee.org/document/4378370](https://ieeexplore.ieee.org/document/4378370)  
**Todo**

- [ ] accept precomputed distance matrix.


* [LSP](#LSP)

    * [new exports.LSP(X, parameters)](#new_LSP_new)

    * [.projection](#DR+projection)

    * [.init(DR, DR_parameters)](#LSP+init)

    * [.transform()](#LSP+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_LSP_new" name="new_LSP_new">#</a> new <code>**exports.LSP**</code>
(X, parameters)

Least Squares Projection.


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.neighbors] <code>Number</code> <code> = Math.max(Math.floor(N / 10), 2)</code> - number of neighbors to consider.
    - [.control_points] <code>Number</code> <code> = Math.ceil(Math.sqrt(N))</code> - number of controlpoints
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*lsP*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#LSP+init" name="LSP+init">#</a> <code>*lsP*.**init**</code>
(DR, DR_parameters)


- DR [<code>DR</code>](#DR) - method used for position control points.
- DR_parameters <code>Object</code> - Object containing parameters for the DR method which projects the control points


<a href="#LSP+transform" name="LSP+transform">#</a> <code>*lsP*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Computes the projection.

**Returns**: <code>Matrix</code> - Returns the projection.  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*lsP*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*lsP*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*lsP*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*lsP*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#LTSA" name="LTSA">#</a> <code>*dimensionality_reduction***LTSA**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://epubs.siam.org/doi/abs/10.1137/S1064827502419154](https://epubs.siam.org/doi/abs/10.1137/S1064827502419154)  

* [LTSA](#LTSA)

    * [new exports.LTSA(X, parameters)](#new_LTSA_new)

    * [.projection](#DR+projection)

    * [.transform()](#LTSA+transform)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_LTSA_new" name="new_LTSA_new">#</a> new <code>**exports.LTSA**</code>
(X, parameters)

Local Tangent Space Alignment


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - .neighbors <code>Number</code> - the number of neighbors [LTSA](#LTSA) should use to project the data.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*ltsA*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#LTSA+transform" name="LTSA+transform">#</a> <code>*ltsA*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimenionality [d](d).


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*ltsA*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*ltsA*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*ltsA*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*ltsA*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#MDS" name="MDS">#</a> <code>*dimensionality_reduction***MDS**</code>


**Extends**: [<code>DR</code>](#DR)  

* [MDS](#MDS)

    * [new exports.MDS(X, parameters)](#new_MDS_new)

    * [.projection](#DR+projection)

    * [.transform()](#MDS+transform)

    * [.stress()](#MDS+stress)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_MDS_new" name="new_MDS_new">#</a> new <code>**exports.MDS**</code>
(X, parameters)

Classical MDS.


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*mdS*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#MDS+transform" name="MDS+transform">#</a> <code>*mdS*.**transform**</code>
()

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimensionality [d](d).


<a href="#MDS+stress" name="MDS+stress">#</a> <code>*mdS*.**stress**</code>
()

**Returns**: <code>Number</code> - - the stress of the projection.  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*mdS*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*mdS*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*mdS*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*mdS*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#PCA" name="PCA">#</a> <code>*dimensionality_reduction***PCA**</code>


**Extends**: [<code>DR</code>](#DR)  

* [PCA](#PCA)

    * [new exports.PCA(X, parameters)](#new_PCA_new)

    * [.projection](#DR+projection)

    * [.transform([A])](#PCA+transform)

    * [.principal_components()](#PCA+principal_components)

    * [.parameter(name, [value])](#DR+parameter)

    * [.generator()](#DR+generator)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_PCA_new" name="new_PCA_new">#</a> new <code>**exports.PCA**</code>
(X, parameters)


- X <code>Matrix</code> | <code>Array.&lt;Array.&lt;Number&gt;&gt;</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.
    - [.eig_args] <code>Number</code> - Parameters for the eigendecomposition algorithm.


<a href="#DR+projection" name="DR+projection">#</a> <code>*pcA*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#PCA+transform" name="PCA+transform">#</a> <code>*pcA*.**transform**</code>
([A])

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimensionality [d](d). If parameter [A](A) is given, then project [A](A) with the principal components of [X](X).


- [A] <code>null</code> | <code>Matrix</code> | <code>Array</code> <code> = </code> - If given, the data to project.

**Returns**: <code>Matrix</code> \| <code>Array</code> - - The projected data.  

<a href="#PCA+principal_components" name="PCA+principal_components">#</a> <code>*pcA*.**principal_components**</code>
()

Computes the [d](d) principal components of Matrix [X](X).


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*pcA*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+generator" name="DR+generator">#</a> <code>*pcA*.**generator**</code>
()

**Overrides**: [<code>generator</code>](#DR+generator)  
Computes the projection.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the dimensionality reduction method.  

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*pcA*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*pcA*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#SAMMON" name="SAMMON">#</a> <code>*dimensionality_reduction***SAMMON**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**: [https://arxiv.org/pdf/2009.01512.pdf](https://arxiv.org/pdf/2009.01512.pdf)  

* [SAMMON](#SAMMON)

    * [new exports.SAMMON(X, parameters)](#new_SAMMON_new)

    * [.projection](#DR+projection)

    * [.transform([max_iter])](#SAMMON+transform)

    * [.generator([max_iter])](#SAMMON+generator)

    * [.parameter(name, [value])](#DR+parameter)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_SAMMON_new" name="new_SAMMON_new">#</a> new <code>**exports.SAMMON**</code>
(X, parameters)

SAMMON's Mapping


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.init] <code>&quot;PCA&quot;</code> | <code>&quot;MDS&quot;</code> | <code>&quot;random&quot;</code> <code> = &quot;random&quot;</code> - Either "PCA" or "MDS", with which SAMMON initialiates the projection. With "random" a random matrix gets used as starting point.
    - [.init_parameters] <code>Object</code> - Parameters for the [init](init)-DR method.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*sammoN*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#SAMMON+transform" name="SAMMON+transform">#</a> <code>*sammoN*.**transform**</code>
([max_iter])

**Overrides**: [<code>transform</code>](#DR+transform)  
Transforms the inputdata [X](X) to dimenionality 2.


- [max_iter] <code>Number</code> <code> = 200</code> - Maximum number of iteration steps.

**Returns**: <code>Matrix</code> \| <code>Array</code> - - The projection of [X](X).  

<a href="#SAMMON+generator" name="SAMMON+generator">#</a> <code>*sammoN*.**generator**</code>
([max_iter])

**Overrides**: [<code>generator</code>](#DR+generator)  
Transforms the inputdata [X](X) to dimenionality 2.


- [max_iter] <code>Number</code> <code> = 200</code> - Maximum number of iteration steps.

**Returns**: <code>Generator</code> - - A generator yielding the intermediate steps of the projection of [X](X).  

<a href="#DR+parameter" name="DR+parameter">#</a> <code>*sammoN*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*sammoN*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*sammoN*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#TriMap" name="TriMap">#</a> <code>*dimensionality_reduction***TriMap**</code>


**Extends**: [<code>DR</code>](#DR)  
**See**

- [https://arxiv.org/pdf/1910.00204v1.pdf](https://arxiv.org/pdf/1910.00204v1.pdf)
- [https://github.com/eamid/trimap](https://github.com/eamid/trimap)


* [TriMap](#TriMap)

    * [new exports.TriMap(X, parameters)](#new_TriMap_new)

    * [.projection](#DR+projection)

    * [.init([pca], [knn])](#TriMap+init)

    * [._generate_triplets(n_inliers, n_outliers, n_random)](#TriMap+_generate_triplets)

    * [._grad(Y)](#TriMap+_grad)

    * [.transform(max_iteration)](#TriMap+transform)

    * [.generator(max_iteration)](#TriMap+generator)

    * [.parameter(name, [value])](#DR+parameter)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_TriMap_new" name="new_TriMap_new">#</a> new <code>**exports.TriMap**</code>
(X, parameters)


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.weight_adj] <code>Number</code> <code> = 500</code> - scaling factor.
    - [.c] <code>Number</code> <code> = 5</code> - number of triplets multiplier.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.tol] <code>Number</code> <code> = 1e-8</code> - -
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*triMap*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#TriMap+init" name="TriMap+init">#</a> <code>*triMap*.**init**</code>
([pca], [knn])


- [pca] <code>Matrix</code> <code> = </code> - Initial Embedding (if null then PCA gets used).
- [knn] <code>KNN</code> <code> = </code> - KNN Object (if null then BallTree gets used).


<a href="#TriMap+_generate_triplets" name="TriMap+_generate_triplets">#</a> <code>*triMap*.**_generate_triplets**</code>
(n_inliers, n_outliers, n_random)

Generates [n_inliers](n_inliers) x [n_outliers](n_outliers) x [n_random](n_random) triplets.


- n_inliers <code>Number</code>
- n_outliers <code>Number</code>
- n_random <code>Number</code>


<a href="#TriMap+_grad" name="TriMap+_grad">#</a> <code>*triMap*.**_grad**</code>
(Y)

Computes the gradient for updating the embedding.


- Y <code>Matrix</code> - The embedding


<a href="#TriMap+transform" name="TriMap+transform">#</a> <code>*triMap*.**transform**</code>
(max_iteration)

**Overrides**: [<code>transform</code>](#DR+transform)  

- max_iteration <code>Number</code> <code> = 400</code>


<a href="#TriMap+generator" name="TriMap+generator">#</a> <code>*triMap*.**generator**</code>
(max_iteration)

**Overrides**: [<code>generator</code>](#DR+generator)  

- max_iteration <code>Number</code> <code> = 800</code>


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*triMap*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*triMap*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*triMap*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#TSNE" name="TSNE">#</a> <code>*dimensionality_reduction***TSNE**</code>


**Extends**: [<code>DR</code>](#DR)  

* [TSNE](#TSNE)

    * [new exports.TSNE(X, parameters)](#new_TSNE_new)

    * [.projection](#DR+projection)

    * [.init(distance_matrix)](#TSNE+init)

    * [.transform([iterations])](#TSNE+transform)

    * [.generator([iterations])](#TSNE+generator)

    * [.parameter(name, [value])](#DR+parameter)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_TSNE_new" name="new_TSNE_new">#</a> new <code>**exports.TSNE**</code>
(X, parameters)


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.perplexity] <code>Number</code> <code> = 50</code> - perplexity.
    - [.epsilon] <code>Number</code> <code> = 10</code> - learning parameter.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code> - the metric which defines the distance between two points.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*tsnE*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#TSNE+init" name="TSNE+init">#</a> <code>*tsnE*.**init**</code>
(distance_matrix)


- distance_matrix <code>Matrix</code> - accepts a precomputed distance matrix


<a href="#TSNE+transform" name="TSNE+transform">#</a> <code>*tsnE*.**transform**</code>
([iterations])

**Overrides**: [<code>transform</code>](#DR+transform)  

- [iterations] <code>Number</code> <code> = 500</code> - number of iterations.


<a href="#TSNE+generator" name="TSNE+generator">#</a> <code>*tsnE*.**generator**</code>
([iterations])

**Overrides**: [<code>generator</code>](#DR+generator)  

- [iterations] <code>Number</code> <code> = 500</code> - number of iterations.


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*tsnE*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*tsnE*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*tsnE*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  

<a href="#UMAP" name="UMAP">#</a> <code>*dimensionality_reduction***UMAP**</code>


**Extends**: [<code>DR</code>](#DR)  

* [UMAP](#UMAP)

    * [new exports.UMAP(X, parameters)](#new_UMAP_new)

    * [.projection](#DR+projection)

    * [.init()](#UMAP+init)

    * [.transform([iterations])](#UMAP+transform)

    * [.generator([iterations])](#UMAP+generator)

    * [.parameter(name, [value])](#DR+parameter)

    * [.check_init()](#DR+check_init)

    * [.transform_async(...args)](#DR+transform_async)



<a href="#new_UMAP_new" name="new_UMAP_new">#</a> new <code>**exports.UMAP**</code>
(X, parameters)


- X <code>Matrix</code> - the high-dimensional data.
- parameters <code>Object</code> - Object containing parameterization of the DR method.
    - [.n_neighbors] <code>Number</code> <code> = 15</code> - size of the local neighborhood.
    - [.local_connectivity] <code>Number</code> <code> = 1</code> - number of nearest neighbors connected in the local neighborhood.
    - [.min_dist] <code>Number</code> <code> = 1</code> - controls how tightly points get packed together.
    - [.d] <code>Number</code> <code> = 2</code> - the dimensionality of the projection.
    - [.metric] <code>function</code> <code> = euclidean</code> - the metric which defines the distance between two points in the high-dimensional space.
    - [._spread] <code>Number</code> <code> = 1</code> - The effective scale of embedded points. (In combination with [parameters.min_dist](parameters.min_dist))
    - [._set_op_mix_ratio] <code>Number</code> <code> = 1</code> - Interpolate between union and intersection.
    - [._repulsion_strength] <code>Number</code> <code> = 1</code> - Weighting applied to negative samples.
    - [._negative_sample_rate] <code>Number</code> <code> = 5</code> - The number of negative samples per positive sample.
    - [._n_epochs] <code>Number</code> <code> = 350</code> - The number of training epochs.
    - [._initial_alpha] <code>Number</code> <code> = 1</code> - The initial learning rate for the optimization.
    - [.seed] <code>Number</code> <code> = 1212</code> - the seed for the random number generator.


<a href="#DR+projection" name="DR+projection">#</a> <code>*umaP*.**projection**</code>


**Overrides**: [<code>projection</code>](#DR+projection)  
**Returns**: <code>Matrix</code> \| <code>Array</code> - Returns the projection.  

<a href="#UMAP+init" name="UMAP+init">#</a> <code>*umaP*.**init**</code>
()

Computes all necessary


<a href="#UMAP+transform" name="UMAP+transform">#</a> <code>*umaP*.**transform**</code>
([iterations])

**Overrides**: [<code>transform</code>](#DR+transform)  

- [iterations] <code>Number</code> <code> = 350</code> - number of iterations.


<a href="#UMAP+generator" name="UMAP+generator">#</a> <code>*umaP*.**generator**</code>
([iterations])

**Overrides**: [<code>generator</code>](#DR+generator)  

- [iterations] <code>Number</code> <code> = 350</code> - number of iterations.


<a href="#DR+parameter" name="DR+parameter">#</a> <code>*umaP*.**parameter**</code>
(name, [value])

**Overrides**: [<code>parameter</code>](#DR+parameter)  
Set and get parameters


- name <code>String</code> - name of the parameter.
- [value] <code>any</code> <code> = </code> - value of the parameter to set.

**Returns**: [<code>DR</code>](#DR) \| <code>any</code> - - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.  
**Example**  
```js
const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
DR.parameter("d"); // returns 3,
DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
```

<a href="#DR+check_init" name="DR+check_init">#</a> <code>*umaP*.**check_init**</code>
()

**Overrides**: [<code>check\_init</code>](#DR+check_init)  
If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.


<a href="#DR+transform_async" name="DR+transform_async">#</a> <code>*umaP*.**transform_async**</code>
(...args)

**Overrides**: [<code>transform\_async</code>](#DR+transform_async)  

- ...args <code>any</code> - Arguments the transform method of the respective DR method takes.

**Returns**: <code>Promise</code> - - A promise yielding the dimensionality reduced dataset.  
