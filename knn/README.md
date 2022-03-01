
<a name="module_knn">#</a> <code>**knn**</code>



* [knn](#module_knn)

    * [BallTree](#BallTree)

        * [new exports.BallTree([elements], [metric])](#new_BallTree_new)

        * [.add(elements)](#BallTree+add)

        * [.search(t, [k])](#BallTree+search)

    * [HNSW](#HNSW)

        * [new exports.HNSW([metric], [heuristic], [m], [ef], [m0], [mL], [seed])](#new_HNSW_new)

        * [.add(elements)](#HNSW+add)

        * [.search(q, K, ef)](#HNSW+search)

        * [.search_iter(q, K, ef)](#HNSW+search_iter)

    * [KNN](#KNN)

        * [new exports.KNN([elements], [metric])](#new_KNN_new)

        * [.search(t, [k])](#KNN+search)

    * [NNDescent](#NNDescent)

        * [new exports.NNDescent([elements], [metric], [K], [rho], [delta], [seed])](#new_NNDescent_new)

        * [.add(elements)](#NNDescent+add)

        * [.search(x, k)](#NNDescent+search)



<a name="BallTree">#</a> <code>*knn***BallTree**</code>


**See**

- [https://en.wikipedia.org/wiki/Ball_tree](https://en.wikipedia.org/wiki/Ball_tree)
- [https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js](https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js)


* [BallTree](#BallTree)

    * [new exports.BallTree([elements], [metric])](#new_BallTree_new)

    * [.add(elements)](#BallTree+add)

    * [.search(t, [k])](#BallTree+search)



<a name="new_BallTree_new">#</a> new <code>**exports.BallTree**</code>
([elements], [metric])

Generates a BallTree with given [elements](elements).


- [elements] <code>Array</code> <code> = </code> - Elements which should be added to the BallTree
- [metric] <code>function</code> <code> = euclidean</code> - metric to use: (a, b) => distance


<a name="BallTree+add">#</a> <code>*ballTree*.**add**</code>
(elements)


- elements <code>Array.&lt;\*&gt;</code> - new elements.


<a name="BallTree+search">#</a> <code>*ballTree*.**search**</code>
(t, [k])


- t <code>\*</code> - query element.
- [k] <code>Number</code> <code> = 5</code> - number of nearest neighbors to return.

**Returns**: <code>Heap</code> - - Heap consists of the [k](k) nearest neighbors.  

<a name="HNSW">#</a> <code>*knn***HNSW**</code>


**See**

- [https://arxiv.org/abs/1603.09320](https://arxiv.org/abs/1603.09320)
- [https://arxiv.org/pdf/1904.02077](https://arxiv.org/pdf/1904.02077)


* [HNSW](#HNSW)

    * [new exports.HNSW([metric], [heuristic], [m], [ef], [m0], [mL], [seed])](#new_HNSW_new)

    * [.add(elements)](#HNSW+add)

    * [.search(q, K, ef)](#HNSW+search)

    * [.search_iter(q, K, ef)](#HNSW+search_iter)



<a name="new_HNSW_new">#</a> new <code>**exports.HNSW**</code>
([metric], [heuristic], [m], [ef], [m0], [mL], [seed])

Hierarchical navigable small world graph. Efficient and robust approximate nearest neighbor search.


- [metric] <code>function</code> <code> = euclidean</code> - metric to use: (a, b) => distance.
- [heuristic] <code>Boolean</code> <code> = true</code> - use heuristics or naive selection.
- [m] <code>Number</code> <code> = 5</code> - max number of connections.
- [ef] <code>Number</code> <code> = 200</code> - size of candidate list.
- [m0] <code>Number</code> <code> = 2 * m</code> - max number of connections for ground layer.
- [mL] <code>Number</code> <code> = 1 / Math.log2(m)</code> - normalization factor for level generation.
- [seed] <code>Number</code> <code> = 1987</code> - seed for random number generator.


<a name="HNSW+add">#</a> <code>*hnsW*.**add**</code>
(elements)


- elements <code>Array.&lt;\*&gt;</code> - new elements.


<a name="HNSW+search">#</a> <code>*hnsW*.**search**</code>
(q, K, ef)


- q <code>\*</code> - query element.
- K <code>\*</code> - number of nearest neighbors to return.
- ef <code>\*</code> <code> = 1</code> - size of the dynamic cnadidate list.

**Returns**: <code>Array</code> - K nearest elements to q.  

<a name="HNSW+search_iter">#</a> <code>*hnsW*.**search_iter**</code>
(q, K, ef)

Iterator for searching the HNSW graphs


- q <code>\*</code> - query element.
- K <code>\*</code> - number of nearest neighbors to return.
- ef <code>\*</code> <code> = 1</code> - size of the dynamic cnadidate list.


<a name="KNN">#</a> <code>*knn***KNN**</code>



* [KNN](#KNN)

    * [new exports.KNN([elements], [metric])](#new_KNN_new)

    * [.search(t, [k])](#KNN+search)



<a name="new_KNN_new">#</a> new <code>**exports.KNN**</code>
([elements], [metric])

Generates a KNN list with given [elements](elements).


- [elements] <code>Array</code> <code> = </code> - Elements which should be added to the KNN list
- [metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code> - metric is either precomputed or a function to use: (a, b) => distance


<a name="KNN+search">#</a> <code>*knN*.**search**</code>
(t, [k])


- t <code>Array</code> | <code>Number</code> - query element or index.
- [k] <code>Number</code> <code> = 5</code> - number of nearest neighbors to return.

**Returns**: <code>Heap</code> - - Heap consists of the [k](k) nearest neighbors.  

<a name="NNDescent">#</a> <code>*knn***NNDescent**</code>


**See**: [http://www.cs.princeton.edu/cass/papers/www11.pdf](http://www.cs.princeton.edu/cass/papers/www11.pdf)  

* [NNDescent](#NNDescent)

    * [new exports.NNDescent([elements], [metric], [K], [rho], [delta], [seed])](#new_NNDescent_new)

    * [.add(elements)](#NNDescent+add)

    * [.search(x, k)](#NNDescent+search)



<a name="new_NNDescent_new">#</a> new <code>**exports.NNDescent**</code>
([elements], [metric], [K], [rho], [delta], [seed])


- [elements] <code>Array.&lt;\*&gt;</code> - called V in paper.
- [metric] <code>function</code> <code> = euclidean</code> - called sigma in paper.
- [K] <code>Number</code> <code> = 10</code> - number of neighbors [search](search) should return.
- [rho] <code>Number</code> <code> = .8</code> - sample rate.
- [delta] <code>Number</code> <code> = 0.0001</code> - precision parameter.
- [seed] <code>Number</code> <code> = 1987</code> - seed for the random number generator.


<a name="NNDescent+add">#</a> <code>*nnDescent*.**add**</code>
(elements)


- elements <code>Array</code>


<a name="NNDescent+search">#</a> <code>*nnDescent*.**search**</code>
(x, k)

**Todo**

- [ ] not implemented yet


- x <code>\*</code>
- k <code>\*</code> <code> = 5</code>

