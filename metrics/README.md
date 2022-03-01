
<a name="module_metrics">#</a> <code>**metrics**</code>



* [metrics](#module_metrics)

    * [canberra(a, b)](#canberra)

    * [chebyshev(a, b)](#chebyshev)

    * [cosine(a, b)](#cosine)

    * [euclidean_squared(a, b)](#euclidean_squared)

    * [euclidean(a, b)](#euclidean)

    * [hamming(a, b)](#hamming)

    * [jaccard(a, b)](#jaccard)

    * [manhattan(a, b)](#manhattan)

    * [sokal_michener(a, b)](#sokal_michener)

    * [yule(a, b)](#yule)



<a name="canberra">#</a> <code>*metrics***canberra**</code>
(a, b)

**See**: [https://en.wikipedia.org/wiki/Canberra_distance](https://en.wikipedia.org/wiki/Canberra_distance)  
Computes the canberra distance between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - The canberra distance between [a](a) and [b](b).  

<a name="chebyshev">#</a> <code>*metrics***chebyshev**</code>
(a, b)

Computes the chebyshev distance (L<sub>∞</sub>) between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the chebyshev distance between [a](a) and [b](b).  

<a name="cosine">#</a> <code>*metrics***cosine**</code>
(a, b)

Computes the cosine distance (not similarity) between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - The cosine distance between [a](a) and [b](b).  
**Example**  
```js
druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4
```

<a name="euclidean_squared">#</a> <code>*metrics***euclidean_squared**</code>
(a, b)

Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the squared euclidean distance between [a](a) and [b](b).  

<a name="euclidean">#</a> <code>*metrics***euclidean**</code>
(a, b)

Computes the euclidean distance (l<sub>2</sub>) between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the euclidean distance between [a](a) and [b](b).  

<a name="hamming">#</a> <code>*metrics***hamming**</code>
(a, b)

Computes the hamming distance between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the hamming distance between [a](a) and [b](b).  

<a name="jaccard">#</a> <code>*metrics***jaccard**</code>
(a, b)

Computes the jaccard distance between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the jaccard distance between [a](a) and [b](b).  

<a name="manhattan">#</a> <code>*metrics***manhattan**</code>
(a, b)

Computes the manhattan distance (l<sub>1</sub>) between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the manhattan distance between [a](a) and [b](b).  

<a name="sokal_michener">#</a> <code>*metrics***sokal_michener**</code>
(a, b)

Computes the Sokal-Michener distance between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the Sokal-Michener distance between [a](a) and [b](b).  

<a name="yule">#</a> <code>*metrics***yule**</code>
(a, b)

Computes the yule distance between [a](a) and [b](b).


- a <code>Array.&lt;Number&gt;</code>
- b <code>Array.&lt;Number&gt;</code>

**Returns**: <code>Number</code> - the yule distance between [a](a) and [b](b).  
