## Modules

<dl>
<dt><a href="#module_linear_algebra">linear_algebra</a></dt>
<dd></dd>
</dl>

## Typedefs

<dl>
<dt><a href="#Eigenpair">Eigenpair</a> : <code><a href="#Eigenpair">Eigenpair</a></code></dt>
<dd></dd>
</dl>


<a name="module_linear_algebra">#</a> <code>**linear_algebra**</code>



* [linear_algebra](#module_linear_algebra)

    * [inner_product(a, b)](#inner_product)

    * [poweriteration_n(data, n, x, beta, max_iter, seed)](#poweriteration_n)

    * [qr_givens(A)](#qr_givens)

    * [qr_householder(A)](#qr_householder)

    * [qr(A)](#qr)

    * [simultaneous_poweriteration(A, k, parameters)](#simultaneous_poweriteration)

    * [svrg(data, x, beta, epoch, m, s, seed)](#svrg)



<a name="inner_product">#</a> <code>*linear_algebra***inner_product**</code>
(a, b)

Computes the inner product between two arrays of the same length.


- a <code>Array</code> | <code>Float64Array</code> - Array a
- b <code>Array</code> | <code>Float64Array</code> - Array b

**Returns**: The inner product between [a](a) and [b](b)  

<a name="poweriteration_n">#</a> <code>*linear_algebra***poweriteration_n**</code>
(data, n, x, beta, max_iter, seed)

Computes the [n](n) biggest Eigenpair of the Matrix [data](data).


- data <code>Matrix</code> - the data matrix
- n <code>int</code> - Number of Eigenvalues / Eigenvectors
- x <code>Matrix</code> - Initial Point as 1 times cols Matrix
- beta <code>number</code> - momentum parameter
- max_iter <code>number</code> - maximum number of iterations
- seed <code>number</code> - seed for the random number generator

**Returns**: [<code>Eigenpair</code>](#Eigenpair) - The [n](n) Eigenpairs.  

<a name="qr_givens">#</a> <code>*linear_algebra***qr_givens**</code>
(A)

**See**: [https://en.wikipedia.org/wiki/QR_decomposition#Using_Givens_rotations](https://en.wikipedia.org/wiki/QR_decomposition#Using_Givens_rotations)  
Computes the QR Decomposition of the Matrix [A](A) with givens rotation.


- A <code>Matrix</code>


<a name="qr_householder">#</a> <code>*linear_algebra***qr_householder**</code>
(A)

**See**

- [https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections](https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections)
- [http://mlwiki.org/index.php/Householder_Transformation](http://mlwiki.org/index.php/Householder_Transformation)

Computes the QR Decomposition of the Matrix [A](A) with householder transformations.


- A <code>Matrix</code>


<a name="qr">#</a> <code>*linear_algebra***qr**</code>
(A)

**See**: [https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process](https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process)  
Computes the QR Decomposition of the Matrix [A](A) using Gram-Schmidt process.


- A <code>Matrix</code>


<a name="simultaneous_poweriteration">#</a> <code>*linear_algebra***simultaneous_poweriteration**</code>
(A, k, parameters)

Computes the [k](k) biggest Eigenvectors and Eigenvalues from Matrix [A](A) with the QR-Algorithm.


- A <code>Matrix</code> - The Matrix
- k <code>Number</code> - The number of eigenvectors and eigenvalues to compute.
- parameters <code>Object</code> - Object containing parameterization of the simultanious poweriteration method.
    - [.max_iterations] <code>Number</code> <code> = 100</code> - The number of maxiumum iterations the algorithm should run.
    - [.seed] <code>Number</code> | <code>Randomizer</code> <code> = 1212</code> - The seed value or a randomizer used in the algorithm.
    - [.qr] <code>function</code> <code> = qr_gramschmidt</code> - The QR technique to use.
    - [.tol] <code>Number</code> <code> = 1e-8</code> - Allowed error for stopping criteria

**Returns**: <code>Object</code> - - The [k](k) biggest eigenvectors and eigenvalues of Matrix [A](A).  

<a name="svrg">#</a> <code>*linear_algebra***svrg**</code>
(data, x, beta, epoch, m, s, seed)

**See**: [https://arxiv.org/abs/1707.02670](https://arxiv.org/abs/1707.02670)  
Computes the eigenvector of [X](X) with an accelerated stochastic power iteration algorithm.


- data <code>Matrix</code> - the data matrix
- x <code>Matrix</code> - Initial Point as 1 times cols Matrix
- beta <code>number</code> - momentum parameter
- epoch <code>number</code> - number of epochs
- m <code>number</code> - epoch length
- s <code>number</code> - mini-batch size
- seed <code>number</code> - seed for the random number generator


<a name="Eigenpair">#</a> <code>**Eigenpair**</code>


**Properties**

| Name | Type | Description |
| --- | --- | --- |
| Eigenvalues | <code>Array</code> | Array of Eigenvalues |
| Eigenvectors | <code>Array.&lt;Array&gt;</code> | Array of Eigenvectors |

