
<a href="#module_clustering" name="module_clustering">#</a> <code>**clustering**</code>



* [clustering](#module_clustering)

    * [Hierarchical_Clustering](#Hierarchical_Clustering)

        * [new exports.Hierarchical_Clustering(matrix, [linkage], [metric])](#new_Hierarchical_Clustering_new)

        * [.get_clusters(value, [type])](#Hierarchical_Clustering+get_clusters)

        * [.init()](#Hierarchical_Clustering+init)

        * [.do()](#Hierarchical_Clustering+do)

    * [KMeans](#KMeans)

        * [new exports.KMeans(matrix, K, [metric], [seed], [init])](#new_KMeans_new)

        * [.get_clusters()](#KMeans+get_clusters)

        * [.init(K)](#KMeans+init)

    * [KMedoids](#KMedoids)

        * [new exports.KMedoids(matrix, K, [max_iter], [metric], [seed])](#new_KMedoids_new)

        * [.get_clusters()](#KMedoids+get_clusters)

        * [._iteration()](#KMedoids+_iteration)

        * [.init(K)](#KMedoids+init)

        * [._get_random_medoids(K)](#KMedoids+_get_random_medoids)

    * [OPTICS](#OPTICS)

        * [new exports.OPTICS(matrix, epsilon, min_points, [metric])](#new_OPTICS_new)

        * [.init()](#OPTICS+init)

        * [.get_clusters()](#OPTICS+get_clusters)

        * [.get_cluster_affirmation()](#OPTICS+get_cluster_affirmation)

    * [XMeans](#XMeans)

        * [new exports.XMeans(matrix, K_max, K_min, [metric], [seed])](#new_XMeans_new)



<a href="#Hierarchical_Clustering" name="Hierarchical_Clustering">#</a> <code>*clustering***Hierarchical_Clustering**</code>


**Todo**

- [ ] needs restructuring.


* [Hierarchical_Clustering](#Hierarchical_Clustering)

    * [new exports.Hierarchical_Clustering(matrix, [linkage], [metric])](#new_Hierarchical_Clustering_new)

    * [.get_clusters(value, [type])](#Hierarchical_Clustering+get_clusters)

    * [.init()](#Hierarchical_Clustering+init)

    * [.do()](#Hierarchical_Clustering+do)



<a href="#new_Hierarchical_Clustering_new" name="new_Hierarchical_Clustering_new">#</a> new <code>**exports.Hierarchical_Clustering**</code>
(matrix, [linkage], [metric])


- matrix <code>Matrix</code> - Data or distance matrix if metric is 'precomputed'
- [linkage] <code>&quot;single&quot;</code> | <code>&quot;complete&quot;</code> | <code>&quot;average&quot;</code> <code> = &quot;complete&quot;</code>
- [metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code>


<a href="#Hierarchical_Clustering+get_clusters" name="Hierarchical_Clustering+get_clusters">#</a> <code>*hierarchical_Clustering*.**get_clusters**</code>
(value, [type])


- value <code>Number</code> - value where to cut the tree.
- [type] <code>&quot;distance&quot;</code> | <code>&quot;depth&quot;</code> <code> = &quot;distance&quot;</code> - type of value.

**Returns**: <code>Array.&lt;Array&gt;</code> - - Array of clusters with the indices of the rows in given [matrix](matrix).  

<a href="#Hierarchical_Clustering+init" name="Hierarchical_Clustering+init">#</a> <code>*hierarchical_Clustering*.**init**</code>
()

computes the tree.


<a href="#Hierarchical_Clustering+do" name="Hierarchical_Clustering+do">#</a> <code>*hierarchical_Clustering*.**do**</code>
()

computes the tree.


<a href="#KMeans" name="KMeans">#</a> <code>*clustering***KMeans**</code>


**Todo**

- [ ] needs restructuring.


* [KMeans](#KMeans)

    * [new exports.KMeans(matrix, K, [metric], [seed], [init])](#new_KMeans_new)

    * [.get_clusters()](#KMeans+get_clusters)

    * [.init(K)](#KMeans+init)



<a href="#new_KMeans_new" name="new_KMeans_new">#</a> new <code>**exports.KMeans**</code>
(matrix, K, [metric], [seed], [init])


- matrix <code>Matrix</code>
- K <code>Numbers</code>
- [metric] <code>function</code> <code> = euclidean</code>
- [seed] <code>Number</code> <code> = 1987</code>
- [init] <code>Boolean</code> <code> = true</code>


<a href="#KMeans+get_clusters" name="KMeans+get_clusters">#</a> <code>*kMeans*.**get_clusters**</code>
()

**Returns**: <code>Array.&lt;Array&gt;</code> - - Array of clusters with the indices of the rows in given [matrix](matrix).  

<a href="#KMeans+init" name="KMeans+init">#</a> <code>*kMeans*.**init**</code>
(K)

Computes [K](K) clusters out of the [matrix](matrix).


- K <code>Number</code> - number of clusters.


<a href="#KMedoids" name="KMedoids">#</a> <code>*clustering***KMedoids**</code>


**See**: [https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16](https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16) Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms  
**Todo**

- [ ] needs restructuring.


* [KMedoids](#KMedoids)

    * [new exports.KMedoids(matrix, K, [max_iter], [metric], [seed])](#new_KMedoids_new)

    * [.get_clusters()](#KMedoids+get_clusters)

    * [._iteration()](#KMedoids+_iteration)

    * [.init(K)](#KMedoids+init)

    * [._get_random_medoids(K)](#KMedoids+_get_random_medoids)



<a href="#new_KMedoids_new" name="new_KMedoids_new">#</a> new <code>**exports.KMedoids**</code>
(matrix, K, [max_iter], [metric], [seed])


- matrix <code>Matrix</code> - data matrix
- K <code>Numbers</code> - number of clusters
- [max_iter] <code>number</code> <code> = </code> - maximum number of iterations. Default is 10 * Math.log10(N)
- [metric] <code>function</code> <code> = euclidean</code> - metric defining the dissimilarity
- [seed] <code>Number</code> <code> = 1212</code> - seed value for random number generator


<a href="#KMedoids+get_clusters" name="KMedoids+get_clusters">#</a> <code>*kMedoids*.**get_clusters**</code>
()

**Returns**: <code>Array.&lt;Array&gt;</code> - - Array of clusters with the indices of the rows in given [matrix](matrix).  

<a href="#KMedoids+_iteration" name="KMedoids+_iteration">#</a> <code>*kMedoids*.**_iteration**</code>
()

Algorithm 2. FastPAM2: SWAP with multiple candidates


<a href="#KMedoids+init" name="KMedoids+init">#</a> <code>*kMedoids*.**init**</code>
(K)

Computes [K](K) clusters out of the [matrix](matrix).


- K <code>Number</code> - number of clusters.


<a href="#KMedoids+_get_random_medoids" name="KMedoids+_get_random_medoids">#</a> <code>*kMedoids*.**_get_random_medoids**</code>
(K)

Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.


- K <code>number</code> - number of clusters


<a href="#OPTICS" name="OPTICS">#</a> <code>*clustering***OPTICS**</code>


**See**

- [https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf](https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf)
- [https://en.wikipedia.org/wiki/OPTICS_algorithm](https://en.wikipedia.org/wiki/OPTICS_algorithm)

**Todo**

- [ ] needs restructuring.


* [OPTICS](#OPTICS)

    * [new exports.OPTICS(matrix, epsilon, min_points, [metric])](#new_OPTICS_new)

    * [.init()](#OPTICS+init)

    * [.get_clusters()](#OPTICS+get_clusters)

    * [.get_cluster_affirmation()](#OPTICS+get_cluster_affirmation)



<a href="#new_OPTICS_new" name="new_OPTICS_new">#</a> new <code>**exports.OPTICS**</code>
(matrix, epsilon, min_points, [metric])

**O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.


- matrix <code>Matrix</code> - the data.
- epsilon <code>Number</code> - the minimum distance which defines whether a point is a neighbor or not.
- min_points <code>Number</code> - the minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
- [metric] <code>function</code> <code> = euclidean</code> - the distance metric which defines the distance between two points of the [matrix](matrix).


<a href="#OPTICS+init" name="OPTICS+init">#</a> <code>*opticS*.**init**</code>
()

Computes the clustering.


<a href="#OPTICS+get_clusters" name="OPTICS+get_clusters">#</a> <code>*opticS*.**get_clusters**</code>
()

Returns an array of clusters.

**Returns**: <code>Array.&lt;Array&gt;</code> - Array of clusters with the indices of the rows in given [matrix](matrix).  

<a href="#OPTICS+get_cluster_affirmation" name="OPTICS+get_cluster_affirmation">#</a> <code>*opticS*.**get_cluster_affirmation**</code>
()

**Returns**: <code>Array</code> - Returns an array, where the ith entry defines the cluster affirmation of the ith point of [matrix](matrix). (-1 stands for outlier)  

<a href="#XMeans" name="XMeans">#</a> <code>*clustering***XMeans**</code>


**See**

- [https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf](https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf)
- [https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py](https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py)
- [https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java](https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java)

**Todo**

- [ ] needs restructuring and repairing!!


<a href="#new_XMeans_new" name="new_XMeans_new">#</a> new <code>**exports.XMeans**</code>
(matrix, K_max, K_min, [metric], [seed])


- matrix <code>Matrix</code>
- K_max <code>Numbers</code> <code> = 10</code>
- K_min <code>Numbers</code> <code> = 2</code>
- [metric] <code>function</code> <code> = euclidean</code>
- [seed] <code>Number</code> <code> = 1987</code>

