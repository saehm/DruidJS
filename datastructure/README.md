
<a name="module_datastructure">#</a> <code>**datastructure**</code>



* [datastructure](#module_datastructure)

    * [DisjointSet](#DisjointSet)

        * [new exports.DisjointSet([elements])](#new_DisjointSet_new)

    * [Heap](#Heap)

        * [new exports.Heap([elements], [accessor], [comparator])](#new_Heap_new)

        * _instance_
            * [.first](#Heap+first)

            * [.length](#Heap+length)

            * [.empty](#Heap+empty)

            * [.push(element)](#Heap+push)

            * [.pop()](#Heap+pop)

            * [.iterate()](#Heap+iterate)

            * [.toArray()](#Heap+toArray)

            * [.data()](#Heap+data)

            * [.raw_data()](#Heap+raw_data)

        * _static_
            * [.heapify(elements, [accessor], [comparator])](#Heap.heapify)



<a name="DisjointSet">#</a> <code>*datastructure***DisjointSet**</code>


**See**: [https://en.wikipedia.org/wiki/Disjoint-set_data_structure](https://en.wikipedia.org/wiki/Disjoint-set_data_structure)  

<a name="new_DisjointSet_new">#</a> new <code>**exports.DisjointSet**</code>
([elements])


- [elements] <code>Array</code> <code> = </code>


<a name="Heap">#</a> <code>*datastructure***Heap**</code>


**See**: [https://en.wikipedia.org/wiki/Binary_heap](https://en.wikipedia.org/wiki/Binary_heap)  

* [Heap](#Heap)

    * [new exports.Heap([elements], [accessor], [comparator])](#new_Heap_new)

    * _instance_
        * [.first](#Heap+first)

        * [.length](#Heap+length)

        * [.empty](#Heap+empty)

        * [.push(element)](#Heap+push)

        * [.pop()](#Heap+pop)

        * [.iterate()](#Heap+iterate)

        * [.toArray()](#Heap+toArray)

        * [.data()](#Heap+data)

        * [.raw_data()](#Heap+raw_data)

    * _static_
        * [.heapify(elements, [accessor], [comparator])](#Heap.heapify)



<a name="new_Heap_new">#</a> new <code>**exports.Heap**</code>
([elements], [accessor], [comparator])

A heap is a datastructure holding its elements in a specific way, so that the top element would be the first entry of an ordered list.


- [elements] <code>Array</code> <code> = </code> - Contains the elements for the Heap. [elements](elements) can be null.
- [accessor] <code>function</code> <code> = (d) &#x3D;&gt; d</code> - Function returns the value of the element.
- [comparator] <code>&quot;min&quot;</code> | <code>&quot;max&quot;</code> | <code>function</code> <code> = &quot;min&quot;</code> - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)


<a name="Heap+first">#</a> <code>*heap*.**first**</code>


Returns the top entry of the heap without removing it.

**Returns**: <code>Object</code> - Object consists of the element and its value (computed by [accessor](accessor)).  

<a name="Heap+length">#</a> <code>*heap*.**length**</code>


The size of the heap.


<a name="Heap+empty">#</a> <code>*heap*.**empty**</code>


Returns false if the the heap has entries, true if the heap has no entries.


<a name="Heap+push">#</a> <code>*heap*.**push**</code>
(element)

Pushes the element to the heap.


- element


<a name="Heap+pop">#</a> <code>*heap*.**pop**</code>
()

Removes and returns the top entry of the heap.

**Returns**: <code>Object</code> - Object consists of the element and its value (computed by [accessor](accessor)).  

<a name="Heap+iterate">#</a> <code>*heap*.**iterate**</code>
()

Yields the raw data


<a name="Heap+toArray">#</a> <code>*heap*.**toArray**</code>
()

Returns the heap as ordered array.

**Returns**: <code>Array</code> - Array consisting the elements ordered by [comparator](comparator).  

<a name="Heap+data">#</a> <code>*heap*.**data**</code>
()

Returns elements of container array.

**Returns**: <code>Array</code> - Array consisting the elements.  

<a name="Heap+raw_data">#</a> <code>*heap*.**raw_data**</code>
()

Returns the container array.

**Returns**: <code>Array</code> - The container array.  

<a name="Heap.heapify">#</a> <code>*Heap*.**heapify**</code>
(elements, [accessor], [comparator])

Creates a Heap from an Array


- elements <code>Array</code> | <code>Set</code> - Contains the elements for the Heap.
- [accessor] <code>function</code> <code> = (d) &#x3D;&gt; d</code> - Function returns the value of the element.
- [comparator] <code>String</code> | <code>function</code> <code> = &quot;min&quot;</code> - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)

