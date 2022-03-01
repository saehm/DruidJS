
<a name="module_utils">#</a> <code>**utils**</code>



* [utils](#module_utils)

    * [max(values)](#max)

    * [min(values)](#min)

    * _instance_
        * [.Randomizer](#Randomizer)

            * [new exports.Randomizer([_seed])](#new_Randomizer_new)

            * _instance_
                * [.seed](#Randomizer+seed)

                * [.random](#Randomizer+random)

                * [.random_int](#Randomizer+random_int)

                * [.choice(A, n)](#Randomizer+choice)

            * _static_
                * [.choice(A, n, seed)](#Randomizer.choice)



<a name="max">#</a> <code>*utils***max**</code>
(values)

Returns maximum in Array [values](values).


- values <code>Array</code>


<a name="min">#</a> <code>*utils***min**</code>
(values)

Returns maximum in Array [values](values).


- values <code>Array</code>


<a name="Randomizer">#</a> <code>*utils*.**Randomizer**</code>


**See**: https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js  

* [.Randomizer](#Randomizer)

    * [new exports.Randomizer([_seed])](#new_Randomizer_new)

    * _instance_
        * [.seed](#Randomizer+seed)

        * [.random](#Randomizer+random)

        * [.random_int](#Randomizer+random_int)

        * [.choice(A, n)](#Randomizer+choice)

    * _static_
        * [.choice(A, n, seed)](#Randomizer.choice)



<a name="new_Randomizer_new">#</a> new <code>**exports.Randomizer**</code>
([_seed])

Mersenne Twister random number generator.


- [_seed] <code>Number</code> <code> = new Date().getTime()</code> - The seed for the random number generator. If <code>_seed == null</code> then the actual time gets used as seed.


<a name="Randomizer+seed">#</a> <code>*randomizer*.**seed**</code>


Returns the seed of the random number generator.

**Returns**: <code>Number</code> - - The seed.  

<a name="Randomizer+random">#</a> <code>*randomizer*.**random**</code>


Returns a float between 0 and 1.

**Returns**: <code>Number</code> - - A random number between [0, 1]  

<a name="Randomizer+random_int">#</a> <code>*randomizer*.**random_int**</code>


Returns an integer between 0 and MAX_INTEGER.

**Returns**: <code>Integer</code> - - A random integer.  

<a name="Randomizer+choice">#</a> <code>*randomizer*.**choice**</code>
(A, n)

Returns samples from an input Matrix or Array.


- A <code>Matrix</code> | <code>Array</code> | <code>Float64Array</code> - The input Matrix or Array.
- n <code>Number</code> - The number of samples.

**Returns**: <code>Array</code> - - A random selection form [A](A) of [n](n) samples.  

<a name="Randomizer.choice">#</a> <code>*Randomizer*.**choice**</code>
(A, n, seed)


- A <code>Matrix</code> | <code>Array</code> | <code>Float64Array</code> - The input Matrix or Array.
- n <code>Number</code> - The number of samples.
- seed <code>Number</code> <code> = 1212</code> - The seed for the random number generator.

**Returns**: <code>Array</code> - - A random selection form [A](A) of [n](n) samples.  
