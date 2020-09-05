## DruidJS — A JavaScript Library for Dimensionality Reduction.

### Resources
- [Documentation](https://saehm.github.io/DruidJS/index.html) 
- [Demo](https://renecutura.eu/druid_demo)

### Installation
If you use npm, install with `npm install @saehrimnir/druidjs`, and use it with
```js
import * as druid from "@saehrimnir/druidjs";
```

Otherwise download the files [here](https://github.com/saehm/DruidJS/releases/latest), or use for instance [unpkg](https://unpkg.com/@saehrimnir/druidjs) this way:

```html
<script src="https://unpkg.com/@saehrimnir/druidjs"></script>
```

### Matrix
DruidJS uses internally the Matrix class for storing data. You can use it by creating a `druid.Matrix` object for instance with the function `from`, in example:
```js
    import * as druid from '@saehrimnir/druidjs';

    let data = [[...], [...], ...];
    let matrix = druid.Matrix.from(data);
```
You can create a `druid.Matrix` object programmatically by:
```js
    let fn = (row, col) => row == col ? 1 : 0;
    let matrix = new druid.Matrix(rows, columns, fn);
```
If `rows == columns`, then `matrix` would be a identity matrix.
A shortcut for a identity matrix is:
```js
    let matrix = new druid.Matrix(rows, columns, "I");
    // or
    let matrix = new druid.Matrix(rows, columnbs, "identity");
```
There are more shortcuts for creating matrices:
```js
    let matrix = new druid.Matrix(3, 3, "zeros"); // matrix would be a 3x3 matrix with zeroes
    let matrix = new druid.Matrix(3, 3, "center"); // matrix would be a 3x3 center matrix;
    let number = 12;
    let matrix = new druid.Matrix(3, 3, number); // matrix would b a 3x3 matrix filled with 'number'
```

If you want to use a `druid.Matrix` object, for instance, with [d3](https://d3js.org), you can use either the `to2dArray` property, the `iterate_rows` generator function, or just use the `druid.Matrix` object as an iterable (works with d3 since version 6).
```js
    let data = await d3.csv("data.csv");
    let matrix = druid.Matrix.from(data);
    d3.selectAll("datapoints").data(matrix.to2dArray)//...
    d3.selectAll("datapoints").data(matrix.iterate_rows)//...
    d3.selectAll("datapoints").data(matrix)//...
```
### DR methods
...