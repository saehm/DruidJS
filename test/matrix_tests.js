var tape = require("tape");
var sci = require('../dist/sci');

tape("Matrix, constructor", function(test) {
    let expected = [[1,2,3], [4,5,6], [7,8,9]]
    let A = new sci.Matrix(expected);
    actual = A.to2dArray;
    test.deepEquals(actual, expected);

    test.end()
})

tape("Matrix, shape", function(test) {
    let A = new sci.Matrix();
    A.entries = [[1,2], [3,4]];
    let expected = [2, 2];
    let actual = A.shape;
    test.deepEquals(actual, expected);

    A.entries = [[1,2,3], [4,5,6], [7,8,9]]
    expected = [3, 3];
    actual = A.shape;
    test.deepEquals(actual, expected);

    test.end()
})

tape("Matrix, row() and col()", function(test) {
    let A = new sci.Matrix();
    A.entries = [[1,2,3], [4,5,6], [7,8,9]]

    expected = [1,4,7]
    actual = A.col(0)
    test.deepEquals(actual, expected);

    expected = [2,5,8]
    actual = A.col(1)
    test.deepEquals(actual, expected);
    test.end()
})

tape("Matrix, entry()", function(test) {
    let A = new sci.Matrix();
    A.entries = [[1,2,3], [4,5,6], [7,8,9]]

    expected = 4
    actual = A.entry(1,0)
    test.deepEquals(actual, expected);

    expected = 9
    actual = A.entry(2,2)
    test.deepEquals(actual, expected);
    test.end()
})

tape("Matrix, to2dArray()", function(test) {
    let A = new sci.Matrix();
    let expected = [[1,2,3], [4,5,6], [7,8,9]];
    A.entries = expected

    actual = A.to2dArray
    test.deepEquals(actual, expected);
    test.end()
})

tape("Matrix, set shape / transpose", function(test) {
    let A = new sci.Matrix([[1,2,3], [4,5,6], [7,8,9]])
    let B = A.transpose();
    let expected = [[1,4,7], [2,5,8], [3,6,9]];
    let actual = B.to2dArray;
    test.deepEquals(actual, expected);

    expected = [[1,2,3], [4,5,6], [7,8,9]];
    actual = A.to2dArray;
    test.deepEquals(actual, expected);

    A = new sci.Matrix();
    A.shape = [3,2];
    expected = [[0,0], [0,0], [0,0]];
    actual = A.to2dArray;
    test.deepEquals(actual, expected);
    test.end()
})

tape("Matrix, dot", function(test) {
    let expected = [[1,2,3], [4,5,6], [7,8,9]];
    let A = new sci.Matrix(expected)
    let B = new sci.Matrix()
    B.shape = [3,3,(i,j) => i === j ? 1 : 0]
    let actual = A.dot(B).to2dArray;
    test.deepEquals(actual, expected);
    test.end()
})
