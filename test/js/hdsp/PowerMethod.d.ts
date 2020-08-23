export declare class PowerMethod {
    static singularValueDecomposition: (U: number[][], K: number) => Decomposition;
    private static product;
    static eigenDecomposition: (U: number[][], K: number) => Decomposition;
    private static selfProd;
    private static scalar;
    private static normalize;
    private static vector;
    private static matrix;
}
export declare class Decomposition {
    values: number[];
    vectors: number[][];
    constructor(values: number[], vectors: number[][]);
}
