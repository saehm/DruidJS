export declare class Utils {
    static array2D: (M: number, N: number) => number[][];
    static array(N: number): number[];
    static fillArray<T>(N: number, value: T): T[];
    static random2D: (M: number, N: number) => number[][];
    static randomArray(N: number): number[];
    static distance: (featureVectors: number[][]) => number[][];
    static euclidean: (u: number[], v: number[]) => number;
    static sub: (u: number[], v: number[]) => number[];
    static add: (u: number[], v: number[]) => number[];
    static mult: (u: number[], v: number[]) => number[];
    static norm: (u: number[]) => number;
    static shuffle: (u: any[]) => any[];
}
