/**
 * @author Daniel Karl, IBM Research
 */
 import {Utils} from "./Utils";
 export class PowerMethod {
 	static singularValueDecomposition = (U: number[][], K: number): Decomposition => {
 		let MMT:number[][] = PowerMethod.selfProd(U);
 		let eigen:Decomposition = PowerMethod.eigenDecomposition(MMT, K);
 		for(let k = 0; k < K; k++) {
 			eigen.values[k] = Math.sqrt(eigen.values[k]);
 		}
 		let eVecs:number[][] = PowerMethod.product(U, eigen.vectors);
		for(let k = 0; k < K; k++) {
			PowerMethod.normalize(eVecs[k]);
		} 		
 		return new Decomposition(eigen.values, eVecs);
 	}

 	private static product = (U:number[][], V:number[][]):number[][] => {
 		let result:number[][] = Utils.array2D(V.length, U[0].length);
 		for(let m = 0; m < result.length; m++) {
 			for(let i = 0; i < U[0].length; i++) {
 				result[m][i] = 0;
 				for(let j = 0; j < U.length; j++) {
 					result[m][i] += U[j][i] * V[m][j];
 				}
 			}
 		}
 		return result;
 	}

 	static eigenDecomposition = (U: number[][], K: number): Decomposition => {
 		let N = U.length;
 		let eVecs:number[][] = PowerMethod.matrix(K, N);
 		let eVals:number[] = PowerMethod.vector(K);
 		let temp:number[][] = Utils.array2D(K, N);
 		let epsilon:number = 0.0000000001;
 		let r:number, sumEVals:number = 0, sumEValsLast:number, fac:number;
 		// PowerMethod.randomize(eVecs);
 		do {
 			for(let k = 0; k < K; k++) {
 				for(let n = 0; n < N; n++) {
 					temp[k][n] = eVecs[k][n];
 					eVecs[k][n] = 0;
 				}
 			}
 			for (let k = 0; k < K; k++) {
		    	for (let i = 0; i < N; i++) {
		        	for (let j = 0; j < N; j++) {
						eVecs[k][j] += U[i][j] * temp[k][i];
		          	}
		        }
		    }
		    for (let k = 0; k < K; k++) {
		    	for (let p = 0; p < k; p++) {
		    		fac = PowerMethod.scalar(eVecs[k], eVecs[p]) / PowerMethod.scalar(eVecs[p], eVecs[p]);
		    		for (let n = 0; n < N; n++) {
		    			eVecs[k][n] -= fac * eVecs[p][n];
		    		}
		    	}
		    }
		    for (let k = 0; k < K; k++) {
				eVals[k] = PowerMethod.normalize(eVecs[k]);
			}
			r = 1;
			sumEValsLast = sumEVals;
			sumEVals = 0;
			for (let k = 0; k < K; k++) {
				r = Math.min(Math.abs(PowerMethod.scalar(eVecs[k], temp[k])), r);
				sumEVals += eVals[k];
			}
 		} while(r < 1 - epsilon && Math.abs(sumEVals - sumEValsLast) > epsilon);
 		return new Decomposition(eVals, eVecs);
 	}

 	private static selfProd = (U: number[][]):number[][] => {
 		let N:number = U.length;
 		let M:number = U[0].length;
 		let result:number[][] = Utils.array2D(N, M);
 		let sum:number;
 		for(let i = 0; i < N; i++) {
 			for(let j = 0; j <= i; j++) {
 				sum = 0;
 				for(let k = 0; k < M; k++) {
 					sum += U[i][k] * U[j][k];
 				}
 				result[i][j] = result[j][i] = sum;
 			}
 		}
 		return result;
 	}

 	private static scalar = (v: number[], u: number[]):number => {
 		let s:number = 0;
 		for(let i = 0; i < v.length; i++) {
 			s+= v[i] * u[i];
 		}
 		return s;
 	}

 	private static normalize = (v: number[]):number => {
 		let norm:number = Math.sqrt(PowerMethod.scalar(v, v));
 		for(let i = 0; i < v.length; i++) {
 			v[i] /= norm;
 		}
 		return norm;
 	}

 	private static vector = (length: number): number[] => {
 		let v:number[] = Utils.array(length);
 		for(let i = 0; i < length; i++) {
 			v[i] = Math.random();
 		}
 		return v;
 	}

 	private static matrix = (i:number, j:number):number[][] => {
 		let m:number[][] = Utils.array2D(i, j);
 		for(let k = 0; k < i; k++) {
 			for(let l = 0; l < j; l++) {
 				m[k][l] = Math.random();
 			}
 		}
 		return m;
 	}
 }

 export class Decomposition {
 	constructor(public values: number[], public vectors: number[][]) {}
 }
