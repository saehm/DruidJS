/**
 * @author Daniel Karl, IBM Research
 */
export class Utils {
	static array2D = (M:number, N:number):number[][] => {
		let result:number[][] = [];
		for(let m = 0; m < M; m++) {
			result.push(Utils.fillArray(N, 0));
		}
		return result;
	};

	static array(N:number):number[] {
		return Utils.fillArray<number>(N, 0);
	};

	static fillArray<T>(N:number, value:T):T[] {
	    let result:T[] = [];
		for(let n = 0; n < N; n++) {
			result.push(value);
		}
		return result;
    };

	static random2D = (M:number, N:number): number[][] => {
		let result:number[][] = [];
		for(let m = 0; m < M; m++) {
			result.push(Utils.randomArray(N));
		}
		return result;
	};

	static randomArray(N:number):number[] {
	    let result:number[] = [];
		for(let n = 0; n < N; n++) {
			result.push(Math.random());
		}
		return result;
    };

	static distance = (featureVectors: number[][]): number[][] => {
		let N = featureVectors.length;
		let result:number[][] = Utils.array2D(N, N);
		for(let n = 0; n < N; n++) {
			for(let m = 0; m <= n; m++) {
				result[n][m] = result[m][n] = Utils.euclidean(featureVectors[n], featureVectors[m]);
			}
		}
		return result;
	};

	static euclidean = (u:number[], v:number[]): number => {
		let sum:number = 0;
		for(let i = 0; i < u.length; i++) {
			sum += Math.pow(u[i] - v[i], 2);
		}
		return Math.sqrt(sum);
	};

	static sub = (u: number[], v:number[]): number[] => {
		let sub:number[] = [];
		for(let i = 0; i < u.length; i++) {
			sub.push(u[i] - v[i]);
		}
		return sub;
	};

	static add = (u: number[], v:number[]): number[] => {
		let add:number[] = [];
		for(let i = 0; i < u.length; i++) {
			add.push(u[i] + v[i]);
		}
		return add;
	};

	static mult = (u: number[], v:number[]): number[] => {
		let mult:number[] = [];
		for(let i = 0; i < u.length; i++) {
			mult.push(u[i] * v[i]);
		}
		return mult;
	};

	static norm = (u: number[]): number => {
		let squaredSum = 0;
		for(let i = 0; i < u.length; i++) {
			squaredSum += u[i] ** 2;
		}
		return Math.sqrt(squaredSum);
	};

	static shuffle = (u:any[]): any[] => {
		return u.sort(() => Math.random() - 0.5);
	};
}
