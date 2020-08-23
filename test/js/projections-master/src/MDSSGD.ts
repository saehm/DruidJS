/**
 * @author Daniel Karl, IBM Research
 *
 * Literature:
 * Zheng, Jonathan X., Samraat Pawar, and Dan FM Goodman. "Graph drawing by stochastic gradient descent."
 * IEEE transactions on visualization and computer graphics 25.9 (2018): 2738-2748.
 */
import {Utils} from "./Utils";

export class MDSSGD {

    /**
	 * Project a feature space into `D` dimensions.
	 * @param featureVectors the feature vectors
	 * @param D the output dimensionality
	 * @return
	 */
    static project = (featureVectors: number[][], D:number): number[][] => {
        if(featureVectors == null) {
			return null;
		}
		let N = featureVectors.length;
		// empty input
		if(N == 0) {
			return [];
		}
		let distances:number[][] = Utils.distance(featureVectors);
		let result:number[][] = Utils.random2D(N, D);

		// add constraints
		let constraints = [];
		let w;
		let w_min = Number.MAX_VALUE;
		let w_max = Number.MIN_VALUE;

		for(let i = 0; i < distances.length; i++) {
		    for(let j = 0; j < i; j++) {
		        w = 1/distances[i][j]**2;
		        w_min = Math.min(w_min, w);
		        w_max = Math.max(w_max, w);
		        constraints.push([i, j, w]);
            }
        }

		// setup annealing schedule
        let num_iter = 60;
		let epsilon = 0.1;
        let eta_max = 1 / w_min;
		let eta_min = epsilon / w_max;
		let lambda = Math.log(eta_min / eta_max) / (num_iter - 1);
		let eta = (t:number) => {
			return eta_max * lambda**t
		};
		let schedule = [];
		for(let i = 0; i < num_iter; i++) {
		    schedule.push(eta(i));
        }

		let wc;
		let pq;
		let mag;
		let r;
		let m;
		for(let c of schedule) {
		    // shuffle order
            constraints = Utils.shuffle(constraints);

            for(let [i, j, w] of constraints) {
                wc = w * c;
                if(wc > 1) {
                    wc = 1;
                }
                pq = Utils.sub(result[i], result[j]);
                mag = Utils.norm(pq);
                r = (distances[i][j] - mag) / 2;
                m = Utils.mult(pq, Utils.fillArray(D, wc * r / mag));

                result[i] = Utils.add(result[i], m);
                result[j] = Utils.sub(result[j], m);
            }
        }
		return result;
    }
}
