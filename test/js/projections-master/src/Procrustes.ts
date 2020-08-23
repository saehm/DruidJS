/**
 * Performs a Procrustes analysis and transforms an input set of points A to optimally match a reference
 * set of points B under pairwise Euclidean distances.
 * 
 * @author Daniel Karl, IBM Research
 * @version 06/18/2018
 */
export class Procrustes
{
	private static toPoint2D = (A:number[][], mirror:boolean):Point2D[] => {
		let points:Point2D[] = [];
		A.forEach((a, index) => {
			if(mirror) {
				points.push(new Point2D(a[0], -a[1]));
			} else {
				points.push(new Point2D(a[0], a[1]));
			}
		});
		return points;
	};

	private static toArray = (A:Point2D[]):number[][] => {
		let points:number[][] = [];
		A.forEach((a, index) => {
			points.push([a.x, a.y]);
		});
		return points;
	};

    /**
     * Compute the mean point among a set of points A.
     *
     * @param {Point2D[]} A the point set
     * @returns {Point2D} mean point
     */
	private static mean = (A:Point2D[]):Point2D => {
	    let mean:Point2D = new Point2D(0,0);
	    A.forEach((a, index) => {
	        mean.x += a.x;
	        mean.y += a.y;
        });
	    mean.x /= A.length;
	    mean.y /= A.length;
	    return mean;
    };

    /**
     * Center a point set A such that the centroid of A results at (0, 0).
     * This method returns the average distance from a point to the centroid.
     *
     * @param {Point2D[]} A the point set to be centered
     * @returns {number} average distance to centroid
     */
    private static center = (A:Point2D[]):number => {
	    let mean:Point2D = Procrustes.mean(A);
	    let scale:number = 0;
	    A.forEach((a, index) => {
	        a.x -=  mean.x;
	        a.y -=  mean.y;
	        scale += a.x * a.x;
	        scale += a.y * a.y;
        });
	    return Math.sqrt(scale / A.length);
    };

    /**
     * Scale coordinates of point set A.
     *
     * @param {Point2D[]} A point set to scale
     * @param {number} scale scaling factor
     */
    private static scale = (A:Point2D[], scale:number):void => {
        A.forEach((a, index) => {
            a.x /= scale;
            a.y /= scale;
        });
    };

    /**
     * Rotates input shape A according to reference shape B. Optionally force rotation angle (in radians).
     *
     * @param {Point2D[]} A input shape A
     * @param {Point2D[]} B reference shape B
     * @param {number} theta Omit to determine optimal rotation, or set to force angle.
     * @returns {number} rotated A
     */
    private static rotate = (A:Point2D[], B:Point2D[], theta?:number):number => {
        if(theta == null) {
            let n: number = 0, d: number = 0;
            for (let i = 0; i < A.length; i++) {
                n += A[i].x * B[i].y - A[i].y * B[i].x;
                d += A[i].x * B[i].x + A[i].y * B[i].y;
            }
            theta = Math.atan(n / d);
        }
        let x = 0, y = 0;
        for(let i = 0; i < A.length; i++) {
            x = Math.cos(theta) * A[i].x - Math.sin(theta) * A[i].y;
            y = Math.sin(theta) * A[i].x + Math.cos(theta) * A[i].y;
            A[i].x = x;
            A[i].y = y;
        }
        return theta;
    };

    /**
     * Performs a Procrustes analysis and transforms an input set of points A to optimally match a reference
     * set of points B under pairwise Euclidean distances.
     *
     * Careful: This method modifies B.
     *
     * @param {number[][]} A the input shape A
     * @param {number[][]} B the reference shape B (will be centered and scaled during this operation)
     * @returns {number[][]} A centered, scaled, and optimally rotated towards B
     */
    public static transform = (A:number[][], B:number[][]):number[][] => {
	    let A_p = Procrustes.toPoint2D(A, false);
	    let A_p_m = Procrustes.toPoint2D(A, true);
	    let B_p = Procrustes.toPoint2D(B, false);
	    Procrustes.scale(A_p, Procrustes.center(A_p));
	    Procrustes.scale(B_p, Procrustes.center(B_p));
	    Procrustes.scale(A_p_m, Procrustes.center(A_p_m));
	    let A_p_r = Procrustes.copy(A_p);
	    let theta = Procrustes.rotate(A_p, B_p);
	    Procrustes.rotate(A_p_r, B_p, theta - Math.PI);
	    let A_p_m_r = Procrustes.copy(A_p_m);
	    theta = Procrustes.rotate(A_p_m, B_p);
	    Procrustes.rotate(A_p_m_r, B_p, Math.PI + theta);
	    let d:number = Procrustes.distance(A_p, B_p);
	    let d_r:number = Procrustes.distance(A_p_r, B_p);
	    let d_m:number = Procrustes.distance(A_p_m, B_p);
	    let d_m_r:number = Procrustes.distance(A_p_m_r, B_p);
	    let min = Math.min(d, Math.min(d_r, Math.min(d_m, d_m_r)));
	    if(d == min) {
			return Procrustes.toArray(A_p);
		}
		else if(d_r == min) {
			return Procrustes.toArray(A_p_r);
		}
		else if(d_m == min) {
	        return Procrustes.toArray(A_p_m);
        }
        else if(d_m_r == min) {
            return Procrustes.toArray(A_p_m_r);
        }
        else {
	        return null;
        }
    };

    /**
     * Copies the input array.
     *
     * @param {Point2D[]} A input array
     * @returns {Point2D[]} copy of A
     */
    private static copy = (A:Point2D[]):Point2D[] => {
	    let copy:Point2D[] = [];
	    for(let a of A) {
	        copy.push(new Point2D(a.x, a.y));
        }
        return copy;
    };

	/**
	 * Computes the procrustes distance between point sets A and B. The lists
	 * must have equal length.
	 * 
	 * @param A
	 *            list of points
	 * @param B
	 *            list of points
	 * @return the least squared error between A and B
	 */
	private static distance = (A:Point2D[], B:Point2D[]):number => {
		let k = A.length;
		let e = 0;
		for (let i = 0; i < k; i++)
		{
			e += Math.sqrt(Math.pow(A[i].x - B[i].x, 2) + Math.pow(A[i].y - B[i].y, 2));
		}
		return e;
	};
}

class Point2D {
	constructor(public x:number, public y:number){}
}
