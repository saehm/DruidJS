import { euclidean } from "../metrics/index";

// https://www.researchgate.net/profile/Hannu_Huhdanpaa/publication/242414606_The_Quickhull_Algorithm_for_Convex_Hulls/links/54e5f6970cf2bff5a4f1e7db.pdf
export default function(X, metric = euclidean) {
    const [n, d] = X.shape;
    let hull = [];

    return 
}

function _find_hull(X) {

}

/*

Input = a set S of n points 
Assume that there are at least 2 points in the input set S of points
QuickHull (S) 
{ 
    // Find convex hull from the set S of n points
    Convex Hull := {} 
    Find left and right most points, say A & B, and add A & B to convex hull 
    Segment AB divides the remaining (n-2) points into 2 groups S1 and S2 
        where S1 are points in S that are on the right side of the oriented line from A to B, 
        and S2 are points in S that are on the right side of the oriented line from B to A 
    FindHull (S1, A, B) 
    FindHull (S2, B, A) 
}
FindHull (Sk, P, Q) 
{ 
    // Find points on convex hull from the set Sk of points 
    // that are on the right side of the oriented line from P to Q
    If Sk has no point, then return. 
    From the given set of points in Sk, find farthest point, say C, from segment PQ 
    Add point C to convex hull at the location between P and Q 
    Three points P, Q, and C partition the remaining points of Sk into 3 subsets: S0, S1, and S2 
        where S0 are points inside triangle PCQ, S1 are points on the right side of the oriented 
        line from  P to C, and S2 are points on the right side of the oriented line from C to Q. 
    FindHull(S1, P, C) 
    FindHull(S2, C, Q) 
}
Output = Convex Hull

*/