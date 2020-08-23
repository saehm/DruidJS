/**
 * @author Daniel Karl, IBM Research
 */
import {Utils} from "./Utils";
export class Scaling {
    /**
     * Normalize the coordinates between 0 and 1 preserving the aspect ratios of dimensions.
     * @param {number[][]} coordinates
     */
    static normalize = (coordinates: number[][]): void => {
        let D:number = coordinates[0].length;
        let min: number[] = Utils.fillArray(D, Number.MAX_VALUE);
        let max: number[] = Utils.fillArray(D, Number.MIN_VALUE);
        let x;
        for (let i = 0; i < coordinates.length; i++)
        {
            for (let d = 0; d < D; d++)
            {
                x = coordinates[i][d];
                if (x < min[d])
                {
                    min[d] = x;
                }
                if (x > max[d])
                {
                    max[d] = x;
                }
            }
        }
        let length, longestLength = 0;
        for (let d = 0; d < D; d++)
        {
            length = max[d] - min[d];
            if (length > longestLength)
            {
                longestLength = length;
            }
        }
        for (let d = 0; d < D; d++)
        {
            for (let i = 0; i < coordinates.length; i++)
            {
                x = coordinates[i][d];
                coordinates[i][d] = (x - min[d]) / longestLength;
            }
        }
    }
}
