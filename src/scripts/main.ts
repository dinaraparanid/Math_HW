import {Matrix} from "./Matrix";

const aMatrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
const bMatrix = new Matrix([[20, 19, 18], [17, 16, 15], [14, 13, 12]]);
const sumMatrix = Matrix.sum(aMatrix, bMatrix);

console.assert(sumMatrix.trace() === aMatrix.trace() + bMatrix.trace());