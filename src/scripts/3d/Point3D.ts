import {Pointable} from "../polymorphism/Pointable";

export class Point3D implements Pointable {
    public x: number;
    public y: number;
    public z: number;

    public constructor(readonly x: number, readonly y: number, readonly z: number) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public toString() {
        return `Point {${this.x}, ${this.y}, ${this.z}}`;
    }
}