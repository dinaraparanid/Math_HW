import {Pointable} from "../polymorphism/Pointable";

export class Point2D implements Pointable {
    public x: number;
    public y: number;

    public constructor(readonly x: number, readonly y: number) {
        this.x = x;
        this.y = y;
    }

    public toString() {
        return `Point2D {${this.x}, ${this.y}}`;
    }
}