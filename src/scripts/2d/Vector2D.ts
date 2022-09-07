import {AbstractVector} from "../polymorphism/AbstractVector";
import {Point2D} from "./Point2D";

export class Vector2D extends AbstractVector<Point2D>{
    public constructor(readonly start: Point2D, readonly end: Point2D) {
        super(start, end);
    }

    public override readonly coordinate: [number, number] =
        [this.end.x - this.start.x, this.end.y - this.start.y];

    public override readonly length = Math.sqrt(
        Math.pow(this.coordinate[0], 2) + Math.pow(this.coordinate[1], 2)
    );
}