import {AbstractVector} from "../polymorphism/AbstractVector";
import {Point3D} from "./Point3D";

export class Vector3D extends AbstractVector<Point3D> {
    public constructor(readonly start: Point3D, readonly end: Point3D) {
        super(start, end);
    }

    public override readonly coordinate: [number, number, number] =
        [this.end.x - this.start.x, this.end.y - this.start.y, this.end.z - this.start.z];

    public override readonly length = Math.sqrt(
        Math.pow(this.coordinate[0], 2) + Math.pow(this.coordinate[1], 2) + Math.pow(this.coordinate[2], 2)
    );
}