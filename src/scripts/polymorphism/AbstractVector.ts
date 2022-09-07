import {Pointable} from "./Pointable";

export abstract class AbstractVector<P extends Pointable> {
    public abstract readonly start: P;
    public abstract readonly end: P;

    public abstract readonly coordinate: Array<number>;
    public abstract readonly length: number;

    public constructor(readonly start: P, readonly end: P) {
        this.start = start;
        this.end = end;
    }
}