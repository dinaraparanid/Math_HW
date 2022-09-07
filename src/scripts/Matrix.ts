export class Matrix {
    private readonly table: Array<Array<number>>;

    public constructor(readonly _table: Array<Array<number>>) {
        this.table = _table;
    }

    public static sum(mtx1: Matrix, mtx2: Matrix): Matrix {
        if (mtx1.table.length !== mtx2.table.length) throw new Error(
            `Height of matrix 1 (${mtx1.table.length}) != height of matrix 2 (${mtx2.table.length})`
        );

        const sum: number[][] = new Array<Array<number>>();

        for (let i = 0; i < mtx1.table.length; i++) {
            const width = mtx1.table[i].length;

            if (width !== mtx2.table[i].length) throw new Error(
                `Width of matrix 1 ($width) != width of matrix 2 (${mtx2.table[i].length}) in column ${i}`
            );

            const row = new Array<number>(width);

            for (let q = 0; q < width; q++)
                row[q] = mtx1.table[i][q] +  mtx2.table[i][q];

            sum.push(row);
        }

        return new Matrix(sum);
    }

    public isSquare() {
        const firstRowLen = this.table.length;
        return this.table.every(row => row.length === firstRowLen);
    }

    public trace(): number | null {
        if (!this.isSquare())
            return null;

        return Array.from(Array(this.table.length).keys())
            .map(ind => this.table[ind][ind])
            .reduce((a, b) => a + b);
    }

    public toString() {
        return `Matrix {${this.table.toString()}}`;
    }
}