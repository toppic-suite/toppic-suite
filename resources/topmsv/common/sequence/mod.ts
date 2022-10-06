class Mod {
  private residue_: string;
  private massShift_: number;
  private name_:string;

  constructor(residue: string, massShift: number, name: string) {
    this.residue_ = residue;
    this.massShift_ = massShift;
    this.name_ = name;
  }
  getResidue(): string {
    return this.residue_;
  }
  getShift(): number {
    return this.massShift_;
  }
  getName(): string {
    return this.name_;
  }
}