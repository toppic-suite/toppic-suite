class TheoMass {
  private mass_: number;
  private pos_: number;
  private ion_ = {} as Ion;

  constructor(mass: number, pos: number) {
    this.mass_ = mass;
    this.pos_ = pos;
  }
  getMass(): number {
    return this.mass_;
  }
  setMass(mass: number) {
    this.mass_ = mass;
  }
  getPos(): number {
    return this.pos_;
  }
  getIon(): Ion {
    return this.ion_;
  }
  setIon(ion: Ion): void {
    this.ion_ = ion;
  }
}