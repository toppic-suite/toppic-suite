class MatchedPeakEnvelopePair{
  private theoMass_: number;
  private peak_: Peak;
  private ion_: Ion;
  private envelope_ : Envelope | null;

  constructor(theoMass: number, monoPeak: Peak, ion: Ion, envelope?: Envelope) {
    this.envelope_ = null;
    this.theoMass_ = theoMass;
    this.peak_ = monoPeak;
    this.ion_ = ion;
    if (envelope) {
      this.envelope_ = envelope; //only one envelope per object
    }  
  }
  getPeak(): Peak {
    return this.peak_;
  }
  getIon(): Ion {
    return this.ion_;
  }
  getEnvelope(): Envelope | null{
    return this.envelope_;
  }
  getTheoMass(): number {
    return this.theoMass_;
  }
  setTheoMass(mass: number): void {
    this.theoMass_ = mass;
  }
  addEnvelope(env: Envelope) {
    this.envelope_ = env;
  }
}