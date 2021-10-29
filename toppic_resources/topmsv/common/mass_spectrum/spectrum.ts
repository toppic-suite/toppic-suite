class Spectrum {
  private id_: string;
  private scanNum_: string;
  private level_: number;
  private peakList_: Peak[];
  private deconvPeakList_ : Peak[] | null;
  private nIon_: Ion[];
  private cIon_: Ion[];
  private precMass_: number;
  private precCharge_: number;
  private precMz_: number;
  private envList_: Envelope[] = [];

  constructor(id: string, scanNum: string, level: number, peakList: Peak[], decovPeakList: Peak[] | null, envList: Envelope[] = [], nIon: Ion[], cIon: Ion[],
    mass: number, charge: number = -1, mz: number = -1) {
    this.id_ = id; 
    this.scanNum_ = scanNum;
    this.level_ = level;
    this.peakList_ = peakList;
    this.deconvPeakList_ = decovPeakList;
    this.envList_ = envList;
    this.nIon_ = nIon;
    this.cIon_ = cIon;
    this.precMass_ = mass;
    this.precCharge_ = charge;
    this.precMz_ = mz;
  }
  getSpectrumId(): string {
    return this.id_;
  }
  getScanNum(): string {
    return this.scanNum_;
  }
  getSpectrumLevel(): number {
    return this.level_;
  }
  getPeaks(): Peak[] {
    return this.peakList_;
  }
  getDeconvPeaks(): Peak[] | null {
    return this.deconvPeakList_;
  }
  getEnvs(): Envelope[] {
    return this.envList_;
  }
  getNTerminalIon(): Ion[] {
    return this.nIon_;
  }
  getCTerminalIon(): Ion[] {
    return this.cIon_;
  }
  getPrecMass(): number {
    return this.precMass_;
  }
  getPrecCharge(): number {
    return this.precCharge_;
  }
  getPrecMz(): number {
    return this.precMz_;
  }
  setPeaks(peaks: Peak[]): void {
    this.peakList_ = peaks;
  }
  setDeconvPeaks(decovPeaks: Peak[]): void {
    this.deconvPeakList_ = decovPeaks;
  }
  setEnvs(envs: Envelope[]): void {
    this.envList_ = envs;
  }
  addNTerminalIon(ion: Ion): void {
    this.nIon_.push(ion);
  }
  addCTerminalIon(ion: Ion): void {
    this.cIon_.push(ion);
  }
}