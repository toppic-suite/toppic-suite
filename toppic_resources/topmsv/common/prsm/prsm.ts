class Prsm {
  private matchedPeaks_: Peak[] = [];
  private matchedIons_: Ion[] = [];
  private matchedPeakEnvelopePair_: MatchedPeakEnvelopePair[] = [];
  private id_: string;
  private proteoform_: Proteoform;
  private ms1Spec_: Spectrum | null;
  private ms2Spec_: Spectrum[] | null;
  private breakPoints_: BreakPoints[];
  private eValue_: number;
  private qValue_: number;
  private fileName_: string;
  private featureInte_: number | undefined;
  private precMass_: number | undefined;
  private fragIonCount_: number | undefined;

  constructor(id: string, proteoform: Proteoform, ms1Spec: Spectrum | null, ms2Spec: Spectrum[] | null, 
    breakPoints: BreakPoints[], matchedPeakEnvelopePair: MatchedPeakEnvelopePair[] = [],  fileName: string = "", eValue: number = -1, qValue: number = -1,
    featureInte?: number, precMass?: number, fragIonCount?: number) {
    this.id_ = id;
    this.proteoform_ = proteoform;
    this.ms1Spec_ = ms1Spec;
    this.ms2Spec_ = ms2Spec;
    this.breakPoints_ = breakPoints;
    this.matchedPeakEnvelopePair_ = matchedPeakEnvelopePair;
    this.fileName_ = fileName;
    this.eValue_ = eValue;
    this.qValue_ = qValue;
    this.fragIonCount_ = fragIonCount;
    if (featureInte) {
      this.featureInte_ = featureInte;
    }
    if (precMass) {
        this.precMass_ = precMass;
    }
  }
  getId(): string {
    return this.id_;
  }
  getFragIonCount(): number | undefined {
    return this.fragIonCount_;
  }
  getfileName(): string {
    return this.fileName_;
  }
  getProteoform(): Proteoform {
    return this.proteoform_;
  }
  getMs1Spectra(): Spectrum | null{
    return this.ms1Spec_;
  }
  getMs2Spectra(): Spectrum[] | null{
    return this.ms2Spec_;
  }
  getMatchedPeakCount(): number {
    return this.getMatchedPeakEnvelopePairs().length;
  }
  getUnexpectedModCount(): number {
    let protObj: Proteoform = this.getProteoform();
    let unexpectedMod: MassShift[] = protObj.getUnknownMassShift();
    return unexpectedMod.length;
  }
  getEValue(): number {
    return this.eValue_;
  }
  getQValue(): number {
    return this.qValue_;
  }
  getFeatureInte(): number | undefined{
    return this.featureInte_;
  }
  getPrecMass(): number | undefined {
    return this.precMass_;
  }
  getBreakPoints(): BreakPoints[] {
    return this.breakPoints_;
  }
  getMatchedPeakEnvelopePairs(): MatchedPeakEnvelopePair[] {
    return this.matchedPeakEnvelopePair_;
  }
  setBreakPoints(breakPoints: BreakPoints[]): void {
    this.breakPoints_ = breakPoints;
  }
  setProteoform(proteoform: Proteoform): void {
    this.proteoform_ = proteoform;
  }
}