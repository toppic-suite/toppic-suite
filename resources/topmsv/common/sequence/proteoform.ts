class Proteoform {
  private massShiftList_: MassShift[] = [] as MassShift[];
  private prefixMasses_: TheoMass[];
  private suffixMasses_: TheoMass[];
  private id_: string;
  private protName_: string;
  private protDesc_: string;
  private seq_: string;
  private seqId_: string;
  private firstPos_: number;
  private lastPos_: number;
  private protMass_: number;

  constructor (id: string, protName: string, protDesc: string = "", seq: string, seqId: string = "", 
    firstPos: number, lastPos: number, 
    protMass: number, massShiftList: MassShift[], fixedPtm: MassShift[], 
    protVarPtm: MassShift[] = [], varPtm: MassShift[] = []) {
    this.id_ = id;
    this.protName_ = protName;
    this.protDesc_ = protDesc;
    this.seq_ = seq;
    this.seqId_ = seqId;
    this.firstPos_ = firstPos;
    this.lastPos_ = lastPos;
    this.protMass_ = protMass;//rename to prot_mass
    if (massShiftList) {
      this.massShiftList_ = massShiftList.concat(fixedPtm, protVarPtm, varPtm);
    }
    this.prefixMasses_ = [];
    this.suffixMasses_ = [];
    this.compPrefixSuffixMasses();
  }
  /**getters */
  getId(): string {
    return this.id_;
  }
  getProtName(): string {
    return this.protName_;
  }
  getProtDesc(): string {
    return this.protDesc_;
  }
  getSeq(): string {
    return this.seq_;
  }
  getSeqId(): string {
    return this.seqId_;
  }
  getFirstPos(): number {
    return this.firstPos_;
  }
  getLastPos(): number {
    return this.lastPos_;
  }
  getMass(): number {
    return this.protMass_;
  }
  getProtVarPtm(): MassShift[] {
    let protVarPtm: MassShift[] = [];
    this.massShiftList_.forEach(shift => {
      if (shift.getType() == ModType.ProteinVariable) {
        protVarPtm.push(shift);
      }
    })
    return protVarPtm;
  }
  getVarPtm(): MassShift[] {
    let varPtm: MassShift[] = [];
    this.massShiftList_.forEach(shift => {
      if (shift.getType() == ModType.Variable) {
        varPtm.push(shift);
      }
    })
    return varPtm;
  }
  getUnknownMassShift(): MassShift[] {
    let massShiftList: MassShift[] = [];
    this.massShiftList_.forEach(shift => {
      if (shift.getType() == ModType.Unexpected) {
        massShiftList.push(shift);
      }
    })
    return massShiftList;
  }
  getUnknownMassShiftAndVarPtm(): MassShift[] {//return unknown mass shift and variable ptm in one array 
    let newSeq: string = this.getSeqForMassGraph();
    let massShiftList: MassShift[] = [];
    let fixedPtmMasses: number[] = [];
    let variablePtmPrefixMasses: number[] = [];
    let variablePtmSuffixMasses: number[] = [];
    let unexpectedPrefixMasses: number[] = [];
    let unexpectedSuffixMasses: number[] = [];
    [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses] = this.compMassShiftMasses();
    
    for (let i = 0; i < newSeq.length; i++) {
      let mass = variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
      if (mass != 0.0) {
        let massShift = new MassShift(i, i+1, mass, "unexpected", mass.toFixed(4));
        massShiftList.push(massShift);
      }
    }
    return massShiftList;
  }
  getFixedPtm(): MassShift[] {
    let fixedPtm: MassShift[] = [];
    this.massShiftList_.forEach(shift => {
      if (shift.getType() == ModType.Fixed) {
        fixedPtm.push(shift);
      }
    })
    return fixedPtm;
  }
  getAllShift(): MassShift[] {
    return this.massShiftList_;
  }
  getSeqForMassGraph(): string {
    //unlike the sequence in prsm graph, remove all residues before/after [ and ]
    let newSeq: string = "";
    for (let j = 0; j < this.seq_.length; j++) {
      if (j > this.lastPos_) {
        return newSeq;
      }
      if (j >= this.firstPos_) {
        newSeq = newSeq + this.seq_[j];
      }
    }
    return newSeq;
  }
  /**compute masses and set to the private members of class*/
  compPrefixSuffixMasses(): void {
    let fixedPtmMasses: number[] = [];
    let variablePtmPrefixMasses: number[] = [];
    let variablePtmSuffixMasses: number[] = [];
    let unexpectedPrefixMasses: number[] = [];
    let unexpectedSuffixMasses: number[] = [];
    let aminoAcidMasses = this.compAminoAcidMasses();
    [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses] = this.compMassShiftMasses();
    let prefixModResidueMasses: number[] = this.compPrefixModResidueMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses);
    this.compPrefixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses);
    this.compSuffixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmSuffixMasses, unexpectedSuffixMasses);
    if (this.protMass_ < 0) { //if it is inspect page
        this.compProteoformMass(prefixModResidueMasses);
    }
  }

  compPrefixModResidueMasses(aminoAcidMasses: number[], fixedPtmMasses: number[], 
    variablePtmPrefixMasses: number[], unexpectedPrefixMasses: number[]): number[] {//compute amino acid mass in prefix
    
    let prefixModResidueMasses: number[] = new Array(this.seq_.length).fill(0);
    for (let i = 0; i < this.seq_.length; i++) {
      let mass: number = aminoAcidMasses[i] + fixedPtmMasses[i]  
        + variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
        prefixModResidueMasses[i] = mass;
    }
    return prefixModResidueMasses;
  }
  compProteoformMass(prefixModResidueMasses: number[]) {//compute proteoform mass
    let newSeq: string = this.getSeqForMassGraph();
    let mass: number = 0;
    if (newSeq) {
      for (let i = 0; i < newSeq.length; i++) {
        mass = mass + prefixModResidueMasses[i];
      }
      mass = mass + getWaterMass(); 
    } 
    this.protMass_ = mass;
  }
  compAminoAcidMasses(): number[] {//compute mass of each amino acid
    let newSeq: string = this.getSeqForMassGraph();
    let residueMasses: number[] = new Array(newSeq.length);
    for (let i = 0; i < newSeq.length; i++) {
      let isotopes = getAminoAcidDistribution(newSeq[i]);
      if (isotopes){
        residueMasses[i] = isotopes[0].mass;
      }
    }
    return residueMasses;
  }
  compMassShiftMasses(): any[] {//compute proteoform mass with ptm shifts added
    let newSeq: string = this.getSeqForMassGraph();
    let fixedPtmMasses: number[] = new Array(newSeq.length).fill(0);
    let variablePtmPrefixMasses: number[] = new Array(newSeq.length).fill(0);
    let variablePtmSuffixMasses: number[] = new Array(newSeq.length).fill(0);
    let unexpectedPrefixMasses: number[] = new Array(newSeq.length).fill(0);
    let unexpectedSuffixMasses: number[] = new Array(newSeq.length).fill(0);
    this.massShiftList_.forEach(massShift => {
      if (massShift.getType() == ModType.Fixed) {
        let pos: number = massShift.getLeftPos();
        fixedPtmMasses[pos-this.firstPos_] = massShift.getShift();
      }
      else if(massShift.getType() == ModType.ProteinVariable) {
        let pos: number = massShift.getLeftPos();
        variablePtmPrefixMasses[pos-this.firstPos_] += massShift.getShift();
        variablePtmSuffixMasses[pos-this.firstPos_ - 1] += massShift.getShift();
      }
      else if(massShift.getType() == ModType.Variable) {
        let leftPos: number = massShift.getLeftPos(); 
        variablePtmPrefixMasses[leftPos-this.firstPos_] += massShift.getShift();
        variablePtmSuffixMasses[leftPos-this.firstPos_] += massShift.getShift();
      }
      else{
        unexpectedPrefixMasses[massShift.getLeftPos() - this.firstPos_] += massShift.getShift();
        unexpectedSuffixMasses[massShift.getLeftPos() - this.firstPos_] += massShift.getShift();  
      }
    })
    return [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses];
  }

  compPrefixMasses(aminoAcidMasses: number[], fixedPtmMasses: number[], 
    variablePtmPrefixMasses: number[], unexpectedPrefixMasses: number[]): void {//compute theoratical prefix mass
      if (this.seq_) {
      //sequence for mass graph
      let newSeq: string = this.getSeqForMassGraph();
      let mass: number = 0; 
      for (let i = 0; i < newSeq.length - 1; i++) {
        mass = mass + aminoAcidMasses[i] + fixedPtmMasses[i]  
        + variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
        let theoMass = new TheoMass(mass, i);
        this.prefixMasses_.push(theoMass);
      }
    } 
  }

  compSuffixMasses(aminoAcidMasses: number[], fixedPtmMasses: number[], variablePtmSuffixMasses: number[], 
    unexpectedSuffixMasses: number[]) {//compute theoratical suffix mass
    if (this.seq_) {
      let newSeq: string = this.getSeqForMassGraph();
      let mass: number = 0;
      for (let i = newSeq.length - 1; i > 0; i--) {
        mass = mass + aminoAcidMasses[i] + fixedPtmMasses[i]  
        + variablePtmSuffixMasses[i] + unexpectedSuffixMasses[i];
        let theoMass = new TheoMass(mass, i);
        this.suffixMasses_.push(theoMass);
      }
    }
  }
  /**other functions that return calculated values */
  getNMasses(nIonType: string): TheoMass[] {//theoretical mass with n-term ion shift added
    let ionMassShift: number | undefined = getIonMassShift(nIonType);
    let massList: TheoMass[] = [];
    massList.push(new TheoMass(0, -1));

    if (ionMassShift === undefined) {
      return massList;
    }
    this.prefixMasses_.forEach(theoMass => {
      let newMass: number = theoMass.getMass() + ionMassShift!;//ionMassShift is already checked for undefined
      let pos: number = theoMass.getPos();
      let ionMassAdjusted: TheoMass = new TheoMass(newMass, pos);
      
      massList.push(ionMassAdjusted);
    })
    massList.push(new TheoMass(this.protMass_, this.prefixMasses_.length));
    return massList;
  }

  getCMasses(cIonType: string): TheoMass[] {//theoretical mass with c-term ion shift added
    let ionMassShift:  number | undefined = getIonMassShift(cIonType);
    //console.log("C mass shift", ionMassShift);
    let massList: TheoMass[] = []; 
    massList.push(new TheoMass(0, -1));

    if (ionMassShift === undefined) {
      return massList;
    }
    this.suffixMasses_.forEach(theoMass => {
      let newMass: number = theoMass.getMass() + ionMassShift!;//ionMassShift is already checked for undefined
      let pos: number = theoMass.getPos();
      let ionMassAdjusted: TheoMass = new TheoMass(newMass, pos);
      
      massList.push(ionMassAdjusted);
    })
    massList.push(new TheoMass(this.protMass_, this.suffixMasses_.length));
    return massList;
  }

}