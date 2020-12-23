class Proteoform {
  sequence = [];
  residueMasses = [];
  fixedPtms=[];
  fixedPtmMasses = [];
  protVarPtms = [];
  variablePtms = [];
  variablePtmPrefixMasses = [];
  variablePtmSuffixMasses = [];
  unexpectedMassShifts=[];
  unexpectedPrefixMasses = [];
  unexpectedSuffixMasses = [];
  prefixModResidueMasses = [];
  prefixMasses = [];
  suffixMasses = [];
  proteoformMass = 0.0;

  constructor(sequence = "", firstPos, fixedPtms =[], protVarPtms = [], variablePtms =[], unexpectedMassShifts=[]) {
    this.sequence = sequence;
    //console.log(this.sequence);
    this.fixedPtms = fixedPtms;
    //console.log(this.fixedPtms);
    this.protVarPtms = protVarPtms;
    //console.log(this.protVarPtms);
    this.variablePtms = variablePtms;
    //console.log(this.variablePtms);
    this.unexpectedMassShifts = unexpectedMassShifts;
    //console.log(this.unexpectedMassShifts);
    this.compResidueMasses();
    this.compFixedPtmMasses(firstPos);
    this.compVariablePtmMasses(firstPos);
    this.compUnexpectedMasses(firstPos);
    this.compPrefixModResidueMasses();
    this.compPrefixMasses();
    this.compSuffixMasses();
    this.compProteoformMass();
  }

  compResidueMasses() {
    this.residueMasses = new Array(this.sequence.length);
    for (let i = 0; i < this.sequence.length; i++) {
      let isotopes = getAminoAcidIsotopes(this.sequence[i]);
      this.residueMasses[i] = isotopes[0].mass;
    }
    //console.log(this.residueMasses);
  }

  compFixedPtmMasses(firstPos) {
    this.fixedPtmMasses = new Array(this.sequence.length).fill(0);
    this.fixedPtms.forEach((element) => {
     //console.log(element);
      for (let i = 0; i < element.posList.length; i++) {
        let pos = element.posList[i].leftPos; 
        //console.log(pos);
        this.fixedPtmMasses[pos-firstPos] = parseFloat(element.mono_mass);
      }
    });
    //console.log(this.fixedPtmMasses);
  }

  compVariablePtmMasses(firstPos) {
    this.variablePtmPrefixMasses = new Array(this.sequence.length).fill(0);
    this.variablePtmSuffixMasses = new Array(this.sequence.length).fill(0);
    this.protVarPtms.forEach((element) => {
      //console.log(element);
      for (let i = 0; i < element.posList.length; i++) {
        let pos = element.posList[i].leftPos; 
        //console.log(pos);
        this.variablePtmPrefixMasses[pos-firstPos] += parseFloat(element.monoMass);
        this.variablePtmSuffixMasses[pos-firstPos] += parseFloat(element.monoMass);
      }
    });
    this.variablePtms.forEach((element) => {
      //console.log(element);
      for (let i = 0; i < element.posList.length; i++) {
        let leftPos = element.posList[i].leftPos; 
        let rightPos = element.posList[i].rightPos; 
        //console.log(pos);
        this.variablePtmPrefixMasses[leftPos-firstPos] += parseFloat(element.monoMass);
        this.variablePtmSuffixMasses[rightPos-firstPos] += parseFloat(element.monoMass);
      }
    });
    //console.log(this.variablePtmPrefixMasses);
    //console.log(this.variablePtmSuffixMasses);
  }

  compUnexpectedMasses(firstPos) {
    this.unexpectedPrefixMasses = new Array(this.sequence.length).fill(0);
    this.unexpectedSuffixMasses = new Array(this.sequence.length).fill(0);
    this.unexpectedMassShifts.forEach((element) => {
      //console.log(element);
      this.unexpectedPrefixMasses[element.leftPos - firstPos] += parseFloat(element.anno);
      this.unexpectedSuffixMasses[element.rightPos - 1 - firstPos] += parseFloat(element.anno);
    });
    //console.log(this.unexpectedMasses);
  }

  compPrefixModResidueMasses() {
    this.prefixModResidueMasses = new Array(this.sequence.length).fill(0);
    for (let i = 0; i < this.sequence.length; i++) {
      let mass = this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmPrefixMasses[i] + this.unexpectedPrefixMasses[i];
      this.prefixModResidueMasses[i] = mass;
    }
  }

  compPrefixMasses() {
    if (this.sequence) {
      let mass = 0;
      for (let i = 0; i < this.sequence.length - 1; i++) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmPrefixMasses[i] + this.unexpectedPrefixMasses[i];
        this.prefixMasses.push(mass);
      }
    } 
  }

  compSuffixMasses() {
    if (this.sequence) {
      let mass = 0;
      for (let i = this.sequence.length - 1; i > 0; i--) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmSuffixMasses[i] + this.unexpectedSuffixMasses[i];
        this.suffixMasses.push(mass);
      }
    }
  }

  compProteoformMass() {
    let mass = 0;
    if (this.sequence) {
      for (let i = 0; i < this.sequence.length; i++) {
        mass = mass + this.prefixModResidueMasses[i];
      }
      mass = mass + getWaterMass(); 
    } 
    this.proteoformMass = mass;
  }

  getNMasses(nIonType) {
    let ionMassShift = getIonMassShift(nIonType);
    //console.log("N mass shift", ionMassShift);
    let massList = []; 
    massList.push(0.0);
    for (let i = 0; i < this.prefixMasses.length; i++) {
      massList.push(this.prefixMasses[i] + ionMassShift);
    }
    massList.push(this.proteoformMass);
    //console.log("massList", massList);
    return massList;
  }

  getCMasses(cIonType) {
    let ionMassShift = getIonMassShift(cIonType);
    //console.log("C mass shift", ionMassShift);
    let massList = []; 
    massList.push(0.0);
    for (let i = 0; i < this.suffixMasses.length; i++) {
      massList.push(this.suffixMasses[i] + ionMassShift);
    }
    massList.push(this.proteoformMass);
    return massList;
  }

  getFixedPtmList() {
    let fixedPtmList = [];
    if (this.fixedPtms.length > 0) {
      //console.log(this.fixedPtms[0]);
      fixedPtmList.push({"name":this.fixedPtms[0].name});
    }
    return fixedPtmList;
  }

  /**
   * Get all the unknwon mass lists form the prsm
   */
  getUnknownMassList() {
    let unknownMassShiftList = [];
    for (let i = 0; i < this.sequence.length; i++) {
      let mass = this.variablePtmPrefixMasses[i] + this.unexpectedPrefixMasses[i];
      if (mass != 0.0) {
        unknownMassShiftList.push({"position":i,"mass":mass})
      }
    }
    return unknownMassShiftList;
  }

}
