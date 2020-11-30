class Proteoform {
  sequence = [];
  residueMasses = [];
  fixedPtms=[];
  fixedPtmMasses = [];
  variablePtms = [];
  variablePtmMasses = [];
  unexpectedMassShifts=[];
  unexpectedPrefixMasses = [];
  unexpectedSuffixMasses = [];
  prefixModResidueMasses = [];
  prefixMasses = [];
  suffixMasses = [];
  proteoformMass = 0.0;

  constructor(sequence = "", firstPos, fixedPtms =[], variablePtms =[], unexpectedMassShifts=[]) {
    this.sequence = sequence;
    this.fixedPtms = fixedPtms;
    this.variablePtms = variablePtms;
    this.unexpectedMassShifts = unexpectedMassShifts;
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
      for (let i = 0; i < element.posList.length; i++) {
        let pos = element.posList[i].pos; 
        //console.log(pos);
        this.fixedPtmMasses[pos-firstPos] = parseFloat(element.mono_mass);
      }
    });
    //console.log(this.fixedPtmMasses);
  }

  compVariablePtmMasses(firstPos) {
    this.variablePtmMasses = new Array(this.sequence.length).fill(0);
    this.variablePtms.forEach((element) => {
      //console.log(element);
      for (let i = 0; i < element.posList.length; i++) {
        let pos = element.posList[i].pos; 
        //console.log(pos);
        this.variablePtmMasses[pos-firstPos] = parseFloat(element.mono_mass);
      }
    });
    //console.log(this.variablePtmMasses);
  }

  compUnexpectedMasses(firstPos) {
    this.unexpectedPrefixMasses = new Array(this.sequence.length).fill(0);
    this.unexpectedSuffixMasses = new Array(this.sequence.length).fill(0);
    this.unexpectedMassShifts.forEach((element) => {
      console.log(element);
      this.unexpectedPrefixMasses[element.leftPos - firstPos] = parseFloat(element.anno);
      this.unexpectedSuffixMasses[element.rightPos - 1 - firstPos] = parseFloat(element.anno);
    });
    //console.log(this.unexpectedMasses);
  }

  compPrefixModResidueMasses() {
    this.prefixModResidueMasses = new Array(this.sequence.length).fill(0);
    for (let i = 0; i < this.sequence.length; i++) {
      let mass = this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmMasses[i] + this.unexpectedPrefixMasses[i];
      this.prefixModResidueMasses[i] = mass;
    }
  }

  compPrefixMasses() {
    if (this.sequence) {
      let mass = 0;
      for (let i = 0; i < this.sequence.length - 1; i++) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmMasses[i] + this.unexpectedPrefixMasses[i];
        this.prefixMasses.push(mass);
      }
    } 
  }

  compSuffixMasses() {
    if (this.sequence) {
      let mass = 0;
      for (let i = this.sequence.length - 1; i > 0; i--) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i]  
        + this.variablePtmMasses[i] + this.unexpectedSuffixMasses[i];
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
}
