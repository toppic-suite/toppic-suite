class Proteoform {
  sequence = [];
  residueMasses = [];
  fixedPtmMasses = [];
  fixedPtmMassShift=[];
  unexpectedMassShift=[];
  unexpectedMasses = [];
  prefixMasses = [];
  suffixMasses = [];

  constructor(sequence = "", fixedPtmMassShift =[], unexpectedMassShift=[]) {
    this.sequence = sequence;
    this.fixedPtmMassShift = fixedPtmMassShift;
    this.unexpectedMassShift = unexpectedMassShift;
    this.getResidueMasses();
    this.getFixedPtmMasses();
    this.getUnexpectedMasses();
    let nTermMassShift = getIonMassShift("B");
    this.prefixMasses = this.getPrefixMassList(nTermMassShift);
    let cTermMassShift = getIonMassShift("Y");
    this.suffixMasses = this.getSuffixMassList(cTermMassShift);
  }

  getResidueMasses() {
    this.residueMasses = new Array(this.sequence.length);
    for (let i = 0; i < this.sequence.length; i++) {
      let isotopes = getAminoAcidIsotopes(this.sequence[i]);
      this.residueMasses[i] = isotopes[0].mass;
    }
  }

  getFixedPtmMasses() {
    this.fixedPtmMasses = new Array(this.sequence.length).fill(0);
    this.fixedPtmMassShift.forEach((element) => {
      this.fixedPtmMasses[element.position] = element.mass;
    });
    return this.fixedPtmMasses;
  }

  getUnexpectedMasses() {
    this.unexpectedMasses = new Array(this.sequence.length).fill(0);
    this.unexpectedMassShift.forEach((element) => {
      this.unexpectedMasses[element.position] = element.mass;
    });
    return this.unexpectedMasses;
  }

  getPrefixMassList(ionMassShift) {
    if (this.sequence) {
      let prefixMassList = []; 
      let mass = 0;
      prefixMassList.push(mass);
      for (let i = 0; i < this.sequence.length; i++) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i] 
          + this.unexpectedMasses[i];
        prefixMassList.push(mass+ionMassShift);
      }
      return prefixMassList;
    } else {
      return [];
    }
  }

  getSuffixMassList(ionMassShift) {
    if (this.sequence) {
      let suffixMassList = []; 
      let mass = 0;
      suffixMassList.push(mass);
      for (let i = this.sequence.length - 1; i >= 0; i--) {
        mass = mass + this.residueMasses[i] + this.fixedPtmMasses[i] + this.unexpectedMasses[i];
        suffixMassList.push(mass + ionMassShift);
      }
      return suffixMassList;
    }else {
      return [];
    }  
  }
}
