"use strict";
class Proteoform {
    constructor(id, protName, protDesc = "", seq, seqId = "", firstPos, lastPos, protMass, massShiftList, fixedPtm, protVarPtm = [], varPtm = []) {
        this.massShiftList_ = [];
        this.id_ = id;
        this.protName_ = protName;
        this.protDesc_ = protDesc;
        this.seq_ = seq;
        this.seqId_ = seqId;
        this.firstPos_ = firstPos;
        this.lastPos_ = lastPos;
        this.protMass_ = protMass; //rename to prot_mass
        if (massShiftList) {
            this.massShiftList_ = massShiftList.concat(fixedPtm, protVarPtm, varPtm);
        }
        this.prefixMasses_ = [];
        this.suffixMasses_ = [];
        this.compPrefixSuffixMasses();
    }
    /**getters */
    getId() {
        return this.id_;
    }
    getProtName() {
        return this.protName_;
    }
    getProtDesc() {
        return this.protDesc_;
    }
    getSeq() {
        return this.seq_;
    }
    getSeqId() {
        return this.seqId_;
    }
    getFirstPos() {
        return this.firstPos_;
    }
    getLastPos() {
        return this.lastPos_;
    }
    getMass() {
        return this.protMass_;
    }
    getProtVarPtm() {
        let protVarPtm = [];
        this.massShiftList_.forEach(shift => {
            if (shift.getType() == ModType.ProteinVariable) {
                protVarPtm.push(shift);
            }
        });
        return protVarPtm;
    }
    getVarPtm() {
        let varPtm = [];
        this.massShiftList_.forEach(shift => {
            if (shift.getType() == ModType.Variable) {
                varPtm.push(shift);
            }
        });
        return varPtm;
    }
    getUnknownMassShift() {
        let massShiftList = [];
        this.massShiftList_.forEach(shift => {
            if (shift.getType() == ModType.Unexpected) {
                massShiftList.push(shift);
            }
        });
        return massShiftList;
    }
    getUnknownMassShiftAndVarPtm() {
        let newSeq = this.getSeqForMassGraph();
        let massShiftList = [];
        let fixedPtmMasses = [];
        let variablePtmPrefixMasses = [];
        let variablePtmSuffixMasses = [];
        let unexpectedPrefixMasses = [];
        let unexpectedSuffixMasses = [];
        [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses] = this.compMassShiftMasses();
        for (let i = 0; i < newSeq.length; i++) {
            let mass = variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
            if (mass != 0.0) {
                let massShift = new MassShift(i, i + 1, mass, "unexpected", mass.toFixed(4));
                massShiftList.push(massShift);
            }
        }
        return massShiftList;
    }
    getFixedPtm() {
        let fixedPtm = [];
        this.massShiftList_.forEach(shift => {
            if (shift.getType() == ModType.Fixed) {
                fixedPtm.push(shift);
            }
        });
        return fixedPtm;
    }
    getAllShift() {
        return this.massShiftList_;
    }
    getSeqForMassGraph() {
        //unlike the sequence in prsm graph, remove all residues before/after [ and ]
        let newSeq = "";
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
    compPrefixSuffixMasses() {
        let fixedPtmMasses = [];
        let variablePtmPrefixMasses = [];
        let variablePtmSuffixMasses = [];
        let unexpectedPrefixMasses = [];
        let unexpectedSuffixMasses = [];
        let aminoAcidMasses = this.compAminoAcidMasses();
        [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses] = this.compMassShiftMasses();
        let prefixModResidueMasses = this.compPrefixModResidueMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses);
        this.compPrefixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses);
        this.compSuffixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmSuffixMasses, unexpectedSuffixMasses);
        if (this.protMass_ < 0) { //if it is inspect page
            this.compProteoformMass(prefixModResidueMasses);
        }
    }
    compPrefixModResidueMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses) {
        let prefixModResidueMasses = new Array(this.seq_.length).fill(0);
        for (let i = 0; i < this.seq_.length; i++) {
            let mass = aminoAcidMasses[i] + fixedPtmMasses[i]
                + variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
            prefixModResidueMasses[i] = mass;
        }
        return prefixModResidueMasses;
    }
    compProteoformMass(prefixModResidueMasses) {
        let newSeq = this.getSeqForMassGraph();
        let mass = 0;
        if (newSeq) {
            for (let i = 0; i < newSeq.length; i++) {
                mass = mass + prefixModResidueMasses[i];
            }
            mass = mass + getWaterMass();
        }
        this.protMass_ = mass;
    }
    compAminoAcidMasses() {
        let newSeq = this.getSeqForMassGraph();
        let residueMasses = new Array(newSeq.length);
        for (let i = 0; i < newSeq.length; i++) {
            let isotopes = getAminoAcidDistribution(newSeq[i]);
            if (isotopes) {
                residueMasses[i] = isotopes[0].mass;
            }
        }
        return residueMasses;
    }
    compMassShiftMasses() {
        let newSeq = this.getSeqForMassGraph();
        let fixedPtmMasses = new Array(newSeq.length).fill(0);
        let variablePtmPrefixMasses = new Array(newSeq.length).fill(0);
        let variablePtmSuffixMasses = new Array(newSeq.length).fill(0);
        let unexpectedPrefixMasses = new Array(newSeq.length).fill(0);
        let unexpectedSuffixMasses = new Array(newSeq.length).fill(0);
        this.massShiftList_.forEach(massShift => {
            if (massShift.getType() == ModType.Fixed) {
                let pos = massShift.getLeftPos();
                fixedPtmMasses[pos - this.firstPos_] = massShift.getShift();
            }
            else if (massShift.getType() == ModType.ProteinVariable) {
                let pos = massShift.getLeftPos();
                variablePtmPrefixMasses[pos - this.firstPos_] += massShift.getShift();
                variablePtmSuffixMasses[pos - this.firstPos_ - 1] += massShift.getShift();
            }
            else if (massShift.getType() == ModType.Variable) {
                let leftPos = massShift.getLeftPos();
                variablePtmPrefixMasses[leftPos - this.firstPos_] += massShift.getShift();
                variablePtmSuffixMasses[leftPos - this.firstPos_] += massShift.getShift();
            }
            else {
                unexpectedPrefixMasses[massShift.getLeftPos() - this.firstPos_] += massShift.getShift();
                unexpectedSuffixMasses[massShift.getLeftPos() - this.firstPos_] += massShift.getShift();
            }
        });
        return [fixedPtmMasses, variablePtmPrefixMasses, variablePtmSuffixMasses, unexpectedPrefixMasses, unexpectedSuffixMasses];
    }
    compPrefixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmPrefixMasses, unexpectedPrefixMasses) {
        if (this.seq_) {
            //sequence for mass graph
            let newSeq = this.getSeqForMassGraph();
            let mass = 0;
            for (let i = 0; i < newSeq.length - 1; i++) {
                mass = mass + aminoAcidMasses[i] + fixedPtmMasses[i]
                    + variablePtmPrefixMasses[i] + unexpectedPrefixMasses[i];
                let theoMass = new TheoMass(mass, i);
                this.prefixMasses_.push(theoMass);
            }
        }
    }
    compSuffixMasses(aminoAcidMasses, fixedPtmMasses, variablePtmSuffixMasses, unexpectedSuffixMasses) {
        if (this.seq_) {
            let newSeq = this.getSeqForMassGraph();
            let mass = 0;
            for (let i = newSeq.length - 1; i > 0; i--) {
                mass = mass + aminoAcidMasses[i] + fixedPtmMasses[i]
                    + variablePtmSuffixMasses[i] + unexpectedSuffixMasses[i];
                let theoMass = new TheoMass(mass, i);
                this.suffixMasses_.push(theoMass);
            }
        }
    }
    /**other functions that return calculated values */
    getNMasses(nIonType) {
        let ionMassShift = getIonMassShift(nIonType);
        let massList = [];
        massList.push(new TheoMass(0, -1));
        if (ionMassShift === undefined) {
            return massList;
        }
        this.prefixMasses_.forEach(theoMass => {
            let newMass = theoMass.getMass() + ionMassShift; //ionMassShift is already checked for undefined
            let pos = theoMass.getPos();
            let ionMassAdjusted = new TheoMass(newMass, pos);
            massList.push(ionMassAdjusted);
        });
        massList.push(new TheoMass(this.protMass_, this.prefixMasses_.length));
        return massList;
    }
    getCMasses(cIonType) {
        let ionMassShift = getIonMassShift(cIonType);
        //console.log("C mass shift", ionMassShift);
        let massList = [];
        massList.push(new TheoMass(0, -1));
        if (ionMassShift === undefined) {
            return massList;
        }
        this.suffixMasses_.forEach(theoMass => {
            let newMass = theoMass.getMass() + ionMassShift; //ionMassShift is already checked for undefined
            let pos = theoMass.getPos();
            let ionMassAdjusted = new TheoMass(newMass, pos);
            massList.push(ionMassAdjusted);
        });
        massList.push(new TheoMass(this.protMass_, this.suffixMasses_.length));
        return massList;
    }
}
