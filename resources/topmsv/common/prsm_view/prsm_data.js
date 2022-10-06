"use strict";
class PrsmViewData {
    constructor() {
        //missing underscore in variable names
        this.residues_ = [];
        this.formFirstPos_ = 0;
        this.formLastPos_ = 0;
        this.breakPoints_ = [];
        this.sequence_ = "";
        this.fixedPtms_ = [];
        this.protVariablePtms_ = [];
        this.variablePtms_ = [];
        this.massShifts_ = [];
        this.annotations_ = [];
        this.proteoform_ = null;
        //information above are retrieved from prsm and proteoform object
        this.rowNum_ = 0;
        this.displayFirstPos_ = 0;
        this.displayLastPos_ = 0;
        // if there is a skipping line at the beginning
        this.showStartSkipped_ = true;
        this.startSkippedInfo_ = "";
        // if there is a skipping line at the end
        this.showEndSkipped_ = true;
        this.endSkippedInfo_ = "";
    }
    //getters and setters
    getResidues() {
        return this.residues_;
    }
    getFormFirstPos() {
        return this.formFirstPos_;
    }
    getFormLastPos() {
        return this.formLastPos_;
    }
    getSequence() {
        return this.sequence_;
    }
    getBreakPoints() {
        return this.breakPoints_;
    }
    getFixedPtms() {
        return this.fixedPtms_;
    }
    ;
    getVariablePtms() {
        return this.variablePtms_;
    }
    getMassShifts() {
        return this.massShifts_;
    }
    getAnnotations() {
        return this.annotations_;
    }
    getProteoform() {
        return this.proteoform_;
    }
    getRowNum() {
        return this.rowNum_;
    }
    getDisplayFirstPos() {
        return this.displayFirstPos_;
    }
    getDisplayLastPos() {
        return this.displayLastPos_;
    }
    getShowStartSkipped() {
        return this.showStartSkipped_;
    }
    getStartSkippedInfo() {
        return this.startSkippedInfo_;
    }
    getShowEndSkipped() {
        return this.showEndSkipped_;
    }
    getEndSkippedInfo() {
        return this.endSkippedInfo_;
    }
    initData(prsm, para) {
        this.parseData(prsm);
        this.updatePara(para);
        this.addColor();
    }
    // form residues from sequence
    formResidues(sequence) {
        let residues = [];
        for (let i = 0; i < sequence.length; i++) {
            let tempObj = {
                position: i,
                acid: sequence.charAt(i).toUpperCase(),
                color: ""
            };
            residues.push(tempObj);
        }
        return residues;
    }
    generateAnnotation(protVarPtm, variablePtm, massShift) {
        let annos = [];
        //annotation in prsm object contains only variable and unknown shifts
        for (let i = 0; i < protVarPtm.length; i++) {
            let temp = { "annoText": protVarPtm[i].getAnnotation(), "leftPos": 0, "rightPos": 0, "type": ModType.ProteinVariable };
            temp.leftPos = protVarPtm[i].getLeftPos();
            temp.rightPos = protVarPtm[i].getRightPos();
            annos.push(temp);
        }
        for (let i = 0; i < variablePtm.length; i++) {
            let temp = { "annoText": variablePtm[i].getAnnotation(), "leftPos": 0, "rightPos": 0, "type": ModType.Variable };
            temp.leftPos = variablePtm[i].getLeftPos();
            temp.rightPos = variablePtm[i].getRightPos();
            annos.push(temp);
        }
        for (let i = 0; i < massShift.length; i++) {
            let temp = { "annoText": FormatUtil.formatFloat(massShift[i].getAnnotation(), "massShift"), "leftPos": massShift[i].getLeftPos(), "rightPos": massShift[i].getRightPos(), "type": ModType.Unexpected };
            annos.push(temp);
        }
        //merge variable ptms with same residue range
        annos.sort(function (a, b) {
            return a.leftPos - b.leftPos;
        });
        let prevAnno = {};
        let prevLeft = -1;
        let prevRight = -1;
        let newAnnos = [];
        annos.forEach((anno, i) => {
            if (anno.leftPos != prevLeft || anno.rightPos != prevRight) {
                if (prevLeft >= 0) {
                    newAnnos.push(prevAnno);
                }
                prevAnno = anno;
            }
            else {
                prevAnno.annoText = prevAnno.annoText + ";" + anno.annoText;
            }
            prevLeft = anno.leftPos;
            prevRight = anno.rightPos;
        });
        newAnnos.push(prevAnno);
        return newAnnos;
    }
    parseData(prsm) {
        let proteoformObj = prsm.getProteoform();
        this.residues_ = this.formResidues(proteoformObj.getSeq());
        this.formFirstPos_ = proteoformObj.getFirstPos();
        this.formLastPos_ = proteoformObj.getLastPos();
        this.breakPoints_ = prsm.getBreakPoints();
        this.fixedPtms_ = proteoformObj.getFixedPtm();
        this.protVariablePtms_ = proteoformObj.getProtVarPtm();
        this.variablePtms_ = proteoformObj.getVarPtm();
        this.massShifts_ = proteoformObj.getUnknownMassShift();
        this.sequence_ = proteoformObj.getSeq();
        this.proteoform_ = proteoformObj;
        this.annotations_ = this.generateAnnotation(this.protVariablePtms_, this.variablePtms_, this.massShifts_);
    }
    updatePara(para) {
        let len = this.residues_.length;
        // Include 5 amino acids before and after the form
        this.displayFirstPos_ = Math.floor((this.formFirstPos_ - 5) / para.getRowLength()) * para.getRowLength();
        if (this.displayFirstPos_ < 0) {
            this.displayFirstPos_ = 0;
        }
        this.displayLastPos_ = Math.ceil((this.formLastPos_ + 6) / para.getRowLength()) * para.getRowLength() - 1;
        //console.log("display last pos ", this.displayLastPos);
        if (this.displayLastPos_ > (len - 1)) {
            this.displayLastPos_ = len - 1;
        }
        this.rowNum_ = Math.ceil((this.displayLastPos_ - this.displayFirstPos_ + 1) / para.getRowLength());
        // skipping line
        this.showStartSkipped_ = false;
        this.showEndSkipped_ = false;
        if (para.getShowSkippedLines()) {
            if (this.displayFirstPos_ !== 0) {
                this.showStartSkipped_ = true;
                this.startSkippedInfo_ = "... " + this.displayFirstPos_
                    + " amino acid residues are skipped at the N-terminus ... ";
            }
            if (this.displayLastPos_ !== len - 1) {
                this.showEndSkipped_ = true;
                this.endSkippedInfo_ = "... " + (len - 1 - this.displayLastPos_)
                    + " amino acid residues are skipped at the C-terminus ... ";
            }
        }
    }
    addColor() {
        for (let i = 0; i < this.residues_.length; i++) {
            let residue = this.residues_[i];
            let pos = residue.position;
            if (pos < this.formFirstPos_ || pos > this.formLastPos_) {
                residue.color = "grey";
            }
            else {
                residue.color = "black";
            }
        }
        for (let i = 0; i < this.fixedPtms_.length; i++) {
            let ptm = this.fixedPtms_[i];
            let pos = ptm.getLeftPos();
            this.residues_[pos].color = "red";
        }
    }
}
