"use strict";
class MassShift {
    constructor(leftPos, rightPos, massShift, type, annotation, ptm = null) {
        this.leftPos_ = leftPos; //for drawing annotation background
        this.rightPos_ = rightPos;
        this.massShift_ = massShift;
        this.type_ = this.setModType(type);
        this.annotation_ = annotation;
        this.ptmList_ = [];
        if (ptm) {
            this.ptmList_.push(ptm); //ignore the possibility of ptm being null
        }
    }
    addNewPtm(ptm) {
        this.ptmList_.push(ptm);
    }
    getLeftPos() {
        return this.leftPos_;
    }
    getRightPos() {
        return this.rightPos_;
    }
    getShift() {
        return this.massShift_;
    }
    getType() {
        return this.type_;
    }
    getAnnotation() {
        return this.annotation_;
    }
    getPtmList() {
        return this.ptmList_;
    }
    setLeftPos(pos) {
        this.leftPos_ = pos;
    }
    setRightPos(pos) {
        this.rightPos_ = pos;
    }
    setModType(type) {
        let modType;
        if (type == "Fixed") {
            modType = ModType.Fixed;
        }
        else if (type == "Protein variable") {
            modType = ModType.ProteinVariable;
        }
        else if (type == "Variable") {
            modType = ModType.Variable;
        }
        else {
            modType = ModType.Unexpected;
        }
        return modType;
    }
    setPtmList(ptmList) {
        this.ptmList_ = ptmList;
    }
}
