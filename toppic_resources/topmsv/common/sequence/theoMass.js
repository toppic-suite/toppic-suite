"use strict";
class TheoMass {
    constructor(mass, pos) {
        this.ion_ = {};
        this.mass_ = mass;
        this.pos_ = pos;
    }
    getMass() {
        return this.mass_;
    }
    setMass(mass) {
        this.mass_ = mass;
    }
    getPos() {
        return this.pos_;
    }
    getIon() {
        return this.ion_;
    }
    setIon(ion) {
        this.ion_ = ion;
    }
}
