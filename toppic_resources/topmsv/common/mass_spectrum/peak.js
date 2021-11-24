"use strict";
class Peak {
    constructor(peakId, pos, monoMz, intensity, monoMass, charge, specId) {
        this.displayLevel_ = -1;
        this.peakId_ = peakId;
        this.pos_ = pos;
        this.monoMass_ = monoMass;
        this.monoMz_ = monoMz;
        this.charge_ = charge;
        this.intensity_ = intensity;
        this.specId_ = specId;
    }
    getId() {
        return this.peakId_;
    }
    getSpecId() {
        return this.specId_;
    }
    getPos() {
        return this.pos_;
    }
    getMonoMass() {
        return this.monoMass_;
    }
    getMonoMz() {
        return this.monoMz_;
    }
    getCharge() {
        return this.charge_;
    }
    getIntensity() {
        return this.intensity_;
    }
    getDisplayLevel() {
        return this.displayLevel_;
    }
    setPos(pos) {
        this.pos_ = pos;
    }
    setDisplayLevel(displayLevel) {
        this.displayLevel_ = displayLevel;
    }
    setMonoMass(mass) {
        this.monoMass_ = mass;
    }
    setMonoMz(mass) {
        this.monoMz_ = mass;
    }
    setIntensity(intensity) {
        this.intensity_ = intensity;
    }
}
