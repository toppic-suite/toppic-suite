"use strict";
class Envelope {
    constructor(monoMass, charge, intensity = -1) {
        this.displayColor_ = "";
        this.displayLevel_ = -1;
        this.peaks_ = [];
        this.monoMass_ = monoMass;
        this.charge_ = charge;
        this.intensity_ = intensity;
    }
    getMonoMass() {
        return this.monoMass_;
    }
    getCharge() {
        return this.charge_;
    }
    getIntensity() {
        return this.intensity_;
    }
    getDisplayColor() {
        return this.displayColor_;
    }
    getDisplayLevel() {
        return this.displayLevel_;
    }
    getPeaks() {
        return this.peaks_;
    }
    setDisplayColor(displayColor) {
        this.displayColor_ = displayColor;
    }
    setDisplayLevel(displayLevel) {
        this.displayLevel_ = displayLevel;
    }
    addPeaks(peak) {
        this.peaks_.push(peak);
    }
}
