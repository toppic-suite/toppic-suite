"use strict";
class Spectrum {
    constructor(id, scanNum, level, peakList, decovPeakList, envList = [], nIon, cIon, mass, charge = -1, mz = -1, minMz, maxMz) {
        this.minMz_ = -1;
        this.maxMz_ = -1;
        this.envList_ = [];
        this.id_ = id;
        this.scanNum_ = scanNum;
        this.level_ = level;
        this.peakList_ = peakList;
        this.deconvPeakList_ = decovPeakList;
        this.envList_ = envList;
        this.nIon_ = nIon;
        this.cIon_ = cIon;
        this.precMass_ = mass;
        this.precCharge_ = charge;
        this.precMz_ = mz;
        if (minMz) {
            this.minMz_ = minMz; // only used for ms1 spectrum. target mz - lower offset
        }
        if (maxMz) {
            this.maxMz_ = maxMz; // only used for ms1 spectrum. target mz + higher offset
        }
    }
    getSpectrumId() {
        return this.id_;
    }
    getScanNum() {
        return this.scanNum_;
    }
    getSpectrumLevel() {
        return this.level_;
    }
    getPeaks() {
        return this.peakList_;
    }
    getDeconvPeaks() {
        return this.deconvPeakList_;
    }
    getEnvs() {
        return this.envList_;
    }
    getNTerminalIon() {
        return this.nIon_;
    }
    getCTerminalIon() {
        return this.cIon_;
    }
    getPrecMass() {
        return this.precMass_;
    }
    getPrecCharge() {
        return this.precCharge_;
    }
    getPrecMz() {
        return this.precMz_;
    }
    getMaxMz() {
        return this.maxMz_;
    }
    getMinMz() {
        return this.minMz_;
    }
    setPeaks(peaks) {
        this.peakList_ = peaks;
    }
    setDeconvPeaks(decovPeaks) {
        this.deconvPeakList_ = decovPeaks;
    }
    setEnvs(envs) {
        this.envList_ = envs;
    }
    addNTerminalIon(ion) {
        this.nIon_.push(ion);
    }
    addCTerminalIon(ion) {
        this.cIon_.push(ion);
    }
}
