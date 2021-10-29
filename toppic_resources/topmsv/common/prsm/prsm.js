"use strict";
class Prsm {
    constructor(id, proteoform, ms1Spec, ms2Spec, breakPoints, matchedPeakEnvelopePair = [], fileName = "", eValue = -1, qValue = -1, featureInte, precMass, fragIonCount) {
        this.matchedPeaks_ = [];
        this.matchedIons_ = [];
        this.matchedPeakEnvelopePair_ = [];
        this.id_ = id;
        this.proteoform_ = proteoform;
        this.ms1Spec_ = ms1Spec;
        this.ms2Spec_ = ms2Spec;
        this.breakPoints_ = breakPoints;
        this.matchedPeakEnvelopePair_ = matchedPeakEnvelopePair;
        this.fileName_ = fileName;
        this.eValue_ = eValue;
        this.qValue_ = qValue;
        this.fragIonCount_ = fragIonCount;
        if (featureInte) {
            this.featureInte_ = featureInte;
        }
        if (precMass) {
            this.precMass_ = precMass;
        }
    }
    getId() {
        return this.id_;
    }
    getFragIonCount() {
        return this.fragIonCount_;
    }
    getfileName() {
        return this.fileName_;
    }
    getProteoform() {
        return this.proteoform_;
    }
    getMs1Spectra() {
        return this.ms1Spec_;
    }
    getMs2Spectra() {
        return this.ms2Spec_;
    }
    getMatchedPeakCount() {
        return this.getMatchedPeakEnvelopePairs().length;
    }
    getUnexpectedModCount() {
        let protObj = this.getProteoform();
        let unexpectedMod = protObj.getUnknownMassShift();
        return unexpectedMod.length;
    }
    getEValue() {
        return this.eValue_;
    }
    getQValue() {
        return this.qValue_;
    }
    getFeatureInte() {
        return this.featureInte_;
    }
    getPrecMass() {
        return this.precMass_;
    }
    getBreakPoints() {
        return this.breakPoints_;
    }
    getMatchedPeakEnvelopePairs() {
        return this.matchedPeakEnvelopePair_;
    }
    setBreakPoints(breakPoints) {
        this.breakPoints_ = breakPoints;
    }
    setProteoform(proteoform) {
        this.proteoform_ = proteoform;
    }
}
