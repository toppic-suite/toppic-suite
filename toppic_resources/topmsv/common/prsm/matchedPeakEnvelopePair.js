"use strict";
class MatchedPeakEnvelopePair {
    constructor(theoMass, monoPeak, ion, envelope) {
        this.envelope_ = null;
        this.theoMass_ = theoMass;
        this.peak_ = monoPeak;
        this.ion_ = ion;
        if (envelope) {
            this.envelope_ = envelope; //only one envelope per object
        }
    }
    getPeak() {
        return this.peak_;
    }
    getIon() {
        return this.ion_;
    }
    getEnvelope() {
        return this.envelope_;
    }
    getTheoMass() {
        return this.theoMass_;
    }
    setTheoMass(mass) {
        this.theoMass_ = mass;
    }
    addEnvelope(env) {
        this.envelope_ = env;
    }
}
