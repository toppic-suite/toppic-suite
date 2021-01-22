class Envelope{
    color_ = "";
    level_ = -1;
    theo_peaks_ = [];

    constructor(monoMass, charge, intensity = -1){
        this.monoMass_ = monoMass;
        this.charge_ = charge;
        this.intensity_ = intensity;
    }
    getMonoMass(){
        return this.monoMass_;
    }
    getCharge(){
        return this.charge_;
    }
    getIntensity(){
        return this.intensity_;
    }
    getColor(){
        return this.color_;
    }
    getLevel(){
        return this.level_;
    }
    getTheoPeaks(){
        return this.theo_peaks_;
    }
    setColor(color){
        this.color_ = color;
    }
    setLevel(level){
        this.level_ = level;
    }
    addTheoPeaks(peak){
        this.theo_peaks_.push(peak);
    }
}