"use strict";
class Ion {
    constructor(id, name, terminal, massShift, massError, ppmError) {
        this.id_ = id;
        this.name_ = name;
        this.terminal_ = terminal; //"N" or "C"
        this.massShift_ = massShift;
        this.massError_ = massError;
        this.ppmError_ = ppmError;
    }
    getId() {
        return this.id_;
    }
    getName() {
        return this.name_;
    }
    getTerminal() {
        return this.terminal_;
    }
    getShift() {
        return this.massShift_;
    }
    getMassError() {
        return this.massError_;
    }
    getPpmError() {
        return this.ppmError_;
    }
}
