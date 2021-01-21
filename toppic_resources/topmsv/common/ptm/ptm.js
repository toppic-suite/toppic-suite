class Ptm{
    constructor(residue, massShift, name){
        this.residue_ = residue;
        this.massShift_ = parseFloat(massShift);
        this.name_ = name;
    }
    getResidue(){
        return this.residue_;
    }
    getShift(){
        return this.massShift_;
    }
    getName(){
        return this.name_;
    }
}