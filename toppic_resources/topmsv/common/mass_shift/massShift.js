class MassShift{
    ptmList_ = [];

    constructor(leftPos, rightPos, massShift, type, annotation, ptm = null){
        this.leftPos_ = parseInt(leftPos);
        this.rightPos_ = parseInt(rightPos);
        this.massShift_ = parseFloat(massShift).toFixed(4);
        this.type_ = type;
        this.annotation_ = annotation;

        /*if (this.type_ == "Fixed"){
            this.annotation_ = "";
        }*/
        if (ptm){
            this.ptmList_.push(ptm);
        }
    }
    addNewPtm(ptm){
        this.ptmList_.push(ptm);
    }
    getLeftPos(){
        return this.leftPos_;
    }
    getRightPos(){
        return this.rightPos_;
    }
    getShift(){
        return this.massShift_;
    }
    getType(){
        return this.type_;
    }
    getAnnotation(){
        return this.annotation_;
    }
    getPtmList(){
        return this.ptmList_;
    }
}