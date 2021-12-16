class PrsmViewData {
  //missing underscore in variable names
  private residues_: Residue[] = []; 
  private formFirstPos_: number = 0;
  private formLastPos_: number = 0;
  private breakPoints_: BreakPoints[] = [];

  private sequence_: string = "";
  private fixedPtms_: MassShift[] = [];
  private protVariablePtms_: MassShift[] = [];
  private variablePtms_: MassShift[] = [];
  private massShifts_: MassShift[] = [];
  private annotations_: Annotation[] = [];
  private proteoform_: Proteoform | null = null;

  //information above are retrieved from prsm and proteoform object

  private rowNum_: number = 0;
  private displayFirstPos_: number = 0;
  private displayLastPos_: number = 0;
  // if there is a skipping line at the beginning
  private showStartSkipped_: boolean = true; 
  private startSkippedInfo_: string = "";
  // if there is a skipping line at the end
  private showEndSkipped_: boolean = true;
  private endSkippedInfo_: string = "";

  //getters and setters
  getResidues(): Residue[] {
    return this.residues_;
  }
  getFormFirstPos(): number {
    return this.formFirstPos_;
  }
  getFormLastPos(): number {
    return this.formLastPos_;
  }
  getSequence(): string {
    return this.sequence_;
  }
  getBreakPoints(): BreakPoints[] {
    return this.breakPoints_;
  }
  getFixedPtms(): MassShift[] {
    return this.fixedPtms_;
  };
  getVariablePtms(): MassShift[] {
    return this.variablePtms_;
  }
  getMassShifts(): MassShift[] {
    return this.massShifts_;
  }
  getAnnotations(): Annotation[] {
    return this.annotations_;
  }
  getProteoform(): Proteoform | null {
    return this.proteoform_;
  }
  getRowNum(): number {
    return this.rowNum_;
  }
  getDisplayFirstPos(): number {
    return this.displayFirstPos_;
  }
  getDisplayLastPos(): number {
    return this.displayLastPos_;
  }
  getShowStartSkipped(): boolean {
    return this.showStartSkipped_;
  }
  getStartSkippedInfo(): string {
    return this.startSkippedInfo_;
  }
  getShowEndSkipped(): boolean {
    return this.showEndSkipped_;
  }
  getEndSkippedInfo(): string {
    return this.endSkippedInfo_;
  }

  initData(prsm: Prsm, para: PrsmPara): void{
    this.parseData(prsm);
    this.updatePara(para);
    this.addColor();
  }
  // form residues from sequence
  formResidues(sequence: string): Residue[] {
    let residues: Residue[] = [];
    for (let i = 0; i < sequence.length; i++) {
      let tempObj: Residue = {
        position: i,
        acid: sequence.charAt(i).toUpperCase(),
        color: ""
      }
      residues.push(tempObj);
    }
    return residues;
  }
  generateAnnotation(protVarPtm: MassShift[], variablePtm: MassShift[], massShift: MassShift[]): Annotation[]{
    let annos: Annotation[] = [];
    //annotation in prsm object contains only variable and unknown shifts
    for (let i = 0; i < protVarPtm.length; i++){
      let temp: Annotation = {"annoText":protVarPtm[i].getAnnotation(), "leftPos":0, "rightPos":0, "type": ModType.ProteinVariable};
      temp.leftPos = protVarPtm[i].getLeftPos();
      temp.rightPos = protVarPtm[i].getRightPos();
      annos.push(temp);
    }
    for (let i = 0; i < variablePtm.length; i++){
      let temp: Annotation= {"annoText":variablePtm[i].getAnnotation(), "leftPos":0, "rightPos":0, "type": ModType.Variable};
      temp.leftPos = variablePtm[i].getLeftPos();
      temp.rightPos = variablePtm[i].getRightPos();
      annos.push(temp);
    }
    for (let i = 0; i < massShift.length; i++){
      let temp: Annotation = {"annoText":FormatUtil.formatFloat(massShift[i].getAnnotation(), "massShift"), "leftPos":massShift[i].getLeftPos(), "rightPos":massShift[i].getRightPos(), "type": ModType.Unexpected};
      annos.push(temp);
    }
    //merge variable ptms with same residue range
    annos.sort(function(a,b) {
      return a.leftPos - b.leftPos;
    })
    let prevAnno: Annotation = {} as Annotation;
    let prevLeft: number = -1;
    let prevRight: number = -1;
    let newAnnos: Annotation[] = [];

    annos.forEach((anno, i) => {
      if (anno.leftPos != prevLeft || anno.rightPos != prevRight) {
        if (prevLeft >= 0) {
          newAnnos.push(prevAnno);
        }
        prevAnno = anno;
      }
      else {
        prevAnno.annoText = prevAnno.annoText + ";" + anno.annoText;
      }
      prevLeft = anno.leftPos;
      prevRight = anno.rightPos;
    })
    newAnnos.push(prevAnno);

    return newAnnos;
  }
  parseData(prsm: Prsm) {
    let proteoformObj = prsm.getProteoform();
    this.residues_ = this.formResidues(proteoformObj.getSeq());
    this.formFirstPos_ = proteoformObj.getFirstPos();
    this.formLastPos_ = proteoformObj.getLastPos();
    this.breakPoints_ = prsm.getBreakPoints();
    this.fixedPtms_ = proteoformObj.getFixedPtm();
    this.protVariablePtms_ = proteoformObj.getProtVarPtm();
    this.variablePtms_ = proteoformObj.getVarPtm();
    this.massShifts_ = proteoformObj.getUnknownMassShift();
    this.sequence_ = proteoformObj.getSeq();
    this.proteoform_ = proteoformObj;
    this.annotations_ = this.generateAnnotation(this.protVariablePtms_, this.variablePtms_, this.massShifts_);
    
  }
 
  updatePara(para: PrsmPara): void {
    let len: number = this.residues_.length; 
    // Include 5 amino acids before and after the form
    this.displayFirstPos_ = Math.floor((this.formFirstPos_ - 5) / para.getRowLength()) * para.getRowLength();
    if (this.displayFirstPos_ < 0) {
      this.displayFirstPos_ = 0;
    }
    this.displayLastPos_ = Math.ceil((this.formLastPos_ + 6) / para.getRowLength()) * para.getRowLength() - 1;
    //console.log("display last pos ", this.displayLastPos);
    if (this.displayLastPos_ > (len -1)) {
      this.displayLastPos_ = len -1;
    }
    this.rowNum_ = Math.ceil((this.displayLastPos_ - this.displayFirstPos_ + 1)/para.getRowLength());
    
    // skipping line
    this.showStartSkipped_ = false;
    this.showEndSkipped_ = false;
    if (para.getShowSkippedLines()) {
      if (this.displayFirstPos_ !== 0) {
        this.showStartSkipped_ = true;
		    this.startSkippedInfo_ = "... "+ this.displayFirstPos_ 
          + " amino acid residues are skipped at the N-terminus ... ";
      }
      if (this.displayLastPos_ !== len - 1) {
        this.showEndSkipped_ = true;
		    this.endSkippedInfo_ =  "... "+ (len - 1 - this.displayLastPos_) 
          + " amino acid residues are skipped at the C-terminus ... ";
      }
    }

  }

  addColor() {
    for (let i = 0; i < this.residues_.length; i++) {
      let residue = this.residues_[i];
      let pos = residue.position;
      if (pos < this.formFirstPos_ || pos > this.formLastPos_) {
        residue.color = "grey";
      }
      else {
        residue.color = "black";
      }
    }
    for (let i = 0; i < this.fixedPtms_.length; i++) {
      let ptm = this.fixedPtms_[i];
      let pos = ptm.getLeftPos();
      this.residues_[pos].color = "red";
    }
  }
}