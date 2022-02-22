class PrsmPara {
  private rowLength_: number = 30 ;
  private blockLength_: number = 10 ;
  private letterWidth_: number = 28;
  private letterSize_: number = 12 ;
  private gapWidth_: number = 20;
  private rowHeight_: number = 40;
  private topMargin_: number = 46;
  private bottomMargin_: number = 10;
  private rightMargin_: number = 50;
  private leftMargin_: number = 50;
  private middleMargin_: number = 40;//extra space is needed when the start of sequence is skipped 
  private extraPadding_: number = 20;//increase svg width by this amount to prevent rightmost side of svg getting cutoff
  private numericalWidth_: number = 20;
  private showNum_: boolean = true;
  private showSkippedLines_: boolean = true;
  private skipLineHeight_: number = 40;
  private fontWidth_: number = 12 ;//12px font width = 9pt
  private fontHeight_: number = 18 ; 
  private svgBackgroundColor_: string = "white";
  private massShiftColor_: string = "#64E9EC";
  private modAnnoYShifts_: number[] = [-15, -30];
  private massShiftColors_: string[] = ["", "#64E9EC", "#ec6f64", "#ec6f64"];
  private isModAllowed_: boolean = false;//allow click to add mods (true in inspect page);

  constructor(isModAllowed: boolean) {
    this.isModAllowed_ = isModAllowed;
  }

  //getters
  getShowNum(): boolean {
    return this.showNum_;
  }
  getSvgBackgroundColor(): string {
    return this.svgBackgroundColor_;
  }
  /*getMassShiftColor(): string {
    return this.massShiftColor_;
  }*/
  getMassShiftColor(index: number): string {
    return this.massShiftColors_[index];
  }
  getMiddleMargin(): number {
    return this.middleMargin_;
  }
  getGapWidth(): number {
    return this.gapWidth_;
  }
  getRowLength(): number {
    return this.rowLength_;
  }
  getRowHeight(): number {
    return this.rowHeight_;
  }
  getLetterWidth(): number {
    return this.letterWidth_;
  }
  getNumericalWidth(): number {
    return this.numericalWidth_;
  }
  getFontWidth(): number {
    return this.fontWidth_;
  }
  getFontHeight(): number {
    return this.fontHeight_;
  }
  getModAnnoYShifts(): number[] {
    return this.modAnnoYShifts_;
  }
  getShowSkippedLines(): boolean {
    return this.showSkippedLines_;
  }
  getIsModAllowed(): boolean {
    return this.isModAllowed_;
  }


  setRowLength(rowLength: number): void {
    this.rowLength_= rowLength;
  }
  setRowHeight(rowHeight: number): void {
    this.rowHeight_ = rowHeight;
  }
  setLetterWidth(letterWidth: number): void {
    this.letterWidth_ = letterWidth;
  }
  setGapWidth(gapWidth: number): void {
    this.gapWidth_ = gapWidth;
  }
  setNumericalWidth(numWidth: number): void {
    this.numericalWidth_ = numWidth;
  }
  setShowSkippedLines(showSkippedLines: boolean): void{
    this.showSkippedLines_ = showSkippedLines;
  }

  setShowNum (isShowNum: boolean) {
    if (isShowNum) {
      this.showNum_ = true;
      this.leftMargin_ = 30;
      this.rightMargin_ = 30;
    }
    else {
      this.showNum_ = false;
      this.leftMargin_ = 10;
      this.rightMargin_ = 10;
    }
  }

  /**
   * Function to get x coordinate based on the position of the acid
   * @param {number} position - Contains position of the Acid
   * @param {number} start_value - Contains position of the first Acid
   */
  getX(pos: number, startPos: number): number{
    let num: number = pos - startPos;
    let posInRow: number = num % this.rowLength_;
    let gapNum: number = Math.floor(posInRow/this.blockLength_);
    let x: number = (posInRow) * this.letterWidth_ + gapNum * this.gapWidth_ + this.leftMargin_;
    if (this.showNum_) {
      x = x + this.numericalWidth_; 
    }
    return x;
  }

  /**
   * Function provides the Y coordinate based on the position of the Acid
   * @param {number} position - Contains position of the Acid
   * @param {number} start_value - Contains position of the first Acid
   */
  getY(pos: number, startPos: number): number {
    let row: number = Math.floor((pos - startPos) / this.rowLength_);
    let y: number  = row * this.rowHeight_ + this.topMargin_; 
    /*if(startPos != 0 && this.showSkippedLines_) {
      y = y + this.skipLineHeight_; 
    }*/
    return y;
  }

  getAACoordinates(pos: number, startPos: number): number[] {
    let x: number = this.getX(pos, startPos);
    let y: number = this.getY(pos, startPos);
    return [x,y];
  }

  getBpCoordinates(pos: number, startPos: number): number[] {
    let x: number = this.getX(pos-1, startPos) + this.letterWidth_/2;
    let y: number = this.getY(pos-1, startPos);
    return [x,y];
  }

  /**
   * Function provides position of the Numbers on the left side of the Acid Sequence
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getLeftNumCoordinates(pos: number, startPos: number): number[] {
    let x: number = this.leftMargin_;
    let y: number = this.getY(pos, startPos);
    return [x,y];
  }

  /**
   * Function provides position of the Numbers on the right side of the Acid Sequence
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getRightNumCoordinates(pos: number, startPos: number): number[] {
    let x: number = this.leftMargin_ + this.numericalWidth_ + (this.rowLength_ - 1 ) * this.letterWidth_;
    //buffer width-anno_width to make left and right numbers symmetrical as left numbers are left aligned 
    x = x + ((this.rowLength_/ this.blockLength_) - 1) * this.gapWidth_ + this.numericalWidth_ + this.fontWidth_; 
    let y: number = this.getY(pos, startPos);
    return [x,y];
  }

  /**
   * Function provides position to write information of skipped amino acids at the top of Sequence SVG
   */
  getSkipStartCoordinates(): number[] {
    let x: number = this.leftMargin_;
    let y: number = this.topMargin_; 
    return [x, y];
  }

  /**
   * Function provides position to write information of skipped amino acids at th bottom of Sequence SVG
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getSkipEndCoordinates(pos: number, startPos: number): number[] {
    let x: number = this.leftMargin_;
    let y: number = this.getY(pos, startPos) ; 
    return [x, y];
  }

  addHeight(): void {
    this.rowHeight_ = this.rowHeight_ * 1.2; 
    this.topMargin_ = 40; 
  }

  getSvgSize(rowNum: number, skipBegin: boolean, skipEnd: boolean): number[] {
    let blockNum: number = this.rowLength_/this.blockLength_ - 1 ;
    let width: number = this.letterWidth_ * (this.rowLength_ - 1) + blockNum * this.gapWidth_ + this.rightMargin_ + this.leftMargin_ + this.extraPadding_; 
    if(this.showNum_) {
      width = width + this.numericalWidth_ * 2;
    }
    let height = this.rowHeight_ * (rowNum -1) + this.letterSize_ + this.bottomMargin_ + this.topMargin_; 
    if (skipBegin) {
      height = height + this.skipLineHeight_;
    }
    if (skipEnd) {
      height = height + this.skipLineHeight_;
    }
    return [width,height];
  }
}	

