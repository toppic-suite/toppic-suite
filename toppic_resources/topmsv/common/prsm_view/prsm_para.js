"use strict";
class PrsmPara {
    constructor(isModAllowed) {
        this.rowLength_ = 30;
        this.blockLength_ = 10;
        this.letterWidth_ = 28;
        this.letterSize_ = 12;
        this.gapWidth_ = 20;
        this.rowHeight_ = 40;
        this.topMargin_ = 46;
        this.bottomMargin_ = 10;
        this.rightMargin_ = 50;
        this.leftMargin_ = 50;
        this.middleMargin_ = 40; //extra space is needed when the start of sequence is skipped 
        this.extraPadding_ = 20; //increase svg width by this amount to prevent rightmost side of svg getting cutoff
        this.numericalWidth_ = 20;
        this.showNum_ = true;
        this.showSkippedLines_ = true;
        this.skipLineHeight_ = 40;
        this.fontWidth_ = 12; //12px font width = 9pt
        this.fontHeight_ = 18;
        this.svgBackgroundColor_ = "white";
        this.massShiftColor_ = "#64E9EC";
        this.modAnnoYShifts_ = [-15, -30];
        this.massShiftColors_ = ["", "#64E9EC", "#ec6f64", "#ec6f64"];
        this.isModAllowed_ = false; //allow click to add mods (true in inspect page);
        this.isModAllowed_ = isModAllowed;
    }
    //getters
    getShowNum() {
        return this.showNum_;
    }
    getSvgBackgroundColor() {
        return this.svgBackgroundColor_;
    }
    /*getMassShiftColor(): string {
      return this.massShiftColor_;
    }*/
    getMassShiftColor(index) {
        return this.massShiftColors_[index];
    }
    getMiddleMargin() {
        return this.middleMargin_;
    }
    getGapWidth() {
        return this.gapWidth_;
    }
    getRowLength() {
        return this.rowLength_;
    }
    getRowHeight() {
        return this.rowHeight_;
    }
    getLetterWidth() {
        return this.letterWidth_;
    }
    getNumericalWidth() {
        return this.numericalWidth_;
    }
    getFontWidth() {
        return this.fontWidth_;
    }
    getFontHeight() {
        return this.fontHeight_;
    }
    getModAnnoYShifts() {
        return this.modAnnoYShifts_;
    }
    getShowSkippedLines() {
        return this.showSkippedLines_;
    }
    getIsModAllowed() {
        return this.isModAllowed_;
    }
    setRowLength(rowLength) {
        this.rowLength_ = rowLength;
    }
    setRowHeight(rowHeight) {
        this.rowHeight_ = rowHeight;
    }
    setLetterWidth(letterWidth) {
        this.letterWidth_ = letterWidth;
    }
    setGapWidth(gapWidth) {
        this.gapWidth_ = gapWidth;
    }
    setNumericalWidth(numWidth) {
        this.numericalWidth_ = numWidth;
    }
    setShowSkippedLines(showSkippedLines) {
        this.showSkippedLines_ = showSkippedLines;
    }
    setShowNum(isShowNum) {
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
    getX(pos, startPos) {
        let num = pos - startPos;
        let posInRow = num % this.rowLength_;
        let gapNum = Math.floor(posInRow / this.blockLength_);
        let x = (posInRow) * this.letterWidth_ + gapNum * this.gapWidth_ + this.leftMargin_;
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
    getY(pos, startPos) {
        let row = Math.floor((pos - startPos) / this.rowLength_);
        let y = row * this.rowHeight_ + this.topMargin_;
        /*if(startPos != 0 && this.showSkippedLines_) {
          y = y + this.skipLineHeight_;
        }*/
        return y;
    }
    getAACoordinates(pos, startPos) {
        let x = this.getX(pos, startPos);
        let y = this.getY(pos, startPos);
        return [x, y];
    }
    getBpCoordinates(pos, startPos) {
        let x = this.getX(pos - 1, startPos) + this.letterWidth_ / 2;
        let y = this.getY(pos - 1, startPos);
        return [x, y];
    }
    /**
     * Function provides position of the Numbers on the left side of the Acid Sequence
     * @param {Integer} position - Provides the position of the left side number
     * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
     */
    getLeftNumCoordinates(pos, startPos) {
        let x = this.leftMargin_;
        let y = this.getY(pos, startPos);
        return [x, y];
    }
    /**
     * Function provides position of the Numbers on the right side of the Acid Sequence
     * @param {Integer} position - Provides the position of the left side number
     * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
     */
    getRightNumCoordinates(pos, startPos) {
        let x = this.leftMargin_ + this.numericalWidth_ + (this.rowLength_ - 1) * this.letterWidth_;
        //buffer width-anno_width to make left and right numbers symmetrical as left numbers are left aligned 
        x = x + ((this.rowLength_ / this.blockLength_) - 1) * this.gapWidth_ + this.numericalWidth_ + this.fontWidth_;
        let y = this.getY(pos, startPos);
        return [x, y];
    }
    /**
     * Function provides position to write information of skipped amino acids at the top of Sequence SVG
     */
    getSkipStartCoordinates() {
        let x = this.leftMargin_;
        let y = this.topMargin_;
        return [x, y];
    }
    /**
     * Function provides position to write information of skipped amino acids at th bottom of Sequence SVG
     * @param {Integer} position - Provides the position of the left side number
     * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
     */
    getSkipEndCoordinates(pos, startPos) {
        let x = this.leftMargin_;
        let y = this.getY(pos, startPos);
        return [x, y];
    }
    addHeight() {
        this.rowHeight_ = this.rowHeight_ * 1.2;
        this.topMargin_ = 40;
    }
    getSvgSize(rowNum, skipBegin, skipEnd) {
        let blockNum = this.rowLength_ / this.blockLength_ - 1;
        let width = this.letterWidth_ * (this.rowLength_ - 1) + blockNum * this.gapWidth_ + this.rightMargin_ + this.leftMargin_ + this.extraPadding_;
        if (this.showNum_) {
            width = width + this.numericalWidth_ * 2;
        }
        let height = this.rowHeight_ * (rowNum - 1) + this.letterSize_ + this.bottomMargin_ + this.topMargin_;
        if (skipBegin) {
            height = height + this.skipLineHeight_;
        }
        if (skipEnd) {
            height = height + this.skipLineHeight_;
        }
        return [width, height];
    }
}
