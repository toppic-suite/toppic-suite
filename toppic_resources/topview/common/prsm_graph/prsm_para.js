class PrsmPara {
  rowLength = 30 ;
  blockLength = 10 ;
  letterWidth = 28;
  letterSize = 12 ;
  gapWidth = 20;
  rowHeight = 40;
  topMargin = 36;
  bottomMargin = 10 ;
  rightMargin = 50;
  leftMargin = 50;
  middleMargin = 40;//extra space is needed when the start of sequence is skipped 
  numericalWidth = 20;
  showNum = true ;
  showSkippedLines = true ;
  skipLineHeight = 40;
  fontWidth = 12 ;//12px font width = 9pt
  fontHeight = 18 ; 
  svgBackgroundColor = "white" ;
  massShiftColor = "#64E9EC";
  modAnnoYShifts = [-15, -30];

  setShowNum = function (isShowNum) {
    if (isShowNum) {
      this.showNum = true;
      this.leftMargin = 30;
      this.rightMargin = 30;
    }
    else {
      this.showNum = false;
      this.leftMargin = 10;
      this.rightMargin = 10;
    }
  }

  /**
   * Function to get x coordinate based on the position of the acid
   * @param {number} position - Contains position of the Acid
   * @param {number} start_value - Contains position of the first Acid
   */
  getX = function (pos,startPos) {
    let num = pos - startPos ;
    let posInRow = num % this.rowLength ;
    let gapNum = parseInt(posInRow/this.blockLength) ;
    let x = (posInRow) * this.letterWidth + gapNum * this.gapWidth + this.leftMargin;
    if (this.showNum) {
      x = x + this.numericalWidth; 
    }
    return x;
  }

  /**
   * Function provides the Y coordinate based on the position of the Acid
   * @param {number} position - Contains position of the Acid
   * @param {number} start_value - Contains position of the first Acid
   */
  getY = function (pos, startPos) {
    let row = parseInt((pos - startPos) / this.rowLength);
    let y  = row * this.rowHeight + this.topMargin; 
    if(startPos != 0 && this.showSkippedLine) {
      y = y + this.skipLineHeight; 
    }
    return y;
  }

  getAACoordinates = function (pos, startPos) {
    let x = this.getX(pos, startPos);
    let y = this.getY(pos, startPos);
    return [x,y];
  }

  getBpCoordinates = function (pos, startPos) {
    let x = this.getX(pos-1, startPos) + this.letterWidth/2;
    let y = this.getY(pos-1, startPos);
    return [x,y];
  }


  /**
   * Function provides position of the Numbers on the left side of the Acid Sequence
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getLeftNumCoordinates = function (pos, startPos) {
    let x = this.leftMargin;
    let y = this.getY(pos, startPos);
    return [x,y];
  }

  /**
   * Function provides position of the Numbers on the right side of the Acid Sequence
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getRightNumCoordinates = function (pos,startPos) {
    let x = this.leftMargin + this.numericalWidth + (this.rowLength - 1 ) * this.letterWidth;
    //buffer width-anno_width to make left and right numbers symmetrical as left numbers are left aligned 
    x = x + ((this.rowLength/ this.blockLength) - 1) * this.gapWidth + this.numericalWidth + this.fontWidth; 
    let y = this.getY(pos, startPos);
    return [x,y];
  }

  /**
   * Function provides position to write information of skipped amino acids at the top of Sequence SVG
   */
  getSkipStartCoordinates = function () {
    let x = this.leftMargin ;
    let y = this.topMargin; 
    return [x, y]
  }

  /**
   * Function provides position to write information of skipped amino acids at th bottom of Sequence SVG
   * @param {Integer} position - Provides the position of the left side number
   * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
   */
  getSkipEndCoordinates = function(pos, startPos) {
    let x = this.leftMargin ;
    let y = this.getY(pos, startPos) ; 
    return [x, y]
  }

  addHeight = function() {
    this.rowHeight = this.rowHeight * 1.2; 
    this.topMargin = 40; 
  }

  getSvgSize = function(rowNum, skipBegin, skipEnd) {
    let blockNum = this.rowLength/this.blockLength - 1 ;
    let width = this.letterWidth * (this.rowLength - 1) + blockNum * this.gapWidth + this.rightMargin + this.leftMargin; 
    if(this.showNum) {
      width = width + this.numericalWidth * 2;
    }
    let height = this.rowHeight * (rowNum -1) + this.letterSize + this.bottomMargin + this.topMargin; 
    if (skipBegin) {
      height = height + this.skipLineHeight;
    }
    if (skipEnd) {
      height = height + this.skipLineHeight;
    }
    return [width,height];
  }
}	

