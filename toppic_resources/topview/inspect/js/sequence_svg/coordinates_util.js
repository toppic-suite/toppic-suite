/**
 * Fixed Parameters to draw the Sequence SVG
 */
function parameters()
{
	this.row_length = 40 ;
	this.block_length = 10 ;
	this.letter_width = 25;
	this.letter_size = 12 ;
	this.gap_width = 20;
	this.row_height = 40;
	this.top_margin = 35;
	this.bottom_margin = 30 ;
	this.right_margin = 50;
	this.left_margin = 40;
	this.numerical_width = 20;
	this.show_num = true ;
	this.show_skipped_lines = true ;
	this.skip_line_height = 40;
	this.font_width = 9 ;//12px font with = 9pt
	this.background_color = "#64E9EC";
	this.svgBackground_color = "white" ;

	return this;
}	
/**
 * Function provides the Y coordinate based on the position of the Acid
 * @param {Object} para - Contains parameters of the graph
 * @param {number} position - Contains position of the Acid
 */
function getY(para, position) 
{
	let row = parseInt((position ) / para.row_length);
	let y  = row * para.row_height + para.top_margin; 
  	return y;
}
/**
 * Function provides with X coordinate
 * @param {Object} para - Contains parameters of the graph
 * @param {number} position - Contains position of the Acid
 */
function getX(para,position)
{
	let position_temp = position ;
	let pos_in_row = position_temp % para.row_length ;
	let gap_num = parseInt(pos_in_row/para.block_length) ;
	let x = (pos_in_row) * para.letter_width + gap_num * para.gap_width + para.left_margin;
	x = x + para.numerical_width;
	return x ;
}
/**
 * Function provides position of the Numbers on the left side of the Acid Sequence
 * @param {Object} para - Contains parameters of the graph
 * @param {number} position - Contains position of the Acid
 */
function calibrateLeftNum(para,position) 
{
	/*console.log("para.left_margin : ", para.left_margin)*/
  let x = para.left_margin ;
  let y = getY(para, position);
  return [x,y];
}
/**
 * Function provides position of the Numbers on the right side of the Acid Sequence
 * @param {object} para - Contains the parameters to draw the SVG
 * @param {Integer} position - Provides the position of the left side number
 */
function calibrateRightNum(para,position) 
{
  let x = para.left_margin + para.numerical_width + (para.row_length - 1 ) * para.letter_width;
  //buffer width-anno_width to make left and right numbers symmetrical as left numbers are left aligned 
  x = x + ((para.row_length/ para.block_length) - 1) * para.gap_width + para.numerical_width + para.font_width; 
  let y = getY(para, position);
  return [x,y];
}
/**
 * Function provides position to write information of skipped amino acids at the top of Sequence SVG
 * @param {object} para - Contains the parameters to draw the SVG
 */
function calibrateSkipStart(para, position, start_value)
{
	x = para.left_margin ;
	y = para.top_margin; 
	return [x, y]
}
/**
 * Function provides position to write information of skipped amino acids at th bottom of Sequence SVG
 * @param {object} para - Contains the parameters to draw the SVG
 * @param {Integer} position - Provides the position of the left side number
 * @param {Integer} start_value - Provides starting number the sequnce after trimming the skipped scids
 */
function calibrateSkipEnd(para, position, start_value)
{
	x = para.left_margin ;
	y = getY(para, position, start_value) ; 
	return [x, y]
}

