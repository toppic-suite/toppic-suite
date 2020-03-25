function parameters()
{
	this.row_length = 30 ;
	this.block_length = 10 ;
	this.letter_width = 28;
	this.letter_size = 12 ;
	this.gap_width = 20;
	this.row_height = 40;
	this.top_margin = 35;
	this.bottom_margin = 10 ;
	this.right_margin = 50;
	this.left_margin = 40;
	this.numerical_width = 20;
	this.show_num = true ;
	this.show_skipped_lines = true ;
	this.skip_line_height = 40;
	this.font_width = 9 ;//12px font with = 9pt
	this.background_color = "#64E9EC";
	this.svgBackground_color = "white" ;
}	

function getY(para, position, start_value) 
{
	let row = parseInt((position - start_value) / para.row_length);
	let y  = row * para.row_height + para.top_margin; 
	if(start_value != 0 && para.show_skipped_lines == true)
	{
		y = y + para.skip_line_height;
	}
  return y;
}
function calibrateCoordinates(para,position,start_value)
{
	let position_temp = position - start_value ;
	let pos_in_row = position_temp % para.row_length ;
	let gap_num = parseInt(pos_in_row/para.block_length) ;
	let x = (pos_in_row) * para.letter_width + gap_num * para.gap_width + para.left_margin;
	if(para.show_num) x = x + para.numerical_width; 
	y = getY(para, position, start_value);
	return [x,y] ;
}
function calibrateLeftNum(para,position,start_value) 
{
  let x = para.left_margin ;
  let y = getY(para, position, start_value);
  return [x,y];
}

function calibrateRightNum(para,position,start_value) 
{
  let x = para.left_margin + para.numerical_width + (para.row_length - 1 ) * para.letter_width;
  //buffer width-anno_width to make left and right numbers symmetrical as left numbers are left aligned 
  x = x + ((para.row_length/ para.block_length) - 1) * para.gap_width + para.numerical_width + para.font_width; 
  let y = getY(para, position, start_value);
  return [x,y];
}

function calibrateSkipStart(para, position, start_value)
{
	x = para.left_margin ;
	y = para.top_margin; 
	return [x, y]
}
function calibrateSkipEnd(para, position, start_value)
{
	x = para.left_margin ;
	y = getY(para, position, start_value) ; 
	return [x, y]
}

