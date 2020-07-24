/**
 * Get the size of the svg based on the no. of rows and row length and other parameter.,.
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * prsm is an json attribute insire prsm_data(Global data variable from prsm.js)
 */
function getSvgSize(parameters,prsm)
{
	let num_of_rows = getNumOfRows(parameters,prsm) ;
	let no_of_blocks = parameters.row_length/parameters.block_length - 1 ;
	let width = parameters.letter_width * (parameters.row_length - 1)+ (no_of_blocks)*parameters.gap_width + parameters.right_margin + parameters.left_margin;
	if(parameters.show_num)
	{
		width = width + parameters.numerical_width * 2 ;
	}
	if(isShiftAnnotationNeeded(parameters,prsm))
	{
		parameters.row_height = parameters.row_height + 0.2*parameters.row_height;
		parameters.top_margin = 45 ;
	}
	let height = parameters.row_height * num_of_rows + parameters.bottom_margin + parameters.top_margin ;
	let first_position,last_position,start_info = null ,end_info = null ;
	[parameters,first_position, last_position,start_info,end_info] = skip_list(parameters,prsm);
	return [width,height];
}
/**
 * Get number of rows of the sequence based on the row length
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * prsm is an json attribute insire prsm_data(Global data variable from prsm.js)
 */
function getNumOfRows(parameters,prsm)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	[parameters,first_position, last_position,start_info,end_info] = skip_list(parameters,prsm);
	let new_sequence_length = last_position - first_position ;
	let num_of_rows = parseInt(new_sequence_length/parameters.row_length) ;
	let skip_acid_count = 0 ;
	if((start_info != null || end_info != null)&& parameters.show_skipped_lines)
	{
		skip_acid_count = skip_acid_count + 1;
	}
	num_of_rows = num_of_rows + skip_acid_count;
	return num_of_rows ;
}
/**
 * draw the sequence on to svg
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function buildSvg(parameters,prsm,id)
{
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	let seqlength = parseInt(prsm.annotated_protein.annotation.protein_length) ;
	// Accomodate large numerical number at the start and end of the sequence without getting trimmed by adding extra margin values
	if(seqlength >= 1000)
	{
		parameters.left_margin = parameters.left_margin + parameters.font_width ;
		parameters.right_margin = parameters.right_margin + parameters.font_width ;
	}
	let width,height ;
	[width,height] = getSvgSize(parameters,prsm) ;
	// create a group under svg with svgId_g
	let id_temp = id + "_g" ;
	let svgContainer = d3.select("#"+id).attr("width",width)
									.attr("height",height)
									.attr("font-family","'FreeMono',Miltonian,monospace")
									.attr("font-size","16px")
									.style("fill", parameters.svgBackground_color)
	svgContainer = svgContainer.append("g")
								.attr("id",id_temp)
								.attr("class",id_temp);
	text = 	svgContainer.selectAll("text");
	let first_position,last_position,start_info = null ,end_info = null ;
	//	Get the new first and last position based on the amount of acids are needed to be skipped 
	[parameters,first_position, last_position,start_info,end_info] = skip_list(parameters,prsm);
	
	prsm.annotated_protein.annotation.residue.forEach(function(input,index){
		if(parseInt(input.position) >= first_position && parseInt(input.position) < last_position )
		{
			let x,y ;
			//	Get the x and y coordinates of the acid position	
			[x,y] = calibrateCoordinates(parameters,parseInt(input.position),first_position) ;
			text.data(input.acid)
				.enter()
				.append("text")
				.attr("id",function(d,i){
					return "id_"+id_temp+"_"+input.position ;
				})
				.attr("x", function(d,i){ 
					return x ;
				})
				.attr("y", function(d,i){ 
					return y ;
				})
				.text(function(d,i){
					return d ;
				})
				.style("fill", function(d,i){
					if(parseInt(input.position) < first_residue_position || parseInt(input.position) > last_residue_position)
					{
						return "grey" ;
					}
					else
					{
						return "black" ;
					}
				})
		}
	});
	return [parameters,id] ;
}

/**
 * Get the terminated/skipped acid information on to the svg
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function skippedAcidNotification(parameters,prsm,id)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	//	Get the new first and last position based on the amount of acids are needed to be skipped
	[parameters,first_position, last_position,start_info,end_info] = skip_list(parameters,prsm);
	let svgContainer = d3.select("#"+id+"_g") ;
	let x,y ;
	if(parameters.show_skipped_lines)
	{
		if(!(start_info == null))
		{
			//	Get the coordinates to write the skip information at the start of acid
			[x,y] = calibrateSkipStart(parameters) ;
			svgContainer.append("text")
				.attr("x", x)
				.attr("y", y)
				.style("fill","black")
				.text(start_info);
		}
		if(end_info != null)
		{
			//	Get the coordinates to write the skip information at the end of acid
			[x,y] = calibrateSkipEnd(parameters,last_position,first_position) ;
			svgContainer.append("text")
			.attr("x", x)
			.attr("y", y)
			.style("fill","black")
			.text(end_info) ;
		}
	}
}
/**
 * Put the numerical positions at the start and end of each row of the sequence
 * @param {*} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function getNumValues(para,prsm,id)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	//	Get the new first and last position based on the amount of acids are needed to be skipped
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	let svgContainer = d3.select("#"+id+"_g") ;
	prsm.annotated_protein.annotation.residue.forEach(function(input,index){
		if(parseInt(input.position) >= first_position && parseInt(input.position) < last_position )
		{
			let x,y ;
			l_position_temp = input.position ;
			//	write the numerical values only at the start position and end position(this is in the
			// 	form of 29,59 etc.,. as the data starts with 0 as 1st element)
			if(parseInt(l_position_temp)%(para.row_length) ==  0 || parseInt(l_position_temp)%(para.row_length)  ==  (para.row_length-1) 
									|| parseInt(l_position_temp) == (last_position-1))
			{
				let id_temp ;
				position = parseInt(input.position) +1;
				if(parseInt(input.position)%para.row_length ==  0)
				{
					//	Get the coordinates of left numerical
					[x,y] = calibrateLeftNum(para,parseInt(l_position_temp),first_position) ;
					x = x ;
					id_temp = "left_align" ;
				}
				else
				{
					//	Get the coordinates of right numerical
					[x,y] = calibrateRightNum(para,parseInt(l_position_temp),first_position) ;
					id_temp = "right_align" ;
				}
				svgContainer.append("text")
					.attr("id", id_temp)
					.attr("x",x)
					.attr("y",y)
					.text(function(d,i){
						return position ;
					})
					.style("text-anchor",function(d,i){
						//	Align the left numerical towards left side
						if(id_temp == "left_align")
						{
							return "end" ;
						}
						return null ;
					})
					.style("fill", "black");
				
				if(parseInt(l_position_temp) == (last_position-1) && last_position%(para.row_length) ==  1)
				{
					[x,y] = calibrateRightNum(para,parseInt(l_position_temp),first_position) ;
					id_temp = "right_align" ;
					svgContainer.append("text")
					.attr("id", id_temp)
					.attr("x",x)
					.attr("y",y)
					.text(function(d,i){
						return position ;
					})
					.style("fill", "black");
				}
			}
		}
	})
}
/**
 * Draw annotations
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function annotations(para,prsm,id)
{
	/*	Get annotation position	and information */
	var annotations = json2CleavagePositions(prsm) ;
	/* 	Drawing Annotation on to screen using D3 polyline	*/
	annotations.forEach(function(annotation,index){
		let l_charge = getIonCharge(annotations,annotation.position);
		
		if(annotation.exist_n_ion == "1" && annotation.exist_c_ion == "1")
		{
			drawAnnotation_YB(para,prsm,annotation,l_charge,id) ;
		}
		else
		{
			if(annotation.exist_n_ion == "1")
			{
				drawAnnotation_B(para,prsm,annotation,l_charge,id) ;
			} 
			if(annotation.exist_c_ion == "1")
			{
				drawAnnotation_Y(para,prsm,annotation,l_charge,id) ;
			}
		}
	})
}
/**
 * Invoke drawAnnotation method to draw the annotation when the ion type is B
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {Object} annotation - Json object with information of the position
 * @param {String} l_charge - Contains the Charge to be displayed on the annotation
 * @param {String} id - Contians id of the SVG tag from html.
 */
function drawAnnotation_B(para,prsm,annotation,l_charge,id)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	let x,y ;
	[x,y] = calibrateCoordinates(para,parseInt(annotation.position)-1,first_position);
	x = x + (para.letter_width/2) ;
	// Setting polyline coordinates
	let coordinates = (x-2)+","+(y-13)+ " " +(x+4)+","+ (y-11)+" "+(x+4)+","+(y+2);
	
	drawAnnotation(annotation,l_charge,id,coordinates,x,y);
}
/**
 * Invoke drawAnnotation method to draw the annotation when the ion type is Y
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {Object} annotation - Json object with information of the position
 * @param {String} l_charge - Contains the Charge to be displayed on the annotation
 * @param {String} id - Contians id of the SVG tag from html.
 */
function drawAnnotation_Y(para,prsm,annotation,l_charge,id)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	let x,y ;
	[x,y] = calibrateCoordinates(para,parseInt(annotation.position)-1,first_position);
	x = x + (para.letter_width/2) ;
	// Setting polyline coordinates
	let coordinates = (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(annotation,l_charge,id,coordinates,x,y);
	
}
/**
 * invoke drawAnnotation method to draw the annotation when the ion type is Y and B
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {Object} annotation - Json object with information of the position
 * @param {String} l_charge - Contains the Charge to be displayed on the annotation
 * @param {String} id - Contians id of the SVG tag from html.
 */
function drawAnnotation_YB(para,prsm,annotation,l_charge,id)
{
	let first_position,last_position,start_info = null ,end_info = null ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	let x,y ;
	[x,y] = calibrateCoordinates(para,parseInt(annotation.position)-1,first_position);
	x = x + (para.letter_width/2) ;
	// Setting polyline coordinates
	let coordinates =  (x-2)+","+(y-13)+ " " + (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(annotation,l_charge,id,coordinates,x,y);
}
/**
 * Function to draw the annotations based on annotation type and coordinates
 * @param {Object} annotation - Json object with information of the position
 * @param {String} l_charge - Contains the Charge to be displayed on the annotation
 * @param {String} id - Contians id of the SVG tag from html.
 * @param {String} coordinates - Contains String with coordinates to draw a polyline in a apecific format understandable by D3
 * @param {Integer} x - Contains x coordinate to draw a tooltip at specific position on hover of annotation
 * @param {Integer} y - Contains y coordinate to draw a tooltip at specific position on hover of annotation
 */
function drawAnnotation(annotation,l_charge,id,coordinates,x,y)
{
	let svgContainer = d3.select("#"+id+"_g");
	svgContainer.append("polyline")
				.attr("points", coordinates)
				.style("fill", "none")
				.style("stroke", "1e90ff")
				.style("stroke-width", "1");	
		/*	Rectangle to have flexible on click and on mouse actions	*/
	svgContainer.append("rect")
				.attr("id","annoTooltip")
				.attr("x", x)
				.attr("y", y-14)
				.attr("width", 13)
				.attr("height", 23)
				.style("opacity", 0)
/* 				.attr("cursor", "pointer")
				.on("click",function(){
					if(id == "l_svg")
					{
						input = annotation.ion_position;
						showIonPeaks(input);
					}
				})
				.on("mouseover", function(){
					appendTooltip(l_charge);
				})
				.on("mouseout", function(d){
					removeToolTip();	
				}); */
}
/**
 * Function to add tooltip to the polylines on mouseOver
 * @param {String} charge - Contains consolidated charge to display on tooltip
 */
function appendTooltip(charge)
{
	var div = d3.select("body").append("div")	
								.attr("class", "tooltip")				
								.style("opacity", 0); 
		div.transition()		
			.duration(10)		
			.style("opacity", .9);
		div.html(charge)	
		.style("left", (d3.event.pageX)  + "px")		
		.style("top", (d3.event.pageY - 28)+ "px") ;
}
/**
 * Function to remove tooltip to the polylines on mouseOver
 */
function removeToolTip()
{
	d3.selectAll(".tooltip").remove();
}
/**
 * Function to show the notification text at the top and bottom of the SVG sequence SVG of skipped amino acids
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 */
function skip_list(para,prsm)
{
	let l_afirst_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let l_alast_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;	
	let l_asequence_length = parseInt(prsm.annotated_protein.annotation.protein_length) ;
	let new_first_position = 0 ;
	let initial_skip_count = 0 ;
	let start_info = null ;
	if(l_afirst_residue_position > (parseInt(para.row_length) + l_afirst_residue_position%parseInt(para.row_length) ) )
	{
		new_first_position = l_afirst_residue_position - ((l_afirst_residue_position%para.row_length)+para.row_length );
		initial_skip_count = new_first_position ;
		start_info = "... "+initial_skip_count + " amino acid residues are skipped at the N-terminus ... ";
	}
	let new_last_position = l_alast_residue_position + 1 ;
	let final_skip_count = 0 ;
	let end_info = null ;
	if(l_alast_residue_position+(para.row_length - (l_alast_residue_position%para.row_length) +para.row_length) < l_asequence_length)
	{
		new_last_position = l_alast_residue_position+( para.row_length - (l_alast_residue_position%para.row_length) + para.row_length)  ;
		end_skip_count = l_asequence_length-new_last_position ;
		end_info = "... "+end_skip_count + " amino acid residues are skipped at the C-terminus ... ";
	}
	else if(l_alast_residue_position+1 < l_asequence_length)
	{
		new_last_position = l_asequence_length ;
	}
	return[para,new_first_position,new_last_position,start_info,end_info];
}
 /**
  * Draw the annotations to show the start and end position of the sequence
  * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
  * @param {Object} prsm - Contains the complete information of prsm. 
  * @param {String} id - Contians id of the SVG tag from html.
  */
function drawAnnoOfStartEndPosition(para,prsm,id)
{
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	let sequence_length = parseInt(prsm.annotated_protein.annotation.protein_length) ;
	let first_position,last_position,start_info = null ,end_info = null ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	let svgContainer = d3.select("#"+id+"_g");
	if(first_residue_position != 0)
	{
		let x,y ;
		[x,y] = calibrateCoordinates(para, first_residue_position-1, first_position );
		x = x + (para.letter_width/2) ;
		let first_position_coordinates = (x)+","+(y+2)+ " " +(x+5)+","+ (y+2)+" "+(x+5)+","+(y-12)+ " "+(x) + ","+(y-12);
		svgContainer.append("polyline")
					.attr("class","none")
					.attr("points", first_position_coordinates)
					.style("fill", "none")
					.style("stroke", "red")
					.style("stroke-width", "1.3") ;
	}
	if((sequence_length != last_residue_position +1 ) )
	{
		let x,y ;
		[x,y] = calibrateCoordinates(para, last_residue_position, first_position );
		x = x + (para.letter_width/2) ;
		let first_position_coordinates = (x+7)+","+(y-12)+ " " +(x+2)+","+ (y-12)+" "+(x+2)+","+(y+2)+ " "+(x+7) + ","+(y+2);
		svgContainer.append("polyline")
					.attr("class","none")
					.attr("points", first_position_coordinates)
					.style("fill", "none")
					.style("stroke", "red")
					.style("stroke-width", "1.3") ;
	}
}
 /**
  * Color the background of the occurence
  * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
  * @param {Object} prsm - Contains the complete information of prsm. 
  * @param {String} id - Contians id of the SVG tag from html.
  */
function massShiftBackgroundColor(para,prsm,id)
{
	let Background_color = json2BackgroundColorArray(prsm);
	let non_data_indicator = false ;
	let first_position,last_position,start_info = null ,end_info = null ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	Background_color.forEach(function(input,index){
		let left_position = parseInt(input.left_position) ;
		let right_position = parseInt(input.right_position) ;
		let annotation = input.anno ;
		let isShiftNeeded = shiftAnnotation(para,prsm,index) ;
		while(left_position < right_position)
		{
			let leftPosition, rightPosition ;
			[leftPosition,rightPosition] = getRightPosition(para,left_position,right_position) ;
			rightPosition = rightPosition - 1 ;
			let x,y;
			[x,y] = calibrateCoordinates(para, leftPosition ,first_position);
			let x1,y1;
			[x1,y1] = calibrateCoordinates(para, rightPosition ,first_position);
			let width = x1-x;
			rect_Backgroundcolor(x,y,id,width,para);
			MassShift(x,y,id,annotation,isShiftNeeded);
			left_position = parseInt(rightPosition) + 1 ;
			annotation = "";
			//break ;
		}
	})
	
}
/**
 * Get the right end position to color the background color when there is a change of row
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Integer} leftPosition - Contains start position of acid to add color to the background
 * @param {Integer} rightPosition - Contains end position of the acid to add color to the background
 */
function getRightPosition(para,leftPosition, rightPosition)
{
	if( parseInt(rightPosition/para.row_length) > parseInt(leftPosition/para.row_length) )
	{
		rightPosition = (parseInt(leftPosition/para.row_length)+1)*para.row_length ;
	}
	
	return [leftPosition, rightPosition] ;
}
/**
 * Modifying the color array to match the array of letters and positions
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function addColorToFixedPtms(prsm,id){
	let known_Change = json2FixedPtmOccurence(prsm);
	known_Change.forEach(function(position,i){
		l_class_id ="#id_"+id+"_g_"+position ;
		d3.select(l_class_id)
			.style("fill", "red")
	})
}
/**
 * Get the charge of the Ion by consilidating all the prefix and the charges at that position
 * @param {Array} l_annotation_array - Contains charges at all the positions of the acids
 * @param {Integer} position - Current position at which the chage has to be displayed on tooltip
 */
function getIonCharge(l_annotation_array,position){
	let l_charge = "";
	for(let j=0;j<l_annotation_array.length;j++){
		if(position == l_annotation_array[j].position )
		{
			 l_charge = l_charge + l_annotation_array[j].ion_type+l_annotation_array[j].ion_display_position
			 						+" "+l_annotation_array[j].peak_charge+"+ " ;
		}
	}
	return l_charge ;
}

/**
 * Code to color the background of a occurence acids
 * @param {Integer} x - Contains x coordinate of the start postion to add background color
 * @param {Integer} y - Contains y coordinate of the end position to add background color
 * @param {String} id - Contains id of the SVG tag from html.
 * @param {Integer} width - Contains width to whuich backgroud color has to be added
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 */
function rect_Backgroundcolor(x,y,id,width,para){
	/*	font-size 16px is equal to 12pt	*/
	let font_width = 12 ;
	/*	to draw the rect color uniformly */
	let font_height = 15 ;
	let svgContainer = d3.select("#"+id+"_g");
	svgContainer.append("rect")
					.attr("x", x)
					.attr("y", y-font_height)
					.attr("width", width+font_width)
					.attr("height", 20)
					.attr("dy", "0em")
					.style("fill", para.background_color)
					.style("fill-opacity", ".4")
					.style("stroke-width", "1.5px");
}
/**
 * MassShift value at the top of the acids
 * @param {Integer} x - Contains x coordinate at which the mass shift is added on top of acid 
 * @param {Integer} y - Contains y coordinate at which the mass shift is added on top of acid
 * @param {String} id - Contains id of the SVG tag from html.
 * @param {String} value - Contains Value to be displyed on the acid
 * @param {Boolean} isShift - Contains true if mass shifts adjacent to each other are overlapping
 */
function MassShift(x,y,id,value,isShift)
{
	let dy = -1.2 ;
	// If Adjacent mass shifts are overlapping, an additional shift is provided to not overlap
	if(isShift)
	{
		dy = 1.7*dy ;
	}
	dy = dy+"em" ;
	let svgContainer = d3.select("#"+id+"_g");
	svgContainer.append("text")
				   .attr("x", x)
				   .attr("y", y)
				   .attr("dy", dy) 
				   .text(function(){
					   return value ;
				   })
					.attr("fill","black")
				   .attr("font-size","15px");
}
/**
 * shift the position of the mass shift when overlapping one another
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {Integer} index - Current index form the list of annotations
 */
function shiftAnnotation(para,prsm,index)
{
	let isshiftNeeded = false ;
	let bgColorAndMassShift = json2BackgroundColorArray(prsm) ; 
	/*	font-size 16px is equal to 12pt*/
	let font_width = 12 ;
	let first_position, last_position,start_info,end_info ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	for(let i = 1; i<bgColorAndMassShift.length;i++)
	{
		if(i == index)
		{
			let x1,y1 ;
			[x1,y1] = calibrateCoordinates(para,parseInt(bgColorAndMassShift[i].left_position),first_position) ;
			
			let x2,y2 ;
			[x2,y2] = calibrateCoordinates(para,parseInt(bgColorAndMassShift[i-1].left_position),first_position) ;
			/* subtract -2 for CSS and alignment purpose*/
			if((Math.abs(x2-x1) + font_width )< bgColorAndMassShift[i-1].anno.length*(font_width-2) )
			{
				isshiftNeeded = true ;
			}
		}
	}
	return isshiftNeeded ;
}
/**
 * Check if over shifting is need when annotations collide with each other
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm.  
 */
function isShiftAnnotationNeeded(para,prsm)
{
	let isshiftNeeded = false ;
	let bgColorAndMassShift = json2BackgroundColorArray(prsm) ; 
	/*	font-size 16px is equal to 12pt*/
	let font_width = 12 ;
	let first_position, last_position,start_info,end_info ;
	[para,first_position, last_position,start_info,end_info] = skip_list(para,prsm);
	for(let i = 1; i<bgColorAndMassShift.length;i++)
	{
			let x1,y1 ;
			[x1,y1] = calibrateCoordinates(para,parseInt(bgColorAndMassShift[i].left_position),first_position) ;
			
			let x2,y2 ;
			[x2,y2] = calibrateCoordinates(para,parseInt(bgColorAndMassShift[i-1].left_position),first_position) ;
			/* subtract -2 for CSS and alignment purpose*/
			if((Math.abs(x2-x1) + font_width )< bgColorAndMassShift[i-1].anno.length*(font_width-2) )
			{
				isshiftNeeded = true ;
				break ;
			}
	}
	return isshiftNeeded ;
}
