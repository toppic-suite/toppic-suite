/**
 * Get the size of the svg based on the no. of rows and row length and other parameter.,.
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {String} seq - Contains protein Sequence 
 */
function getSvgSize(parameters,seq)
{
	let seqLen = seq.length;
	let num_of_rows = parseInt(seqLen/parameters.row_length) ;
	let no_of_blocks = parameters.row_length/parameters.block_length - 1 ;
	let width = parameters.letter_width * (parameters.row_length - 1)+ (no_of_blocks)*parameters.gap_width + parameters.right_margin + parameters.left_margin + 2*parameters.numerical_width;
	let height = parameters.row_height * num_of_rows + parameters.bottom_margin + parameters.top_margin ;
	return [width,height];
}
/**
 * Draw the sequence on to svg
 * @param {Object} parameters - Contains parameters of width, letter space etc., to draw SVG
 * @param {String} seq - Contains the protein sequence
 * @param {String} id -Contians id of the SVG tag from html.
 * @param {Array} massShiftList - Contains list of all the mass shifts
 * @param {Array} monoMassList - Contains Mono Mass list data
 */
function buildSvg(parameters,seq,id,massShiftList,monoMassList)
{
	let massShiftListLen = massShiftList.length ;
	d3.selectAll(".residue").remove();
	let width,height ;
	[width,height] = getSvgSize(parameters,seq) ;
	/*create a group under svg with svgId_g*/
	let id_temp = id + "_g" ;
	
	let svgContainer = d3.select("#"+id)
						.attr("width", width)
						.attr("height", height)
						.attr("font-family","'FreeMono',Miltonian,monospace")
						.attr("font-size","16px")
						.style("fill", parameters.svgBackground_color)
	svgContainer = svgContainer.append("g")
								.attr("id",id_temp)
								.attr("class",id_temp);
	text = 	svgContainer.selectAll("text");
	let x,y ;
	/*	Get the x and y coordinates of the acid position	*/
	text.data(seq)
		.enter()
		.append("text")
		.attr("class","residue")
		.attr("x", function(d,i){ 
			x = getX(parameters,i) ;
			return x ;
		})
		.attr("y", function(d,i){ 
			y = getY(parameters,i) ;
			return y ;
		})
		.text(function(d,i){
			for(let k=0;k<massShiftListLen;k++)
			{
				if( i == massShiftList[k].position && massShiftList[k].mass != 0)
				{
					MassShift(this,massShiftList[k].mass,i);
					break;
				}
				else if(i == massShiftList[k].position && massShiftList[k].mass == 0)
				{
					MassShift(this,massShiftList[k].mass,i);
				}
			}
			return d ;
		})
		.on("mouseover",function(d,i){
			d3.select(this).style("cursor","pointer")
			let id = "massshift_" + i;
			d3.select("#"+id).attr("font-size","18px");
		})
		.on("mouseout",function(d,i){
			d3.select(this).style("cursor","default")
			let id = "massshift_" + i;
			d3.select("#"+id).attr("font-size","11px");
		})
		.on("click",function(d,i){
			handleOnClick(d,i,id,seq,massShiftList,monoMassList)
		})
		.style("fill",function(d,i){
			for(let k=0;k<massShiftListLen;k++)
			{
				if( i == massShiftList[k].position )
				{
					MassShift(this,massShiftList[k].mass,i);
					//return "red" ;
					break;
				}
				else if(i == massShiftList[k].position)
				{
					MassShift(this,massShiftList[k].mass,i);
				}
			}
			return "black" ;
		})
	return parameters;
}
/**
 * Handles on click actions. 
 * on click of any amino acid, provides a box to enter mass shift. 
 * On click of ok, will re calculate and redraws entire page.
 * @param {Char} d Current Amino Acid
 * @param {Integer} i Index or position of the amino acid
 * @param {String} id Contains Id of the SVG on which the sequence is drawn
 * @param {String} seq Sequence of the amino acid
 * @param {Array} massShiftList List of all the amino acids
 * @param {Array} monoMassList List of Mono Mass data
 */
function handleOnClick(d,i,id,seq,massShiftList,monoMassList){
	
	d3.selectAll("#tooltip_pop").remove() ;
	var div = d3.select("body").append("div")
	   .attr("class", "tooltip")
	   .attr("id","tooltip_pop")
	   .style("opacity", 1);
	let colorsDropdown = addColorsToDropdown();
	div.transition()
     .duration(200)
     .style("opacity", .9);
	div.html(
			'<input list="browsers" name="myBrowser" type="text" id= "mass_shift" />'+
			colorsDropdown +
			'<button id="ok" style = "none" type="button">ok</button>'
	        )
	        .style("left", (d3.event.pageX - 30) + "px")             
	        .style("top", (d3.event.pageY - 45) + "px");
	
	let shiftPosition = i;
	
	d3.select("#ok").on("click",function(){
	let massShiftVal =  document.getElementById("mass_shift").value ;
	let tooltipcolor =  document.getElementById("tooltip_color").value ;
		d3.select("#tooltip_pop").remove() ;
		massShiftVal = parseFloat(massShiftVal);
		if(!isNaN(massShiftVal)){
			let errorType = $(".error_dropdown").val();
			let errorVal = parseFloat($("#errorval").val().trim());
			let executionObj = new SeqOfExecution();
			executionObj.onClickMassShift = {mass:massShiftVal,position:shiftPosition,color:tooltipcolor}
			executionObj.onClickSequenceOfExecution(errorType,errorVal);
		}
	});
}

/**
 * Put the numerical positions at the start and end of each row of the sequence
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} seq Amino Acid Sequence
 * @param {String} id Id of the SVG from html
 */
function getNumValues(para,seq,id)
{
	//remove all the numbers if exist
	d3.selectAll(".numbers").remove();
	let svgContainer = d3.select("#"+id+"_g") ;
	let len = seq.length;
	for(let i = 0;i < len ;i++ )
	{
		let x,y ;
		let l_position_temp = i ;
		/*	write the numerical values only at the start position and end position(this is in the
		 * 	form of 29,59 etc.,. as the data starts with 0 as 1st element) 							*/
		if(l_position_temp%(para.row_length) ==  0 || l_position_temp%(para.row_length)  ==  (para.row_length-1) 
								|| l_position_temp == len-1)
		{
			let id_temp ;
			let position = i+1 ;
			if(i%para.row_length ==  0)
			{
				/*	Get the coordinates of left numerical	*/
				[x,y] = calibrateLeftNum(para,l_position_temp) ;
				x = x ;
				id_temp = "left_align" ;
			}
			else
			{
				/*	Get the coordinates of right numerical	*/
				[x,y] = calibrateRightNum(para,l_position_temp) ;
				id_temp = "right_align" ;
			}
			svgContainer.append("text")
				.attr("class","numbers")
				.attr("id", id_temp)
				.attr("x",x)
				.attr("y",y)
				.text(function(d,i){
					return position ;
				})
				.style("text-anchor",function(d,i){
					/*	Align the left numerical towards left side */
					if(id_temp == "left_align") return "end" ;
					return null ;
				})
				.style("fill", "black");
			
			// When there exist only one number in the last row	
			if(l_position_temp == len-1 && len%(para.row_length) ==  1)
			{
				[x,y] = calibrateRightNum(para,l_position_temp) ;
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
}
/**
 * Draw annotations
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Array} matchedPeaks Contains Matched List
 * @param {String} id Contains id of the SVG tag from html to draw sequence 
 */
function annotations(para,matchedPeaks,id)
{
	// remove all existing polylines with polyline id
	d3.selectAll("#annoTooltip").remove();
	d3.selectAll("#polyline").remove();
	matchedPeaks.forEach(function(matchedPeak,index){
		let position = matchedPeak.position ;
		let charge = getIonCharge(matchedPeaks,position);
		
		let utilFunctionsObj = new utilFunctions();
		let tempIon = matchedPeak.ion[0].toLowerCase();
		let ionType = utilFunctionsObj.getTerminus(tempIon);
		
		if("NTERMINUS" == ionType)
		{
			drawAnnotation_B(para,position,charge,id) ;
		}
		else
		{
			drawAnnotation_Y(para,position,charge,id) ;
		}
	})
}
/**
 * Draw the annotation when the ion type is B
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Integer} position Position of the amino acid
 * @param {String} charge Contains the charge at the current position to display on hover of the annotation
 * @param {String} id Contains id of the SVG tag from html to draw sequence 
 */
function drawAnnotation_B(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2) ;
	let coordinates = (x-2)+","+(y-13)+ " " +(x+4)+","+ (y-11)+" "+(x+4)+","+(y+2);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
/**
 * Draw the annotation when the ion type is Y
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Integer} position Position of the amino acid
 * @param {String} charge Contains the charge at the current position to display on hover of the annotation
 * @param {String} id Contains id of the SVG tag from html to draw sequence 
 */
function drawAnnotation_Y(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2);
	let coordinates = (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
/**
 * generate cooordinates to dray Y or B annotation
 * @param {Object} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Integer} position Position of the amino acid
 * @param {String} charge Contains the charge at the current position to display on hover of the annotation
 * @param {String} id Contains id of the SVG tag from html to draw sequence 
 */
function drawAnnotation_YB(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2) ;
	let coordinates =  (x-2)+","+(y-13)+ " " + (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
/**
 * Draw Annotations
 * @param {Integer} position Position of the amino acid
 * @param {String} charge Contains the charge at the current position to display on hover of the annotation
 * @param {String} id Contains id of the SVG tag from html to draw sequence 
 * @param {String} coordinates Contains coordinates to draw the annotaion
 * @param {Integer} x Contains x coordinate of start point to draw the annotation
 * @param {Integer} y Contains y coordinate of start point to draw the annotation
 */
function drawAnnotation(position,charge,id,coordinates,x,y)
{
	let svgContainer = d3.select("#"+id+"_g");
	let l_polyline = svgContainer.append("polyline")
							.attr("id","polyline")
							.attr("points", coordinates)
							.style("fill", "none")
							.style("stroke", "1e90ff")
							.style("stroke-width", 1 );	
	//	Rectangle to have flexible on click and on mouse actions
	svgContainer.append("rect")
				.attr("id","annoTooltip")
				.attr("x", x)
				.attr("y", y-14)
				.attr("width", 13)
				.attr("height", 23)
				.style("opacity", 0)
				.attr("cursor", "pointer")
				.on("click",function(){
					showIonPeaks(position);
				})
				.on("mouseover", function(){
					appendTooltip(charge);
				})
				.on("mouseout", function(d){
					removeToolTip();	
				});
}
/**
 * Append Charge to the tool tip
 * @param {String} charge Contains the charge at the current position to display on hover of the annotation
 */
function appendTooltip(charge)
{
	var div = d3.select("body").append("div")	
								.attr("class", "tooltip annotation_tooltip")				
								.style("opacity", 0); 
		div.transition()		
			.duration(10)		
			.style("opacity", .9);
		div.html(charge)	
		.style("left", (d3.event.pageX)  + "px")		
		.style("top", (d3.event.pageY - 28)+ "px") ;
}
/**
 * Remove tootltip on remove of the mouse from annotation
 */
function removeToolTip()
{
	d3.selectAll(".annotation_tooltip").remove();
}
/**
 * Get the charge of the Ion
 * @param {Array} matchedPeaks Array with charge data at all the matched positions
 * @param {Integer} position Position at which the charge data has to be consolidated to show on hover of an annotation
 */
function getIonCharge(matchedPeaks,position){
	let l_charge = "";
	for(let j=0;j<matchedPeaks.length;j++){
		if(position == matchedPeaks[j].position )
		{
			 l_charge = l_charge + matchedPeaks[j].ion
			 						+" "+matchedPeaks[j].charge+"+ " ;
		}
	}
	return l_charge ;
}

/**
 * MassShift value at the top of the acids
 * @param {Object} thisElem Current element of the amino acid
 * @param {Float} MassShift Value of the mass shift
 * @param {Integer} position Position of the amino acid
 */
function MassShift(thisElem,MassShift,position)
{
	let id = "massshift_" + position;
	d3.select("#"+id).remove();
	let dy = -1.2 ;
	let x = $(thisElem).attr("x");
	let y = $(thisElem).attr("y");
	dy = dy+"em" ;
	let svgContainer = d3.select("#seqsvg");
	svgContainer.append("text")
				.attr("class","massshift_class")
				.attr("id",id)
			  	.attr("x", x)
			   	.attr("y", y)
			   	.attr("dy", dy) 
			   	.text(function(){
					if(MassShift != 0)
					{
						return MassShift ;
					}
					return "";
			   	})
				//.attr("text-anchor","middle")
				.attr("fill","black")
			   .attr("font-size","11px");
}
/**
 * Provide a drop down to add different background colors to the mass shifted elements
 */
function addColorsToDropdown(){
	let colors = ["white","green","yellow","blue","red"]
	let startStatement = "<select id=\"tooltip_color\" style=\"background-color:"+colors[0]+"\">";
	let endStatement = "</select>";
	let stringyfyingHTML = startStatement;
	let len = colors.length;
	for(let i=0;i<len;i++)
	{
		stringyfyingHTML += drawRectagleWithColors(colors[i],i);
	}
	stringyfyingHTML = stringyfyingHTML + endStatement;
	console.log("stringyfyingHTML : ", stringyfyingHTML);
	return stringyfyingHTML;
}
/**
 * Add rectagular block with selected color
 * @param {String} color Backgroung color added tot he amino acid
 * @param {Integer} index Position at which the color needed to be added
 */
function drawRectagleWithColors(color,index){
	let option = "<option value=\""+color +"\"style=\"width:80px;height:5px;border:1px solid #000;background-color:"+color+"\">";
	if(index == 0)
	{
		option = "<option value=\""+color +"\"style=\"width:80px;height:5px;border:1px solid #000;background-color:"+color+"\""+"selected"+">";
	}

	//let div = "<div style=\"width:80px;height:5px;border:1px solid #000;background-color:"+color+"\"></div>";
	let finalOption = option+color+"</option>"; 
	return finalOption;
}