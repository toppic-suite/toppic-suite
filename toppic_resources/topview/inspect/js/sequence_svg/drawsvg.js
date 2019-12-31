/*	Get the size of the svg based on the no. of rows and row length and other parameter.,.*/
function getSvgSize(parameters,seq)
{
	let seqLen = seq.length;
	let num_of_rows = parseInt(seqLen/parameters.row_length) ;
	let no_of_blocks = parameters.row_length/parameters.block_length - 1 ;
	let width = parameters.letter_width * (parameters.row_length - 1)+ (no_of_blocks)*parameters.gap_width + parameters.right_margin + parameters.left_margin + 2*parameters.numerical_width;
	let height = parameters.row_height * num_of_rows + parameters.bottom_margin + parameters.top_margin ;
	return [width,height];
}
/* get the sequence on to svg */
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
							//.style("fill","red")
		})
		.on("mouseout",function(d,i){
			d3.select(this).style("cursor","default")
			let id = "massshift_" + i;
			d3.select("#"+id).attr("font-size","11px");
							//.style("fill","black")
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
					return "red" ;
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
function handleOnClick(d,i,id,seq,massShiftList,monoMassList){
	
	d3.selectAll("#tooltip_pop").remove() ;
	var div = d3.select("body").append("div")
	   .attr("class", "tooltip")
	   .attr("id","tooltip_pop")
	   .style("opacity", 1);
	
	div.transition()
     .duration(200)
     .style("opacity", .9);
	div.html(
	        '<input list="browsers" name="myBrowser" type="text" id= "mass_shift" />'+
			'<button id="ok" style = "none" type="button">ok</button>'
	        )
	        .style("left", (d3.event.pageX - 30) + "px")             
	        .style("top", (d3.event.pageY - 45) + "px");
	
	//i starts with 0
	let shiftPosition = i;
	
	d3.select("#ok").on("click",function(){
	let massShiftVal =  document.getElementById("mass_shift").value ;
		d3.select("#tooltip_pop").remove() ;
		massShiftVal = parseFloat(massShiftVal);
		if(!isNaN(massShiftVal)){
			let errorType = $(".error_dropdown").val();
			let errorVal = parseFloat($("#errorval").val().trim());
			let executionObj = new SeqOfExecution();
			executionObj.onClickMassShift = {mass:massShiftVal,position:shiftPosition}
			executionObj.onClickSequenceOfExecution(errorType,errorVal);
		}
	});
}

/* Put the numerical positions at the start and end of each row of the sequence	*/
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
				//x = x + para.anno_width ;
				//position =  position ; 
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
/* Draw annotations*/
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
/*	Draw the annotation when the ion type is B	*/
function drawAnnotation_B(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2) ;
	let coordinates = (x-2)+","+(y-13)+ " " +(x+4)+","+ (y-11)+" "+(x+4)+","+(y+2);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
/* Draw the annotation when the ion type is Y*/
function drawAnnotation_Y(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2);
	let coordinates = (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
function drawAnnotation_YB(para,position,charge,id)
{
	x = getX(para,position-1);
	y = getY(para,position-1);
	x = x + (para.letter_width/2) ;
	let coordinates =  (x-2)+","+(y-13)+ " " + (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
	drawAnnotation(position,charge,id,coordinates,x,y);
}
function drawAnnotation(position,charge,id,coordinates,x,y)
{
	let svgContainer = d3.select("#"+id+"_g");
	let l_polyline = svgContainer.append("polyline")
							.attr("id","polyline")
							.attr("points", coordinates)
							.style("fill", "none")
							.style("stroke", "1e90ff")
							.style("stroke-width", 1 );	
		/*	Rectangle to have flexible on click and on mouse actions	*/
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
function removeToolTip()
{
	d3.selectAll(".tooltip").remove();
}

/* Get the charge of the Ion */
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

/* MassShift value at the top of the acids */
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