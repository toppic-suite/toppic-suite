// maximum number of peaks and envelopes points to be renderedon to graph per each bin(range)
const circlesPerRange = 200;
const peaksPerRange = 200;
/**
 * Function draws the graph, binds zoom and drag function to the graph
 * @param {String} svgId - SVG id on which the graph needed to be drawn
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {list} peakData - contains peakList and envelope list 
 * @param {list} ionData - Contains data with mass and ACID name to plot on the graph
 */
SpectrumGraph = function(svgId,spectrumParameters,peakData, ionData){
	this.svg = d3.select("body").select(svgId);
  	this.id = svgId;
  	this.para = spectrumParameters;
	this.data = peakData;
	this.ionData = ionData;
	let graph = this;
	let tempid = svgId.split("#")[1];
	this.redraw = function(mono_mz,graphFeatures){
		this.para = compSpectrumParameters(this.data.peak_list, this.data.envelope_list, mono_mz);
		this.para.graphFeatures = graphFeatures;
		spectrumParameters = drawSpectrum(this.id, this.para, this.data,this.ionData);
		//Copying as a new variable than referencing. Referencing will change the properties of parent if child properties are changes
		correspondingSpecParams_g[tempid] = jQuery.extend(true, {}, spectrumParameters);
	}

	this.zoomed = function () {
		let transform = d3.event.transform;
		let distance = transform.x - spectrumParameters.specX;
		let ratio = transform.k / spectrumParameters.specScale;
		spectrumParameters.specX = transform.x;
		spectrumParameters.specScale = transform.k;
		let mousePos = d3.mouse(this);
		if(ratio == 1) {
			spectrumParameters.drag(distance);
		}
		else{
			spectrumParameters.zoom(mousePos[0], mousePos[1], ratio);
		}
		spectrumParameters = drawSpectrum(svgId, spectrumParameters, peakData, ionData);
		//Copying as a new variable than referencing. Referencing will change the properties of parent if child properties are changes
		correspondingSpecParams_g[tempid] = jQuery.extend(true, {}, spectrumParameters);
	}

	this.zoom = d3.zoom()
		.on("zoom", this.zoomed);
	this.svg.attr("viewBox", "0 0 "+ spectrumParameters.graphFeatures.svgWidth+" "+ spectrumParameters.graphFeatures.svgHeight)
					.attr("width", "100%")
					.attr("height", "100%")
					.call(this.zoom);
	this.svg.call(this.zoom.transform, d3.zoomIdentity);

	spectrumParameters = drawSpectrum(svgId, spectrumParameters, peakData, ionData); 
	//Copying as a new variable than referencing. Referencing will change the properties of parent if child properties are changes
	correspondingSpecParams_g[tempid] = jQuery.extend(true, {}, spectrumParameters);

	return graph;
}
/**
 * Function that draws ticks on x-axis and y-axis
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
drawTicks = function(svg,spectrumParameters){
	// Creating a group under svg node with id 'ticks' under which ticks are drawn 
	this.addXTicks = svg.append("g").attr("id","ticks");
	for(let i=0; i <= spectrumParameters.xTicks ; i++)
	{
		// Get the default tick width and calculate the actual tick width position based on the current minMz value on the xaxis
		let tickWidth = spectrumParameters.getTickWidth();
		if(tickWidth < 1 && tickWidth != 0)
		{
			tickWidth = (i*tickWidth + spectrumParameters.minMz) - parseFloat((i*tickWidth + spectrumParameters.minMz)%tickWidth) ;
		}
		else if(tickWidth != 0)
		{
			tickWidth = i*tickWidth + spectrumParameters.minMz - (i*tickWidth + spectrumParameters.minMz)%tickWidth ;
		}
		// get the x position of the tick 
		x = spectrumParameters.getPeakXPos(tickWidth);
		// Below condition helps the ticks to be to the right of the y - axis 
		if(x >= spectrumParameters.padding.left && 
				x <= (spectrumParameters.svgWidth - spectrumParameters.padding.right))
		{
			this.addXTicks.append("line")
							.attr("x1",x)
							.attr("y1",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
							.attr("x2",x)
							.attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom + spectrumParameters.ticklength)
							.attr("stroke","black")
							.attr("stroke-width","1")
		}
	}
	this.addYTicks = svg.append("g").attr("id","ticks")
									.attr("class","ticks");
	for(let i=0; i <= spectrumParameters.yTicks ; i++)
	{
		// Get the default tick height and calculate the actual tick height position
		let tickHeight = spectrumParameters.getTickHeight();
		tickHeight = i*tickHeight * spectrumParameters.dataMaxInte /100;
		tickHeight = parseFloat(spectrumParameters.getPeakYPos(tickHeight)) ;
		let y =  tickHeight;
		if(!isNaN(y) && y >= spectrumParameters.padding.head)//y >= spectrumParameters.padding.head helps the ticks to be in the length of Y axis
		{
			this.addYTicks.append("line")
						.attr("x1",spectrumParameters.padding.left)
						.attr("y1",y)
						.attr("x2",spectrumParameters.padding.left - spectrumParameters.ticklength)
						.attr("y2",y)
						.attr("stroke","black")
						.attr("stroke-width","1")
		}
	}
}
/**
 * Function to draw x-axis and y-axis
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
drawAxis = function(svg,spectrumParameters){
	//Draw x-axis
	this.xAxis = svg.append("g").attr("id", "xaxis").append("line")
					.attr("x1",spectrumParameters.padding.left)
					.attr("y1",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
					.attr("x2",spectrumParameters.specWidth+spectrumParameters.padding.left)
					.attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
					.attr("stroke","black")
					.attr("stroke-width","2")
	// Draw y-axis
	this.yAxis = svg.append("g").attr("id", "yaxis").append("line")
					.attr("x1",spectrumParameters.padding.left)
					.attr("y1",spectrumParameters.padding.head)
					.attr("x2",spectrumParameters.padding.left)
					.attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
					.attr("stroke","black")
					.attr("stroke-width","2")
}
/**
 * Function to add tick numbers on x and y axis
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
addDatatoAxis = function(svg,spectrumParameters){
	let maxMz = spectrumParameters.maxMz;
	let minMz = spectrumParameters.minMz ;
	// Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
	this.xAxisData = svg.append("g")
						.attr("id", "xAxisPoints");
						
	for(let i = 0 ; i <=spectrumParameters.xTicks ; i++)
	{
		// Get the default tick width and calculate the actual tick width position based on the current minMz value on the xaxis
		let tickWidth = spectrumParameters.getTickWidth();
		if(tickWidth < 1 && tickWidth != 0)
		{
			tickWidth = i*tickWidth + spectrumParameters.minMz - parseFloat((i*tickWidth + spectrumParameters.minMz)%tickWidth) ;
		}
		else if(tickWidth != 0)
		{
			tickWidth = i*tickWidth + spectrumParameters.minMz - (i*tickWidth + spectrumParameters.minMz)%tickWidth ;
		}
		x = spectrumParameters.getPeakXPos(tickWidth);
		let l_tickWidth = tickWidth;
		if(x >= spectrumParameters.padding.left && 
				x <= (spectrumParameters.svgWidth - spectrumParameters.padding.right))
		{
			this.xAxisData.append("text").attr("id","xtext").attr("x",x)
						.attr("y",(spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.padding.bottom + 20))// Dividing with 1.6 to set the position of the numbers under the ticks appropriately
						.attr("text-anchor","middle")
						.text(function(){
							// conditions to show more decimal values as we zoom in further and limit decimals when zoomed back
							if(maxMz - minMz <=0.0001) return parseFloat(l_tickWidth).toFixed(6);
							else if(maxMz - minMz <=0.001) return parseFloat(l_tickWidth).toFixed(5);
							else if(maxMz - minMz <=0.01) return parseFloat(l_tickWidth).toFixed(4);
							else if(maxMz - minMz <=0.1) return parseFloat(l_tickWidth).toFixed(3);
							else if(maxMz - minMz <=1) return parseFloat(l_tickWidth).toFixed(2);
							else if(maxMz - minMz <= 3) return parseFloat(l_tickWidth).toFixed(2)
							else if(maxMz - minMz <= 5) return parseFloat(l_tickWidth).toFixed(1)
							return parseInt(l_tickWidth)
						})
						.style("font-size","14px")
		}
	}
	// Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
	this.yAxisData = svg.append("g")
							.attr("id", "yAxisPoints");
	for(let i = 0 ; i <= spectrumParameters.yTicks ; i++)
	{
		let tickHeight = 0;
		// Get the default tick height and calculate the actual tick height position
		tickHeight = spectrumParameters.getTickHeight();
		let data = i*tickHeight ;
		if(data <= 1 && data != 0) data = data.toFixed(1);
		tickHeight = i*tickHeight * spectrumParameters.dataMaxInte /100;
		tickHeight = parseFloat(spectrumParameters.getPeakYPos(tickHeight)) ;

		let y =  tickHeight;
		if(!isNaN(y) && y >= spectrumParameters.padding.head)
		{
			this.yAxisData.append("text").attr("class","ytext").attr("x",spectrumParameters.padding.left - spectrumParameters.ticklength)
						.attr("y",y)
						.attr("text-anchor","end")
						.attr("alignment-baseline","middle")
						.text(data + "%")
						.style("font-size","14px")
		}
	}
}
/**
 * Function to add backGround color to the spectrum graph for MS1 spectrum at precursor mz
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
addBackGround = function(svg,spectrumParameters){
	let svg_temp = svg.append("g")
					.attr("id", "svg_bgColor");

	if(!((spectrumParameters.graphFeatures.bgMinMz < spectrumParameters.minMz && spectrumParameters.graphFeatures.bgMaxMz < spectrumParameters.minMz) 
		|| (spectrumParameters.graphFeatures.bgMinMz > spectrumParameters.maxMz && spectrumParameters.graphFeatures.bgMaxMz > spectrumParameters.maxMz)))
	{
		let x = spectrumParameters.getPeakXPos(spectrumParameters.graphFeatures.bgMinMz);
		let x2 = spectrumParameters.getPeakXPos(spectrumParameters.graphFeatures.bgMaxMz);
		if(spectrumParameters.graphFeatures.bgMinMz < spectrumParameters.minMz)
		{
			x = spectrumParameters.getPeakXPos(spectrumParameters.minMz);
		}
		if(spectrumParameters.graphFeatures.bgMaxMz > spectrumParameters.maxMz)
		{
			x2 = spectrumParameters.getPeakXPos(spectrumParameters.maxMz);
		}
		svg_temp.append("rect")
			.attr("x", x)
			.attr("y", spectrumParameters.padding.head)
			.attr("width", x2-x)
			.attr("height", function(){
				let y1 = spectrumParameters.svgHeight - spectrumParameters.padding.bottom;
				let y2 = spectrumParameters.padding.head;
				return y1-y2;
			})
			.style("fill", spectrumParameters.graphFeatures.bgColor)
			.style("fill-opacity", ".4")
			.style("stroke-width", "1.5px");
	}
}
/**
 * Function to draw peak lines on the graph
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
drawPeaks = function(svg,spectrumParameters,peakdata){
	let peaks = svg.append("g")
    				.attr("id", "peaks");
	  var len = peakdata.peak_list.length;
	  // limits provide current count of number of peaks drawn on graph per bin(range) 
	  // so that we can limit tha peak count to peaksPerRange count
	  let limits=[0,0,0,0,0,0,0,0];
	  for(let i =0;i<len;i++)
	  {
		let peak = peakdata.peak_list[i];
		let inLimit = false;
		for(let j=0;j<(spectrumParameters.ranges.length-1);j++)
		{
			if(peak.mz > spectrumParameters.ranges[j] && peak.mz < spectrumParameters.ranges[j+1])
			{
				limits[j] = limits[j]+1;
				if(limits[j] <= peaksPerRange) inLimit = true;
				break;
			}
		}
		if(peak.mz >= spectrumParameters.minMz && peak.mz <= spectrumParameters.maxMz && inLimit == true)
		{
			peaks.append("line")
		  .attr("x1",function(d,i){
			  return spectrumParameters.getPeakXPos(peak.mz);
			  })
		  .attr("y1",function(d,i){
				let y = spectrumParameters.getPeakYPos(peak.intensity);
				if(y<=spectrumParameters.padding.head) return spectrumParameters.padding.head ;
				else return y ;
			  })
		  .attr("x2",function(d,i){
			  	return spectrumParameters.getPeakXPos(peak.mz);
			  })
		  .attr("y2",spectrumParameters.svgHeight - spectrumParameters.padding.bottom )
		  .attr("stroke","black")
		  .attr("stroke-width","2")
		  .on("mouseover",function(d,i){
					onMouseOverPeak(this,peak,spectrumParameters);
				})
			.on("mouseout",function(d,i){
				onPeakMouseOut(this);
			});
		}
	  }
}
/**
 * Function to add circles for the envelope data
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
addCircles = function(svg,spectrumParameters,peakData){
	let circles = svg.append("g").attr("id", "circles");
	let minPercentage = 0.0;
	let maxIntensity = spectrumParameters.dataMaxInte ;
	// limits provide current count of number of peaks drawn on graph per bin(range) 
	// so that we can limit tha peak count to circlesPerRange count
	let limits=[0,0,0,0,0,0,0,0];
	peakData.envelope_list.forEach(function(envelope_list){
		envelope_list.env_peaks.forEach(function(env_peaks){
			//Show only envelopes with minimum of 0.5% 
			let percentInte = env_peaks.intensity/maxIntensity * 100 ;
			let inLimit = false;
			for(let i=0;i<(spectrumParameters.ranges.length-1);i++)
			{
				if(env_peaks.mz > spectrumParameters.ranges[i] && env_peaks.mz < spectrumParameters.ranges[i+1])
				{
					limits[i] = limits[i]+1;
					if(limits[i] <= circlesPerRange) inLimit = true;
					break;
				}
			}
			// Condition keeps the circles to the right of the y axis till the end of the x-axis
			if(env_peaks.mz > spectrumParameters.minMz && env_peaks.mz <= spectrumParameters.maxMz && percentInte >= minPercentage && inLimit == true)
			{
				circles.append("circle")
				.attr("id","circles")
				.attr("cx",function(d,i){
					return spectrumParameters.getPeakXPos(env_peaks.mz);
				})
				.attr("cy",function(d,i){
					let cy = spectrumParameters.getPeakYPos(env_peaks.intensity);
					if(cy < spectrumParameters.padding.head) return spectrumParameters.padding.head ;
					else return cy ;
				})
				.attr("r",function(d,i){
					return spectrumParameters.getCircleSize();
				})
				.style("fill","white")
				.style("opacity", "0.6")
				.style("stroke",envelope_list.color)
				.style("stroke-width","2")
				.on("mouseover",function(d,i){
					onMouseOverCircle(this,envelope_list,spectrumParameters);
				})
				.on("mouseout",function(d,i){
					onCircleMouseOut(this);
				});
			}
		})
	})
}
/**
 * Function to add IONS at the top of the peaks for each cluster of envelopes
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
drawIons = function(svg,spectrumParameters,ionData){
	let ions = svg.append("g").attr("id", "graph_ions");
	// Get the default tick width and calculate the actual tick width position based on the current minMz value on the xaxis
	let tickWidth = spectrumParameters.getTickWidth();
	ionData.forEach((element)=>{
		if(tickWidth <= spectrumParameters.graphFeatures.tickWidthThreshholdval)
		{
			element.forEach((innerElement)=>{
				placeIonOnGraph_innerFunc(innerElement);
			})
		}else{
			let innerElement = element[0];
			placeIonOnGraph_innerFunc(innerElement);
		}
	})

	// Inner function to draw ions on the graph at respective position
	function placeIonOnGraph_innerFunc(innerElement){
		if(innerElement.mz > spectrumParameters.minMz && innerElement.mz <= spectrumParameters.maxMz)
		{
			ions.append("text")
			.attr("id","graph_matched_ions")
			.attr("x",(spectrumParameters.getPeakXPos((innerElement.mz))-spectrumParameters.graphFeatures.adjustableIonPosition))
			.attr("y",function(d,i){
				let y = spectrumParameters.getPeakYPos(innerElement.intensity + (0.1*innerElement.intensity));// Adding 10% to get the Ions on the max Intensity Peak
				y = y - spectrumParameters.graphFeatures.fixedHeightOfIonAboveThePeak;
				if(y <= spectrumParameters.padding.head) return spectrumParameters.padding.head ;
				else return y ;
				})
			.style("fill","black")
			.style("opacity", "0.6")
			.style("stroke-width","2")
			.text(innerElement.ion);
		}
	}
}
/**
 * Function to add Acid names at the top of the graph divided by | symbol
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
drawSequence = function(svg,spectrumParameters){
	let seqSvg = svg.append("g").attr("id", "graph_sequence");
	let sequenceData = spectrumParameters.graphFeatures.prefixSequenceData;
	let x,y,text;
	// Draw | at 0 for prefix mass list
	x = spectrumParameters.getPeakXPos((0));
	y = spectrumParameters.padding.head-40;
	text = "|";
	// Add "|" At the start of the prefix sequence
	if(0 > spectrumParameters.minMz && 0 <= spectrumParameters.maxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	sequenceData.forEach(function(element,index,prefixSequenceData){
		if(element.mass > spectrumParameters.minMz && element.mass <= spectrumParameters.maxMz)
		{
			let x1 = spectrumParameters.getPeakXPos(0);
			if(index != 0) x1 = spectrumParameters.getPeakXPos(prefixSequenceData[index-1].mass);
			x2 = spectrumParameters.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = spectrumParameters.padding.head-40;
			drawAcids_innerFunc(seqSvg,x,y,element.acid);
			drawAcids_innerFunc(seqSvg,x2,y,"|");
		}
	})

	sequenceData = spectrumParameters.graphFeatures.suffixSequeceData;
	// Draw | at water mass for suffix mass list
	let massOfWater = 18.010564683704;
	x = spectrumParameters.getPeakXPos(massOfWater);// Mass of water=18.010564683704
	y = spectrumParameters.padding.head-20;
	text = "|";
	if(massOfWater > spectrumParameters.minMz && massOfWater <= spectrumParameters.maxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	sequenceData.forEach(function(element,index,suffixSequeceData){
		if(element.mass > spectrumParameters.minMz && element.mass <= spectrumParameters.maxMz)
		{
			let x1 = spectrumParameters.getPeakXPos(18.010564683704);// Mass of water=18.010564683704
			if(index != 0) x1 = spectrumParameters.getPeakXPos(suffixSequeceData[index-1].mass);
			x2 = spectrumParameters.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = spectrumParameters.padding.head-20;
			drawAcids_innerFunc(seqSvg,x,y,element.acid);
			drawAcids_innerFunc(seqSvg,x2,y,"|");
		}
	})
	function drawAcids_innerFunc(svgId,x,y,text){
		svgId.append("text")
			.attr("id","")
			.attr("x",x)
			.attr("y",y)
			.style("fill","black")
			.style("opacity", "0.6")
			.style("stroke-width","2")
			.text(text);
	}
}
/**
 * Add Error Plot to the MonoMass Spectrum
 */
addErrorPlot = function(svg, spectrumParameters){
	//Draw x-axis
	this.xAxis = svg.append("g").attr("id", "xaxis_errorplot").append("line")
					.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot/2 - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("x2",spectrumParameters.graphFeatures.specWidth + spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot/2 - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("stroke","black")
					.style("stroke-dasharray", ("5, 3"))
					.attr("stroke-width","2")
	// Draw y-axis
	this.yAxis = svg.append("g").attr("id", "yaxis_errorplot").append("line")
					.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("x2",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("stroke","black")
					.attr("stroke-width","2")
}
drawErrorYticks = function(svg, spectrumParameters){

	let addYTicks = svg.append("g").attr("id","yErrorTicks")
									.attr("class","yErrorTicks");
	let tempTick = spectrumParameters.graphFeatures.errorThreshHoldVal/spectrumParameters.errorYticks;
	//Draw tick at 0th position
	let y = spectrumParameters.getErrorYPos(0);
	inner_drawYTicks(y);
	inner_addErrorYTickValues(0,y);
	// Draw positive ticks above error x axis
	for(let i=1;i<=spectrumParameters.errorYticks;i++)
	{
		y = spectrumParameters.getErrorYPos(i*tempTick);
		inner_drawYTicks(y);
		inner_addErrorYTickValues(i*tempTick,y);
	}
	//Draw negative ticks below error x axis
	for(let i=1;i<=spectrumParameters.errorYticks;i++)
	{
		y = spectrumParameters.getErrorYPos(-(i*tempTick));
		inner_drawYTicks(y);
		inner_addErrorYTickValues(-(i*tempTick),y);
	}
	function inner_drawYTicks(y){
		if(!isNaN(y) && y >= spectrumParameters.padding.head)//y >= spectrumParameters.padding.head helps the ticks to be in the length of Y axis
		{
			addYTicks.append("line")
						.attr("x1",spectrumParameters.padding.left)
						.attr("y1",y)
						.attr("x2",spectrumParameters.padding.left - spectrumParameters.ticklength)
						.attr("y2",y)
						.attr("stroke","black")
						.attr("stroke-width","1")
		}
	}
	function inner_addErrorYTickValues(data,y){
		if(!isNaN(y) && y >= spectrumParameters.padding.head)
		{
			addYTicks.append("text").attr("class","ytext").attr("x",spectrumParameters.padding.left - spectrumParameters.ticklength)
						.attr("y",y)
						.attr("text-anchor","end")
						.attr("alignment-baseline","middle")
						.text(data)
						.style("font-size","14px")
		}
	}
}
drawErrorPoints = function(svg, spectrumParameters){
	let circles = svg.append("g").attr("id", "error_circles");
	spectrumParameters.graphFeatures.errorListData.forEach((element)=>{
		if(parseFloat(element.theoretical_mass) > spectrumParameters.minMz && parseFloat(element.theoretical_mass) <= spectrumParameters.maxMz){
			circles.append("circle")
			.attr("class","error_circles")
			.attr("cx",function(d,i){
				return spectrumParameters.getPeakXPos(parseFloat(element.theoretical_mass));
			})
			.attr("cy",function(d,i){
				let cy = spectrumParameters.getErrorYPos(parseFloat(element.mass_error));
				if(cy < spectrumParameters.padding.head) return spectrumParameters.padding.head ;
				else return cy ;
			})
			.attr("r",function(d,i){
				return 3;
			})
			.style("fill","black")
			.style("opacity", "1")
			//.style("stroke",envelope_list.color)
			.style("stroke-width","2");
		}
	})
}
/**
 * Function to add labels on x and y axis
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
addLabels = function(svg, spectrumParameters){

	svg.append("text").attr("id","label")
						.attr("transform","translate(" + (spectrumParameters.svgWidth/2) + "," + (spectrumParameters.svgHeight-spectrumParameters.graphFeatures.padding.bottom +spectrumParameters.graphFeatures.adjustableHeightVal) + ")")
					.attr("fill","black")
					    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
					    .attr("font-size","16px")
					    .text("m/z");
	svg.append("text").attr("id","label")
					.attr("transform", "translate("+ spectrumParameters.padding.left/3 +","+(spectrumParameters.svgHeight/2+spectrumParameters.labelAdjustVal)+")rotate(-90)")
					.attr("fill","black")
					    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
					    .attr("font-size","16px")
					    .text("Intensity");
}
/**
 * Function to show the data of Mass and Intensity on mouse over of peaks
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
onMouseOverPeak = function(this_element,peak,spectrumParameters)
{
	let x = spectrumParameters.getPeakXPos(peak.mz);
	let y = spectrumParameters.getPeakYPos(peak.intensity);
	intensity =" I:"+ parseFloat(peak.intensity).toFixed(3);
	mz = "M:"+parseFloat(peak.mz).toFixed(3);
	y = y - spectrumParameters.mouseOverPadding.head ;
	if(y<=spectrumParameters.mouseOverPadding.head)
	{
		y = spectrumParameters.mouseOverPadding.head;
	}
	d3.select(this_element).style("stroke","red")
							.style("stroke-width","2");
	
	let tooltipData = mz + "<br>" + intensity ;						
	/*	Rectangle to have flexible on click and on mouse actions	*/
	var div = d3.select("body").append("div")
							.attr("id", "MyTextMZIN")	
							.attr("class", "tooltip")	
							
	div.transition().duration(30)	
					.style("opacity", 2);
	div.html(tooltipData).style("left", (d3.event.pageX + 12)  + "px")		
				.style("top", (d3.event.pageY - 28)+ "px")
				.style("fill", "black");
}
/**
 * Function to show the data of Mass and Intensity on mouse over of circles
 * @param{node} svg -  is a html node on which the graph is being ploted
 */
onMouseOverCircle = function(this_element,envelope_list,spectrumParameters)
{
	let x = parseFloat(d3.select(this_element).attr("cx"));
	let y = parseFloat(d3.select(this_element).attr("cy")) ;
	let mass = "Mass:"+envelope_list.mono_mass.toFixed(2);
	let charge = "Charge:"+ envelope_list.charge ;
	y = y - spectrumParameters.mouseOverPadding.head ;
	if(y<=spectrumParameters.mouseOverPadding.head)
	{
		y = spectrumParameters.mouseOverPadding.head;
	}
	
	let tooltipData = mass + "<br>" + charge ;						
	/*	Rectangle to have flexible on click and on mouse actions	*/
	var div = d3.select("body").append("div")
							.attr("id", "MyTextMassCharge")	
							.attr("class", "tooltip")
							
	div.transition().duration(30)	
					.style("opacity", 2);
	div.html(tooltipData).style("left", (d3.event.pageX + 12)  + "px")		
				.style("top", (d3.event.pageY - 28)+ "px")
				.style("fill", "black");
}
/**
 * Function to reset to the original on mouse out of peaks
 */
onPeakMouseOut = function(this_element)
{
	onMouseOut();
	d3.select(this_element).style("stroke","black");
}
/**
 * Function to reset to the original on mouse out of circle
 */
onCircleMouseOut= function(){
	onMouseOut();
}
/**
 * Function to remove the tooltips on mouseout
 */
onMouseOut = function(){
	d3.selectAll("#MyTextMZIN").remove();
	d3.selectAll("#MyTextMassCharge").remove();
}
/**
 * Function gets invokes whenever zoom or drag is invoked and redraws the graph whenever there is zoom or draw
 * This function invokes all the functions that draw the graph
 * @param {string} svgId - Contains the Id from html on which the graph should be drawn
 * @param {object} spectrumParameters - Contains the parameters to help draw the graph
 * @param {list} peakData - Contains the data of Peaks and Envelopes to draws lines and circles on the graph
 * @param {list} ionData - Contains Ions to draw upon the peaks
 */
function drawSpectrum(svgId, spectrumParameters, peakData,ionData){
	let svg = d3.select("body").select(svgId);
	// Removes all the elements under SVG group 'svgGroup' everytime there this function is called
	svg.selectAll("#svgGroup").remove();
	// Create a group under which all the fucntions of the graph will be appended
	svg = svg.append("g").attr("id","svgGroup");
	
	/*call onMouseOut everytime to fix onHover bug adding multiple data when mouseover and zoomed up*/
	onMouseOut();
	drawTicks(svg, spectrumParameters);
	drawAxis(svg,spectrumParameters);
	addDatatoAxis(svg,spectrumParameters);
	addLabels(svg, spectrumParameters);
	// Condition if background color is needed to be added for specific graph
	if(spectrumParameters.graphFeatures.isAddbgColor)
	{
		addBackGround(svg, spectrumParameters);
	}
	drawPeaks(svg, spectrumParameters, peakData);
	if(spectrumParameters.graphFeatures.showCircles && peakData.envelope_list != null)
	{
		addCircles(svg,spectrumParameters,peakData);
	}
	if(spectrumParameters.graphFeatures.showIons && ionData != null)
	{
		drawIons(svg,spectrumParameters,ionData);
	}
	if(spectrumParameters.graphFeatures.showSequene)
	{
		drawSequence(svg,spectrumParameters);
	}
	if(spectrumParameters.graphFeatures.addErrorPlot)
	{
		addErrorPlot(svg,spectrumParameters);
		drawErrorYticks(svg,spectrumParameters);
		drawErrorPoints(svg,spectrumParameters);
	}
	return spectrumParameters;
}
