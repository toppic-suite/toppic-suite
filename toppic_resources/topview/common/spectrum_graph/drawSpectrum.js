/**
 * Function gets invokes whenever zoom or drag is invoked and redraws the graph whenever there is zoom or draw
 * This function invokes all the functions that draw the graph
 * @param {string} svgId - Contains the Id from html on which the graph should be drawn
 * @param {object} spectrumParameters - Contains the parameters to help draw the graph
 * @param {list} peakData - Contains the data of Peaks and Envelopes to draws lines and circles on the graph
 * @param {list} ionData - Contains Ions to draw upon the peaks
 */
drawSpectrum  = function(svgId, spectrumParameters, peaks, envPeaks, ions) {
  let svg = d3.select("body").select("#"+svgId);
  //svg.attr("width", spectrumParameters.svgWidth).attr("height", spectrumParameters.svgHeight);
  // Removes all the elements under SVG group 'svgGroup' everytime there this function is called
  svg.selectAll("#svgGroup").remove();
  // Create a group under which all the fucntions of the graph will be appended
  svg = svg.append("g").attr("id","svgGroup");

  /*call onMouseOut everytime to fix onHover bug adding multiple data when mouseover and zoomed up*/
  onMouseOut();
  drawTicks(svg, spectrumParameters);
  drawAxis(svg,spectrumParameters);
  addDatatoAxis(svg,spectrumParameters);
  if (spectrumParameters.showHighlight) 
  {
    addHighlight(svg, spectrumParameters);
  }
  drawPeaks(svg, spectrumParameters, peaks);
  if (spectrumParameters.showEnvelopes) {
    drawEnvelopes(svg, spectrumParameters, envPeaks);
  }
  if (spectrumParameters.showIons) 
  {
    drawIons(svg, spectrumParameters, ions);
  }
}

/**
 * @function onPeakMouseOut
 * @description Function to reset to the original on mouse out of peaks
 * @param {Node} this_element - is a html node. On mouse over generates tooltip based on the current peak
 */
onPeakMouseOut = function(this_element)
{
  this.onMouseOut();
  d3.select(this_element).style("stroke","black");
}
/**
 * @function onCircleMouseOut
 * @description Function to reset to the original on mouse out of peaks
 */
onCircleMouseOut = function(){
  this.onMouseOut();
}
/**
 * @function onMouseOut
 * @description Function to remove the tooltips on mouseout
 */
onMouseOut = function(){
  d3.selectAll("#MyTextMZIN").remove();
  d3.selectAll("#MyTextMassCharge").remove();
}

/**
 * @function onMouseOverPeak
 * @description Function to show the data of Mass and Intensity on mouse over of peaks
 * @param {Node} this_element -  is a html node. On mouse over generates tooltip based on the current peak
 * @param {Object} - Contains mz and intensity value of the current peak
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
onMouseOverPeak = function(this_element,peak,spectrumParameters)
{
  let x = spectrumParameters.getPeakXPos(peak.mz);
  let y = spectrumParameters.getPeakYPos(peak.intensity);
  intensity =" inte:"+ parseFloat(peak.intensity).toFixed(3);
  mz = "m/z:"+parseFloat(peak.mz).toFixed(3);
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
 * @function onMouseOverCircle
 * @description Function to show the data of Mass and Intensity on mouse over of circles
 * @param {Node} this_element - is a html node. On mouse over generates tooltip based on the current peak
 * @param {Array} envelope_list - Contains Envelope List
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
onMouseOverCircle = function(this_element,envelope, peak, spectrumParameters)
{
  let x = parseFloat(d3.select(this_element).attr("cx"));
  let y = parseFloat(d3.select(this_element).attr("cy")) ;
  let mz = "m/z:"+peak.mz.toFixed(3);
  let inte = "inte:"+peak.intensity.toFixed(2);
  let mass = "mass:"+envelope.mono_mass.toFixed(3);
  let charge = "charge:"+ envelope.charge ;
  y = y - spectrumParameters.mouseOverPadding.head ;
  if(y<=spectrumParameters.mouseOverPadding.head)
  {
    y = spectrumParameters.mouseOverPadding.head;
  }

  let tooltipData = mz + "<br>" + inte + "<br>" + mass + "<br>" + charge ;
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
 * @function drawTicks
 * @description Function that draws ticks on x-axis and y-axis
 * @param{HTMLBaseElement} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
drawTicks = function(svg,spectrumParameters){
  // Creating a group under svg node with id 'ticks' under which ticks are drawn 
  let addXTicks = svg.append("g").attr("id","x_ticks").attr("class", "ticks");
  let xTickPosList = spectrumParameters.getXTickPosList();
  for(let i=0; i < xTickPosList.length ; i++)
  {
    let tickMz = xTickPosList[i];
    // get the x position of the tick 
    x = spectrumParameters.getPeakXPos(tickMz);
    // Below condition helps the ticks to be to the right of the y - axis 
    if(x >= spectrumParameters.padding.left && 
      x <= (spectrumParameters.svgWidth - spectrumParameters.padding.right))
    {
      addXTicks.append("line")
        .attr("x1",x)
        .attr("y1",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
        .attr("x2",x)
        .attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom + spectrumParameters.tickLength)
        .attr("stroke","black")
        .attr("stroke-width","1")
    }
  }
  addYTicks = svg.append("g").attr("id","y_ticks").attr("class","ticks");
  let tickHeight = spectrumParameters.getTickHeight();
  for(let i=0; i <= spectrumParameters.yTickNum ; i++)
  {
    // Get the default tick height and calculate the actual tick height position
    tickPos = i*tickHeight; 
    //* spectrumParameters.dataMaxInte /100;
    let y = parseFloat(spectrumParameters.getPeakYPos(tickHeight)) ;
    if(!isNaN(y) && y >= spectrumParameters.padding.head)//y >= spectrumParameters.padding.head helps the ticks to be in the length of Y axis
    {
      addYTicks.append("line")
        .attr("x1",spectrumParameters.padding.left)
        .attr("y1",y)
        .attr("x2",spectrumParameters.padding.left - spectrumParameters.tickLength)
        .attr("y2",y)
        .attr("stroke","black")
        .attr("stroke-width","1")
    }
  }
}

/**
 * @function drawAxis
 * @description Function to draw x-axis and y-axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
drawAxis = function(svg,spectrumParameters){
  //Draw x-axis
  xAxis = svg.append("g").attr("id", "xaxis").append("line")
    .attr("x1",spectrumParameters.padding.left)
    .attr("y1",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
    .attr("x2",spectrumParameters.specWidth+spectrumParameters.padding.left)
    .attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
  // Draw y-axis
  yAxis = svg.append("g").attr("id", "yaxis").append("line")
    .attr("x1",spectrumParameters.padding.left)
    .attr("y1",spectrumParameters.padding.head)
    .attr("x2",spectrumParameters.padding.left)
    .attr("y2",spectrumParameters.svgHeight -spectrumParameters.padding.bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
}

/**
 * @function addDatatoAxis
 * @description Function to add tick numbers on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
addDatatoAxis = function(svg,spectrumParameters){
  let maxMz = spectrumParameters.winMaxMz;
  let minMz = spectrumParameters.winMinMz ;
  // Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
  xAxisData = svg.append("g")
    .attr("id", "xAxisPoints");
  let xTickPosList = spectrumParameters.getXTickPosList();
  for(let i = 0 ; i < xTickPosList.length ; i++)
  {
    tickMz = xTickPosList[i];
    x = spectrumParameters.getPeakXPos(tickMz);
    if(x >= spectrumParameters.padding.left && 
      x <= (spectrumParameters.svgWidth - spectrumParameters.padding.right))
    {
      xAxisData.append("text").attr("id","xtext").attr("x",x)
        .attr("y",(spectrumParameters.svgHeight - spectrumParameters.padding.bottom + 20))// Dividing with 1.6 to set the position of the numbers under the ticks appropriately
        .attr("text-anchor","middle")
        .text(function(){
          // conditions to show more decimal values as we zoom in further and limit decimals when zoomed back
          if(maxMz - minMz <=0.0001) return parseFloat(tickMz).toFixed(6);
          else if(maxMz - minMz <=0.001) return parseFloat(tickMz).toFixed(5);
          else if(maxMz - minMz <=0.01) return parseFloat(tickMz).toFixed(4);
          else if(maxMz - minMz <=0.1) return parseFloat(tickMz).toFixed(3);
          else if(maxMz - minMz <=1) return parseFloat(tickMz).toFixed(2);
          else if(maxMz - minMz <= 3) return parseFloat(tickMz).toFixed(2)
          else if(maxMz - minMz <= 5) return parseFloat(tickMz).toFixed(1)
          return parseInt(tickMz)
        })
        .style("font-size","14px")
    }
  }
  // Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
  this.yAxisData = svg.append("g")
    .attr("id", "yAxisPoints");
  for(let i = 0 ; i <= spectrumParameters.yTickNum; i++)
  {
    let tickHeight = 0;
    // Get the default tick height and calculate the actual tick height position
    tickHeight = spectrumParameters.getTickHeight();
    let data = i*tickHeight ;
    if(data <= 1 && data != 0) data = data.toFixed(1);
    tickInt = i*tickHeight * spectrumParameters.dataMaxInte /100;
    y = parseFloat(spectrumParameters.getPeakYPos(tickInt)) ;

    if(!isNaN(y) && y >= spectrumParameters.padding.head)
    {
      this.yAxisData.append("text").attr("class","ytext").attr("x",spectrumParameters.padding.left - spectrumParameters.tickLength)
        .attr("y",y)
        .attr("text-anchor","end")
        .attr("alignment-baseline","middle")
        .text(data + "%")
        .style("font-size","14px")
    }
  }
}

/**
 * @function addHighlight
 * @description Function to add backGround color to the spectrum graph for MS1 spectrum at precursor mz
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
addHighlight = function(svg,spectrumParameters){
  let svg_temp = svg.append("g")
    .attr("id", "svg_bgColor");
  if(!((spectrumParameters.hlMinMz < spectrumParameters.winMinMz 
    && spectrumParameters.hlMaxMz < spectrumParameters.winMinMz) 
    || (spectrumParameters.hlMinMz > spectrumParameters.winMaxMz 
      && spectrumParameters.hlMaxMz > spectrumParameters.winMaxMz)))
  {
    let x = spectrumParameters.getPeakXPos(spectrumParameters.hlMinMz);
    let x2 = spectrumParameters.getPeakXPos(spectrumParameters.hlMaxMz);
    if(spectrumParameters.hlMinMz < spectrumParameters.winMinMz)
    {
      x = spectrumParameters.getPeakXPos(spectrumParameters.winMinMz);
    }
    if(spectrumParameters.hlMaxMz > spectrumParameters.winMaxMz)
    {
      x2 = spectrumParameters.getPeakXPos(spectrumParameters.winMaxMz);
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
      .style("fill", spectrumParameters.hlColor)
      .style("fill-opacity", ".4")
      .style("stroke-width", "1.5px");
  }
}

/**
 * @function drawPeaks
 * @description Function to draw peak lines on the graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
drawPeaks = function(svg,spectrumParameters,peakList){
  let peaks = svg.append("g")
    .attr("id", "peaks");
  var len = peakList.length;
  // limits provide current count of number of peaks drawn on graph per bin(range) 
  // so that we can limit tha peak count to peaksPerRange count
  let limits = new Array(spectrumParameters.binNum).fill(0);
  let binWidth = spectrumParameters.getBinWidth();
  for(let i =0;i<len;i++)
  {
    let peak = peakList[i];
    if(peak.mz >= spectrumParameters.winMinMz && peak.mz < spectrumParameters.winMaxMz)
    {
      let binIndex = Math.floor((peak.mz - spectrumParameters.winMinMz)/binWidth); 
      if (binIndex < spectrumParameters.binNum) 
      {
        limits[binIndex] = limits[binIndex]+1;
        if (limits[binIndex] <= spectrumParameters.peakNumPerBin) {
          peaks.append("line")
            .attr("x1",function(){
              return spectrumParameters.getPeakXPos(peak.mz);
            })
            .attr("y1",function(){
              let y = spectrumParameters.getPeakYPos(peak.intensity);
              if(y<=spectrumParameters.padding.head) return spectrumParameters.padding.head ;
              else return y ;
            })
            .attr("x2",function(){
              return spectrumParameters.getPeakXPos(peak.mz);
            })
            .attr("y2",spectrumParameters.svgHeight - spectrumParameters.padding.bottom )
            .attr("stroke","black")
            .attr("stroke-width","2")
            .on("mouseover",function(){
              onMouseOverPeak(this,peak,spectrumParameters);
            })
            .on("mouseout",function(){
              onPeakMouseOut(this);
            });
        }
      }
    }
  }
}

/**
 * @function drawEnvelopes
 * @description Function to add circles for the envelope data
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
drawEnvelopes = function(svg,spectrumParameters,envPeakList) {
  let circles = svg.append("g").attr("id", "circles");
  let minPercentage = 0.0;
  let maxIntensity = spectrumParameters.dataMaxInte ;
  // limits provide current count of number of peaks drawn on graph per bin(range)
  // so that we can limit tha peak count to circlesPerRange count
  let limits = new Array(spectrumParameters.binNum).fill(0);
  let binWidth = spectrumParameters.getBinWidth();
  for (let i = 0; i < envPeakList.length; i++) {
    let peak = envPeakList[i]; 
    let env = peak.env; 
    let color = env.color;
    //Show only envelopes with minimum of 0.5%
    let percentInte = peak.intensity/maxIntensity * 100 ;
    if(peak.mz >= spectrumParameters.winMinMz && peak.mz < spectrumParameters.winMaxMz && percentInte >= minPercentage) 
    { 
      let binIndex = Math.floor((peak.mz - spectrumParameters.winMinMz)/binWidth); 
      if (binIndex < spectrumParameters.binNum) 
      {
        limits[binIndex] = limits[binIndex]+1;
        if (limits[binIndex] <= spectrumParameters.peakNumPerBin) 
        {
          circles.append("circle")
            .attr("id","circles")
            .attr("cx",function(){
              return spectrumParameters.getPeakXPos(peak.mz);
            })
            .attr("cy",function(){
              let cy = spectrumParameters.getPeakYPos(peak.intensity);
              if(cy < spectrumParameters.padding.head) return spectrumParameters.padding.head;
              else return cy ;
            })
            .attr("r",function(){
              return spectrumParameters.getCircleSize();
            })
            .style("fill","white")
            .style("opacity", "0.6")
            .style("stroke",color)
            .style("stroke-width","2")
            .on("mouseover",function(){
              onMouseOverCircle(this,env,peak,spectrumParameters);
            })
            .on("mouseout",function(){
              onCircleMouseOut(this);
            });
        }
      }
    }
  }
}

/**
 * @function drawIons
 * @description Function to add IONS at the top of the peaks for each cluster of envelopes
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} ionData - Contians Ion list to display on the graph
 */
drawIons = function(svg,spectrumParameters,ionData){
  let ions = svg.append("g").attr("id", "graph_ions");
  // Get the default tick width and calculate the actual tick width position based on the current minMz value on the xaxis
  let tickWidth = spectrumParameters.getTickWidth();
  ionData.forEach((element)=>{
    if(tickWidth <= spectrumParameters.tickWidthThreshhold)
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
        .attr("x",(spectrumParameters.getPeakXPos((innerElement.mz))-spectrumParameters.adjustableIonPosition))
        .attr("y",function(){
          let y = spectrumParameters.getPeakYPos(innerElement.intensity + (0.1*innerElement.intensity));// Adding 10% to get the Ions on the max Intensity Peak
          y = y - spectrumParameters.fixedHeightOfIonAboveThePeak;
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
 * @function drawSequence
 * @description Draw Sequence on spectrum graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
drawSequence = function(svg,spectrumParameters, sequence){
	let seqSvg = svg.append("g").attr("id", "graph_sequence");
	let x,y,text;
	// Draw | at 0 for prefix mass list
	x = spectrumParameters.getPeakXPos((0));
	y = spectrumParameters.padding.head-35;
	text = "|";
	// Add "|" At the start of the prefix sequence
	if(0 >= spectrumParameters.minMz && 0 <= spectrumParameters.maxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	sequence.forEach(function(element,index,sequenceData){
		if(element.mass >= spectrumParameters.minMz && element.mass <= spectrumParameters.maxMz)
		{
			let x1 = spectrumParameters.getPeakXPos(0);
			if(index != 0) x1 = spectrumParameters.getPeakXPos(prefixSequenceData[index-1].mass);
			x2 = spectrumParameters.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = spectrumParameters.padding.head-35;
			drawAcids_innerFunc(seqSvg,x,y,element.acid);
			drawAcids_innerFunc(seqSvg,x2,y,"|");
		}
	})

	sequenceData = spectrumParameters.graphFeatures.suffixSequeceData;
	// Draw | at water mass for suffix mass list
	let massOfWater = 18.010564683704;
	x = spectrumParameters.getPeakXPos(massOfWater);// Mass of water=18.010564683704
	y = spectrumParameters.padding.head-15;
	text = "|";
	if(massOfWater >= spectrumParameters.minMz && massOfWater <= spectrumParameters.maxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	sequenceData.forEach(function(element,index,suffixSequeceData){
		if(element.mass >= spectrumParameters.minMz && element.mass <= spectrumParameters.maxMz)
		{
			let x1 = spectrumParameters.getPeakXPos(18.010564683704);// Mass of water=18.010564683704
			if(index != 0) x1 = spectrumParameters.getPeakXPos(suffixSequeceData[index-1].mass);
			x2 = spectrumParameters.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = spectrumParameters.padding.head-15;
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
*/
/**
 * @function addErrorPlot
 * @description Add Error Plot to the MonoMass Spectrum
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
addErrorPlot = function(svg, spectrumParameters){
	//Draw x-axis
	this.xAxis = svg.append("g").attr("id", "xaxis_errorplot").append("line")
					.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot/2 - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("x2",spectrumParameters.graphFeatures.specWidth + spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot/2 - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("stroke","black")
					.style("stroke-dasharray", ("5, 3"))
					.attr("stroke-width","1.5")
	// Draw y-axis
	this.yAxis = svg.append("g").attr("id", "yaxis_errorplot").append("line")
					.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("x2",spectrumParameters.graphFeatures.errorplot_padding.left)
					.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot - spectrumParameters.graphFeatures.errorplot_padding.bottom)
					.attr("stroke","black")
					.attr("stroke-width","1")
}
*/
/**
 * @function addErrorBlock
 * @description Draw Error plot
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
addErrorBlock = function(svg, spectrumParameters){
	let rectBlock = svg.append("g").attr("id", "rect_errorplot");
	rectBlock.append("line")
			.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
			.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("x2",spectrumParameters.graphFeatures.specWidth + spectrumParameters.graphFeatures.errorplot_padding.left)
			.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
	rectBlock.append("line")
			.attr("x1",spectrumParameters.graphFeatures.errorplot_padding.left)
			.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("x2",spectrumParameters.graphFeatures.specWidth + spectrumParameters.graphFeatures.errorplot_padding.left)
			.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
	rectBlock.append("line")
			.attr("x1",spectrumParameters.graphFeatures.svgWidth - spectrumParameters.graphFeatures.errorplot_padding.right)
			.attr("y1",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("x2",spectrumParameters.graphFeatures.svgWidth - spectrumParameters.graphFeatures.errorplot_padding.right)
			.attr("y2",spectrumParameters.graphFeatures.svgHeight - spectrumParameters.graphFeatures.heightForErrorPlot - spectrumParameters.graphFeatures.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
}
*/
/**
 * @function drawErrorYticks
 * @description Draw Error plot y ticks
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
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
*/
/**
 * @function drawErrorPoints
 * @description Draw Error points on the error graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
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
*/
/**
 * @function addLabels
 * @description Function to add labels on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
addLabels = function(svg, spectrumParameters){

	svg.append("text").attr("id","label")
						// -5 is added simply as buffer to place m/z on top of error rect plot
						.attr("transform","translate(" + (spectrumParameters.svgWidth/2) + "," + (spectrumParameters.svgHeight-spectrumParameters.graphFeatures.padding.bottom +spectrumParameters.graphFeatures.adjustableHeightVal - 5) + ")")
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
*/
