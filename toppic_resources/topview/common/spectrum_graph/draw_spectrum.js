/**
 * Function gets invokes whenever zoom or drag is invoked and redraws the graph whenever there is zoom or draw
 * This function invokes all the functions that draw the graph
 * @param {string} svgId - Contains the Id from html on which the graph should be drawn
 * @param {object} para - Contains the parameters to help draw the graph
 * @param {list} peaks - Contains the data of Peaks to draw lines on the graph
 * @param {list} envPeaks - Contains the data of Envelope peaks to draw circles on the graph
 * @param {list} ions - Contains Ions to draw upon the peaks
 */
drawSpectrum  = function(svgId, para, peaks, envPeaks, ions, prefixSequence, suffixSequence, errorList) {
  let svg = d3.select("body").select("#"+svgId);
  // svg.attr("width", para.svgWidth).attr("height", para.svgHeight);
  // Removes all the elements under SVG group 'svgGroup' everytime there this function is called
  svg.selectAll("#svgGroup").remove();
  // Create a group under which all the fucntions of the graph will be appended
  svg = svg.append("g").attr("id","svgGroup");

  /*call onMouseOut everytime to fix onHover bug adding multiple data when mouseover and zoomed up*/
  onMouseOut();
  drawTicks(svg, para);
  drawAxis(svg,para);
  addDatatoAxis(svg,para);
  if (para.showHighlight) {
    addHighlight(svg, para);
  }
  drawPeaks(svg, para, peaks);
  if (para.showEnvelopes) {
    drawEnvelopes(svg, para, envPeaks);
  }
  if (para.showIons) {
    drawIons(svg, para, ions);
  }
  if (para.showSequence) {
    drawSequence(svg, para, prefixSequence, suffixSequence);
  }
  if(para.showErrorPlots)
	{
		addErrorPlot(svg, para);
		drawErrorYticks(svg,para);
		drawErrorPoints(svg,para, errorList);
		addErrorBlock(svg,para);
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
onMouseOverPeak = function(this_element,peak,para)
{
  let x = para.getPeakXPos(peak.mz);
  let y = para.getPeakYPos(peak.intensity);
  intensity =" inte:"+ parseFloat(peak.intensity).toFixed(3);
  mz = "m/z:"+parseFloat(peak.mz).toFixed(3);
  y = y - para.mouseOverPadding.head ;
  if(y<=para.mouseOverPadding.head)
  {
    y = para.mouseOverPadding.head;
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
onMouseOverCircle = function(this_element,envelope, peak, para)
{
  let x = parseFloat(d3.select(this_element).attr("cx"));
  let y = parseFloat(d3.select(this_element).attr("cy")) ;
  let mz = "m/z:"+peak.mz.toFixed(3);
  let inte = "inte:"+peak.intensity.toFixed(2);
  let mass = "mass:"+envelope.mono_mass.toFixed(3);
  let charge = "charge:"+ envelope.charge ;
  y = y - para.mouseOverPadding.head ;
  if(y<=para.mouseOverPadding.head)
  {
    y = para.mouseOverPadding.head;
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
drawTicks = function(svg,para){
  // Creating a group under svg node with id 'ticks' under which ticks are drawn 
  let addXTicks = svg.append("g").attr("id","x_ticks").attr("class", "ticks");
  let xTickPosList = para.getXTickPosList();
  for(let i=0; i < xTickPosList.length ; i++)
  {
    let tickMz = xTickPosList[i];
    // get the x position of the tick 
    x = para.getPeakXPos(tickMz);
    // Below condition helps the ticks to be to the right of the y - axis 
    if(x >= para.padding.left && 
      x <= (para.svgWidth - para.padding.right))
    {
      addXTicks.append("line")
        .attr("x1",x)
        .attr("y1",para.svgHeight -para.padding.bottom)
        .attr("x2",x)
        .attr("y2",para.svgHeight -para.padding.bottom + para.tickLength)
        .attr("stroke","black")
        .attr("stroke-width","1")
    }
  }
  addYTicks = svg.append("g").attr("id","y_ticks").attr("class","ticks");
  let tickHeight = para.getTickHeight();
  for(let i=0; i <= para.yTickNum ; i++)
  {
    // Get the default tick height and calculate the actual tick height position
    tickPos = i*tickHeight; 
    //* para.dataMaxInte /100;
    let y = parseFloat(para.getPeakYPos(tickHeight)) ;
    if(!isNaN(y) && y >= para.padding.head)//y >= para.padding.head helps the ticks to be in the length of Y axis
    {
      addYTicks.append("line")
        .attr("x1",para.padding.left)
        .attr("y1",y)
        .attr("x2",para.padding.left - para.tickLength)
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
drawAxis = function(svg,para){
  //Draw x-axis
  xAxis = svg.append("g").attr("id", "xaxis").append("line")
    .attr("x1",para.padding.left)
    .attr("y1",para.svgHeight -para.padding.bottom)
    .attr("x2",para.specWidth+para.padding.left)
    .attr("y2",para.svgHeight -para.padding.bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
  // Draw y-axis
  yAxis = svg.append("g").attr("id", "yaxis").append("line")
    .attr("x1",para.padding.left)
    .attr("y1",para.padding.head)
    .attr("x2",para.padding.left)
    .attr("y2",para.svgHeight -para.padding.bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
}

/**
 * @function addDatatoAxis
 * @description Function to add tick numbers on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
addDatatoAxis = function(svg,para){
  let maxMz = para.winMaxMz;
  let minMz = para.winMinMz ;
  // Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
  xAxisData = svg.append("g")
    .attr("id", "xAxisPoints");
  let xTickPosList = para.getXTickPosList();
  for(let i = 0 ; i < xTickPosList.length ; i++)
  {
    tickMz = xTickPosList[i];
    x = para.getPeakXPos(tickMz);
    if(x >= para.padding.left && 
      x <= (para.svgWidth - para.padding.right))
    {
      xAxisData.append("text").attr("id","xtext").attr("x",x)
        .attr("y",(para.svgHeight - para.padding.bottom + 20))// Dividing with 1.6 to set the position of the numbers under the ticks appropriately
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
  for(let i = 0 ; i <= para.yTickNum; i++)
  {
    let tickHeight = 0;
    // Get the default tick height and calculate the actual tick height position
    tickHeight = para.getTickHeight();
    let data = i*tickHeight ;
    if(data <= 1 && data != 0) data = data.toFixed(1);
    tickInt = i*tickHeight * para.dataMaxInte /100;
    y = parseFloat(para.getPeakYPos(tickInt)) ;

    if(!isNaN(y) && y >= para.padding.head)
    {
      this.yAxisData.append("text").attr("class","ytext").attr("x",para.padding.left - para.tickLength)
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
addHighlight = function(svg,para){
  let svg_temp = svg.append("g")
    .attr("id", "svg_bgColor");
  if(!((para.hlMinMz < para.winMinMz 
    && para.hlMaxMz < para.winMinMz) 
    || (para.hlMinMz > para.winMaxMz 
      && para.hlMaxMz > para.winMaxMz)))
  {
    let x = para.getPeakXPos(para.hlMinMz);
    let x2 = para.getPeakXPos(para.hlMaxMz);
    if(para.hlMinMz < para.winMinMz)
    {
      x = para.getPeakXPos(para.winMinMz);
    }
    if(para.hlMaxMz > para.winMaxMz)
    {
      x2 = para.getPeakXPos(para.winMaxMz);
    }
    //console.log(para.winMaxMz, x, x2);
    svg_temp.append("rect")
      .attr("x", x)
      .attr("y", para.padding.head)
      .attr("width", x2-x)
      .attr("height", function(){
        let y1 = para.svgHeight - para.padding.bottom;
        let y2 = para.padding.head;
        return y1-y2;
      })
      .style("fill", para.hlColor)
      .style("fill-opacity", ".4")
      .style("stroke-width", "1.5px");
  }
}

/**
 * @function drawPeaks
 * @description Function to draw peak lines on the graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
drawPeaks = function(svg,para,peakList){
  let peaks = svg.append("g")
    .attr("id", "peaks");
  var len = peakList.length;
  // limits provide current count of number of peaks drawn on graph per bin(range) 
  // so that we can limit tha peak count to peaksPerRange count
  let limits = new Array(para.binNum).fill(0);
  let binWidth = para.getBinWidth();
  for(let i =0;i<len;i++)
  {
    let peak = peakList[i];
    if(peak.mz >= para.winMinMz && peak.mz < para.winMaxMz)
    {
      let binIndex = Math.floor((peak.mz - para.winMinMz)/binWidth); 
      if (binIndex < para.binNum) 
      {
        limits[binIndex] = limits[binIndex]+1;
        if (limits[binIndex] <= para.peakNumPerBin) {
          peaks.append("line")
            .attr("x1",function(){
              return para.getPeakXPos(peak.mz);
            })
            .attr("y1",function(){
              let y = para.getPeakYPos(peak.intensity);
              if(y<=para.padding.head) return para.padding.head ;
              else return y ;
            })
            .attr("x2",function(){
              return para.getPeakXPos(peak.mz);
            })
            .attr("y2",para.svgHeight - para.padding.bottom )
            .attr("stroke","black")
            .attr("stroke-width","2")
            .on("mouseover",function(){
              onMouseOverPeak(this,peak,para);
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
drawEnvelopes = function(svg,para,envPeakList) {
  let circles = svg.append("g").attr("id", "circles");
  let minPercentage = 0.0;
  let maxIntensity = para.dataMaxInte ;
  // limits provide current count of number of peaks drawn on graph per bin(range)
  // so that we can limit tha peak count to circlesPerRange count
  let limits = new Array(para.binNum).fill(0);
  let binWidth = para.getBinWidth();
  for (let i = 0; i < envPeakList.length; i++) {
    let peak = envPeakList[i]; 
    let env = peak.env; 
    //console.log(env);
    let color = env.color;
    //Show only envelopes with minimum of 0.5%
    let percentInte = peak.intensity/maxIntensity * 100 ;
    if(peak.mz >= para.winMinMz && peak.mz < para.winMaxMz && percentInte >= minPercentage) 
    { 
      let binIndex = Math.floor((peak.mz - para.winMinMz)/binWidth); 
      if (binIndex < para.binNum) 
      {
        limits[binIndex] = limits[binIndex]+1;
        if (limits[binIndex] <= para.peakNumPerBin) 
        {
          circles.append("circle")
            .attr("id","circles")
            .attr("cx",function(){
              return para.getPeakXPos(peak.mz);
            })
            .attr("cy",function(){
              let cy = para.getPeakYPos(peak.intensity);
              if(cy < para.padding.head) return para.padding.head;
              else return cy ;
            })
            .attr("r",function(){
              return para.getCircleSize();
            })
            .style("fill","white")
            .style("opacity", "0.8")
            .style("stroke",color)
            .style("stroke-width","2")
            .on("mouseover",function(){
              onMouseOverCircle(this,env,peak,para);
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
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} ionData - Contians Ion list to display on the graph
 */
drawIons = function(svg,para,ions){
  let ionGroup = svg.append("g").attr("id", "graph_ions");
  // console.log(ions);
  for (let i = 0; i < ions.length; i++) {
    let ion = ions[i];
    let x = ion.mz;
    if(x >= para.winMinMz && x <= para.winMaxMz) {
      let xPos = para.getPeakXPos(x) + para.ionXShift;
      let yPos = para.getPeakYPos(ion.intensity) + para.ionYShift;
      // let color = ion.env.color;
      if (typeof ion.color !== "undefined") {
        let color = ion.color;
        ionGroup.append("text")
        .attr("id","graph_matched_ions")
        .attr("x", xPos)
        .attr("y", yPos) 
        .style("fill", color)
        .style("opacity", "0.8")
        .style("stroke-width","2")
        .text(ion.text);
      } else {
        ionGroup.append("text")
        .attr("id","graph_matched_ions")
        .attr("x", xPos)
        .attr("y", yPos) 
        .style("opacity", "0.8")
        .style("stroke-width","2")
        .text(ion.text);
      }
      
    }
  }
}

/**
 * @function drawSequence
 * @description Draw Sequence on spectrum graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
drawSequence = function(svg, para, prefixSequence, suffixSequence){
	let seqSvg = svg.append("g").attr("id", "graph_sequence");
  let x,y,text;
	// Draw | at 0 for prefix mass list
	x = para.getPeakXPos((0));
	y = para.padding.head-35;
	text = "|";
	// Add "|" At the start of the prefix sequence
	if(0 >= para.winMinMz && 0 <= para.winMaxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	prefixSequence.forEach(function(element,index,prefixSequence){
		if(element.mass >= para.winMinMz && element.mass <= para.winMaxMz)
		{
			let x1 = para.getPeakXPos(0);
			if(index !== 0) x1 = para.getPeakXPos(prefixSequence[index-1].mass);
			x2 = para.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = para.padding.head-35;
			drawAcids_innerFunc(seqSvg,x,y,element.acid);
			drawAcids_innerFunc(seqSvg,x2,y,"|");
		}
	})

	// Draw | at water mass for suffix mass list
	let massOfWater = 18.010564683704;
	x = para.getPeakXPos(massOfWater);// Mass of water=18.010564683704
	y = para.padding.head-15;
	text = "|";
	if(massOfWater >= para.winMinMz && massOfWater <= para.winMaxMz)
	{
		drawAcids_innerFunc(seqSvg,x,y,text);
	}
	suffixSequence.forEach(function(element,index,suffixSequence){
		if(element.mass >= para.winMinMz && element.mass <= para.winMaxMz)
		{
			let x1 = para.getPeakXPos(18.010564683704);// Mass of water=18.010564683704
			if(index != 0) x1 = para.getPeakXPos(suffixSequence[index-1].mass);
			x2 = para.getPeakXPos((element.mass));
			x = (x1+x2)/2;
			y = para.padding.head-15;
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
 * @function addErrorPlot
 * @description Add Error Plot to the MonoMass Spectrum
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */

addErrorPlot = function(svg, para){
	//Draw x-axis
	this.xAxis = svg.append("g").attr("id", "xaxis_errorplot").append("line")
					.attr("x1",para.errorplot_padding.left)
					.attr("y1",para.svgHeight - para.heightForErrorPlot/2 - para.errorplot_padding.bottom)
					.attr("x2",para.specWidth + para.errorplot_padding.left)
					.attr("y2",para.svgHeight - para.heightForErrorPlot/2 - para.errorplot_padding.bottom)
					.attr("stroke","black")
					.style("stroke-dasharray", ("5, 3"))
					.attr("stroke-width","1.5")
	// Draw y-axis
	this.yAxis = svg.append("g").attr("id", "yaxis_errorplot").append("line")
					.attr("x1",para.errorplot_padding.left)
					.attr("y1",para.svgHeight - para.errorplot_padding.bottom)
					.attr("x2",para.errorplot_padding.left)
					.attr("y2",para.svgHeight - para.heightForErrorPlot - para.errorplot_padding.bottom)
					.attr("stroke","black")
					.attr("stroke-width","1")
}

/**
 * @function addErrorBlock
 * @description Draw Error plot
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */

addErrorBlock = function(svg, para){
	let rectBlock = svg.append("g").attr("id", "rect_errorplot");
	rectBlock.append("line")
			.attr("x1",para.errorplot_padding.left)
			.attr("y1",para.svgHeight - para.heightForErrorPlot - para.errorplot_padding.bottom)
			.attr("x2",para.specWidth + para.errorplot_padding.left)
			.attr("y2",para.svgHeight - para.heightForErrorPlot - para.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
	rectBlock.append("line")
			.attr("x1",para.errorplot_padding.left)
			.attr("y1",para.svgHeight - para.errorplot_padding.bottom)
			.attr("x2",para.specWidth + para.errorplot_padding.left)
			.attr("y2",para.svgHeight - para.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
	rectBlock.append("line")
			.attr("x1",para.svgWidth - para.errorplot_padding.right)
			.attr("y1",para.svgHeight - para.errorplot_padding.bottom)
			.attr("x2",para.svgWidth - para.errorplot_padding.right)
			.attr("y2",para.svgHeight - para.heightForErrorPlot - para.errorplot_padding.bottom)
			.attr("stroke","black")
			.attr("stroke-width","1")
}

/**
 * @function drawErrorYticks
 * @description Draw Error plot y ticks
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */

drawErrorYticks = function(svg, para){

	let addYTicks = svg.append("g").attr("id","yErrorTicks")
									.attr("class","yErrorTicks");
	let tempTick = para.errorThreshHoldVal/para.errorYticks;
	//Draw tick at 0th position
	let y = para.getErrorYPos(0);
	inner_drawYTicks(y);
	inner_addErrorYTickValues(0,y);
	// Draw positive ticks above error x axis
	for(let i=1;i<=para.errorYticks;i++)
	{
		y = para.getErrorYPos(i*tempTick);
		inner_drawYTicks(y);
		inner_addErrorYTickValues(i*tempTick,y);
	}
	//Draw negative ticks below error x axis
	for(let i=1;i<=para.errorYticks;i++)
	{
		y = para.getErrorYPos(-(i*tempTick));
		inner_drawYTicks(y);
		inner_addErrorYTickValues(-(i*tempTick),y);
	}
	function inner_drawYTicks(y){
		if(!isNaN(y) && y >= para.padding.head)//y >= para.padding.head helps the ticks to be in the length of Y axis
		{
			addYTicks.append("line")
						.attr("x1",para.padding.left)
						.attr("y1",y)
						.attr("x2",para.padding.left - para.tickLength)
						.attr("y2",y)
						.attr("stroke","black")
						.attr("stroke-width","1")
		}
	}
	function inner_addErrorYTickValues(data,y){
		if(!isNaN(y) && y >= para.padding.head)
		{
			addYTicks.append("text").attr("class","ytext").attr("x",para.padding.left - para.tickLength)
						.attr("y",y)
						.attr("text-anchor","end")
						.attr("alignment-baseline","middle")
						.text(data)
						.style("font-size","14px")
		}
	}
}

/**
 * @function drawErrorPoints
 * @description Draw Error points on the error graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */

drawErrorPoints = function(svg, para, errorList){
	let circles = svg.append("g").attr("id", "error_circles");
	errorList.forEach((element)=>{
		if(parseFloat(element.theoretical_mass) > para.winMinMz && parseFloat(element.theoretical_mass) <= para.winMaxMz){
			circles.append("circle")
			.attr("class","error_circles")
			.attr("cx",function(d,i){
				return para.getPeakXPos(parseFloat(element.theoretical_mass));
			})
			.attr("cy",function(d,i){
				let cy = para.getErrorYPos(parseFloat(element.mass_error));
				if(cy < para.padding.head) return para.padding.head ;
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
 * @function addLabels
 * @description Function to add labels on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
/*
addLabels = function(svg, para){

	svg.append("text").attr("id","label")
						// -5 is added simply as buffer to place m/z on top of error rect plot
						.attr("transform","translate(" + (para.svgWidth/2) + "," + (para.svgHeight-para.graphFeatures.padding.bottom +para.graphFeatures.adjustableHeightVal - 5) + ")")
					.attr("fill","black")
					    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
					    .attr("font-size","16px")
					    .text("m/z");
	svg.append("text").attr("id","label")
					.attr("transform", "translate("+ para.padding.left/3 +","+(para.svgHeight/2+para.labelAdjustVal)+")rotate(-90)")
					.attr("fill","black")
					    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
					    .attr("font-size","16px")
					    .text("Intensity");
}
*/
