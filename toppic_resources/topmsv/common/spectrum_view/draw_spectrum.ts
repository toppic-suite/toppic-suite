/**
 * Function gets invokes whenever zoom or drag is invoked and redraws the graph whenever there is zoom or draw
 * This function invokes all the functions that draw the graph
 * @param {string} svgId - Contains the Id from html on which the graph should be drawn
 * @param {object} para - Contains the parameters to help draw the graph
 * @param {list} peaks - Contains the data of Peaks to draw lines on the graph
 * @param {list} envPeaks - Contains the data of Envelope peaks to draw circles on the graph
 * @param {list} ions - Contains Ions to draw upon the peaks
 */
 function drawBasicSpectrum(svgId: string, para: SpectrumViewParameters, peaks: Peak[], 
  ions: MatchedIon[] | null) {
  let svg = d3.select("body").select("#"+svgId);
  // svg.attr("width", para.svgWidth).attr("height", para.svgHeight);
  // Removes all the elements under SVG group 'svgGroup' everytime there this function is called
  svg.selectAll("#svgGroup").remove();
  // Create a group under which all the fucntions of the graph will be appended
  //@ts-ignore
  svg = svg.append("g").attr("id","svgGroup");

  //call onMouseOut everytime to fix onHover bug adding multiple data when mouseover and zoomed up
  onMouseOut();
  drawTicks(svg, para);
  drawAxis(svg,para);
  addDatatoAxis(svg,para);
  addLabels(svg,para);
  if (para.getShowHighlight()) {
    addHighlight(svg, para);
  }
  drawPeaks(svg, para, peaks);
  if (para.getShowIons() && ions) {//only run when ion information is not undefined
    drawIons(svg, para, ions!);
  }
}

function drawRawSpectrum(svgId: string, para: SpectrumViewParameters, envPeaks: Envelope[]) {
  let svg = d3.select("body").select("#"+svgId).select("#svgGroup");
  if (para.getShowEnvelopes() && envPeaks != null) {
    drawEnvelopes(svg, para, envPeaks);
  }
}

function drawMonoMassSpectrum(svgId: string, para: SpectrumViewParameters, proteoform: Proteoform | null, nMasses: TheoMass[], cMasses: TheoMass[], 
  ions: MatchedIon[]) {
  let svg = d3.select("body").select("#"+svgId).select("#svgGroup");

  //check whether to draw the annotation lines
  if ($("#checkbox-anno-line").length) {//if the tab exists
    if ($("#checkbox-anno-line").is(":checked")) {
      para.setShowLines(true);
    }
    else {
      para.setShowLines(false);
    }
  }

  //updateViewBox(svgId, para.getSVGWidth(), para.getSVGHeight());
  drawSequence(svg, para, proteoform, nMasses, cMasses, ions);
  if (para.getShowError()) {
    addErrorPlot(svg, para);
    addErrorBox(svg, para);
    drawErrorYTicks(svg, para);
    drawErrorPoints(svg, para, ions);  
  }
}

/**
 * @function onPeakMouseOut
 * @description Function to reset to the original on mouse out of peaks
 * @param {Node} this_element - is a html node. 
 * On mouse over generates tooltip based on the current peak
 */
function onPeakMouseOut(this_element: any): void {
  onMouseOut();
  d3.select(this_element).style("stroke","black");
}
/**
 * @function onCircleMouseOut
 * @description Function to reset to the original on mouse out of peaks
 */
function onCircleMouseOut(): void{
  onMouseOut();
}

function onFragmentMassAndIonTypeMouseOut(this_element: any) {
  onMouseOut();
  d3.select(this_element).style("stroke","black")
    .style("stroke-width","1");
}

/**
 * @function onMouseOut
 * @description Function to remove the tooltips on mouseout
 */
function onMouseOut(){
  d3.selectAll("#MyTextMZIN").remove();
  d3.selectAll("#MyTextMassCharge").remove();
  d3.selectAll("#MyTextFragmentMass").remove();
}

/**
 * @function onMouseOverPeak
 * @description Function to show the data of Mass and Intensity on mouse over of peaks
 * @param {Node} this_element -  is a html node. On mouse over generates tooltip based on the current peak
 * @param {Object} - Contains mz and intensity value of the current peak
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
 function onMouseOverPeak(this_element: any,peak: Peak, para: SpectrumViewParameters) {
  let intensity: string = " inte:"+ peak.getIntensity().toFixed(3);
  let pos: string = peak.getPos().toFixed(3);
  if (para.getIsMonoMassGraph()) {
    pos = "mass:" + pos;
  }
  else {
    pos = "m/z:"+ pos;
  }

  d3.select(this_element).style("stroke","red")
    .style("stroke-width","2");
  let tooltipData = pos + "<br>" + intensity ;
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
function onMouseOverCircle(this_element: any, envelope: Envelope, peak: Peak) {
  let mz: string = "m/z:"+peak.getPos().toFixed(3);
  let inte: string = "inte:"+peak.getIntensity().toFixed(2);
  let mass: string = "mass:"+envelope.getMonoMass().toFixed(3);
  let charge: string = "charge:"+ envelope.getCharge() ;
  let tooltipData: string = mz + "<br>" + inte + "<br>" + mass + "<br>" + charge ;
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
 * @function onMouseOverFragmentMassAndIonType
 * @description Function to show the theoretical mass and matched ion type on mouse over of peaks
 */
 function onMouseOverFragmentMassAndIonType(this_element: any, mass: number, ionData: string | null) {
  let pos: string = mass.toFixed(3);
  let tooltipData: string = "mass: " + pos + ", " + "ion type: " + ionData;
  if (ionData == null){
    tooltipData = "mass: " + pos;
  }

  d3.select(this_element).style("stroke","red")
    .style("stroke-width","2");

  /*	Rectangle to have flexible on click and on mouse actions	*/
  var div = d3.select("body").append("div")
    .attr("id", "MyTextFragmentMass")
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
function drawTicks(svg: any, para: SpectrumViewParameters){
  // Creating a group under svg node with id 'ticks' under which ticks are drawn 
  let addXTicks = svg.append("g").attr("id","x_ticks").attr("class", "ticks");
  let xTickPosList: number[] = para.getXTickPosList();
  for(let i=0; i < xTickPosList.length ; i++) {
    let tickMz: number = xTickPosList[i];
    // get the x position of the tick 
    let x: number = para.getPeakXPos(tickMz);
    // Below condition helps the ticks to be to the right of the y - axis 
    if(x >= para.getPadding().left && 
      x <= (para.getSVGWidth() - para.getPadding().right))
    {
      addXTicks.append("line")
        .attr("x1",x)
        .attr("y1",para.getSVGHeight() -para.getPadding().bottom)
        .attr("x2",x)
        .attr("y2",para.getSVGHeight() -para.getPadding().bottom + para.getTickLength())
        .attr("stroke","black")
        .attr("stroke-width","1")
    }
  }
  let addYTicks = svg.append("g").attr("id","y_ticks").attr("class","ticks");
  let tickHeight: number = para.getTickHeight();
  for(let i=0; i <= para.getYTickNum() ; i++)
  {
    // Get the default tick height and calculate the actual tick height position
    let tickPos: number = i*tickHeight; 
    //* para.dataMaxInte /100;
    let y: number = para.getPeakYPos(tickHeight) ;
    if(!isNaN(y) && y >= para.getPadding().head)//y >= para.padding.head helps the ticks to be in the length of Y axis
    {
      addYTicks.append("line")
        .attr("x1",para.getPadding().left)
        .attr("y1",y)
        .attr("x2",para.getPadding().left - para.getTickLength())
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
function drawAxis(svg: any, para: SpectrumViewParameters){
  //Draw x-axis
  let xAxis = svg.append("g").attr("id", "xaxis").append("line")
    .attr("x1",para.getPadding().left)
    .attr("y1",para.getSVGHeight() -para.getPadding().bottom)
    .attr("x2",para.getSVGWidth() +para.getPadding().left)
    .attr("y2",para.getSVGHeight() -para.getPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
  // Draw y-axis
  let yAxis = svg.append("g").attr("id", "yaxis").append("line")
    .attr("x1",para.getPadding().left)
    .attr("y1",para.getPadding().head)
    .attr("x2",para.getPadding().left)
    .attr("y2",para.getSVGHeight() -para.getPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","2")
}

/**
 * @function addDatatoAxis
 * @description Function to add tick numbers on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function addDatatoAxis(svg: any, para: SpectrumViewParameters){
  let maxMz: number = para.getWinMaxMz();
  let minMz: number = para.getWinMinMz();
  // Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
  let xAxisData = svg.append("g")
    .attr("id", "xAxisPoints");
  let xTickPosList: number[] = para.getXTickPosList();
  let lastXTickPos: number = -1;
  for(let i = 0 ; i < xTickPosList.length ; i++)
  {
    let tickMz: number = xTickPosList[i];
    let x: number = para.getPeakXPos(tickMz);
    if(x >= para.getPadding().left && 
      x <= (para.getSVGWidth() - para.getPadding().right))
    {
      lastXTickPos = x;
      xAxisData.append("text").attr("id","xtext").attr("x",x)
        .attr("y",(para.getSVGHeight() - para.getPadding().bottom + 20))// Dividing with 1.6 to set the position of the numbers under the ticks appropriately
        .attr("text-anchor","middle")
        .text(function(){
          // conditions to show more decimal values as we zoom in further and limit decimals when zoomed back
          if(maxMz - minMz <=0.0001) return tickMz.toFixed(6);
          else if(maxMz - minMz <=0.001) return tickMz.toFixed(5);
          else if(maxMz - minMz <=0.01) return tickMz.toFixed(4);
          else if(maxMz - minMz <=0.1) return tickMz.toFixed(3);
          else if(maxMz - minMz <=1) return tickMz.toFixed(2);
          else if(maxMz - minMz <= 3) return tickMz.toFixed(2)
          else if(maxMz - minMz <= 5) return tickMz.toFixed(1)
          return tickMz;
        })
        .style("font-size","14px")
    }
  }
  // Creating a group wih id 'axisPoints' under which the code to add tick numbers is added  
  let yAxisData = svg.append("g")
    .attr("id", "yAxisPoints");
  for(let i = 0 ; i <= para.getYTickNum(); i++)
  {
    let tickHeight: number = 0;
    // Get the default tick height and calculate the actual tick height position
    tickHeight = para.getTickHeight();
    let data: number = i*tickHeight ;

    if (data <= 1 && data != 0) {
      data = parseFloat(data.toFixed(3));
    }

    let tickInt: number = i*tickHeight * para.getDataMaxInte() /100;
    let y: number = para.getPeakYPos(tickInt);
    if(!isNaN(y) && y >= para.getPadding().head)
    {
      yAxisData.append("text").attr("class","ytext").attr("x",para.getPadding().left - para.getTickLength())
        .attr("y",y)
        .attr("text-anchor","end")
        .attr("alignment-baseline","middle")
        .text(data + "%")
        .style("font-size","14px")
    }
  }
}

/**
 * @function addLabels
 * @description Function to add labels on x and y axis
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function addLabels(svg: any, para: SpectrumViewParameters){
  let text: string = "M/Z";
  if (para.getIsMonoMassGraph()) {
    text = "Mass";
  }
  svg.append("text").attr("id","label")
  // -5 is added simply as buffer to place m/z on top of error rect plot
    .attr("transform","translate(" + (para.getSVGWidth()-40) + "," 
      + (para.getSVGHeight() - para.getPadding().bottom + 20) + ")")
    .attr("fill","black")
    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
    .attr("font-size","16px")
    .text(text);
  svg.append("text").attr("id","label")
    .attr("transform", "translate("+ para.getPadding().left/3 +","+(para.getSVGHeight()/2)+")rotate(-90)")
    .attr("fill","black")
    .attr("font-family","Helvetica Neue,Helvetica,Arial,sans-serif")
    .attr("font-size","16px")
    .text("Intensity");
}

/**
 * @function addHighlight
 * @description Function to add backGround color to the spectrum graph for MS1 spectrum at precursor mz
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function addHighlight(svg: any, para: SpectrumViewParameters){
  let svg_temp = svg.append("g")
    .attr("id", "svg_bgColor");
  if(!((para.getHlMinMz() < para.getWinMinMz() 
    && para.getHlMaxMz() < para.getWinMinMz()) 
    || (para.getHlMinMz() > para.getWinMaxMz() 
      && para.getHlMaxMz() > para.getWinMaxMz())))
  {
    let x : number= para.getPeakXPos(para.getHlMinMz());
    let x2: number = para.getPeakXPos(para.getHlMaxMz());
    if(para.getHlMinMz() < para.getWinMinMz())
    {
      x = para.getPeakXPos(para.getWinMinMz());
    }
    if(para.getHlMaxMz() > para.getWinMaxMz())
    {
      x2 = para.getPeakXPos(para.getWinMaxMz());
    }
    //console.log(para.winMaxMz, x, x2);
    svg_temp.append("rect")
      .attr("x", x)
      .attr("y", para.getPadding().head)
      .attr("width", x2-x)
      .attr("height", function(){
        let y1 = para.getSVGHeight() - para.getPadding().bottom;
        let y2 = para.getPadding().head;
        return y1-y2;
      })
      .style("fill", para.getHlColor())
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
function drawPeaks(svg: any, para: SpectrumViewParameters, peakList: Peak[]){
  let peaks = svg.append("g")
    .attr("id", "peaks");
  var len: number = peakList.length;
  // limits provide current count of number of peaks drawn on graph per bin(range) 
  // so that we can limit tha peak count to peaksPerRange count
  let ratio: number = (para.getWinMaxMz() - para.getWinMinMz()) / (para.getDataMaxMz() - para.getDataMinMz());
  ratio = Math.min(1, ratio);
  let spectrumData = new SpectrumFunction();

  for(let i =0;i<len;i++)
  {
    let peak: Peak = peakList[i];

    if(peak.getPos() >= para.getWinMinMz() && peak.getPos()  < para.getWinMaxMz())
    {
      if (peak.getDisplayLevel() / (spectrumData.getMzLevel().length - 3)>= ratio || ratio <= 0.2){
        peaks.append("line")
        .attr("x1",function(){
          return para.getPeakXPos(peak.getPos() );
        })
        .attr("y1",function(){
          let y = para.getPeakYPos(peak.getIntensity());
          if(y<=para.getPadding().head) return para.getPadding().head ;
          else return y ;
        })
        .attr("x2",function(){
          return para.getPeakXPos(peak.getPos());
        })
        .attr("y2",para.getSVGHeight() - para.getPadding().bottom )
        .attr("stroke","black")
        .attr("stroke-width","2")
        .on("mouseover",function(){
          //@ts-ignore - allow using this to pass interacted html node
          onMouseOverPeak(this,peak,para);
        })
        .on("mouseout",function(){
          //@ts-ignore
          onPeakMouseOut(this);
        });
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
function drawEnvelopes(svg: any, para: SpectrumViewParameters,envList: Envelope[]) {
  let circles = svg.append("g").attr("id", "circles");
  let minPercentage: number = 0.0;
  let maxIntensity: number = para.getDataMaxInte();
  let spectrumData = new SpectrumFunction();
  // limits provide current count of number of peaks drawn on graph per bin(range)
  // so that we can limit tha peak count to circlesPerRange count
  let ratio = (para.getWinMaxMz() - para.getWinMinMz()) / (para.getDataMaxMz() - para.getDataMinMz());
  ratio = Math.min(1, ratio);

  envList.forEach(env => {
    let peaks = env.getPeaks(); 
    let color = env.getDisplayColor();
    
    if(peaks[0].getPos() >= para.getWinMinMz() && peaks[0].getPos() < para.getWinMaxMz()) 
    { 
      if (env.getDisplayLevel() / (spectrumData.getMzLevel().length - 3) >= ratio || ratio <= 0.2){
        //display envelopes based on level, but when the ratio falls below threshold, show all envs in the range
        peaks.forEach(peak => {
          let percentInte = peak.getIntensity()/maxIntensity * 100 ;
          if (percentInte >= minPercentage){//Show only envelopes with minimum of 0.5%
            circles.append("circle")
            .attr("id","circles")
            .attr("cx",function(){
              return para.getPeakXPos(peak.getPos());
            })
            .attr("cy",function(){
              let cy = para.getPeakYPos(peak.getIntensity());
              if(cy < para.getPadding().head) return para.getPadding().head;
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
              //@ts-ignore
              onMouseOverCircle(this,env,peak);
            })
            .on("mouseout",function(){
              //@ts-ignore
              onCircleMouseOut(this);
            });
          }
        })
      }
    }
  })
}
/*function drawEnvelopes(svg,para,envPeakList) {
  let circles = svg.append("g").attr("id", "circles");
  let minPercentage = 0.0;
  let maxIntensity = para.dataMaxInte ;
  let spectrumData = new SpectrumFunction();

  // limits provide current count of number of peaks drawn on graph per bin(range)
  // so that we can limit tha peak count to circlesPerRange count
  let ratio = (para.winMaxMz - para.winMinMz) / (para.dataMaxMz - para.dataMinMz);
  ratio = Math.min(1, ratio);

  for (let i = 0; i < envPeakList.length; i++) {
    let peak = envPeakList[i]; 
    let env = peak.env; 
    console.log(peak);
    let color = env.getDisplayColor();
    //Show only envelopes with minimum of 0.5%
    let percentInte = peak.intensity/maxIntensity * 100 ;
    if(peak.mz >= para.winMinMz && peak.mz < para.winMaxMz && percentInte >= minPercentage) 
    { 
      if (env.level / (spectrumData.mzLevel.length - 3) >= ratio || ratio <= 0.2){
        //display envelopes based on level, but when the ratio falls below threshold, show all envs in the range
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
          onMouseOverCircle(this,env,peak);
        })
        .on("mouseout",function(){
          onCircleMouseOut(this);
        });
      }
    }
  }
}*/

/**
 * @function drawIons
 * @description Function to add IONS at the top of the peaks for each cluster of envelopes
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} ionData - Contians Ion list to display on the graph
 */
function drawIons(svg: any, para: SpectrumViewParameters,
  ions: MatchedIon[]){
  let ionGroup = svg.append("g").attr("id", "graph_ions");
  
  ions.sort(function(x,y){
    return d3.ascending(x.mz, y.mz);
  })
  //console.log(ions);

  for (let i = 0; i < ions.length; i++) {
    let ion = ions[i];
    let x: number = ion.mz;
    let xPos: number = para.getPeakXPos(x) + para.getIonXShift();
    let yPos: number = para.getPeakYPos(ion.intensity) + para.getIonYShift() + 7;
    //console.log("yPos", yPos);
   // console.log("para.getPeakYPos(ion.intensity)",para.getPeakYPos(ion.intensity))
    
    if(x >= para.getWinMinMz() && x <= para.getWinMaxMz()) {
      if (para.getIsMonoMassGraph()){
        //if mass graph, to avoid overlapping with annotation with the sequence
        if (yPos < 50){
          continue;
        }
      }
      let color: string = "black";
      if (typeof ion.env !== "undefined") {
        color = ion.env.getDisplayColor();
      }
      if (i > 0){
        //check if this ion is annotating the same peak as the previous ion
        //if so, its yPos should be adjusted
        if (ion.mz == ions[i-1].mz){
          yPos = yPos + 15;
        }
      }
      if (yPos >= para.getPadding().head) {
        ionGroup.append("text")
        .attr("id","graph_matched_ions")
        .attr("x", xPos)
        .attr("y", yPos) 
        .style("fill", color)
        .style("opacity", "0.8")
        .style("stroke-width","2")
        .text(ion.text);
      }
    } /*else {
      ionGroup.append("text")
      .attr("id","graph_matched_ions")
      .attr("x", xPos)
      .attr("y", yPos) 
      .style("opacity", "0.8")
      .style("stroke-width","2")
      .text(ion.text);
    }*/
  }
}

/**
 * @function drawSequence
 * @description Draw Sequence on spectrum graph
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function drawSequence(svg: any, para: SpectrumViewParameters, proteoform: Proteoform | null, 
  nMasses: TheoMass[], cMasses: TheoMass[], ions: MatchedIon[]){
	if (!proteoform){
    console.error("ERROR: proteoform object is empty!");
    return;
  }
  let x: number;
  let y: number;
	// Draw | at 0 for prefix mass list
	x = para.getPeakXPos(0);
	y = 15;
  let seq: string = proteoform.getSeq();
  seq = seq.slice(proteoform.getFirstPos(), proteoform.getLastPos() + 1);
  let prevMass: number = -1.0;
  let residues: string = "";
  for (let i = 0; i < nMasses.length; i++) {
    let curMass: number = nMasses[i].getMass();
    if (i > 0) {
      residues = residues + seq[i-1];
    }
    if (curMass > prevMass) {
      if (curMass >= para.getWinMinMz() - 10 && curMass <= para.getWinMaxMz()) {
        let x = para.getPeakXPos(curMass);
          interAddLineAndAnno(svg,x,y, curMass, i, ions, "n");
      }
      if (i > 0) {
        let mz = (prevMass + curMass)/2
        if (mz >= para.getWinMinMz() && mz <= para.getWinMaxMz()) {
          let x = para.getPeakXPos(mz) - 5;
          interAddAminoAcid(svg, x,y+12,residues);
        }
        residues = "";
      }
      prevMass = curMass;
    }
  }

	x = para.getPeakXPos(0);
  y = 35;
  prevMass = -1.0;
  residues = "";
  for (let i = 0; i < cMasses.length; i++) {
    let curMass = cMasses[i].getMass();
    if (i > 0) {
      residues = residues + seq[seq.length - i];
    }
    if (curMass > prevMass) {
      if (curMass >= para.getWinMinMz() - 10 && curMass <= para.getWinMaxMz()) {
        let x = para.getPeakXPos(curMass);
          interAddLineAndAnno(svg,x,y, curMass, i, ions, "c");
      }
      if (i > 0) {
        let mz = (prevMass + curMass)/2
        if (mz >= para.getWinMinMz() && mz <= para.getWinMaxMz()) {
          let x = para.getPeakXPos(mz) - 5;
          interAddAminoAcid(svg, x,y+12,residues);
        }
        residues = "";
      }
      prevMass = curMass;
    }
  }
  function getMatchedIon(residuePos: number, ions: MatchedIon[], 
  terminal: string): MatchedIon | null {
    for (let i = 0; i < ions.length; i++){
      if (!ions[i].pos) {
        console.error("ERROR: ion information is empty");
        return null;
      }
      if (residuePos == parseInt(ions[i].pos!)){
        let text: string = ions[i].text.split(ions[i].pos!)[0]
        if ((terminal == "n" && (text[0] == "A" || text[0] == "B" || text[0] == "C")) ||
            (terminal == "c" && (text[0] == "X" || text[0] == "Y" || text[0] == "Z"))){
              //return text;
              return ions[i];
        }
      }
    }
    return null;
  }
 
  function interAddLineAndAnno(svg: any, x: number, y: number, mass: number, residuePos: number, 
    ions: MatchedIon[], 
    terminal: string) {
    let lineGroup = svg.append("g").attr("id", "seq_peak_line_anno");
    let seqGroup = svg.append("g").attr("id", "graph_sequence");
    let ionData: MatchedIon | null = getMatchedIon(residuePos, ions, terminal);

    let lineYPos: number = 180;

    if (terminal == "c"){
      lineYPos = 160;
    }

    seqGroup.append("line")
      .attr("x1",x)
      .attr("y1",y)
      .attr("x2",x)
      .attr("y2",y+15)
      .attr("stroke","black")
      .attr("stroke-width","1")
      .on("mouseover",function(){
        //@ts-ignore
        onMouseOverFragmentMassAndIonType(this, mass, ionData.text);
      })
      .on("mouseout",function(){
        //@ts-ignore
        onFragmentMassAndIonTypeMouseOut(this);
      });
      if (ionData && para.getShowLines()){//if matched peak exists, draw a dotted line
        let lineYPosEnd: number = para.getPeakYPos(ionData.intensity) + para.getIonYShift();
        if (lineYPosEnd - 8 >= para.getPadding().head) {
          lineGroup.append("line")
          .attr("x1",x)
          .attr("y1",para.getPadding().head)
          .attr("x2",x)
          .attr("y2",lineYPosEnd - 8)
          .attr("stroke","black")
          .attr("stroke-width","1")
          .style("stroke-dasharray", ("5, 6"))
        }
      }
  }
 
	function interAddAminoAcid(svg: any, x: number, y: number, text: string){
    let svgGroup = svg.select("#seq_peak_line_anno");
		svgGroup.append("text")
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
 function addErrorPlot(svg: any, para: SpectrumViewParameters){
  //Draw x-axis
  let xAxis = svg.append("g").attr("id", "xaxis_errorplot").append("line")
    .attr("x1",para.getErrorPlotPadding().left)
    .attr("y1", para.getSVGHeight() - para.getErrorPlotHeight() / 2 - para.getErrorPlotPadding().bottom)
    .attr("x2", para.getSVGWidth() - para.getErrorPlotPadding().right)
    .attr("y2", para.getSVGHeight() - para.getErrorPlotHeight() / 2 - para.getErrorPlotPadding().bottom)
    .attr("stroke","black")
    .style("stroke-dasharray", ("5, 3"))
    .attr("stroke-width","1.5")
  // Draw y-axis
  let yAxis = svg.append("g").attr("id", "yaxis_errorplot").append("line")
    .attr("x1",para.getErrorPlotPadding().left)
    .attr("y1",para.getSVGHeight() - para.getErrorPlotPadding().bottom)
    .attr("x2",para.getErrorPlotPadding().left)
    .attr("y2",para.getSVGHeight() - para.getErrorPlotHeight() - para.getErrorPlotPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","1")
}

/**
 * @function addErrorBox
 * @description Draw Error plot
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function addErrorBox(svg: any, para: SpectrumViewParameters){
  let rectBlock = svg.append("g").attr("id", "rect_error_plot");
  rectBlock.append("line")
    .attr("x1",para.getErrorPlotPadding().left)
    .attr("y1", para.getSVGHeight() - para.getErrorPlotHeight() - para.getErrorPlotPadding().bottom)
    .attr("x2", para.getSVGWidth() - para.getErrorPlotPadding().right)
    .attr("y2", para.getSVGHeight() - para.getErrorPlotHeight() - para.getErrorPlotPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","1")
  rectBlock.append("line")
    .attr("x1",para.getErrorPlotPadding().left)
    .attr("y1", para.getSVGHeight() - para.getErrorPlotPadding().bottom)
    .attr("x2", para.getSVGWidth() - para.getErrorPlotPadding().right)
    .attr("y2", para.getSVGHeight() - para.getErrorPlotPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","1")
  rectBlock.append("line")
    .attr("x1",para.getSVGWidth() - para.getErrorPlotPadding().right)
    .attr("y1",para.getSVGHeight() - para.getErrorPlotPadding().bottom)
    .attr("x2",para.getSVGWidth() - para.getErrorPlotPadding().right)
    .attr("y2",para.getSVGHeight() - para.getErrorPlotHeight() - para.getErrorPlotPadding().bottom)
    .attr("stroke","black")
    .attr("stroke-width","1")
}

/**
 * @function drawErrorYticks
 * @description Draw Error plot y ticks
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 */
function drawErrorYTicks(svg: any, para: SpectrumViewParameters){

	let addYTicks = svg.append("g").attr("id","yErrorTicks")
									.attr("class","yErrorTicks");
  let tickSize: number = para.getErrorThreshold()/para.getErrorYTickNum();
	// Draw ticks
	for(let i=-para.getErrorYTickNum();i<=para.getErrorYTickNum();i+=2) {
		let y: number = para.getErrorYPos(i*tickSize);
		innerDrawYTick(y);
		innerAddErrorYTickValue(i*tickSize,y);
	}
	function innerDrawYTick(y: number){
    //y >= para.padding.head helps the ticks to be in the length of Y axis
    if(!isNaN(y) && y >= para.getPadding().head) {
			addYTicks.append("line")
						.attr("x1",para.getPadding().left)
						.attr("y1",y)
						.attr("x2",para.getPadding().left - para.getTickLength())
						.attr("y2",y)
						.attr("stroke","black")
						.attr("stroke-width","1")
		}
	}
  function innerAddErrorYTickValue(data: number, y: number) {
    if(!isNaN(y) && y >= para.getPadding().head) {
      addYTicks.append("text").attr("class","ytext").attr("x",para.getPadding().left - para.getTickLength())
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
function drawErrorPoints(svg: any, para: SpectrumViewParameters, 
  ionList: MatchedIon[]){
	let circles = svg.append("g").attr("id", "error_circles");
  ionList.forEach((element)=>{
    let mass: number = element.mz;
    if(mass > para.getWinMinMz() && mass <= para.getWinMaxMz()){
      let cx: number = para.getPeakXPos(mass);
      let cy: number = para.getErrorYPos(element.error);
      circles.append("circle")
        .attr("class","error_circles")
        .attr("cx",cx)
        .attr("cy",cy)
        .attr("r", 3)
        .style("fill","black")
        .style("opacity", "1")
        .style("stroke-width","2");
    }
  })
}

function updateViewBox(svgId: string, width: number, height: number) {
  let svg = d3.select("body").select("#"+svgId);
  svg.attr("viewBox", "0 0 "+ width +" "+ height);
}
