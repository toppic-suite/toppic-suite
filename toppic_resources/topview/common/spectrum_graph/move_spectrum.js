/**
 * Function gets invoked whenever drag is invoked. It redraws the spectrum but keeps old peaks and envelope data.
 * This function invokes all the functions that draw the graph
 * @param {string} svgId - Contains the Id from html on which the graph should be drawn
 * @param {object} para - Contains the parameters to help draw the graph
 * @param {list} peaks - Contains the data of Peaks to draw lines on the graph
 * @param {list} envPeaks - Contains the data of Envelope peaks to draw circles on the graph
 * @param {list} ions - Contains Ions to draw upon the peaks
 */
function moveSpectrum(svgId, para, peaks, envPeaks, proteoform, ions) {
  let svg = d3.select("body").select("#"+svgId);
  // svg.attr("width", para.svgWidth).attr("height", para.svgHeight);
  // Removes all the elements under SVG group 'svgGroup' everytime there this function is called
  svg.selectAll("#svgGroup").remove();
  // Create a group under which all the fucntions of the graph will be appended
  svg = svg.append("g").attr("id","svgGroup");

  //call onMouseOut everytime to fix onHover bug adding multiple data when mouseover and zoomed up
  onMouseOut();
  drawTicks(svg, para);
  drawAxis(svg,para);
  addDatatoAxis(svg,para);
  addLabels(svg,para);
  if (para.showHighlight) {
    addHighlight(svg, para);
  }
  movePeaks(svg, para, peaks);
  if (para.showEnvelopes) {
    moveEnvelopes(svg, para, envPeaks);
  }
  if (para.showIons) {
    drawIons(svg, para, ions);
  }
  if (para.isMonoMassGraph) {
    updateViewBox(svgId, para.svgWidth, para.svgHeight);
    drawSequence(svg, para, proteoform);
    addErrorPlot(svg, para);
    addErrorBox(svg, para);
    drawErrorYTicks(svg, para);
    drawErrorPoints(svg, para, ions);
  }
}
/**
 * @function movePeaks
 * @description Function to redraw peak lines on the graph when dragging
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
function movePeaks(svg,para,peakList){
  let peaks = svg.append("g").attr("id", "peaks");
  var len = peakList.length;
  // limits provide current count of number of peaks drawn on graph per bin(range) 
  // so that we can limit tha peak count to peaksPerRange count
  let limits = new Array(para.binNum).fill(0);
  let binWidth = para.getBinWidth();

  let minBin = Math.floor(para.winMinMz/binWidth);

  for(let i =0;i<len;i++)
  {
    let peak = peakList[i];
    if(peak.mz >= para.winMinMz && peak.mz < para.winMaxMz)
    {
      let minMz = para.winMinMz;
      if (minMz < 0){
        minMz = 0;
      }
      let binIndex = Math.floor((peak.mz)/binWidth); 

      if (binIndex <= minBin + para.binNum) 
      {
        limits[binIndex - minBin - 1] = limits[binIndex - minBin - 1]+1;
        if (limits[binIndex - minBin - 1] <= para.peakNumPerBin) {
         // console.log("peak drawing")
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
 * @function moveEnvelopes
 * @description Function to add circles for the envelope data when dragging
 * @param {Node} svg -  is a html node on which the graph is being ploted
 * @param {object} para - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakdata - Contians both peak list and envelopelist
 */
function moveEnvelopes(svg, para, envPeakList) {
  let circles = svg.append("g").attr("id", "circles");
  let minPercentage = 0.0;
  let maxIntensity = para.dataMaxInte ;
  // limits provide current count of number of peaks drawn on graph per bin(range)
  // so that we can limit tha peak count to circlesPerRange count
  let limits = new Array(para.binNum).fill(0);
  let binWidth = para.getBinWidth();

  let minBin = Math.floor(para.winMinMz/binWidth);

  for (let i = 0; i < envPeakList.length; i++) {
    let peak = envPeakList[i]; 
    let env = peak.env; 
    //console.log(env);
    let color = env.color;
    //Show only envelopes with minimum of 0.5%
    let percentInte = peak.intensity/maxIntensity * 100 ;
    
    if(peak.mz >= para.winMinMz && peak.mz < para.winMaxMz && percentInte >= minPercentage) 
    { 
      let binIndex = Math.floor((peak.mz)/binWidth); 
      
      if (binIndex <= minBin + para.binNum)//so that this evaluation is looking at the current view range of graph
      {
        limits[binIndex - minBin - 1] = limits[binIndex - minBin - 1]+1;
        //binIndex - minBin -1 because limits array starts from index 0. 
        if (limits[binIndex - minBin - 1] <= para.peakNumPerBin) 
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