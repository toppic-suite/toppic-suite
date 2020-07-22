/**
 * @function SpectrumGraph
 * @description Function draws the graph, binds zoom and drag function to the graph
 * @param {String} svgId - SVG id on which the graph needed to be drawn
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakData - contains peakList and envelope list 
 * @param {Array} ionData - Contains data with mass and ACID name to plot on the graph
 */
SpectrumGraph = function(svgId, peakList, envList, ionList){
  // parameters for zoom
  this.transformX = 0;
  this.transformScale = 1.0;

  this.id = svgId;
  this.para = new SpectrumParameters();
  this.peakList = peakList;
  this.para.initParameters(peakList);
  this.peakList.sort(function(x,y){
    return y.intensity - x.intensity; 
  });
  this.envList = envList;
  this.ionList = ionList;
  $("#" + svgId).data("graph", this);

  this.redraw = function(){
    drawSpectrum(this.id, this.para, this.peakList, this.envList, this.ionList);
  }

  this.zoomed = function () {
    let transform = d3.event.transform;
    let graph = $("#"+ this.id).data("graph");
    let distance = transform.x - graph.transformX;
    let ratio = transform.k / graph.transformScale;
    graph.transformX = transform.x;
    graph.transformScale = transform.k;
    let mousePos = d3.mouse(this);
    if(ratio == 1) {
      graph.para.drag(distance);
    }
    else{
      graph.para.zoom(mousePos[0], mousePos[1], ratio);
    }
    graph.redraw(); 
  }

  this.zoom = d3.zoom()
    .on("zoom", this.zoomed);

  // add zoom function
  this.svg = d3.select("body").select("#"+svgId);
  this.svg.attr("viewBox", "0 0 "+ this.para.svgWidth+" "+ this.para.svgHeight)
    .attr("width", "100%")
    .attr("height", "100%")
    .call(this.zoom);
  this.redraw();
}
