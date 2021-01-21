/**
 * @function SpectrumGraph
 * @description Function draws the graph, binds zoom and drag function to the graph
 * @param {String} svgId - SVG id on which the graph needed to be drawn
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakData - contains peakList and envelope list 
 * @param {Array} ionData - Contains data with mass and ACID name to plot on the graph
 */
class SpectrumGraph {
  // parameters for zoom
  transformX = 0;
  transformScale = 1.0;

  constructor(svgId, peakList, sequenceLength = -1) {
    this.id = svgId;
    this.para = new SpectrumParameters();
    this.para.seqLength = sequenceLength;
    this.peakList = peakList;
    this.para.initParameters(peakList);
    this.peakList.sort(function(x,y){
      return y.getIntensity() - x.getIntensity(); 
    });
    $("#" + svgId).data("graph", this);
    // add zoom function
    this.svg = d3.select("body").select("#"+svgId);
    this.svg.attr("viewBox", "0 0 "+ this.para.svgWidth+" "+ this.para.svgHeight)
      .attr("width", "100%")
      .attr("height", "100%")
      .call(this.zoom);
  }


  addRawSpectrumAnno(envList, ionList){
    this.envList = envList;
    this.para.addColorToEnvelopes(envList);
    this.envPeakList = this.getEnvPeakList(this.envList);
    this.ionList = ionList; 
  }

  addMonoMassSpectrumAnno(ionList, proteoform, nIonType, cIonType){
    this.ionList = ionList; 
    this.proteoform = proteoform;
    this.nIon = nIonType;
    this.cIon = cIonType;
    //console.log(nIonType, cIonType);
    this.nMassList = proteoform.getNMasses(nIonType); 
    //console.log(proteoform);
    //console.log(this.nMassList);
    this.cMassList = proteoform.getCMasses(cIonType);
    //console.log(this.cMassList);
  }

  redraw = function(monoMz){
    if (this.para.isMonoMassGraph && monoMz) {
      this.para.updateMassRange(monoMz);
    } else if(monoMz) {
      this.para.updateMzRange(monoMz);
    }
    drawBasicSpectrum(this.id, this.para, this.peakList, this.ionList);

    if (this.para.isMonoMassGraph) {
      drawMonoMassSpectrum(this.id, this.para, this.proteoform, this.nMassList, this.cMassList, this.ionList);
    }
    else {
      drawRawSpectrum(this.id, this.para, this.envPeakList);
    }
  }

  zoomed = function () {
    let transform = d3.event.transform;
    let graph = $("#"+ this.id).data("graph");
    let distance = transform.x - graph.transformX;
    let ratio = transform.k / graph.transformScale;
    graph.transformX = transform.x;
    graph.transformScale = transform.k;
    let mousePos = d3.mouse(this);
    if (ratio == 1) {
      graph.para.drag(distance);
    }
    graph.para.zoom(mousePos[0], mousePos[1], ratio);
    graph.redraw(); 
  }

  zoom = d3.zoom()
    .on("zoom", this.zoomed)

  getEnvPeakList = function(envList) {
    if (!envList || envList.length === 0 || typeof envList[0].env_peaks === "undefined") {return [];}
    let envPeakList = [];
    for (let i = 0; i < envList.length; i++) {
      let env = envList[i];
      for (let j = 0; j < env.env_peaks.length; j++) {
        let peak = env.env_peaks[j];
        peak.env = env;
        envPeakList.push(peak);
      }
    }
    envPeakList.sort(function(x,y){
      return y.intensity - x.intensity;
    });
    return envPeakList;
  }
}
