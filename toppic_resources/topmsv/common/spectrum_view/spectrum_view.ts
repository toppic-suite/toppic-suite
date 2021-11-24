/**
 * @function SpectrumView
 * @description Function draws the graph, binds zoom and drag function to the graph
 * @param {String} svgId - SVG id on which the graph needed to be drawn
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakData - contains peakList and envelope list 
 * @param {Array} ionData - Contains data with mass and ACID name to plot on the graph
 */
class SpectrumView {
  // parameters for zoom
  private transformX_: number = 0;
  private transformScale_: number = 1.0;
  private id_: string;
  private para_: SpectrumViewParameters;
  private peakList_: Peak[];
  private envList_: Envelope[] = [];
  private ionList_: MatchedIon[] | null = null;
  private svg_: any;
  private proteoform_: Proteoform | null  = null;
  private nIon_: string = "";
  private cIon_: string = "";
  private nMassList_: TheoMass[] = []; 
  private cMassList_: TheoMass[] = [];
  private centerPos_: number = -1;//center m/z or mono mass of the current view

  constructor(svgId: string, peakList: Peak[], sequenceLength: number = -1) {
    this.id_ = svgId;
    this.para_ = new SpectrumViewParameters();
    this.para_.setSeqLength(sequenceLength);
    this.peakList_ = peakList;
    this.para_.initParameters(peakList);
    this.peakList_.sort(function(x,y){
      return y.getIntensity() - x.getIntensity(); 
    });
    $("#" + svgId).data("graph", this);
    // add zoom function
    this.svg_ = d3.select("body").select("#"+svgId);
    this.svg_.attr("viewBox", "0 0 "+ this.para_.getSVGWidth()+" "+ this.para_.getSVGHeight())
      .attr("width", "100%")
      .attr("height", "100%")
      //@ts-ignore
      .call(this.zoom.bind(this));
  }
  //getter
  getPeakList(): Peak[]{
    return this.peakList_;
  }
  getEnvList(): Envelope[]{
    return this.envList_;
  }
  getIonList(): MatchedIon[] | null{
    return this.ionList_;
  }
  getProteoform(): Proteoform | null{
    return this.proteoform_;
  }
  getNIon(): string{
    return this.nIon_;
  }
  getCIon(): string{
    return this.cIon_;
  }
  getNMassList(): TheoMass[] {
    return this.nMassList_;
  }
  getCMassList(): TheoMass[] {
    return this.cMassList_;
  }
  getPara(): SpectrumViewParameters{
    return this.para_;
  }
  getTransformX(): number{
    return this.transformX_;
  }
  getTransformScale(): number{
    return this.transformScale_;
  }
  getSvgId(): string{
    return this.id_;
  }
  getCenterPos(): number {
    return this.centerPos_;
  }
  setTransformX(transformX: number): void{
    this.transformX_ = transformX;
  }
  setTransformScale(transformScale: number): void{
    this.transformScale_ = transformScale;
  }
  setCenterPos(newCenter: number): void {
    this.centerPos_ = newCenter;
  }
  addRawSpectrumAnno(envList: Envelope[], ionList: MatchedIon[] | null): void{
    this.envList_ = envList;
    this.para_.addColorToEnvelopes(envList);
    //this.envPeakList = this.getEnvPeakList(this.envList);
    if (!ionList){
      console.error("ERROR: invalid input for spectrum graph");
      return;
    }
    this.ionList_ = ionList; 
  }

  addMonoMassSpectrumAnno(ionList: MatchedIon[] | null, 
  proteoform: Proteoform, nIonType: string, cIonType: string): void{
    this.ionList_ = ionList; 
    this.proteoform_ = proteoform;
    this.nIon_ = nIonType;
    this.cIon_ = cIonType;
    this.nMassList_ = proteoform.getNMasses(nIonType); 
    this.cMassList_ = proteoform.getCMasses(cIonType);
  }

  redraw(monoMz?: number): void {
    if (this.para_.getIsMonoMassGraph() && monoMz) {
      this.para_.updateMassRange(monoMz);
    } else if(monoMz) {
      this.para_.updateMzRange(monoMz);
    }
    this.setCenterPos(this.para_.getWinCenterMz());
    drawBasicSpectrum(this.id_, this.para_, this.peakList_, this.ionList_);
    if (this.para_.getIsMonoMassGraph() && this.ionList_) {
      drawMonoMassSpectrum(this.id_, this.para_, this.proteoform_, this.nMassList_, this.cMassList_, this.ionList_);
    }
    else {
      //drawRawSpectrum(this.id, this.para, this.envPeakList);
      drawRawSpectrum(this.id_, this.para_, this.envList_);
    }
  }

  zoomed(): void {
    let transform = d3.event.transform;
    let graph = $("#"+ this.getSvgId()).data("graph");
    let svg = document.getElementById(this.getSvgId());
    if (svg) {
      let distance = transform.x - graph.getTransformX();
      let ratio = transform.k / graph.getTransformScale();
      graph.setTransformX(transform.x);
      graph.setTransformScale(transform.k);
      let mousePos = d3.mouse(svg);
      if (ratio == 1) {
        graph.getPara().drag(distance);
      }
  
      graph.getPara().zoom(mousePos[0], mousePos[1], ratio);
      graph.redraw(); 
    }
    
  }

  zoom = d3.zoom()
    .on("zoom", this.zoomed.bind(this))
}
