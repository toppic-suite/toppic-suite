/**	@function SpectrumViewParameters
 * @description Get data from global variable spectrum_data and utilities to manupulate
 * the data
 */

 class SpectrumViewParameters {
  // Ratio between average and monoisopotic mass
  private avgToMonoRatio_: number = 1.000684;

  // SVG size
  private svgWidth_: number = 910;
  private svgHeight_: number = 270;
	// SVG padding 
	private padding_: Padding = {left:70, right:50, head:10, bottom:50};
  // spectrum size
	private specWidth_: number = this.svgWidth_ - this.padding_.left - this.padding_.right;
  private specHeight_: number = this.svgHeight_ - this.padding_.head - this.padding_.bottom;

  // M/z range of visuable window
  private winMinMz_: number = 0;
  private winMaxMz_: number = 2000;
  private winCenterMz_: number = 1000;

  //minimum possible m/z after zooming/dragging, to prevent dragging/zooming into negative m/z value
  //if new m/z is less than this value, it is reset to this value
  private minPossibleMz_: number = -100;
  private maxPossibleMzMargin_: number = 300;

  // M/z range of peaks
  private dataMinMz_: number = 0;
  private dataMaxMz_: number = 2000;

  // M/z range, color of highlighted part.
  private showHighlight_: boolean = false;
  private hlMinMz_: number = 0;
  private hlMaxMz_: number = 0;
  private hlColor_: string = "gray";

  // Max intensity of visuable window
  private winMaxInte_: number = 30000;

  // Intensity range of peaks
  private dataMaxInte_: number = 30000;
  private dataMinInte_: number = 0;
  // add a margin so that the visuable intensity range is [0, dataMaxInte * inteMargin]
  private inteMargin_: number = 1.2;

  // scale m/z to x coordinate
  private xScale_: number = 0.35;
  // scale intensity to y coordinate
  private yScale_: number = 0.005;

  // Numbers of ticks
  private xTickNum_: number = 10;
  private yTickNum_: number = 5 ;
  private tickLength_: number = 7 ;
  // Tick width list used in the function getTickWidth
  private tickWidthList_: number[] = [10000,8000,6000,5000,4000,3000,2000,1000,800,700,600,500,450,400,350,300,250,200,150,100,50,20,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,0.00001,0.000005,0.000001];
  // Tick height list used in the function getTickHeight
  private tickHeightList_: number[] = [50,40,30,25,20,15,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001];

  //Limiting the peaks and envelopes to 4000 using 20 bins
  private binNum_: number = 20;
  private peakNumPerBin_: number = 50;
  //Padding for mouse over peak floatings.
  private mouseOverPadding_: {"head": number, "middle": number} = {head:20,middle:14};

  // Envelope circle size: min and max radius	
  private showEnvelopes_: boolean = true;
  private defaultRadius_: number = 0.05;
  private minRadius_: number = 2;
  private maxRadius_: number = 5;
  //	Colors for the envelope circles	
  private envColorList_: string[] = ["red","darkorange","blue"];

  // Parameters related to annoated ions
  private showIons_: boolean = true;
  private ionXShift_: number = -5;
  private ionYShift_: number = -15;

  // Mono mass graph
  private showError_: boolean = true;
  private showLines_: boolean = true;
  private isMonoMassGraph_: boolean = false;
  private errorPlotPadding_: Padding = {left:70, right:50, head:10, bottom:10};
  private errorPlotHeight_: number = 40;
  private errorThreshold_: number = 0.2;
  private errorYTickNum_: number = 2;
  
  //restrict x zoom
  private isXZoomAllowed_ = true;

  //sequence length
  //for determining max m/z window based on seq length in mass graph
  private seqLength_: number = -1; 

  constructor() {
  }

  //getters and setters
  getErrorYTickNum(): number{
    return this.errorYTickNum_;
  }
  getErrorThreshold(): number{
    return this.errorThreshold_;
  }
  getErrorPlotPadding(): Padding{
    return this.errorPlotPadding_;
  }
  getErrorPlotHeight(): number{
    return this.errorPlotHeight_;
  }
  getIonXShift(): number{
    return this.ionXShift_;
  }
  getIonYShift(): number{
    return this.ionYShift_;
  }
  getDataMaxMz(): number{
    return this.dataMaxMz_;
  }
  getDataMinMz(): number{
    return this.dataMinMz_;
  }
  getHlColor(): string{
    return this.hlColor_;
  }
  getHlMinMz(): number{
    return this.hlMinMz_;
  }
  getHlMaxMz(): number{
    return this.hlMaxMz_;
  }
  getDataMaxInte(): number{
    return this.dataMaxInte_;
  }
  getWinMaxMz(): number{
    return this.winMaxMz_;
  }
  getWinMinMz(): number{
    return this.winMinMz_;
  }
  getWinMaxInte(): number {
    return this.winMaxInte_;
  }
  getWinCenterMz(): number {
    return this.winCenterMz_;
  }
  getTickLength(): number{
    return this.tickLength_; 
  }
  getYTickNum(): number{
    return this.yTickNum_;
  }
  getSVGWidth(): number{
    return this.svgWidth_; 
  }
  getSVGHeight(): number{
    return this.svgHeight_; 
  }
  getIsMonoMassGraph(): boolean{
    return this.isMonoMassGraph_; 
  }
  getShowHighlight(): boolean{
    return this.showHighlight_;
  }
  getShowEnvelopes(): boolean{
    return this.showEnvelopes_;
  }
  getShowIons(): boolean{
    return this.showIons_;
  }
  getShowError(): boolean{
    return this.showError_;
  }
  getShowLines(): boolean{
    return this.showLines_;
  }
  getSpecHeight(): number {
    return this.specHeight_;
  }

  getPadding(): Padding {
    return this.padding_;
  }
  getIsXZoomAllowed(): boolean {
    return this.isXZoomAllowed_;
  }
  setSeqLength(seqLength: number): void{
    this.seqLength_ = seqLength; 
  }
  setShowEnvelopes_(showEnvelopes: boolean): void{
    this.showEnvelopes_ = showEnvelopes;
  }
  setShowIons(showIons: boolean): void{
    this.showIons_ = showIons;
  }
  setShowError(showError: boolean): void{
    this.showError_ = showError;
  }
  setShowLines(showLines: boolean) : void {
    this.showLines_ = showLines;
  }
  setIsXZoomAllowed_(allowXZoom: boolean): void {
    this.isXZoomAllowed_ = allowXZoom;
  }
  setSVGHeight(newHeight: number): void {
    this.svgHeight_ = newHeight;
  }
  setPadding(left: number, right: number, head: number, bottom: number): void {
    this.padding_.left = left;
    this.padding_.right = right;
    this.padding_.head = head;
    this.padding_.bottom = bottom;
  }
  setSpecHeight(height: number): void {
    this.specHeight_ = height;
  }

  /**
   * @function getTickWidth
   * @description Function Provides width between each tick when zoomed in and out or dragged
   */
  getTickWidth(): number{
    let tempDiff: number = this.winMaxMz_ - this.winMinMz_;
    let tickWidth: number = Math.floor(this.tickWidthList_[0]);
    for(let i = 0; i < this.tickWidthList_.length; i++)
    {
      if(tempDiff/this.xTickNum_ <= Math.floor(this.tickWidthList_[i]) && 
         tempDiff/this.xTickNum_ > Math.floor(this.tickWidthList_[i+1]))
      {
        tickWidth = Math.floor(this.tickWidthList_[i]);
        break ;
      }
    }
    return tickWidth;
  }

  getXTickPosList(): number[] {
    let posList: number[] = new Array(this.xTickNum_ + 1);
    let tickWidth: number = this.getTickWidth();
    for(let i=0; i <= this.xTickNum_; i++)
    {
      // calculate the actual tick position based on the current minMz value on the xaxis
      let tickMz: number = 0;
      if(tickWidth < 1 && tickWidth != 0)
      {
        tickMz = (i*tickWidth + this.winMinMz_) - Math.floor((i*tickWidth + this.winMinMz_)%tickWidth) ;
      }
      else if(tickWidth != 0)
      {
        tickMz = i*tickWidth + this.winMinMz_ - (i*tickWidth + this.winMinMz_)%tickWidth ;
      }
      posList[i] = tickMz;
    }
    return posList;
  }

  /**
   * @function getTickHeight
   * @description Function Provides height between each tick when zoomed in and out or dragged
   */
  getTickHeight(): number{
    let tickheight: number = Math.floor(this.tickHeightList_[0]) ;
    let maxIntPercent: number = this.winMaxInte_/this.dataMaxInte_ * 100;
		for(let i = 0; i < this.tickHeightList_.length; i++)
    {
			if(maxIntPercent/this.yTickNum_ <= this.tickHeightList_[i] 
        && maxIntPercent/this.yTickNum_ > this.tickHeightList_[i+1])
			{
				tickheight = this.tickHeightList_[i];
				return tickheight;
			}
    }
	  return -1;
  }

  /**
   * @function getPeakXPos
   * @description Function provides the x coordinate for the mass
   */
  getPeakXPos(mz: number): number {
    let peakX: number = (mz - this.winMinMz_) * this.xScale_ + this.padding_.left;
    return peakX;
  }
  /**
   * @function getPeakYPos
   * @description Function provides the y coordinate for the intensity
   */
  getPeakYPos(intensity: number): number {
    let peakY: number = this.svgHeight_ - intensity * this.yScale_ - this.padding_.bottom;
    return peakY;
  }

  /**
   * @function getErrorYPos
   * @description Function provides the y coordinate for the error val on the error plot
   */
  getErrorYPos(errorVal: number): number {
    // Multiply with 2 as the coordinates has to be both positive and negative
    let yErrorScale: number = this.errorPlotHeight_/(this.errorThreshold_*2);
    let pos: number = this.svgHeight_ - (errorVal * yErrorScale) 
      - this.errorPlotPadding_.bottom - this.errorPlotHeight_/2;
    return pos;
  }

  /**
   * @function getBinWidth
   * @description Function to compute bin width
   **/
  getBinWidth(): number {
    let width: number = (this.winMaxMz_ - this.winMinMz_)/this.binNum_;
    return width;
  }

  /**
   * @function getCircleSize
   * @description Function provides the radius of the circles drawn on the graph as zoomed in and out
   */
  getCircleSize(): number {
    let radius: number = this.defaultRadius_ * this.xScale_;
    if (radius < this.minRadius_) {
      radius = this.minRadius_;
    }
    if (radius > this.maxRadius_) {
      radius = this.maxRadius_;
    }
    return radius;
  }

  /**
   * Function to set spectrum perameters based on the data
   * @param {Array} peakList - contains the list of data with mz and intensity used to draw lines on the graph 
   */
  compDataRanges(peakList: Peak[]): number[]{
    let minMz: number = 0;
    let maxMz: number = 2000;
    let maxInte: number = 100;
    if (peakList != null && peakList.length > 0) {
      // Sort by mz
      peakList.sort(function(x,y){
        return x.getPos() - y.getPos();
      });
      let listSize: number = peakList.length;
      maxMz = Math.floor(peakList[listSize-1].getPos());

      // Sort by intensity
      peakList.sort(function(x,y){
        return x.getIntensity() - y.getIntensity();
      });
      maxInte = Math.floor(peakList[listSize-1].getIntensity());
    }
    return [minMz, maxMz, maxInte];
  }

  /**
   * @function updateScale
   * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
   * initializing xScale, yScale.
   */
  updateScale(winMinMz: number, winMaxMz: number, winMaxInte: number): void {
    this.winMinMz_ = winMinMz;
    this.winMaxMz_ = winMaxMz;
    if(winMinMz == this.dataMinMz_ && winMaxMz == this.dataMaxMz_)
    {
      this.winMinMz_ = 0;
      this.winMaxMz_ = 1.1 * this.winMaxMz_;
    }
    this.winCenterMz_ = (this.winMinMz_ + this.winMaxMz_)/2.0;
    this.xScale_ = this.specWidth_/(this.winMaxMz_ - this.winMinMz_);

    this.winMaxInte_ = winMaxInte;
    this.yScale_ = this.specHeight_/this.winMaxInte_;
  }

  /**
   * @function initParameters
   * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
   * initializing xScale, yScale.
   */
  initParameters(peakList: Peak[]): void {
    let [dataMinMz, dataMaxMz, dataMaxInte] = this.compDataRanges(peakList);
    this.dataMinMz_ = dataMinMz;
    this.dataMaxMz_ = dataMaxMz + (0.10 * dataMaxMz);
    this.dataMaxInte_ = dataMaxInte;
    // add 1/4th of max intensity to keep the max point at 3/4 of the y axis*
    this.updateScale(this.dataMinMz_, this.dataMaxMz_, this.dataMaxInte_ * this.inteMargin_);
  }

  /**
   * @function drag
   * @description 
   * Function provides minMz and maxMz based on the amount of drag done
   */
  drag(distX: number): void {
    let mzDist: number = distX / this.xScale_;
    
    this.winMinMz_ = this.winMinMz_ - mzDist; 
    this.winMaxMz_ = this.winMaxMz_ - mzDist;
    this.winCenterMz_ = this.winCenterMz_ - mzDist;
    
    //allow drag up to -50 m/z (this.minPossibleMz) to give some padding 
    if (this.winMinMz_ - mzDist < this.minPossibleMz_){
        let minMaxDiff: number = this.winMaxMz_ - this.winMinMz_;
        let centerDiff: number = this.winCenterMz_ - this.winMinMz_;

        this.winMinMz_ = this.minPossibleMz_;
        this.winMaxMz_ = this.winMinMz_ + minMaxDiff;
        this.winCenterMz_ = this.winMinMz_ + centerDiff;
    }
    if (this.winMaxMz_ > this.dataMaxMz_ + this.maxPossibleMzMargin_) {
      let minMaxDiff: number = this.winMaxMz_ - this.winMinMz_;
      this.winMaxMz_ = this.dataMaxMz_ + this.maxPossibleMzMargin_;
      this.winMinMz_ = this.winMaxMz_ - minMaxDiff;
      let centerDiff: number = this.winMaxMz_ - this.winCenterMz_;
      this.winCenterMz_ = this.winMinMz_ + centerDiff;
    }
  }

  /**
   * @function xZoom
   * @description Function provides with current xScale, current minMz and MaxMz based on the zoom on x-axis.
   * Function also calls setLimita which helps in drawing limited number of peaks and circles per eachbin/range of mz values.
   */
  xZoom(mouseSvgX: number, ratio: number): void {
    if (!this.isXZoomAllowed_) {
      return;
    }
    let oriValues: {"min": number, "max": number, "center": number, "xScale": number} = {} as {"min": number, "max": number, "center": number, "xScale": number}; //so that the view range can be restored when the view shouldn't be zoomed
    oriValues.min = this.winMinMz_;
    oriValues.max = this.winMaxMz_;
    oriValues.center = this.winCenterMz_;
    oriValues.xScale = this.xScale_;

    let mouseSpecX: number = mouseSvgX - this.padding_.left;
    this.winCenterMz_ =  mouseSpecX/this.xScale_ + this.winMinMz_;
    /*self is a global variable of datasource object containing all the data needed to use when zoomed*/
    this.xScale_ = this.xScale_ * ratio ; 
    this.winMinMz_ = this.winCenterMz_ - mouseSpecX / this.xScale_; 
    this.winMaxMz_ = this.winCenterMz_ + (this.specWidth_ - mouseSpecX) / this.xScale_;
    //console.log(this.winMaxMz_, this.dataMaxMz_ + 500)
    if (this.winMinMz_ < this.minPossibleMz_){//prevent zooming out into negative mass
      this.winMinMz_ = this.minPossibleMz_;
    }
    if (this.winMaxMz_ > this.dataMaxMz_ + this.maxPossibleMzMargin_){
      this.winMaxMz_ = this.dataMaxMz_ + this.maxPossibleMzMargin_;
    }
    if (this.winCenterMz_ > this.winMaxMz_) {
      this.winMinMz_ = oriValues.min;
      this.winMaxMz_ = oriValues.max;
      this.winCenterMz_ = oriValues.center;
      this.xScale_ = oriValues.xScale;
    }
  }
  /**
   * @function yZoom
   * @description Function provides with current yScale, current max Intensity based on the zoom on y-axis
   */
  yZoom(ratio: number): void {
    //Reducing zoom factor to smoothenup and remove gliches
    if(ratio > 1 ) ratio = 1.4;
    else if(ratio < 1) ratio = 0.9;
    //restrict zooming in when current max intensity is smaller than 0.01% of the max intensity of entire data
    if ((ratio > 1.0 && (this.winMaxInte_ >= this.dataMaxInte_ * 0.01 / 100)) 
    || (ratio < 1.0 && (this.winMaxInte_ <= this.dataMaxInte_ * this.inteMargin_))) {
      this.yScale_ = this.yScale_ * ratio;
      this.winMaxInte_ = this.specHeight_ / this.yScale_;
    }
  }

  /**
   * @function zoom
   * @description 
   * Function to invoke respective zoom functionality(zoom on x or y) based on position of X, Y 
   * It fixes amount of zoom based on zooming in or out 
   */
  zoom(mouseSvgX: number, mouseSvgY: number, ratio: number): void {
    if(ratio > 1 ) ratio = 1.4; // Zooming in and fixing ration to 1.4 (fixed values based on testing the smooting of zoom)
    else if(ratio < 1) ratio = 0.9; // Zooming out and fixing ration to 0.9 (fixed values based on testing the smooting of zoom)
    if (mouseSvgY > this.svgHeight_ - this.padding_.bottom) {
      this.xZoom(mouseSvgX, ratio);
    }
    else {
      this.yZoom(ratio);
    }
  }

  /**
   * @function addColorToEnvelopes
   * @description 
   * Add color to envelopes.
   */
  addColorToEnvelopes(envList: Envelope[]): void{
    if(!envList || envList.length === 0 || typeof envList[0].getPeaks() === "undefined") return;
    envList.sort(function(x,y){
      return (x.getPeaks()[0].getPos() - y.getPeaks()[0].getPos());
    })
    let colorNum: number = this.envColorList_.length; 
    for (let i = 0; i < envList.length; i++) 
    {
      envList[i].setDisplayColor(this.envColorList_[i%colorNum]);
    }
  }

  /**
   * @function setHighlight
   * @description 
   * set highlight region for MS1 precursor envelope
   */
  setHighlight(ms1Spec: Spectrum): void {
    this.showHighlight_ = true;
    this.hlMinMz_ = ms1Spec.getMinMz();
    this.hlMaxMz_ = ms1Spec.getMaxMz();
    //console.log(precMonoMz, this.hlMinMz, this.hlMaxMz);
  }

  setDefaultPadding(): void {
    this.padding_.head = 20;
    this.padding_.bottom = 50;
  }

  setMonoMassPaddding(): void {
    this.padding_.head = 60;
    this.padding_.bottom = 75;
  }

  setMonoMassGraph(isMonoMass: boolean): void {
    this.isMonoMassGraph_ = isMonoMass;
    if (isMonoMass) {
      this.setMonoMassPaddding();
    }
    else {
      this.setDefaultPadding();
    }
    this.specHeight_ = this.svgHeight_ - this.padding_.head - this.padding_.bottom;
    this.updateScale(this.winMinMz_, this.winMaxMz_, this.winMaxInte_);
  }

  /**
   * @function updataMzRange
   * @description 
   */
  updateMzRange(monoMz: number): void {
    let centerMz: number = monoMz * this.avgToMonoRatio_;
    this.winMinMz_ = centerMz - 3;
    this.winMaxMz_ = centerMz + 3;
    this.updateScale(this.winMinMz_, this.winMaxMz_, this.winMaxInte_);
  }

  updateMassRange(mass: number): void {
    let centerMass: number = mass;
    this.winMinMz_ = centerMass - 3;
    this.winMaxMz_ = centerMass + 3;
    this.updateScale(this.winMinMz_, this.winMaxMz_, this.winMaxInte_);
  }
}
