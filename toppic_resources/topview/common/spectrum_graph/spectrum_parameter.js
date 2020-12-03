/**	@function SpectrumParameters
 * @description Get data from global variable spectrum_data and utilities to manupulate
 * the data
 */
class SpectrumParameters {
  // Ratio between average and monoisopotic mass
  avgToMonoRatio = 1.000684;

  // SVG size
  svgWidth = 910;
  svgHeight = 270;
	// SVG padding 
	padding = {left:70, right:50, head:10, bottom:50};
  // spectrum size
	specWidth = this.svgWidth - this.padding.left - this.padding.right;
  specHeight = this.svgHeight - this.padding.head - this.padding.bottom;

  // M/z range of visuable window
  winMinMz = 0 ;
  winMaxMz = 2000;
  winCenterMz = 1000;

  // M/z range of peaks
  dataMinMz = 0;
  dataMaxMz = 2000;

  // M/z range, color of highlighted part.
  showHighlight = false;
  hlMinMz = 0;
  hlMaxMz = 0;
  hlColor = "gray";

  // Max intensity of visuable window
  winMaxInte = 30000;

  // Intensity range of peaks
  dataMaxInte = 30000;
  dataMinInte = 0;
  // add a margin so that the visuable intensity range is [0, dataMaxInte * inteMargin]
  inteMargin = 1.2;

  // scale m/z to x coordinate
  xScale = 0.35;
  // scale intensity to y coordinate
  yScale = 0.005;

  // Numbers of ticks
  xTickNum = 10;
  yTickNum = 5 ;
  tickLength = 7 ;
  // Tick width list used in the function getTickWidth
  tickWidthList = [10000,8000,6000,5000,4000,3000,2000,1000,800,700,600,500,450,400,350,300,250,200,150,100,50,20,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,0.00001,0.000005,0.000001];
  // Tick height list used in the function getTickHeight
  tickHeightList = [50,40,30,25,20,15,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001];

  //Limiting the peaks and envelopes to 4000 using 20 bins
  binNum = 20;
  peakNumPerBin = 50;
  //Padding for mouse over peak floatings.
  mouseOverPadding = {head:20,middle:14};

  // Envelope circle size: min and max radius	
  showEnvelopes = true;
  defaultRadius = 0.05;
  minRadius = 2;
  maxRadius = 5;
  //	Colors for the envelope circles	
  envColorList = ["red","darkorange","blue"];

  // Parameters related to annoated ions
  showIons = true;
  ionXShift = -5;
  ionYShift = -15;

  // Mono mass graph
  isMonoMassGraph = false;
  errorPlotPadding = {left:70, right:50, head:10, bottom:10};
  errorPlotHeight = 40;
  errorThreshold = 0.2;
  errorYTickNum = 2;
  
  constructor() {
  }

  /**
   * @function getTickWidth
   * @description Function Provides width between each tick when zoomed in and out or dragged
   */
  getTickWidth = function(){
    let tempDiff = this.winMaxMz - this.winMinMz;
    let tickWidth = parseInt(this.tickWidthList[0]) ;
    for(let i = 0; i < this.tickWidthList.length; i++)
    {
      if(tempDiff/this.xTickNum <= parseFloat(this.tickWidthList[i]) && 
         tempDiff/this.xTickNum > parseFloat(this.tickWidthList[i+1]))
      {
        tickWidth = parseFloat(this.tickWidthList[i]);
        break ;
      }
    }
    return 	tickWidth ;
  }

  getXTickPosList = function() {
    let posList = new Array(this.xTickNum + 1);
    let tickWidth = this.getTickWidth();
    for(let i=0; i <= this.xTickNum ; i++)
    {
      // calculate the actual tick position based on the current minMz value on the xaxis
      let tickMz = 0;
      if(tickWidth < 1 && tickWidth != 0)
      {
        tickMz = (i*tickWidth + this.winMinMz) - parseFloat((i*tickWidth + this.winMinMz)%tickWidth) ;
      }
      else if(tickWidth != 0)
      {
        tickMz = i*tickWidth + this.winMinMz - (i*tickWidth + this.winMinMz)%tickWidth ;
      }
      posList[i] = tickMz;
    }
    return posList;
  }

  /**
   * @function getTickHeight
   * @description Function Provides height between each tick when zoomed in and out or dragged
   */
  getTickHeight = function(){
    let tickheight = parseInt(this.tickHeightList[0]) ;
    let maxIntPercent = this.winMaxInte/this.dataMaxInte * 100;
		for(let i = 0; i < this.tickHeightList.length; i++)
		{
			if(maxIntPercent/this.yTickNum <= parseFloat(this.tickHeightList[i]) 
        && maxIntPercent/this.yTickNum > parseFloat(this.tickHeightList[i+1]))
			{
				tickheight = parseFloat(this.tickHeightList[i]);
				break ;
			}
    }
	  return tickheight ;
  }

  /**
   * @function getPeakXPos
   * @description Function provides the x coordinate for the mass
   */
  getPeakXPos = function (mz) {
    let peakX = (mz - this.winMinMz) * this.xScale + this.padding.left;
    return peakX;
  }
  /**
   * @function getPeakYPos
   * @description Function provides the y coordinate for the intensity
   */
  getPeakYPos = function (intensity) {
    let peakY = this.svgHeight - intensity * this.yScale - this.padding.bottom;
    return peakY;
  }

  /**
   * @function getErrorYPos
   * @description Function provides the y coordinate for the error val on the error plot
   */
  getErrorYPos = function(errorVal) {
    // Multiply with 2 as the coordinates has to be both positive and negative
    let yErrorScale = this.errorPlotHeight/(this.errorThreshold*2);
    let pos = this.svgHeight - (errorVal * yErrorScale) 
      - this.errorPlotPadding.bottom - this.errorPlotHeight/2;
    return pos;
  }

  /**
   * @function getBinWidth
   * @description Function to compute bin width
   **/
  getBinWidth = function() {
    let width = (this.winMaxMz - this.winMinMz)/this.binNum;
    return width;
  }

  /**
   * @function getCircleSize
   * @description Function provides the radius of the circles drawn on the graph as zoomed in and out
   */
  getCircleSize = function() {
    let radius = this.defaultRadius * this.xScale;
    if (radius < this.minRadius) {
      radius = this.minRadius;
    }
    if (radius > this.maxRadius) {
      radius = this.maxRadius;
    }
    return radius;
  }

  /**
   * Function to set spectrum perameters based on the data
   * @param {Array} peakList - contains the list of data with mz and intensity used to draw lines on the graph 
   */
  compDataRanges = function(peakList){
    let minMz = 0;
    let maxMz = 2000;
    let maxInte = 100;
    if (peakList != null && peakList.length > 0) {
      // Sort by mz
      peakList.sort(function(x,y){
        return x.mz - y.mz;
      });
      let listSize = peakList.length;
      maxMz = parseFloat(peakList[listSize-1].mz);

      // Sort by intensity
      peakList.sort(function(x,y){
        return x.intensity - y.intensity;
      });
      maxInte = parseFloat(peakList[listSize-1].intensity);
    }
    return [minMz, maxMz, maxInte];
  }

  /**
   * @function updateScale
   * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
   * initializing xScale, yScale.
   */
  updateScale = function(winMinMz, winMaxMz, winMaxInte) {
    this.winMinMz = winMinMz;
    this.winMaxMz = winMaxMz;
    if(winMinMz == this.dataMinMz && winMaxMz == this.dataMaxMz)
    {
      this.winMinMz = 0;
      this.winMaxMz = 1.1 * this.winMaxMz;
    }
    this.winCenterMz = (this.winMinMz + this.winMaxMz)/2.0;
    this.xScale = this.specWidth/(this.winMaxMz - this.winMinMz);

    this.winMaxInte = winMaxInte;
    this.yScale = this.specHeight/this.winMaxInte;
  }

  /**
   * @function initParameters
   * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
   * initializing xScale, yScale.
   */
  initParameters = function(peakList) {
    let [dataMinMz, dataMaxMz, dataMaxInte] = this.compDataRanges(peakList);
    this.dataMinMz = dataMinMz;
    this.dataMaxMz = dataMaxMz + (0.10 * dataMaxMz);
    this.dataMaxInte = dataMaxInte;
    // add 1/4th of max intensity to keep the max point at 3/4 of the y axis*
    this.updateScale(this.dataMinMz, this.dataMaxMz, this.dataMaxInte * this.inteMargin);
  }

  /**
   * @function drag
   * @description 
   * Function provides minMz and maxMz based on the amount of drag done
   */
  drag = function(distX) {
    let mzDist = distX / this.xScale;
    let minMz = -50;
    //allow drag up to -50 m/z to give some padding 

    if (this.winMinMz - mzDist < minMz){
        let minMaxDiff = this.winMaxMz - this.winMinMz;
        let centerDiff = this.winCenterMz - this.winMinMz;

        this.winMinMz = minMz;
        this.winMaxMz = this.winMinMz + minMaxDiff;
        this.winCenterMz = this.winMinMz + centerDiff;
    }
    else{
      this.winMinMz = this.winMinMz - mzDist; 
      this.winMaxMz = this.winMaxMz - mzDist;
      this.winCenterMz = this.winCenterMz - mzDist;
    }
  }

  /**
   * @function xZoom
   * @description Function provides with current xScale, current minMz and MaxMz based on the zoom on x-axis.
   * Function also calls setLimita which helps in drawing limited number of peaks and circles per eachbin/range of mz values.
   */
  xZoom = function (mouseSvgX, ratio) {
   if ((ratio > 1.0) || ((this.winMaxMz - this.winMinMz) < this.dataMaxMz) ) {
      let mouseSpecX = mouseSvgX - this.padding.left;
      this.winCenterMz =  mouseSpecX/this.xScale + this.winMinMz;
      /*self is a global variable of datasource object containing all the data needed to use when zoomed*/
      this.xScale = this.xScale * ratio ; 
      this.winMinMz = this.winCenterMz - mouseSpecX / this.xScale; 
      this.winMaxMz = this.winCenterMz + (this.specWidth - mouseSpecX) / this.xScale;
    }
  }
  /**
   * @function yZoom
   * @description Function provides with current yScale, current max Intensity based on the zoom on y-axis
   */
  yZoom = function (ratio) {
    //Reducing zoom factor to smoothenup and remove gliches
    if(ratio > 1 ) ratio = 1.4;
    else if(ratio < 1) ratio = 0.9;
    if ((ratio > 1.0 && this.winMaxInte >= this.dataMinInte * this.inteMargin) 
      || (ratio < 1.0 && this.winMaxInte <= this.dataMaxInte * this.inteMargin)) {
      this.yScale = this.yScale * ratio;
      this.winMaxInte = this.specHeight / this.yScale;
    }
  }

  /**
   * @function zoom
   * @description 
   * Function to invoke respective zoom functionality(zoom on x or y) based on position of X, Y 
   * It fixes amount of zoom based on zooming in or out 
   */
  zoom = function(mouseSvgX, mouseSvgY, ratio) {
    if(ratio > 1 ) ratio = 1.4; // Zooming in and fixing ration to 1.4 (fixed values based on testing the smooting of zoom)
    else if(ratio < 1) ratio = 0.9; // Zooming out and fixing ration to 0.9 (fixed values based on testing the smooting of zoom)
    if (mouseSvgY > this.svgHeight - this.padding.bottom) {
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
  addColorToEnvelopes = function(envList){
    if(!envList || envList.length === 0 || typeof envList[0].env_peaks === "undefined") return;
    envList.sort(function(x,y){
      return (x.env_peaks[0].mz - y.env_peaks[0].mz);
    })
    let colorNum = this.envColorList.length; 
    for (let i = 0; i < envList.length; i++) 
    {
      envList[i].color = this.envColorList[i%colorNum];
    }
  }

  /**
   * @function setHighlight
   * @description 
   * set highlight region for MS1 precursor envelope
   */
  setHighlight = function(precMonoMz, charge) {
    this.showHighlight = true;
    let monoMz = parseFloat(precMonoMz);
    let centerMz = monoMz * this.avgToMonoRatio;
    let dist = centerMz - monoMz + 0.02;
    this.hlMinMz = centerMz - dist; 
    this.hlMaxMz = centerMz + dist;
    //console.log(precMonoMz, this.hlMinMz, this.hlMaxMz);
  }

  setMonoMassGraph(isMonoMass) {
    this.isMonoMassGraph = isMonoMass;
    if (isMonoMass) {
      this.padding.head = 60;
      this.padding.bottom = 75;
    }
    else {
      this.padding.head = 20;
      this.padding.bottom = 50;
    }
    this.specHeight = this.svgHeight - this.padding.head - this.padding.bottom;
    this.updateScale(this.winMinMz, this.winMaxMz, this.winMaxInte);
  }

  /**
   * @function updataMzRange
   * @description 
   */
  updateMzRange = function(monoMz) {
    let centerMz = parseFloat(monoMz) * this.avgToMonoRatio;
    this.winMinMz = centerMz - 3;
    this.winMaxMz = centerMz + 3;
    this.updateScale(this.winMinMz, this.winMaxMz, this.winMaxInte);
  }

  updateMassRange = function(mass) {
    let centerMass = parseFloat(mass);
    this.winMinMz = centerMass - 3;
    this.winMaxMz = centerMass + 3;
    this.updateScale(this.winMinMz, this.winMaxMz, this.winMaxInte);
  }
}
