/**	@function SpectrumParameters
 * @description Get data from global variable spectrum_data and utilities to manupulate
 * the data
 */
SpectrumParameters = function() {
  // SVG size
  this.svgWidth = 910;
  this.svgHeight = 220;
	// SVG padding 
	this.padding = {left:70, right:20, head:10, bottom:50};
  // spectrum size
	this.specWidth = this.svgWidth - this.padding.left - this.padding.right;
  this.specHeight = this.svgHeight - this.padding.head - this.padding.bottom;

  // M/z range of visuable window
  this.winMinMz = 0 ;
  this.winMaxMz = 2000;
  this.winCenterMz = 1000;

  // M/z range of peaks
  this.dataMinMz = 0;
  this.dataMaxMz = 2000;

  // M/z range, color of highlighted part.
  this.showHighlight = false;
  this.hlMinMz = 600;
  this.hlMaxMz = 700;
  this.hlColor = "orange";

  // Max intensity of visuable window
  this.winMaxInte = 30000;

  // Intensity range of peaks
  this.dataMaxInte = 30000;
  this.dataMinInte = 0;
  // add a margin so that the visuable intensity range is [0, dataMaxInte * inteMargin]
  this.inteMargin = 1.2;

  // scale m/z to x coordinate
  this.xScale = 0.35;
  // scale intensity to y coordinate
  this.yScale = 0.005;

  // Numbers of ticks
  this.xTickNum = 10;
  this.yTickNum = 5 ;
  this.tickLength = 7 ;
  // Tick width list used in the function getTickWidth
  this.tickWidthList = [10000,8000,6000,5000,4000,3000,2000,1000,800,700,600,500,450,400,350,300,250,200,150,100,50,20,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,0.00001,0.000005,0.000001];
  // Tick height list used in the function getTickHeight
  this.tickHeightList = [50,40,30,25,20,15,10,5,3,2,1,0.5,0.2,0.1,0.05,0.01,0.005,0.001];

  //Limiting the peaks and envelopes to 4000 using 20 bins
  this.binNum = 20;
  this.peakNumPerBin = 200;
  //Padding for mouse over peak floatings.
  this.mouseOverPadding = {head:20,middle:14};

  // Envelope circle size: min and max radius	
  this.showEnvelopes = true;
  this.defaultRadius = 0.05;
  this.minRadius = 2;
  this.maxRadius = 5;
  //	Colors for the envelope circles	
  this.envColorList = ["red","orange","blue","green"];

  // Parameters related to annoated ions
  this.showIons = false;
  this.tickWidthThreshhold = 0.5;
  this.adjustableIonPosition = 4;
  this.fixedHeightOfIonAboveThePeak = 10;

  /**
   * @function getTickWidth
   * @description Function Provides width between each tick when zoomed in and out or dragged
   */
  this.getTickWidth = function(){
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

  this.getXTickPosList = function() {
    posList = new Array(this.xTickNum + 1);
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
  this.getTickHeight = function(){
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
  this.getPeakXPos = function (mz) {
    let peakX = (mz - this.winMinMz) * this.xScale + this.padding.left;
    return peakX;
  }
  /**
   * @function getPeakYPos
   * @description Function provides the y coordinate for the intensity
   */
  this.getPeakYPos = function (intensity) {
    let peakY = this.svgHeight - intensity * this.yScale - this.padding.bottom;
    return peakY;
  }

  /**
   * @function getBinWidth
   * @description Function to compute bin width
   **/
  this.getBinWidth = function() {
    let width = (this.winMaxMz - this.winMinMz)/this.binNum;
    return width;
  }

  /**
   * @function getCircleSize
   * @description Function provides the radius of the circles drawn on the graph as zoomed in and out
   */
  this.getCircleSize = function() {
    radius = this.defaultRadius * this.xScale;
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
  this.compDataRanges = function(peakList){
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
  this.updateScale = function(winMinMz, winMaxMz, winMaxInte) {
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
  this.initParameters = function(peakList) {
    [dataMinMz, dataMaxMz, dataMaxInte] = this.compDataRanges(peakList);
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
  this.drag = function(distX) {
    let mzDist = distX / this.xScale;
    this.winMinMz = this.winMinMz - mzDist; 
    this.winMaxMz = this.winMaxMz - mzDist;
    this.winCenterMz = this.winCenterMz - mzDist;
  }

  /**
   * @function xZoom
   * @description Function provides with current xScale, current minMz and MaxMz based on the zoom on x-axis.
   * Function also calls setLimita which helps in drawing limited number of peaks and circles per eachbin/range of mz values.
   */
  this.xZoom = function (mouseSvgX, ratio) {
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
  this.yZoom = function (ratio) {
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
  this.zoom = function(mouseSvgX, mouseSvgY, ratio) {
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
  this.addColorToEnvelopes = function(envList){
    envList.sort(function(x,y){
      return (x.env_peaks[0].mz - y.env_peaks[0].mz);
    })
    let colorNum = this.envColorList.length; 
    for (let i = 0; i < envList.length; i++) 
    {
      envList[i].color = this.envColorList[i%colorNum];
    }
  }
}
