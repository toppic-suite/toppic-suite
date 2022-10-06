"use strict";
/**	@function SpectrumViewParameters
 * @description Get data from global variable spectrum_data and utilities to manupulate
 * the data
 */
class SpectrumViewParameters {
    constructor() {
        // Ratio between average and monoisopotic mass
        this.avgToMonoRatio_ = 1.000684;
        // SVG size
        this.svgWidth_ = 910;
        this.svgHeight_ = 270;
        // SVG padding 
        this.padding_ = { left: 70, right: 50, head: 10, bottom: 50 };
        // spectrum size
        this.specWidth_ = this.svgWidth_ - this.padding_.left - this.padding_.right;
        this.specHeight_ = this.svgHeight_ - this.padding_.head - this.padding_.bottom;
        // M/z range of visuable window
        this.winMinMz_ = 0;
        this.winMaxMz_ = 2000;
        this.winCenterMz_ = 1000;
        //minimum possible m/z after zooming/dragging, to prevent dragging/zooming into negative m/z value
        //if new m/z is less than this value, it is reset to this value
        this.minPossibleMz_ = -100;
        this.maxPossibleMzMargin_ = 300;
        // M/z range of peaks
        this.dataMinMz_ = 0;
        this.dataMaxMz_ = 2000;
        // M/z range, color of highlighted part.
        this.showHighlight_ = false;
        this.hlMinMz_ = 0;
        this.hlMaxMz_ = 0;
        this.hlColor_ = "gray";
        // Max intensity of visuable window
        this.winMaxInte_ = 30000;
        // Intensity range of peaks
        this.dataMaxInte_ = 30000;
        this.dataMinInte_ = 0;
        // add a margin so that the visuable intensity range is [0, dataMaxInte * inteMargin]
        this.inteMargin_ = 1.2;
        // scale m/z to x coordinate
        this.xScale_ = 0.35;
        // scale intensity to y coordinate
        this.yScale_ = 0.005;
        // Numbers of ticks
        this.xTickNum_ = 10;
        this.yTickNum_ = 5;
        this.tickLength_ = 7;
        // Tick width list used in the function getTickWidth
        this.tickWidthList_ = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1000, 800, 700, 600, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 20, 10, 5, 3, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001];
        // Tick height list used in the function getTickHeight
        this.tickHeightList_ = [50, 40, 30, 25, 20, 15, 10, 5, 3, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001];
        //Limiting the peaks and envelopes to 4000 using 20 bins
        this.binNum_ = 20;
        this.peakNumPerBin_ = 50;
        //Padding for mouse over peak floatings.
        this.mouseOverPadding_ = { head: 20, middle: 14 };
        // Envelope circle size: min and max radius	
        this.showEnvelopes_ = true;
        this.defaultRadius_ = 0.05;
        this.minRadius_ = 2;
        this.maxRadius_ = 5;
        //	Colors for the envelope circles	
        this.envColorList_ = ["red", "darkorange", "blue"];
        // Parameters related to annoated ions
        this.showIons_ = true;
        this.ionXShift_ = -5;
        this.ionYShift_ = -15;
        // Mono mass graph
        this.showError_ = true;
        this.showLines_ = true;
        this.isMonoMassGraph_ = false;
        this.errorPlotPadding_ = { left: 70, right: 50, head: 10, bottom: 10 };
        this.errorPlotHeight_ = 40;
        this.errorThreshold_ = 0.2;
        this.errorYTickNum_ = 2;
        //restrict x zoom
        this.isXZoomAllowed_ = true;
        //sequence length
        //for determining max m/z window based on seq length in mass graph
        this.seqLength_ = -1;
    }
    //getters and setters
    getErrorYTickNum() {
        return this.errorYTickNum_;
    }
    getErrorThreshold() {
        return this.errorThreshold_;
    }
    getErrorPlotPadding() {
        return this.errorPlotPadding_;
    }
    getErrorPlotHeight() {
        return this.errorPlotHeight_;
    }
    getIonXShift() {
        return this.ionXShift_;
    }
    getIonYShift() {
        return this.ionYShift_;
    }
    getDataMaxMz() {
        return this.dataMaxMz_;
    }
    getDataMinMz() {
        return this.dataMinMz_;
    }
    getHlColor() {
        return this.hlColor_;
    }
    getHlMinMz() {
        return this.hlMinMz_;
    }
    getHlMaxMz() {
        return this.hlMaxMz_;
    }
    getDataMaxInte() {
        return this.dataMaxInte_;
    }
    getWinMaxMz() {
        return this.winMaxMz_;
    }
    getWinMinMz() {
        return this.winMinMz_;
    }
    getWinMaxInte() {
        return this.winMaxInte_;
    }
    getWinCenterMz() {
        return this.winCenterMz_;
    }
    getTickLength() {
        return this.tickLength_;
    }
    getYTickNum() {
        return this.yTickNum_;
    }
    getSVGWidth() {
        return this.svgWidth_;
    }
    getSVGHeight() {
        return this.svgHeight_;
    }
    getIsMonoMassGraph() {
        return this.isMonoMassGraph_;
    }
    getShowHighlight() {
        return this.showHighlight_;
    }
    getShowEnvelopes() {
        return this.showEnvelopes_;
    }
    getShowIons() {
        return this.showIons_;
    }
    getShowError() {
        return this.showError_;
    }
    getShowLines() {
        return this.showLines_;
    }
    getSpecHeight() {
        return this.specHeight_;
    }
    getPadding() {
        return this.padding_;
    }
    getIsXZoomAllowed() {
        return this.isXZoomAllowed_;
    }
    setSeqLength(seqLength) {
        this.seqLength_ = seqLength;
    }
    setShowEnvelopes_(showEnvelopes) {
        this.showEnvelopes_ = showEnvelopes;
    }
    setShowIons(showIons) {
        this.showIons_ = showIons;
    }
    setShowError(showError) {
        this.showError_ = showError;
    }
    setShowLines(showLines) {
        this.showLines_ = showLines;
    }
    setIsXZoomAllowed_(allowXZoom) {
        this.isXZoomAllowed_ = allowXZoom;
    }
    setSVGHeight(newHeight) {
        this.svgHeight_ = newHeight;
    }
    setPadding(left, right, head, bottom) {
        this.padding_.left = left;
        this.padding_.right = right;
        this.padding_.head = head;
        this.padding_.bottom = bottom;
    }
    setSpecHeight(height) {
        this.specHeight_ = height;
    }
    /**
     * @function getTickWidth
     * @description Function Provides width between each tick when zoomed in and out or dragged
     */
    getTickWidth() {
        let tempDiff = this.winMaxMz_ - this.winMinMz_;
        let tickWidth = Math.floor(this.tickWidthList_[0]);
        for (let i = 0; i < this.tickWidthList_.length; i++) {
            if (tempDiff / this.xTickNum_ <= Math.floor(this.tickWidthList_[i]) &&
                tempDiff / this.xTickNum_ > Math.floor(this.tickWidthList_[i + 1])) {
                tickWidth = Math.floor(this.tickWidthList_[i]);
                break;
            }
        }
        return tickWidth;
    }
    getXTickPosList() {
        let posList = new Array(this.xTickNum_ + 1);
        let tickWidth = this.getTickWidth();
        for (let i = 0; i <= this.xTickNum_; i++) {
            // calculate the actual tick position based on the current minMz value on the xaxis
            let tickMz = 0;
            if (tickWidth < 1 && tickWidth != 0) {
                tickMz = (i * tickWidth + this.winMinMz_) - Math.floor((i * tickWidth + this.winMinMz_) % tickWidth);
            }
            else if (tickWidth != 0) {
                tickMz = i * tickWidth + this.winMinMz_ - (i * tickWidth + this.winMinMz_) % tickWidth;
            }
            posList[i] = tickMz;
        }
        return posList;
    }
    /**
     * @function getTickHeight
     * @description Function Provides height between each tick when zoomed in and out or dragged
     */
    getTickHeight() {
        let tickheight = Math.floor(this.tickHeightList_[0]);
        let maxIntPercent = this.winMaxInte_ / this.dataMaxInte_ * 100;
        for (let i = 0; i < this.tickHeightList_.length; i++) {
            if (maxIntPercent / this.yTickNum_ <= this.tickHeightList_[i]
                && maxIntPercent / this.yTickNum_ > this.tickHeightList_[i + 1]) {
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
    getPeakXPos(mz) {
        let peakX = (mz - this.winMinMz_) * this.xScale_ + this.padding_.left;
        return peakX;
    }
    /**
     * @function getPeakYPos
     * @description Function provides the y coordinate for the intensity
     */
    getPeakYPos(intensity) {
        let peakY = this.svgHeight_ - intensity * this.yScale_ - this.padding_.bottom;
        return peakY;
    }
    /**
     * @function getErrorYPos
     * @description Function provides the y coordinate for the error val on the error plot
     */
    getErrorYPos(errorVal) {
        // Multiply with 2 as the coordinates has to be both positive and negative
        let yErrorScale = this.errorPlotHeight_ / (this.errorThreshold_ * 2);
        let pos = this.svgHeight_ - (errorVal * yErrorScale)
            - this.errorPlotPadding_.bottom - this.errorPlotHeight_ / 2;
        return pos;
    }
    /**
     * @function getBinWidth
     * @description Function to compute bin width
     **/
    getBinWidth() {
        let width = (this.winMaxMz_ - this.winMinMz_) / this.binNum_;
        return width;
    }
    /**
     * @function getCircleSize
     * @description Function provides the radius of the circles drawn on the graph as zoomed in and out
     */
    getCircleSize() {
        let radius = this.defaultRadius_ * this.xScale_;
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
    compDataRanges(peakList) {
        let minMz = 0;
        let maxMz = 2000;
        let maxInte = 100;
        if (peakList != null && peakList.length > 0) {
            // Sort by mz
            peakList.sort(function (x, y) {
                return x.getPos() - y.getPos();
            });
            let listSize = peakList.length;
            maxMz = Math.floor(peakList[listSize - 1].getPos());
            // Sort by intensity
            peakList.sort(function (x, y) {
                return x.getIntensity() - y.getIntensity();
            });
            maxInte = Math.floor(peakList[listSize - 1].getIntensity());
        }
        return [minMz, maxMz, maxInte];
    }
    /**
     * @function updateScale
     * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
     * initializing xScale, yScale.
     */
    updateScale(winMinMz, winMaxMz, winMaxInte) {
        this.winMinMz_ = winMinMz;
        this.winMaxMz_ = winMaxMz;
        if (winMinMz == this.dataMinMz_ && winMaxMz == this.dataMaxMz_) {
            this.winMinMz_ = 0;
            this.winMaxMz_ = 1.1 * this.winMaxMz_;
        }
        this.winCenterMz_ = (this.winMinMz_ + this.winMaxMz_) / 2.0;
        this.xScale_ = this.specWidth_ / (this.winMaxMz_ - this.winMinMz_);
        this.winMaxInte_ = winMaxInte;
        this.yScale_ = this.specHeight_ / this.winMaxInte_;
    }
    /**
     * @function initParameters
     * @description Initializing the spectrum Parameters with the data from the peak list and envilopelist.
     * initializing xScale, yScale.
     */
    initParameters(peakList) {
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
    drag(distX) {
        let mzDist = distX / this.xScale_;
        this.winMinMz_ = this.winMinMz_ - mzDist;
        this.winMaxMz_ = this.winMaxMz_ - mzDist;
        this.winCenterMz_ = this.winCenterMz_ - mzDist;
        //allow drag up to -50 m/z (this.minPossibleMz) to give some padding 
        if (this.winMinMz_ - mzDist < this.minPossibleMz_) {
            let minMaxDiff = this.winMaxMz_ - this.winMinMz_;
            let centerDiff = this.winCenterMz_ - this.winMinMz_;
            this.winMinMz_ = this.minPossibleMz_;
            this.winMaxMz_ = this.winMinMz_ + minMaxDiff;
            this.winCenterMz_ = this.winMinMz_ + centerDiff;
        }
        if (this.winMaxMz_ > this.dataMaxMz_ + this.maxPossibleMzMargin_) {
            let minMaxDiff = this.winMaxMz_ - this.winMinMz_;
            this.winMaxMz_ = this.dataMaxMz_ + this.maxPossibleMzMargin_;
            this.winMinMz_ = this.winMaxMz_ - minMaxDiff;
            let centerDiff = this.winMaxMz_ - this.winCenterMz_;
            this.winCenterMz_ = this.winMinMz_ + centerDiff;
        }
    }
    /**
     * @function xZoom
     * @description Function provides with current xScale, current minMz and MaxMz based on the zoom on x-axis.
     * Function also calls setLimita which helps in drawing limited number of peaks and circles per eachbin/range of mz values.
     */
    xZoom(mouseSvgX, ratio) {
        if (!this.isXZoomAllowed_) {
            return;
        }
        let oriValues = {}; //so that the view range can be restored when the view shouldn't be zoomed
        oriValues.min = this.winMinMz_;
        oriValues.max = this.winMaxMz_;
        oriValues.center = this.winCenterMz_;
        oriValues.xScale = this.xScale_;
        let mouseSpecX = mouseSvgX - this.padding_.left;
        this.winCenterMz_ = mouseSpecX / this.xScale_ + this.winMinMz_;
        /*self is a global variable of datasource object containing all the data needed to use when zoomed*/
        this.xScale_ = this.xScale_ * ratio;
        this.winMinMz_ = this.winCenterMz_ - mouseSpecX / this.xScale_;
        this.winMaxMz_ = this.winCenterMz_ + (this.specWidth_ - mouseSpecX) / this.xScale_;
        //console.log(this.winMaxMz_, this.dataMaxMz_ + 500)
        if (this.winMinMz_ < this.minPossibleMz_) { //prevent zooming out into negative mass
            this.winMinMz_ = this.minPossibleMz_;
        }
        if (this.winMaxMz_ > this.dataMaxMz_ + this.maxPossibleMzMargin_) {
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
    yZoom(ratio) {
        //Reducing zoom factor to smoothenup and remove gliches
        if (ratio > 1)
            ratio = 1.4;
        else if (ratio < 1)
            ratio = 0.9;
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
    zoom(mouseSvgX, mouseSvgY, ratio) {
        if (ratio > 1)
            ratio = 1.4; // Zooming in and fixing ration to 1.4 (fixed values based on testing the smooting of zoom)
        else if (ratio < 1)
            ratio = 0.9; // Zooming out and fixing ration to 0.9 (fixed values based on testing the smooting of zoom)
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
    addColorToEnvelopes(envList) {
        if (!envList || envList.length === 0 || typeof envList[0].getPeaks() === "undefined")
            return;
        envList.sort(function (x, y) {
            return (x.getPeaks()[0].getPos() - y.getPeaks()[0].getPos());
        });
        let colorNum = this.envColorList_.length;
        for (let i = 0; i < envList.length; i++) {
            envList[i].setDisplayColor(this.envColorList_[i % colorNum]);
        }
    }
    /**
     * @function setHighlight
     * @description
     * set highlight region for MS1 precursor envelope
     */
    setHighlight(ms1Spec) {
        this.showHighlight_ = true;
        this.hlMinMz_ = ms1Spec.getMinMz();
        this.hlMaxMz_ = ms1Spec.getMaxMz();
        //console.log(precMonoMz, this.hlMinMz, this.hlMaxMz);
    }
    setDefaultPadding() {
        this.padding_.head = 20;
        this.padding_.bottom = 50;
    }
    setMonoMassPaddding() {
        this.padding_.head = 60;
        this.padding_.bottom = 75;
    }
    setMonoMassGraph(isMonoMass) {
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
    updateMzRange(monoMz) {
        let centerMz = monoMz * this.avgToMonoRatio_;
        this.winMinMz_ = centerMz - 3;
        this.winMaxMz_ = centerMz + 3;
        this.updateScale(this.winMinMz_, this.winMaxMz_, this.winMaxInte_);
    }
    updateMassRange(mass) {
        let centerMass = mass;
        this.winMinMz_ = centerMass - 3;
        this.winMaxMz_ = centerMass + 3;
        this.updateScale(this.winMinMz_, this.winMaxMz_, this.winMaxInte_);
    }
}
