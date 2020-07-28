/** Class with variables to set graph properties */
class GraphFeatures{
    constructor(){
        this.showCircles = true;
        this.showIons = true;
        this.showSequence = false;
        this.isAddbgColor = false;
        this.addErrorPlot = false;
        this.fixedWidthOfBgColorForMs1 = 2;
        this.ratio = 1.000684;
        this.tickWidthThreshholdval = 0.5;
        this.bgMinMz = 0;
        this.bgMaxMz = 0;
        this.bgColor = "orange";
        this.prefixSequenceData = [];
        this.suffixSequeceData = [];
        this.errorListData = [];
        this.adjustableIonPosition = 4;
        this.svgWidth = 910;
        this.svgHeight = 220;
        this.padding = {left:70, right:20, head:10, bottom:50};
        this.errorplot_padding = {left:70, right:20, head:10, bottom:10};
        this.adjustableHeightVal = 40;
        this.fixedHeightOfIonAboveThePeak = 10;
        this.heightForErrorPlot = 60;
        this.errorThreshHoldVal = 0.2;
        this.specWidth = this.svgWidth - this.padding.left - this.padding.right;
        this.specHeight = this.svgHeight - this.padding.head - this.padding.bottom;
    }
}