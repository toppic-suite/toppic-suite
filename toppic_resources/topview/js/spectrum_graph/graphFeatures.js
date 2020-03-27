class GraphFeatures{
    constructor(){
        this.showCircles = true;
        this.showIons = true;
        this.showSequene = false;
        this.isAddbgColor = false;
        this.fixedWidthOfBgColorForMs1 = 2;
        this.ratio = 1.000684;
        this.tickWidthThreshholdval = 0.5;
        this.bgMinMz = 0;
        this.bgMaxMz = 0;
        this.bgColor = "orange";
        this.prefixSequenceData = [];
        this.suffixSequeceData = [];
        this.adjustableIonPosition = 4;
        this.svgWidth = 910;
        this.svgHeight = 220;
        this.padding = {left:70, right:20, head:10, bottom:50};
        this.adjustableHeightVal = 60;
        this.fixedHeightOfIonAboveThePeak = 10;
        this.specWidth = this.svgWidth - this.padding.left - this.padding.right;
        this.specHeight = this.svgHeight - this.padding.head - this.padding.bottom;
    }
}