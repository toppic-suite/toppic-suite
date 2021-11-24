"use strict";
/**
 * @function SpectrumView
 * @description Function draws the graph, binds zoom and drag function to the graph
 * @param {String} svgId - SVG id on which the graph needed to be drawn
 * @param {object} spectrumParameters - Contains the parameters like height, width etc.,. tht helps to draw the graph
 * @param {Array} peakData - contains peakList and envelope list
 * @param {Array} ionData - Contains data with mass and ACID name to plot on the graph
 */
class SpectrumView {
    constructor(svgId, peakList, sequenceLength = -1) {
        // parameters for zoom
        this.transformX_ = 0;
        this.transformScale_ = 1.0;
        this.envList_ = [];
        this.ionList_ = null;
        this.proteoform_ = null;
        this.nIon_ = "";
        this.cIon_ = "";
        this.nMassList_ = [];
        this.cMassList_ = [];
        this.centerPos_ = -1; //center m/z or mono mass of the current view
        this.zoom = d3.zoom()
            .on("zoom", this.zoomed.bind(this));
        this.id_ = svgId;
        this.para_ = new SpectrumViewParameters();
        this.para_.setSeqLength(sequenceLength);
        this.peakList_ = peakList;
        this.para_.initParameters(peakList);
        this.peakList_.sort(function (x, y) {
            return y.getIntensity() - x.getIntensity();
        });
        $("#" + svgId).data("graph", this);
        // add zoom function
        this.svg_ = d3.select("body").select("#" + svgId);
        this.svg_.attr("viewBox", "0 0 " + this.para_.getSVGWidth() + " " + this.para_.getSVGHeight())
            .attr("width", "100%")
            .attr("height", "100%")
            //@ts-ignore
            .call(this.zoom.bind(this));
    }
    //getter
    getPeakList() {
        return this.peakList_;
    }
    getEnvList() {
        return this.envList_;
    }
    getIonList() {
        return this.ionList_;
    }
    getProteoform() {
        return this.proteoform_;
    }
    getNIon() {
        return this.nIon_;
    }
    getCIon() {
        return this.cIon_;
    }
    getNMassList() {
        return this.nMassList_;
    }
    getCMassList() {
        return this.cMassList_;
    }
    getPara() {
        return this.para_;
    }
    getTransformX() {
        return this.transformX_;
    }
    getTransformScale() {
        return this.transformScale_;
    }
    getSvgId() {
        return this.id_;
    }
    getCenterPos() {
        return this.centerPos_;
    }
    setTransformX(transformX) {
        this.transformX_ = transformX;
    }
    setTransformScale(transformScale) {
        this.transformScale_ = transformScale;
    }
    setCenterPos(newCenter) {
        this.centerPos_ = newCenter;
    }
    addRawSpectrumAnno(envList, ionList) {
        this.envList_ = envList;
        this.para_.addColorToEnvelopes(envList);
        //this.envPeakList = this.getEnvPeakList(this.envList);
        if (!ionList) {
            console.error("ERROR: invalid input for spectrum graph");
            return;
        }
        this.ionList_ = ionList;
    }
    addMonoMassSpectrumAnno(ionList, proteoform, nIonType, cIonType) {
        this.ionList_ = ionList;
        this.proteoform_ = proteoform;
        this.nIon_ = nIonType;
        this.cIon_ = cIonType;
        this.nMassList_ = proteoform.getNMasses(nIonType);
        this.cMassList_ = proteoform.getCMasses(cIonType);
    }
    redraw(monoMz) {
        if (this.para_.getIsMonoMassGraph() && monoMz) {
            this.para_.updateMassRange(monoMz);
        }
        else if (monoMz) {
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
    zoomed() {
        let transform = d3.event.transform;
        let graph = $("#" + this.getSvgId()).data("graph");
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
}
