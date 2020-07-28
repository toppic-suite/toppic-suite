/**
 * Generates Corresponding graph on to the html by calling addSpectrum from invokespectrum.js
 * @param {Array} peakDataList - Contains Peak list data
 * @param {Array} distributionList - Contains envelope data
 * @param {Float} prec_mz - contains mass to which the graph has to be zoomed
 */
function generateCorrespondingGraph(peakDataList,distributionList,prec_mz){
    let ionData = null;
    let graphFeatures = new GraphFeatures();
    ms2_graph = addSpectrum("spectrum",peakDataList,distributionList,prec_mz,ionData,graphFeatures);
}
/**
 * Generating Mono Mass Graph with calculated masses
 * @param {Array} monoMassList - Contians the mono Mass list data
 * @param {Float} prec_mz - Contains mass to which the graph has to be zoomed
 */
function generateMonoMassGraph(monoMassList,prec_mz){
    let ionData = [];
    //Draw MonoMass Graph
    const monoMassList_mz = monoMassList.map((element)=>{
        let tempElement = {
            "mz":element.mass,
            "intensity":element.intensity,
            "charge": element.charge,
            "peakId": element.peakId,
            "ion": element.ion,
            "position": element.position,
            "massError": element.massError,
            "thMass": element.thMass,
            "PPMerror": element.PPMerror,
            "matchedInd": element.matchedInd
        };
        let tempIonData = {"mz":element.mass,"intensity":element.intensity,"ion": element.ion};
        ionData.push(tempIonData);
        return tempElement;
    });
    let graphFeatures = new GraphFeatures();
    graphFeatures.showSequence = true;
    let prefixSequenceData = [];
    let suffixSequeceData = [];
    completeCalData.prefixmasslist.forEach((eachElem,index)=>{
        eachElem["acid"] = completeCalData.sequence[index];
        prefixSequenceData.push(eachElem);
    });
    let len = completeCalData.sequence.length;
    completeCalData.suffixmasslist.forEach((eachElem,index)=>{
        eachElem["acid"] = completeCalData.sequence[len-index-1];
        suffixSequeceData.push(eachElem);
    })
    graphFeatures.addErrorPlot = true;
    graphFeatures.prefixSequenceData = prefixSequenceData;
    graphFeatures.suffixSequeceData = suffixSequeceData;
    graphFeatures.svgHeight = graphFeatures.svgHeight + graphFeatures.adjustableHeightVal + graphFeatures.heightForErrorPlot;
    graphFeatures.padding.head = graphFeatures.padding.head + graphFeatures.adjustableHeightVal;
    graphFeatures.padding.bottom = graphFeatures.padding.bottom + graphFeatures.heightForErrorPlot;
    graphFeatures.adjustableIonPosition = 10;
    let errorListData = completeCalData.matchedandunmatcheddata.map((element)=>{
        let errorListObj = {"theoretical_mass":element.thMass,"mass_error":element.massError};
        return errorListObj;
    })
    graphFeatures.errorListData = errorListData;
    graphFeatures.errorThreshHoldVal = getAbsoluteMaxValfromList(graphFeatures.errorListData);
    
    addSpectrum("monoMassGraph",monoMassList_mz,null,prec_mz,ionData,graphFeatures);
}
/**
 * Function return the max error value from the list
 * @param {Array} errorDataList - Contains list of error data to plot on the error graph
 */
function getAbsoluteMaxValfromList(errorDataList){
    let max = 0;
    errorDataList.forEach((element)=>{
        let val = Math.abs(element.mass_error);
        if(max < val) max = val; 
    })
    //Getting the round off fraction value
    max = max * 100;
    max = Math.ceil(max)/100;
    return max;
}