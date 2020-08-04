/**
 * Available functions in the file are 
 *      {generateCorrespondingGraph}, {graphOnClickActions}
 * Generating multiple graphs with different scan Id.
 * 1.   Once the data is loaded from local files from the function {promiseLoadDataJS} from multiscan.js,
 *      it sets the data into global variables ms2_ScansWithData and ms1_ScansWithData.
 * 2.   Once the data is set to the global variables, then it creates navigation elements by calling 
 *      {createMs1NavEements}, {createMs2NavEements}. 
 * 3.   Once navigation elements are generated, {graphOnClickActions} function is called to set the 
 *      on click actions on the navigation buttons of the spectrum graph
 */

/**
 * 
 * @param {Array} current_data - Contains a list with ion data if it is generating Mono Mass spectrum, 
 * contians peaklist and enveloplist if generating spectrums of MS1/MS2
 * @param {String} id - Contains the id of spectrum svg with _ followed by scan number (monoMassSvg_872/ms2svg_870)
 * @param {Float} prec_mz - Contains the mass to which the graph needs to be zoomed
 * @param {Int} specId - Contains the spectrum number
 */
function generateCorrespondingGraph(current_data,id,prec_mz,specId){
    // Gets the svg id of the spectrum
    let startOfId = id.split("_")[0];
    // console.log("current_data:", current_data);
    let graphFeatures = new GraphFeatures();
    if(startOfId == "monoMassSvg")
    {
        let sequenceObj = new Sequence(prsm_data);
        let sequence = sequenceObj.getSequence();

        let ionMassShiftObj = new IonMassShift("B");
        let ionShift = ionMassShiftObj.getIonTypeMass();

        let massShiftListObj = new MassShiftList(prsm_data);
        let massShiftList = massShiftListObj.getMassShiftList();

        let cpsmObj = new CalcPrefixSuffixMassList(sequence, massShiftList);
        let prefixMassList = cpsmObj.getPrefixMassList(ionShift);

        ionMassShiftObj.ionType = "Y";
        ionShift = ionMassShiftObj.getIonTypeMass();

        let suffixMassList = cpsmObj.getSuffixMassList(ionShift);

        // Setting the graphFeatures object with all the features needed for the graph
        graphFeatures.showSequence = true;
        graphFeatures.addErrorPlot = true;
        graphFeatures.prefixSequenceData = prefixMassList;
        graphFeatures.suffixSequeceData = suffixMassList;
        graphFeatures.svgHeight = graphFeatures.svgHeight + graphFeatures.adjustableHeightVal + graphFeatures.heightForErrorPlot;
        graphFeatures.padding.head = graphFeatures.padding.head + graphFeatures.adjustableHeightVal;
        graphFeatures.padding.bottom = graphFeatures.padding.bottom + graphFeatures.heightForErrorPlot;
        graphFeatures.adjustableIonPosition = 10; // Random tested value for alignment
        // Gets the data list with mass error to plot in the monomass spectrum 
        let prsmDataUtilObj = new PrsmDataUtil(prsm_data);
        graphFeatures.errorListData = prsmDataUtilObj.json2ErrorDataList(); 
        // console.log(graphFeatures.errorListData);
        // Gets the absolute max and minimum value for upper bound and lower bound of y axis to draw the error plot
        graphFeatures.errorThreshHoldVal = getAbsoluteMaxValfromList(graphFeatures.errorListData);
        // Invoking spectrum function to draw the spectrum
        spectrumgraph = new addSpectrum(id, current_data, null, prec_mz, current_data,graphFeatures);
    }
}

/**
 * This generates spectrum for each spec Id
 * @param {String} divId - Contains Id of the div tag under which the monomass graphs are drawn
 * @param {String} svgId - Contains id as "monoMassSvg_" to which scan Id is added
 * @param {String} className - Contains class name to which the corresponding svg graphs are drawn 
 * @param {Array} data - contains the data of all the scans available
 */
function createMultipleSvgs(divId,svgId,className,data){
    let div = document.getElementById(divId); 
    data.forEach(function(element,i){
        let id = svgId+element.scanId;
        let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
        svg.setAttribute("id",id);
        svg.setAttribute("class",className);
        svg.style.backgroundColor = "#F8F8F8"; 
        if(i != 0)
        {
            svg.style.display = "none"; 
        }
        div.appendChild(svg);
        generateCorrespondingGraph(element.value,id,null,element.specId);
    });
}

/**
 * Function to hide graphs of all other spectrums other than current showing spectrum
 * @param {String} id - Contains Id of the SVG tag to be shown
 * @param {String} className - className of the spectrums to be hidden 
 */
function showCorrespondingGraph(id,className)
{
    $(className).hide();
    document.getElementById(id).style = "block";
}
/**
 * Function return the absolute max error value from the list
 * @param {Array} errorDataList - Contains the list of mass errors
 */
function getAbsoluteMaxValfromList(errorDataList){
    let max = 0;
    errorDataList.forEach((element,index)=>{
        let val = Math.abs(element.mass_error);
        if(max < val) max = val; 
    })
    //Getting the round off fraction value
    max = max * 100;
    max = Math.ceil(max)/100;
    console.log("max : ", max);
    return max;
}
