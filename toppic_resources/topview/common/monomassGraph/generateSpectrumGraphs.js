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
// function svgIds()
// {
//     this.ms2SvgId = "ms2svg_";
//     this.ms2popUpSvgId = "ms2svg_popup_";
//     return this;
// }
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
    let graphFeatures = new GraphFeatures();
    if(startOfId == "monoMassSvg")
    {
        let calculatePrefixAndSuffixMassObj = new CalculatePrefixAndSuffixMass();
        let massShift_in = calculatePrefixAndSuffixMassObj.getIonTypeMass("B");
        let seq = calculatePrefixAndSuffixMassObj.getSequence(prsm_data);
        let massShiftList = calculatePrefixAndSuffixMassObj.getUnknownMassList();
        let prefixMassList = calculatePrefixAndSuffixMassObj.getPrefixMassList(seq,massShiftList,massShift_in);
        massShift_in = calculatePrefixAndSuffixMassObj.getIonTypeMass("Y");
        let suffixMassList = calculatePrefixAndSuffixMassObj.getSuffixMassList(seq,massShiftList,massShift_in);
        // Setting the graphFeatures object with all the features needed for the graph
        graphFeatures.showSequene = true;
        graphFeatures.addErrorPlot = true;
        graphFeatures.prefixSequenceData = prefixMassList;
        graphFeatures.suffixSequeceData = suffixMassList;
        graphFeatures.svgHeight = graphFeatures.svgHeight + graphFeatures.adjustableHeightVal + graphFeatures.heightForErrorPlot;
        graphFeatures.padding.head = graphFeatures.padding.head + graphFeatures.adjustableHeightVal;
        graphFeatures.padding.bottom = graphFeatures.padding.bottom + graphFeatures.heightForErrorPlot;
        graphFeatures.adjustableIonPosition = 10; // Random tested value for alignment
        // Gets the data list with mass error to plot in the monomass spectrum 
        graphFeatures.errorListData = json2ErrorDataList(prsm_data.prsm); 
        console.log(graphFeatures.errorListData);
        // Gets the absolute max and minimum value for upper bound and lower bound of y axis to draw the error plot
        graphFeatures.errorThreshHoldVal = getAbsoluteMaxValfromList(graphFeatures.errorListData);
        // Invoking spectrum function to draw the spectrum
        spectrumgraph = new addSpectrum(id, current_data, null, prec_mz, current_data,graphFeatures);
    }   
    else{
        let peak_data = new PeakData();
        // Get the formatted peak data list  
        let peak_list = peak_data.getPeakData(current_data);
        // Get formatted envelope data list
        let envelope_list = peak_data.getEnvelopeData(current_data);
        // Check if MS1 spectrum should be drawn
        if(id != "popupspectrum")
        {
            let ionData = peak_data.getIonData(prsm_data,specId,current_data);
            ms2_graph = new addSpectrum(id, peak_list, envelope_list, prec_mz, ionData, graphFeatures);
        }
        else{
            graphFeatures.isAddbgColor = true;
            // Get precursor mass to zoom to that point on launch of MS1 spectrum
            let precursor_mass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
            graphFeatures.bgMinMz = (graphFeatures.ratio * prec_mz) - graphFeatures.fixedWidthOfBgColorForMs1;
            graphFeatures.bgMaxMz = (graphFeatures.ratio * prec_mz) + graphFeatures.fixedWidthOfBgColorForMs1;
            spectrumgraph = new addSpectrum(id, peak_list, envelope_list, prec_mz, null,graphFeatures);
        }
    }
}

/**
 * Function to redraw the spectrum on popup window on click of download button for respected scan number
 * @param {Array} current_data - contains peak list and envelope list data of the corresponding spectrum(Scan Number) 
 * @param {String} id - Contains the Id of the SVG tag to which the spectrum need to be drawn
 * @param {Object} spectrumParameters - Contains spectrum parameters to apply same to the pop window spectrum
 * @param {Int} specId - Contains the Spec Id of the corresponding graph
 */
function reDrawWithSpecParams(current_data,id,spectrumParameters,specId)
{
    let peak_data = new PeakData();
    let peak_list = peak_data.getPeakData(current_data);
    let envelope_list = peak_data.getEnvelopeData(current_data);
    let ionData = peak_data.getIonData(prsm_data,specId,current_data);
    peak_list.sort(function(x,y){
		return d3.descending(x.intensity, y.intensity);
	})
    envelope_list = sortEnvelopes(envelope_list) ;
    let svgId = "#"+id;
    let peakData = {peak_list:peak_list, envelope_list:envelope_list};
    let graphFeatures = new GraphFeatures();
    new SpectrumGraph(svgId,spectrumParameters,peakData, ionData,graphFeatures);
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
 * Generating Jquery Onclick actions
 */
function graphOnClickActions(){

    // On click of download button on the main page of the prsm
    $("#graph_download").click(function(){
        // Get the Scan number of the Current spectrum graph showing on the screen
        let currentGraphNode = $('.ms2_scanIds.active')[0].text;
        let scanId = currentGraphNode.split(" ")[1];
        // Get the data tot he correcponding Scan Id
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanId);
        let id = "ms2svg_"+scanId;
        // Get SpectrumParametes and copy to a new variable. This way of copying doesn't pass reference rather creates new copy
        // correspondingSpecParams_g is a global variable contians the spec parameters of each spectrum 
        let specparams = jQuery.extend(true, {}, correspondingSpecParams_g[id]);
        let graphFeatures = new GraphFeatures();
        document.getElementsByName("show_envelops")[0].checked = graphFeatures.showCircles;
        document.getElementsByName("show_ions")[0].checked = graphFeatures.showIons;
        // Redraw the graph with correspondong specparameters
        reDrawWithSpecParams(current_data,"popup_ms2_spectrum",specparams,specId);
        // This allows pop up window to be moved
		$("#ms2spectrumpop").draggable({
			appendTo: "body"
		});
    })

    // ms2_scanIds is the Id of the nav tabs for multiple navs.
    // On Click shows corresponding graph by hiding others.
	$(".ms2_scanIds").click(function(){
        let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(ms2_ScansWithData,value);
        id = "ms2svg_"+value;
        // Hide all the graphs except the one clicked
        showCorrespondingGraph(id,".ms2_svg_graph_class");
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
    
    // Hide all othe graphs of monomass spectrum other than the one clicked
    $(".monoMass_scanIds").click(function(){
		let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(monoMassDataList,value);
        id = "monoMassSvg_"+value;
        // Hide all the graphs except the one clicked
        showCorrespondingGraph(id,".monoMass_svg_graph_class");
		$("#monoMass_nav .active").removeClass("active");
   		$(this).addClass("active");
    })
    
    // On click of mono mass mz, zoom all the graph to the corresponding point
    $(".peakRows").click(function() {
		let parent_id  = $(this).parent().parent().prop('id');
		let CurrentScanVal = document.getElementById(parent_id).firstChild.innerHTML;
		/*	get Mono M/z value till 3 decimal values	*/
		let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
		let [currentData,specId] = getCurrentData(ms2_ScansWithData,CurrentScanVal);
        let id = "ms2svg_"+CurrentScanVal;
        showCorrespondingGraph(id,".ms2_svg_graph_class");
        generateCorrespondingGraph(currentData,id,peak_value,specId);
        id = "monoMassSvg_"+CurrentScanVal;
        showCorrespondingGraph(id,".monoMass_svg_graph_class");
        [currentData,specId] = getCurrentData(monoMassDataList,CurrentScanVal);
        let CurrentMonoMassVal = $("#"+parent_id+" .row_monoMass").html();
        generateCorrespondingGraph(currentData,id,parseFloat(CurrentMonoMassVal),specId);
        activateCurrentnavbar("ms2_graph_nav",CurrentScanVal);
        activateCurrentnavbar("monoMass_nav",CurrentScanVal);
        
		showSpectrun();
    });
    // By changing the graph features, on click of redraw invoke this and generate the pop spectrun repectively to be downloaded
    $("#ms2_popup_redraw").click(function(){
        let currentGraphNode = $('.ms2_scanIds.active')[0].text;
        let scanId = currentGraphNode.split(" ")[1];
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanId);
        let id = "popup_ms2_spectrum";
        // Copying as a new variable than referencing. Referencing will change the properties of parent if child properties are changes
        // However, this is a shallow copying, we need to do this for all needed objects in side specparams
        let specparams = jQuery.extend(true, {}, correspondingSpecParams_g[id]);
        specparams.graphFeatures = jQuery.extend(true, {}, correspondingSpecParams_g[id].graphFeatures);
        specparams.graphFeatures.showCircles = document.getElementsByName("show_envelops")[0].checked ;
        specparams.graphFeatures.showIons = document.getElementsByName("show_ions")[0].checked ;
        reDrawWithSpecParams(current_data,"popup_ms2_spectrum",specparams,specId);
    })
    
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
