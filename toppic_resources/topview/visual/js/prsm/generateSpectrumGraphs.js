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
function svgIds()
{
    this.ms2SvgId = "ms2svg_";
    this.ms2popUpSvgId = "ms2svg_popup_";
    return this;
}
function generateCorrespondingGraph(current_data,id,prec_mz,specId){
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
        graphFeatures.showSequene = true;
        graphFeatures.addErrorPlot = true;
        graphFeatures.prefixSequenceData = prefixMassList;
        graphFeatures.suffixSequeceData = suffixMassList;
        graphFeatures.svgHeight = graphFeatures.svgHeight + graphFeatures.adjustableHeightVal + graphFeatures.heightForErrorPlot;
        graphFeatures.padding.head = graphFeatures.padding.head + graphFeatures.adjustableHeightVal;
        graphFeatures.padding.bottom = graphFeatures.padding.bottom + graphFeatures.heightForErrorPlot;
        graphFeatures.adjustableIonPosition = 10;
        graphFeatures.errorListData = json2ErrorDataList(prsm_data.prsm);
        console.log(graphFeatures.errorListData);
        graphFeatures.errorThreshHoldVal = getAbsoluteMaxValfromList(graphFeatures.errorListData);
        // Current Data itself contains Peak data
        console.log("current_data : ", current_data);
        console.log("graphFeatures : ", graphFeatures);
        spectrumgraph = new addSpectrum(id, current_data, null, prec_mz, current_data,graphFeatures);
    }   
    else{
        let peak_data = new PeakData();
        let peak_list = peak_data.getPeakData(current_data);
        let envelope_list = peak_data.getEnvelopeData(current_data);
        if(id != "popupspectrum")
        {
            let ionData = peak_data.getIonData(prsm_data,specId,current_data);
            ms2_graph = new addSpectrum(id, peak_list, envelope_list, prec_mz, ionData, graphFeatures);
        }
        else{
            graphFeatures.isAddbgColor = true;
            let precursor_mass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
            graphFeatures.bgMinMz = (graphFeatures.ratio * prec_mz) - graphFeatures.fixedWidthOfBgColorForMs1;
            graphFeatures.bgMaxMz = (graphFeatures.ratio * prec_mz) + graphFeatures.fixedWidthOfBgColorForMs1;
            spectrumgraph = new addSpectrum(id, peak_list, envelope_list, prec_mz, null,graphFeatures);
        }
    }
}
function getMinMzFromSpectrumData(envelope_list,precursor_mass)
{ 
    let env_peaks = [];
    let minMz = 0;
    let maxMz = 0;
    envelope_list.forEach((element)=>{
        if(parseFloat(element.mono_mass).toFixed(3) == parseFloat(precursor_mass).toFixed(3))
        {
            env_peaks = element.env_peaks;
        }
    });
    if(env_peaks.length != 0)
    {
        env_peaks.sort(function(a, b){
            return a.mz-b.mz;
        })
    }
    minMz = env_peaks[0].mz;
    maxMz = env_peaks[env_peaks.length-1].mz;
    return[minMz,maxMz];
}
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
function createMultipleSvgs(divId,svgId,className,current_data){
    let div = document.getElementById(divId); 
    current_data.forEach(function(element,i){
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

function createMultiplePopUpSvgs(current_data)
{
    let div = document.getElementById("ms2_graph_popup_nav_div"); 
    current_data.forEach(function(element,i){
        let id = "ms2svg_popup_"+element.scanId;
        let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
        svg.setAttribute("id",id);
        svg.setAttribute("class","ms2_svg_graph_popup_class");
        svg.style.backgroundColor = "#F8F8F8"; 
        if(i != 0)
        {
            svg.style.display = "none"; 
        }
        div.appendChild(svg);
        generateCorrespondingGraph(element.value,id,null,element.specId);
    });
}
function createDummySpecParamsList(current_data){
    let spectrumParams_global = [];
    let specParams = {};
    let tempSpecParams = {};
    current_data.forEach(function(element,i){
        tempSpecParams = {"scanId":element.scanId,"specParams":specParams};
        spectrumParams_global.push(tempSpecParams);
    })
    return spectrumParams_global;
}
function showCorrespondingGraph(id,className)
{
    $(className).hide();
    document.getElementById(id).style = "block";
}
function graphOnClickActions(){

    $("#graph_download").click(function(){
        let currentGraphNode = $('.ms2_scanIds.active')[0].text;
        let scanId = currentGraphNode.split(" ")[1];
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanId);
        let id = "ms2svg_"+scanId;
        let specparams = jQuery.extend(true, {}, correspondingSpecParams_g[id]);//correspondingSpecParams_g[id];
        let graphFeatures = new GraphFeatures();
        document.getElementsByName("show_envelops")[0].checked = graphFeatures.showCircles;
        document.getElementsByName("show_ions")[0].checked = graphFeatures.showIons;
        reDrawWithSpecParams(current_data,"popup_ms2_spectrum",specparams,specId);
		$("#ms2spectrumpop").draggable({
			appendTo: "body"
		});
    })

    //ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
        let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(ms2_ScansWithData,value);
        id = "ms2svg_"+value;
        showCorrespondingGraph(id,".ms2_svg_graph_class");
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	
    $(".monoMass_scanIds").click(function(){
		let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(monoMassDataList,value);
        id = "monoMassSvg_"+value;
        showCorrespondingGraph(id,".monoMass_svg_graph_class");
		//generateCorrespondingGraph(currentData,"monoMassSvg_",null,specId);
		$("#monoMass_nav .active").removeClass("active");
   		$(this).addClass("active");
    })
    
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

    $("#ms2_popup_redraw").click(function(){
        let currentGraphNode = $('.ms2_scanIds.active')[0].text;
        let scanId = currentGraphNode.split(" ")[1];
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanId);
        let id = "popup_ms2_spectrum";
        //Copying as a new variable than referencing. Referencing will change the properties of parent if child properties are changes
        //However, this is a shallow copying, we need to do this for all needed objects in side specparams
        let specparams = jQuery.extend(true, {}, correspondingSpecParams_g[id]);
        specparams.graphFeatures = jQuery.extend(true, {}, correspondingSpecParams_g[id].graphFeatures);
        specparams.graphFeatures.showCircles = document.getElementsByName("show_envelops")[0].checked ;
        specparams.graphFeatures.showIons = document.getElementsByName("show_ions")[0].checked ;
        reDrawWithSpecParams(current_data,"popup_ms2_spectrum",specparams,specId);
    })
    
}
/**
 * Function return the max error value from the list
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
