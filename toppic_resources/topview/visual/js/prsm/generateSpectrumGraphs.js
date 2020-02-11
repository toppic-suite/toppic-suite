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
 * 4.   
 */
function svgIds()
{
    this.ms2SvgId = "ms2svg_";
    this.ms2popUpSvgId = "ms2svg_popup_";
    return this;
}
function generateCorrespondingGraph(current_data,id,prec_mz,specId){
    console.log("current_data : ", current_data );
    let peak_data = new PeakData();
    let peak_list = peak_data.getPeakData(current_data);
    let envelope_list = peak_data.getEnvelopeData(current_data);
    let ionData = peak_data.getIonData(prsm_data,specId);
    let graphParams = graphOptions();
    graphParams.scanId = current_data.scanId;
    if(id != "popupspectrum")
    {
        console.log("In graph");
        ms2_graph = addSpectrum(id, peak_list, envelope_list, prec_mz, ionData, graphParams);
    }
    else{
        spectrumgraph = addSpectrum(id, peak_list, envelope_list, prec_mz, null , graphParams);
    }
}
function createMultipleSvgs(current_data)
{
    //<svg id = "ms2svg" style = "background-color:#F8F8F8;display: none;"></svg>
    console.log("current_data : ", current_data);
    let div = document.getElementById("ms2svg_div"); 
    current_data.forEach(function(element,i){
        let id = "ms2svg_"+element.scanId;
        let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
        svg.setAttribute("id",id);
        svg.setAttribute("class","ms2_svg_graph_class");
        svg.style.backgroundColor = "#F8F8F8"; 
        if(i != 0)
        {
            svg.style.display = "none"; 
        }
        div.appendChild(svg);
        console.log("id : ",id);
        console.log("element.specId : ",element.specId);
        generateCorrespondingGraph(element.value,id,null,element.specId);
    });
}
function createMultiplePopUpSvgs(current_data)
{
    //<svg id = "ms2svg" style = "background-color:#F8F8F8;display: none;"></svg>
    console.log("current_data : ", current_data);
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
        console.log("id : ",id);
        console.log("element.specId : ",element.specId);
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
function showCorrespondingGraph(scanId)
{
    id = "ms2svg_"+scanId;
    $(".ms2_svg_graph_class").hide();
    document.getElementById(id).style = "block";
}


function graphOnClickActions(){

    $("#graph_download").click(function(){
		//$("#ms2_graph_popup_nav li").remove();
        //let scanIdList = prsm_data.prsm.ms.ms_header.scans.split(" ");
        let currentGraphNode = $('.ms2_scanIds.active')[0].text;
        let scanId = currentGraphNode.split(" ")[1];
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanId);
        let id = "ms2_popup_scanIds_"+scanId;
        let element = document.getElementById(id);
        console.log("element : ", element);
        $(".ms2_popup_scanIds.active").removeClass("active");
        element.classList.add("active");  
        
		$("#ms2spectrumpop").draggable({
			appendTo: "body"
		});
    })

    //ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
        let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(ms2_ScansWithData,value);
        console.log("value : ", value);
        id = "ms2svg_"+value;
        showCorrespondingGraph(value);
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	$(".ms1_scanIds").click(function(){
		let value = this.getAttribute('value');
		let [currentData,specId] = getCurrentData(ms1_ScansWithData,value);
		generateCorrespondingGraph(currentData,"popupspectrum",null,specId);
		$("#ms1_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
    })
    
    $(".peakRows").click(function() {
		let parent_id  = $(this).parent().parent().prop('id');
		let CurrentScanVal = document.getElementById(parent_id).firstChild.innerHTML;
		/*	get Mono M/z value till 3 decimal values	*/
		let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
		console.log("ms2_ScansWithData in buttons: ", ms2_ScansWithData);
		console.log("CurrentScanVal : ", CurrentScanVal);
		let [currentData,specId] = getCurrentData(ms2_ScansWithData,CurrentScanVal);
		console.log("currentData : ", currentData);
        let id = "ms2svg_"+CurrentScanVal;
        showCorrespondingGraph(CurrentScanVal);
		generateCorrespondingGraph(currentData,id,peak_value,specId);
		activateCurrentnavbar("ms2_graph_nav",CurrentScanVal)
		showSpectrun();
    });
    
}
