
class MultiScan{
    constructor(){
    }
    //ms1_ms2_json this give ms1_json folder or ms2_json folder
    promiseLoadDataJS(specidList,scanIdList,ms1_ms2_json){
        //let timeoutVal = 2000;
        let len = specidList.length;
        let scanListWithData = [];
        let count = 0;
        //get the list of lists with scan Id and value
        specidList.forEach(function(element,i){
            let temp_filename = "../../topfd/"+ms1_ms2_json+"/spectrum"+element+".js";
            let temp_script= document.createElement('script');
            temp_script.src = temp_filename;
            document.head.appendChild(temp_script);
            temp_script.onload = function(){
                let temp_scanid ;
                let temp_value ;
                if(ms1_ms2_json == "ms1_json")
                {
                    temp_scanid = ms1_data.scan;
                    temp_value = ms1_data;
                }else{
                    temp_scanid = ms2_data.scan;
                    temp_value = ms2_data;
                }
                let temp_ms2DataObj = {"scanid":temp_scanid,"specId":element,"value":temp_value};
                scanListWithData.push(temp_ms2DataObj);
                count++;
                if(count == len)
                {
                    let MultiScanObj = new MultiScan();
                    let current_data ;
                    let specId;
                    [current_data,specId] = getCurrentData(scanListWithData,scanIdList[0]);
                    
                    if(ms1_ms2_json == "ms1_json")
                    {
                        ms1_ScansWithData = scanListWithData;
                        MultiScanObj.createMs1NavEements(scanIdList,"popupspectrum");
                        let prec_mz = prsm_data.prsm.ms.ms_header.precursor_mz;
                        document.getElementById("popup_dataLoading").remove(); 
                        generateCorrespondingGraph(current_data,"popupspectrum",prec_mz,specId);
                    }
                    else{
                        //Setting data to a global variable
                        ms2_ScansWithData = scanListWithData;
                        MultiScanObj.createMs2NavEements(scanIdList,"ms2svg");
                        document.getElementById("dataLoading").remove();
                        generateCorrespondingGraph(current_data,"ms2svg",null,specId);
                    }
                    scanbuttons();
                }
            }
        });
    }
    getUniqueScanIdList(MultiScanList){
        let uniqueList = [];
        let uniqueIdSet = new Set();
        MultiScanList.forEach(function(element){
            uniqueIdSet.add(element);
        });
        //spread operator converts array to list
        uniqueList = [...uniqueIdSet];
        return uniqueList ;
    }
    
    createMs2NavEements(scanidList,id){
        let _ul = document.getElementById("ms2_graph_nav");
        scanidList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let a = document.createElement("a");
            a.setAttribute("class","nav-link ms2_scanIds");
            if(i == 0)
            {
                a.setAttribute("class","nav-link ms2_scanIds active");
            }
            a.setAttribute("href","#!");
            a.setAttribute("value",element);
            a.innerHTML = "Scan "+ element;
            li.appendChild(a);
            _ul.appendChild(li);
         })
    }
    createMs1NavEements(scanidList,id){
        let _ul = document.getElementById("ms1_graph_nav");
        scanidList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let a = document.createElement("a");
            a.setAttribute("class","nav-link ms1_scanIds");
            if(i == 0)
            {
                a.setAttribute("class","nav-link ms1_scanIds active");
            }
            a.setAttribute("href","#!");
            a.setAttribute("value",element);
            a.innerHTML = "Scan "+ element;
            li.appendChild(a);
            _ul.appendChild(li);
         })
    }
} 
function generateCorrespondingGraph(current_data,id,prec_mz,specId){
    let peak_data = new PeakData();
    let peak_list = peak_data.getPeakData(current_data);
    let envelope_list = peak_data.getEnvelopeData(current_data);
    let ionData = peak_data.getIonData(prsm_data,specId);
    let graphParams = graphOptions();
    if(id == "ms2svg")
    {
        ms2_graph = addSpectrum(id, peak_list, envelope_list, prec_mz, ionData, graphParams);
    }
    else{
        spectrumgraph = addSpectrum(id, peak_list, envelope_list, prec_mz, ionData, graphParams);
    }
}
function getCurrentData(dataList,scanid){
    let current_data;
    let len = dataList.length;
    let specId;
    for(let i=0;i<len;i++)
    {
        if(dataList[i].scanid == parseInt(scanid))
        {
            current_data = dataList[i].value;
            specId = dataList[i].specId;
            break;
        }
    }
    return [current_data,specId];
}
function activateCurrentnavbar(id,currentValue){
    let childs = $("#"+id).children();
    let len = childs.length;
    for(let i=0;i<len;i++)
    {
        let grandchildValue = $(childs[i]).children().attr('value');
        if(grandchildValue == currentValue)
        {
            $("#"+id+" .active").removeClass("active");
            $(childs[i]).children().addClass("active");
            break;
        }
    }
}
function scanbuttons(){
	//ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
        let value = this.getAttribute('value');
        let [currentData,specId] = getCurrentData(ms2_ScansWithData,value);
		generateCorrespondingGraph(currentData,"ms2svg",null,specId);
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
}