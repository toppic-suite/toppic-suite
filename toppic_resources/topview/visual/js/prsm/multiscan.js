
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
                let temp_ms2DataObj = {"key":temp_scanid,"value":temp_value};
                scanListWithData.push(temp_ms2DataObj);
                count++;
                if(count == len)
                {
                    let MultiScanObj = new MultiScan();
                    let current_data = getCurrentData(scanListWithData,scanIdList[0]);
                    
                    if(ms1_ms2_json == "ms1_json")
                    {
                        ms1_ScansWithData = scanListWithData;
                        MultiScanObj.createMs1NavEements(scanIdList,"popupspectrum");
                        console.log("scanListWithData[0] : ", scanListWithData[0]);
                        let prec_mz = prsm_data.prsm.ms.ms_header.precursor_mz;
                        document.getElementById("popup_dataLoading").remove(); 
                        generateCorrespondingGraph(current_data,"popupspectrum",prec_mz);
                    }
                    else{
                        //Setting data to a global variable
                        ms2_ScansWithData = scanListWithData;
                        MultiScanObj.createMs2NavEements(scanIdList,"ms2svg");
                        console.log("scanListWithData[0] : ", scanListWithData[0]);
                        document.getElementById("dataLoading").remove();
                        generateCorrespondingGraph(current_data,"ms2svg",null);
                    }
                    scanbuttons(scanListWithData);
                }
            }
        });
    }
    getUniqueScanIdList(MultiScanList){
        console.log("MultiScanList 1", MultiScanList );
        let uniqueList = [];
        let uniqueIdSet = new Set();
        MultiScanList.forEach(function(element){
            console.log("element : ", element);
            uniqueIdSet.add(element);
        });
        //spread operator converts array to list
        uniqueList = [...uniqueIdSet];
        return uniqueList ;
    }
    
    createMs2NavEements(scanidList,id){
        let _ul = document.getElementById("ms2_graph_nav");
        console.log("scanidList : ", scanidList);
        scanidList.forEach(function(element,i){
            console.log("i", i);
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
        console.log("self",self);
        let _ul = document.getElementById("ms1_graph_nav");
        console.log("ul : ", _ul);
        //console.log("nav scanidList : ", scanidList);
        scanidList.forEach(function(element,i){
            console.log();
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
function generateCorrespondingGraph(current_data,id,prec_mz){
    console.log("current_data : ", current_data);
    let peak_data = new PeakData();
    let peak_list = peak_data.getPeakData(current_data);
    let envelope_list = peak_data.getEnvelopeData(current_data);
    if(id == "ms2svg")
    {
        ms2_graph = addSpectrum(id, peak_list, envelope_list, prec_mz);
    }
    else{
        spectrumgraph = addSpectrum(id, peak_list, envelope_list, prec_mz);
    }
}
function getCurrentData(dataList,scanid){
    console.log("dataList : ",dataList, "scanid : ", scanid);
    let current_data;
    let len = dataList.length;
    for(let i=0;i<len;i++)
    {
        if(dataList[i].key == parseInt(scanid))
        {
            current_data = dataList[i].value;
            break;
        }
    }
    return current_data;
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
function scanbuttons(scanListWithData){
	//ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms2_ScansWithData.length;
		let currentData = getCurrentData(scanListWithData,value);
		generateCorrespondingGraph(currentData,"ms2svg",null);
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	$(".ms1_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms1_ScansWithData.length;
		let currentData = getCurrentData(scanListWithData,value);
		generateCorrespondingGraph(currentData,"popupspectrum",null);
		$("#ms1_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
}