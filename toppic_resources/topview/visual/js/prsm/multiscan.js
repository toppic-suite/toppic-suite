
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
                let temp_ms2DataObj = {"scanId":temp_scanid,"specId":element,"value":temp_value};
                scanListWithData.push(temp_ms2DataObj);
                count++;
                // This comparison help to call the generate graph function after loading the data completely
                // As data loading is an asynchronous process, we need to wait till all the data is loaded to execute the below functions
                if(count == len)
                {
                    let MultiScanObj = new MultiScan();
                    let current_data ;
                    let specId;
                    if(ms1_ms2_json == "ms1_json")
                    {

                        ms1_ScansWithData = scanListWithData;
                        MultiScanObj.createMs1NavEements(scanIdList,"popupspectrum");
                        let prec_mz = prsm_data.prsm.ms.ms_header.precursor_mz;
                        document.getElementById("popup_dataLoading").remove(); 
                        generateCorrespondingGraph(ms1_ScansWithData[0].value,"popupspectrum",prec_mz,specId);
                    }
                    else{
                        [current_data,specId] = getCurrentData(scanListWithData,scanIdList[0]);
                        //Setting data to a global variable
                        ms2_ScansWithData = scanListWithData;
                        MultiScanObj.createMs2NavEements(scanIdList,"ms2_graph_nav");
                        document.getElementById("dataLoading").remove();
                        createMultipleSvgs("ms2svg_div","ms2svg_","ms2_svg_graph_class",ms2_ScansWithData);
                    }
                    graphOnClickActions();
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
        let _ul = document.getElementById(id);
        scanidList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let li_id = id+"_"+element;
            li.setAttribute("id",li_id);
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
    createMs2PopUpNavEements(scanidList,id){
        let _ul = document.getElementById(id);
        scanidList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let li_id = id+"_"+element;
            li.setAttribute("id",li_id);
            let a = document.createElement("a");
            a.setAttribute("class","nav-link ms2_popup_scanIds");
            let a_id = "ms2_popup_scanIds_"+element;
            a.setAttribute("id",a_id);
            if(i == 0)
            {
                a.setAttribute("class","nav-link ms2_popup_scanIds active");
            }
            a.setAttribute("href","#!");
            a.setAttribute("value",element);
            a.innerHTML = "Scan "+ element;
            li.appendChild(a);
            _ul.appendChild(li);
         })
    }
    createMs1NavEements(element,id){
        let _ul = document.getElementById("ms1_graph_nav");
        let li = document.createElement("li");
        li.setAttribute("class","nav-item");
        let a = document.createElement("a");
        a.setAttribute("class","nav-link ms1_scanIds active");
        a.setAttribute("href","#!");
        a.setAttribute("value",element);
        a.innerHTML = "Scan "+ element;
        li.appendChild(a);
        _ul.appendChild(li);
    }
    createMonoMassNavEements(scanidList,id){
        let _ul = document.getElementById(id);
        scanidList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let li_id = id+"_"+element;
            li.setAttribute("id",li_id);
            let a = document.createElement("a");
            a.setAttribute("class","nav-link monoMass_scanIds");
            if(i == 0)
            {
                a.setAttribute("class","nav-link monoMass_scanIds active");
            }
            a.setAttribute("href","#!");
            a.setAttribute("value",element);
            a.innerHTML = "Scan "+ element;
            li.appendChild(a);
            _ul.appendChild(li);
         })
    }
} 
function getCurrentData(dataList,scanId){
    let current_data;
    let len = dataList.length;
    let specId;
    for(let i=0;i<len;i++)
    {
        if(dataList[i].scanId == parseInt(scanId))
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