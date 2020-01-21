
class MultiScan{

    constructor(){
    }
    //ms1_ms2_json this give ms1_json folder or ms2_json folder
    promiseLoadDataJS(scanidList,ms1_ms2_json){
        let timeout = 700;
        if(scanidList.length == 1 ) timeout = 200;
        return new Promise((resolve,reject) => {
            let scanListWithData = [];
            //get the list of lists with scan Id and value
            scanidList.forEach(function(element,i){
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
                    scanListWithData[i] = temp_ms2DataObj;
                }
            });
            setTimeout(function(){
                this.scanListWithData = scanListWithData ;
                resolve(scanListWithData);
            },timeout)
        })
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
        console.log("uniqueIdSet : "+uniqueList);
        return uniqueList ;
    }
    
    createMs2NavEements(scanidList,id){
        console.log("self",self);
        let _ul = document.getElementById("ms2_graph_nav");
        console.log("ul : ", _ul);
        console.log("nav scanidList : ", scanidList);
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
            a.setAttribute("value",element.key);
            a.innerHTML = "Scan "+ element.key;
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
            a.setAttribute("value",element.key);
            a.innerHTML = "Scan "+ element.key;
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