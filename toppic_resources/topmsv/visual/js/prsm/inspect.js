function onclickTopView(e){
  let title = document.getElementById("Protein-Spectrum-Match-Id-SpecId").innerHTML;
  //let specID = title.split("#").pop();
  let specID = e.currentTarget.getAttribute('specid');
  let folderName = "../../topfd";
  folderName = folderName.split("&")[0];
  let script= document.createElement('script');
  let body = document.getElementsByTagName('body')[0];
  let fileName = folderName+"/ms2_json/spectrum"+specID+".js";

  let massAndIntensityList = [];
  let protVarPtmList = [];
  let variablePtmList = [];

  script.src = fileName;
  body.append(script);
  script.onload = function(){
    let peakAndIntensityList = getDataFromPRSMtoSpectralView(ms2_data);
    let ionList = [];
    ionList.push(ms2_data.n_ion_type);
    ionList.push(ms2_data.c_ion_type);
    //console.log(ionList);
    massAndIntensityList = getMassAndIntensityData(specID);
    //console.log(prsmGraph.data.proteoform);
    let sequence = prsmGraph.data.proteoform.sequence; 
    let fixedPtmList = prsmGraph.data.proteoform.getFixedPtmList(); 
    //console.log(fixedPtmList);
    let unknownMassShiftList = prsmGraph.data.proteoform.getUnknownMassList();
    //console.log(unknownMassShiftList);
    let precursorMass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
    // Stores all the data in the variables respectively
    window.localStorage.setItem('peakAndIntensityList',  JSON.stringify(peakAndIntensityList));
    window.localStorage.setItem('massAndIntensityList',  JSON.stringify(massAndIntensityList));
    window.localStorage.setItem('ionType', ionList);
    window.localStorage.setItem('sequence',  JSON.stringify(sequence));
    window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));
    window.localStorage.setItem('protVarPtmsList', JSON.stringify(protVarPtmList));
    window.localStorage.setItem('variablePtmsList', JSON.stringify(variablePtmList));
    window.localStorage.setItem('unknownMassShiftList', JSON.stringify(unknownMassShiftList));
    window.localStorage.setItem('precursorMass', JSON.stringify(precursorMass));
    window.open("../inspect/spectrum.html"); 
  }    
}

/**
 * Get the peaklist from respective spectrum.js to set the data for inspect page
 * @param {object} ms2_data - json object with complete data spectrum for corresponding scan Id
 */
function getDataFromPRSMtoSpectralView(ms2_data){
    let peakAndIntensity = [];
    ms2_data.peaks.forEach(function(eachrow){
        let tempObj = eachrow.mz + " " + eachrow.intensity;
        peakAndIntensity.push(tempObj);
    })
    return peakAndIntensity;
}
/**
 * Get the masslist from respective prsm.js to set the data for inspect page
 * @param {Integer} specId - Contians spec Id to get the data of corrsponding mass list
 */
function getMassAndIntensityData(specId){
  let massAndIntensityList = [];
  prsm_data.prsm.ms.peaks.peak.forEach(function(eachPeak,i){
    if (eachPeak.spec_id == specId) {
      let tempObj = eachPeak.monoisotopic_mass + " "+eachPeak.intensity+ " "+eachPeak.charge;
      massAndIntensityList.push(tempObj);
    }
  })
  return massAndIntensityList; 
}

/**
 * Create HTML dropdown buttons based on the scan list
 * @param {Array} scanIdList - Contains Scan id numbers
 * @param {Array} specIdList - Contains Spec Id numbers
 */
function setDropDownItemsForInspectButton(scanIdList,specIdList){
    let dropdown_menu = $(".dropdownscanlist .dropdown-menu");
    let len = scanIdList.length;
    for(let i=0; i<len;i++)
    {
        let value = scanIdList[i];
        let specId = specIdList[i];
        let id = "scan_"+ value ;
        let a = document.createElement("a");
        a.setAttribute("class","dropdown-item");
        a.setAttribute("href","#!");
        a.setAttribute("id",id);
        a.setAttribute("value",value);
        a.setAttribute("specid", specId);
        a.innerHTML = "Scan "+value;
        dropdown_menu.append(a);
    }
    
}

/**
 * Onclick function, invoked on click of the inspect scn button
 */
function onClickToInspect(){
    $(".dropdownscanlist .dropdown-item ").click(function(e){
        onclickTopView(e);
    });  
}
