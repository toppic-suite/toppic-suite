/**
 * Function gets invoked when clicked on inspect button.
 * This function stores all the information of a spectrum that is being inspected using local storage.
 * @param {object} scansWithData - This is a json object with data of all the scan Ids
 * @param {Integer} scanId - Contians the Scan numbers
 */
/*function onclickTopView(scansWithData,scanId){
    let currentSpectrumData;
    [currentSpectrumData,specId] = getCurrentData(scansWithData,scanId);
    let peakAndIntensityList = getDataFromPRSMtoSpectralView(currentSpectrumData);
    let massAndIntensityList = getMassAndIntensityData(specId);
    let sequence = getSequence();
    let fixedPtmList = getFixedPTMMassList();
    let unknownMassShiftList = getUnknownMassList();
    let precursorMass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
    // Stores all the data in the variables respectively
    window.localStorage.setItem('peakAndIntensityList',  JSON.stringify(peakAndIntensityList));
    window.localStorage.setItem('massAndIntensityList',  JSON.stringify(massAndIntensityList));
    window.localStorage.setItem('sequence',  JSON.stringify(sequence));
    window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));
    window.localStorage.setItem('unknownMassShiftList', JSON.stringify(unknownMassShiftList));
    window.localStorage.setItem('precursorMass', JSON.stringify(precursorMass));
    window.open("../inspect/spectrum.html");
}*/
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
    let ionArray = [];

	script.src = fileName;
	body.append(script);
	script.onload = function(){
        let peakAndIntensityList = getDataFromPRSMtoSpectralView(ms2_data);
        [massAndIntensityList, ionArray] = getMassAndIntensityData(specID);
        let sequence = getSequence();
        let fixedPtmList = getFixedPTMMassList();
        let variablePtmList = getVariablePTMMassList();
        let unknownMassShiftList = getUnknownMassList();
        let precursorMass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
        // Stores all the data in the variables respectively
        window.localStorage.setItem('peakAndIntensityList',  JSON.stringify(peakAndIntensityList));
        window.localStorage.setItem('massAndIntensityList',  JSON.stringify(massAndIntensityList));
        window.localStorage.setItem('ionType', ionArray);
        window.localStorage.setItem('sequence',  JSON.stringify(sequence));
        window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));
        window.localStorage.setItem('variablePtmList', JSON.stringify(variablePtmList));
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
  //to save time, retrieve ion array from along with peaks and intensity data, instead of iterating the same data again
  let massAndIntensityList = [];
  let ionArray = [];
 
  prsm_data.prsm.ms.peaks.peak.forEach(function(eachPeak,i){
    if (eachPeak.spec_id == specId) {
      let tempObj = eachPeak.monoisotopic_mass + " "+eachPeak.intensity+ " "+eachPeak.charge;
      massAndIntensityList.push(tempObj);
    }
    if (parseInt(eachPeak.matched_ions_num)>0){
        let ion = eachPeak.matched_ions.matched_ion.ion_type;
        if (ionArray.length < 1){
            ionArray.push(ion);
        }			
        else if (ionArray.indexOf(ion) < 0) {
            ionArray.push(ion);
        }
    }
  })
  return [massAndIntensityList, ionArray];
}
/**
 * Function to get the sequence of the protein from prsm.js
 */
function getSequence(){
    let sequence = [];
    let firstposition = prsm_data.prsm.annotated_protein.annotation.first_residue_position;
    let lastposition = prsm_data.prsm.annotated_protein.annotation.last_residue_position;
    prsm_data.prsm.annotated_protein.annotation.residue.forEach(function(eachrow,i){
        if(parseInt(eachrow.position) >= parseInt(firstposition)&&
            parseInt(eachrow.position) <= parseInt(lastposition))
        {
            sequence = sequence+eachrow.acid;
        }
    })
    return sequence;
}
/**
 * Gets all the masslist shifts of the protein variable ptms for prsm
 */
function getVariablePTMMassList()
{
    let variablePTMList = [];
    let l_prsm = prsm_data;
    let firstPos = parseInt(l_prsm.prsm.annotated_protein.annotation.first_residue_position);

    if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('ptm'))
    {
        let ptm = l_prsm.prsm.annotated_protein.annotation.ptm ;
		if(Array.isArray(ptm))
		{
			ptm.forEach(function(ptm, index){
				if(ptm.ptm_type == "Protein variable")
				{
                    let abbrevation = ptm.ptm.abbreviation ;
                    let position = parseInt(ptm.occurence.left_pos) - firstPos;
                    let tempObj = {"name":abbrevation, "position":position, "mono_mass":ptm.ptm.mono_mass};
                    variablePTMList.push(tempObj);
				}
			})
        }
        else
        {
            if(ptm.ptm_type == "Protein variable")
            {
                let abbrevation = ptm.ptm.abbreviation ;
                let position = parseInt(ptm.occurence.left_pos) - firstPos;
                let tempObj = {"name":abbrevation, "position":position, "mono_mass":ptm.ptm.mono_mass};
                variablePTMList.push(tempObj);
            }
        }
    }
    return variablePTMList;
}
/**
 * Gets all the masslist shifts of the Fixed ptms for prsm
 */
function getFixedPTMMassList()
{
    let fixedPTMList = [];
    let l_prsm = prsm_data;
    if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('ptm'))
    {
        let ptm = l_prsm.prsm.annotated_protein.annotation.ptm ;
		if(Array.isArray(ptm))
		{
			ptm.forEach(function(ptm, index){
				if(ptm.ptm_type == "Fixed")
				{
                    let abbrevation = ptm.ptm.abbreviation ;
                    let tempObj = {name:abbrevation};
                    fixedPTMList.push(tempObj);
				}
			})
        }
        else
        {
            if(ptm.ptm_type == "Fixed")
            {
                let abbrevation = ptm.ptm.abbreviation ;
                let tempObj = {name:abbrevation};
                fixedPTMList.push(tempObj);
            }
        }
    }
    return fixedPTMList;
}
/**
 * Get all the unknwon mass lists form the prsm
 */
function getUnknownMassList()
{
    let unknownMassShiftList = [];
    let l_prsm = prsm_data;
    if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
	{
        let mass_shift = l_prsm.prsm.annotated_protein.annotation.mass_shift ;
        //adjust position based on firstPos and lastPos of sequence
        //console.log(l_prsm)
        let firstPos = parseInt(l_prsm.prsm.annotated_protein.annotation.first_residue_position);
        //let lastPos = l_prsm.prsm.annotated_protein.annotation.last_residue_position;

        if(Array.isArray(mass_shift))
        {
            let len = mass_shift.length;
            mass_shift.forEach(function(each_mass_shift, i){
                // Removing -1 as the sequece in inspect elements takes from 0
                let position = parseInt(each_mass_shift.left_position) - firstPos;
                let mass = parseFloat(each_mass_shift.anno);
                unknownMassShiftList.push({"position":position,"mass":mass})
            })
        }
        else if(mass_shift.shift_type == "unexpected")
        {
                // Removing -1 as the sequece in inspect elements takes from 0
            let position = parseInt(mass_shift.left_position) - firstPos;
            let mass = parseFloat(mass_shift.anno);
            unknownMassShiftList.push({"position":position,"mass":mass})
        }
    }
    return unknownMassShiftList;
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
