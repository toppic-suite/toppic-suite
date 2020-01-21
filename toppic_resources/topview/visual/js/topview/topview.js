/**
 * On Click implimentation of topview button
 */
function onclickTopView(){
    window.open("../inspect/spectrum.html");
    //let url = "../Spectrum-Visualization/spectrum.html"+"?"+"prsmid=87/"+"specId=";
    let peakAndIntensityList = getDataFromPRSMtoSpectralView();
    let massAndIntensityList = getMassAndIntensityData();
    let sequence = getSequence();
    let fixedPtmList = getFixedPTMMassList();
    let unknownMassShiftList = getUnknownMassList();
    let precursorMass = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
    console.log("unknownMassShiftList : ", unknownMassShiftList);

    window.localStorage.setItem('peakAndIntensityList',  JSON.stringify(peakAndIntensityList));
    window.localStorage.setItem('massAndIntensityList',  JSON.stringify(massAndIntensityList));
    window.localStorage.setItem('sequence',  JSON.stringify(sequence));
    window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));
    window.localStorage.setItem('unknownMassShiftList', JSON.stringify(unknownMassShiftList));
    window.localStorage.setItem('precursorMass', JSON.stringify(precursorMass));
}

function getDataFromPRSMtoSpectralView(){
    let peakdata = new PeakData();
    let peakAndIntensity = [];
    ms2_data.peaks.forEach(function(eachrow){
        let tempObj = eachrow.mz + " " + eachrow.intensity;
        peakAndIntensity.push(tempObj);
    })
    return peakAndIntensity;
}

function getMassAndIntensityData(){
  let massAndIntensityList = [];
  let ms2_ids = prsm_data.prsm.ms.ms_header.ids;
  let ms2_id_list = ms2_ids.split(" ");
  let ms2_id_1 = ms2_id_list[0];

  prsm_data.prsm.ms.peaks.peak.forEach(function(eachPeak,i){
    if (eachPeak.spec_id == ms2_id_1) {
      let tempObj = eachPeak.monoisotopic_mass + " "+eachPeak.intensity+ " "+eachPeak.charge;
      massAndIntensityList.push(tempObj);
    }
  })
  return massAndIntensityList;
}
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
function getUnknownMassList()
{
    let unknownMassShiftList = [];
    let l_prsm = prsm_data;
    if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
	{
        let mass_shift = l_prsm.prsm.annotated_protein.annotation.mass_shift ;
			if(Array.isArray(mass_shift))
			{
				let len = mass_shift.length;
				mass_shift.forEach(function(each_mass_shift, i){
                    // Removing -1 as the sequece in inspect elements takes from 0
                    let position = parseInt(each_mass_shift.left_position) - 1 ;
                    let mass = parseFloat(each_mass_shift.anno);
					unknownMassShiftList.push({"position":position,"mass":mass})
				})
			}
			else if(mass_shift.shift_type == "unexpected")
			{
                 // Removing -1 as the sequece in inspect elements takes from 0
                let position = parseInt(mass_shift.left_position) - 1;
                let mass = parseFloat(mass_shift.anno);
                unknownMassShiftList.push({"position":position,"mass":mass})
			}
    }
    return unknownMassShiftList;
}
