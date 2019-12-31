/**
 * On Click implimentation of topview button
 */
function onclickTopView(){
    //let url = "../Spectrum-Visualization/spectrum.html"+"?"+"prsmid=87/"+"specId=";
    let peakAndIntensityList = getDataFromPRSMtoSpectralView();
    let massAndIntensityList = getMassAndIntensityData();
    let sequence = getSequence();
    let fixedPtmList = getFixedPTMMassList();

    window.localStorage.setItem('peakAndIntensityList',  JSON.stringify(peakAndIntensityList));
    window.localStorage.setItem('massAndIntensityList',  JSON.stringify(massAndIntensityList));
    window.localStorage.setItem('sequence',  JSON.stringify(sequence));
    window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));

    window.open("../inspect/spectrum.html");
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
    prsm_data.prsm.ms.peaks.peak.forEach(function(eachPeak,i){
        let tempObj = eachPeak.monoisotopic_mass + " "+eachPeak.intensity+ " "+eachPeak.charge;
        massAndIntensityList.push(tempObj);
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
