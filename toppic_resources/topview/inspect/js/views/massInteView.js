/**
 * Set Mass List on to html
 */
function setDataToMassAndIntensity(massAndIntensityList){
    jqueryElements.massData.val(massAndIntensityList);
}
/**
 * getlist of Mass,Intensity and charge from UI
 */
function getMassListFromUI(){
    let spectrumDataList = [];
    // Read data line by line from the mass and intensity box
    var lines = $('#data').val().split('\n');
    for(var i = 0; i < lines.length;i++){
        let massAndInte = lines[i].trim() ;
        if(massAndInte.length !=  0 )
        {
            let spectrumData = {};
            // Get Mass,intensity and charge either by space seperated or tab seperated
            let massInte = massAndInte.split(/[\s]+/);
            if(massInte[0] != undefined && massInte[1] != undefined 
                    && !isNaN(massInte[0]) && !isNaN(massInte[1]))
            {
                spectrumData.mass = parseFloat(massInte[0]);
                spectrumData.intensity = parseFloat(massInte[1]);
                spectrumData.charge = parseFloat(massInte[2]);
            }
            if(!jQuery.isEmptyObject(spectrumData))
            {
                if((spectrumData.mass !== undefined) &&
                    (spectrumData.intensity !== undefined))
                {
                    spectrumDataList.push(spectrumData) ;
                }
            }
        }
    }
    completeCalData.monomasslist = spectrumDataList;
    return spectrumDataList ;
}