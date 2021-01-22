/**
 * Set Peak List on to html
 */
function setDataToPeakAndIntensity(peakAndIntensityList){
    jqueryElements.peakData.val(peakAndIntensityList);
}
/**
 * Function to get data of peaks and intensity from UI
 */
function getPeakListFromUI(){
    let spectrumDataList = [];
    // Read data line by line from peak and intensity box
    var lines = jqueryElements.peakData.val().split('\n');
    for(var i = 0; i < lines.length;i++){
        let peakAndInte = lines[i].trim();
        if(peakAndInte.length !== 0)
        {
            let peak = null;
            let peakInte = peakAndInte.split(/[\s]+/);
            if(peakInte[0] != undefined && peakInte[1] != undefined 
                    && !isNaN(peakInte[0]) && !isNaN(peakInte[1]))
            {
                peak = new Peak(i, parseFloat(peakInte[0]), parseFloat(peakInte[1]))
            }
            if(peak)
            {
                if((peak.getMz() !== undefined) &&
                    (peak.getIntensity() !== undefined))
                {
                    spectrumDataList.push(peak) ;
                }
            }
        }
    }	
    // completeCalData.peakdatalist = spectrumDataList;
    return spectrumDataList ;
}