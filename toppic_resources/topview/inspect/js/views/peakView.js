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
    var lines = $('#peakdata').val().split('\n');
    for(var i = 0; i < lines.length;i++){
        let peakAndInte = lines[i].trim();
        if(peakAndInte.length !=  0 )
        {
            let spectrumData = {};
            let peakInte = peakAndInte.split(/[\s]+/);
            if(peakInte[0] != undefined && peakInte[1] != undefined 
                    && !isNaN(peakInte[0]) && !isNaN(peakInte[1]))
            {
                spectrumData.mz = parseFloat(peakInte[0])///spectrumData.charge ;
                spectrumData.intensity = parseFloat(peakInte[1]);
            }
            if(!jQuery.isEmptyObject(spectrumData))
            {
                if((spectrumData.mz !== undefined) &&
                    (spectrumData.intensity !== undefined))
                {
                    spectrumDataList.push(spectrumData) ;
                }
            }
        }
    }	
    completeCalData.peakdatalist = spectrumDataList;
    return spectrumDataList ;
}