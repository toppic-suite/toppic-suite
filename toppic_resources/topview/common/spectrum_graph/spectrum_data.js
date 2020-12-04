/**	@function SpectrumData
 * @description assign level to the input data to the spectrum graphs 
 */
class SpectrumData{
    mzLevel = [{"interval": 25, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 10, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 5, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 3, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 2, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 1.5, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000},
               {"interval": 1, "intervalNum": 0, "maxPeakIdx":-1, "maxInte": -10000}]//m/z interval
    constructor() {}

/**
 * @function assignLevel
 * @description assigns level to each peak and envelope to control the number of peaks/envelopes shown at once
 * @param {object} peakOrEnvList - fragment ion peaks list or envelope list
 */
 assignLevel(peakOrEnvList){
     //iterating once through the peakList or envList, assign level to each peak/envelope
     //initially, all peaks/envs are assigned to the largest level
     //and if a peak/env is the hightest peak in the given interval (ex: betwen 200 - 300 m/z or 150m/z - 200m/z)
     //the peak's level is updated
     
    peakOrEnvList.sort(function(x, y){
        return d3.ascending(parseFloat(x.mz), parseFloat(y.mz));
    })

    let minMz = parseFloat(peakOrEnvList[0].mz);

    for (let i = 0; i < peakOrEnvList.length; i++){
        let mz = parseFloat(peakOrEnvList[i].mz) - minMz;
        let inte = parseFloat(peakOrEnvList[i].intensity);
        peakOrEnvList[i]["level"] = this.mzLevel.length;//as an initial value, set to the +1 of the largest level

        //for each peak, check with each interval level to see if the peak is the maximum intensity in the range
        for (let k = 0; k < this.mzLevel.length; k++){
            let eachInterval = this.mzLevel[k];
            let idx = Math.floor(mz / eachInterval.interval);

            if (idx > eachInterval.intervalNum){
                if (peakOrEnvList[eachInterval.maxPeakIdx]["level"] > k){
                    //don't update the peak level if it is already assigned a higher level (= if it is a max peak in larger range) 
                    peakOrEnvList[eachInterval.maxPeakIdx]["level"] = k;
                }
                eachInterval.maxInte = -10000;
                eachInterval.intervalNum++;
            }
            if (inte > eachInterval.maxInte){
                eachInterval.maxInte = inte;
                eachInterval.maxPeakIdx = i;
            }
        }
    }
  }
}