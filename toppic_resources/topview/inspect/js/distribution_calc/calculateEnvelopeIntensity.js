/**
 * Code to normalize the Intensity. 
 * Take the average of intensity from the peaks entered by the user.
 * Take the average of the calculated distribution for each Array element in the Array. 
 * Make both of them equal and calculating the rest of the 
 * distribution intensity based on the avg value from the peak list.
 * @param {Array} totDistributionList - Total distribution calculated
 * @param {Array} peakDataList - Peak Data entered by the user
 */
getNormalizedIntensity(totDistributionList,peakDataList)
{
    let len = totDistributionList.length;
    let peakListLen = peakDataList.length;
    let count = 0 ;
    let maxinte = 0;
    let mininte = 100;

    let peakMaxinte = 0;
    let peakMininte = 10000000;

    let maxMz = 0;
    let minMz = 10000000;
    for(let i=0;i<len;i++)
    {
        for(let j=0;j<peakListLen;j++)
        {
            if((totDistributionList[i].mz - peakDataList[j].mz) <= this.toleraceMassDiff
                && (totDistributionList[i].mz - peakDataList[j].mz) >= 0-this.toleraceMassDiff )
            {
                if(maxMz < totDistributionList[i].mz){
                    maxMz = totDistributionList[i].mz;
                }
                if(minMz > totDistributionList[i].mz){
                    minMz = totDistributionList[i].mz;
                }
                count = count + 1;
            }
        }
    }
    maxMz = maxMz + this.toleraceMassDiff;
    minMz = minMz - this.toleraceMassDiff;
    for(let i=0;i<len;i++)
    {
        if(minMz <= totDistributionList[i].mz &&  totDistributionList[i].mz <= maxMz)
        {
            if(maxinte < totDistributionList[i].intensity){
                maxinte = totDistributionList[i].intensity;
            }
            if(mininte > totDistributionList[i].intensity){
                mininte = totDistributionList[i].intensity;
            }
        }
    }
    for(let j=0;j<peakListLen;j++)
    {
        if(peakDataList[j].mz >= minMz && peakDataList[j].mz <= maxMz)
        {
            if(peakMaxinte < peakDataList[j].intensity){
                peakMaxinte = peakDataList[j].intensity;
            }
            if(peakMininte > peakDataList[j].intensity){
                peakMininte = peakDataList[j].intensity;
            }
        }
    }
    if(count !=0 )
    {
        let avg ;
        let distributionAvgInte;

        avg = (peakMaxinte + peakMininte)/2;
        distributionAvgInte = (maxinte+mininte)/2;
        for(let i=0;i<len;i++)
        {
            totDistributionList[i].intensity = (avg * totDistributionList[i].intensity)/distributionAvgInte ;
        }
    }
    return totDistributionList ;
}