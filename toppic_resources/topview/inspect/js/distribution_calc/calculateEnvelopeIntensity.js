/**
 * Code to normalize the Intensity. 
 * Take the average of intensity from the peaks entered by the user.
 * Take the average of the calculated distribution for each Array element in the Array. 
 * Make both of them equal and calculating the rest of the 
 * distribution intensity based on the avg value from the peak list.
 * @param {Array} totDistributionList - Total distribution calculated
 * @param {Array} peakDataList - Peak Data entered by the user
 */
function getNormalizedIntensity(totDistributionList,peakDataList)
{
    let len = totDistributionList.length;
    let peakListLen = peakDataList.length;
    let intensity = 0;
    let count = 0 ;
    let distributionInte = 0;
    for(let i=0;i<len;i++)
    {
        for(let j=0;j<peakListLen;j++)
        {
            if(Math.abs(totDistributionList[i].mz - peakDataList[j].mz) <= this.toleraceMassDiff )
            {
                intensity = intensity + peakDataList[j].intensity ;
                distributionInte = distributionInte + totDistributionList[i].intensity;
                count = count + 1;
            }
        }
    }
    let avg = intensity/count ;
    let distributionAvgInte = distributionInte/count;
    for(let i=0;i<len;i++)
    {
        totDistributionList[i].intensity = (avg * totDistributionList[i].intensity)/distributionAvgInte ;
    }
    return totDistributionList ;
}