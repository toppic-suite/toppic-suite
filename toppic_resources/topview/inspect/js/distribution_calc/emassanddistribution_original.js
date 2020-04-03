/**
 * Class to calculate Emass and Distribution for the given input
 */
class CalculateEmassAndDistribution{

	constructor(){
		this.PREFIX = "PREFIX";
		this.SUFFIX = "SUFFIX";
		this.minintensity = 0.000001;
		this.protonMass = 1.00727647  ;
		this.intensityTolerance = 1 ;
		this.toleraceMassDiff = 0.02;
	}
	/**
	 * Function to calculate the emass distrubution fo the given sequence
	 * @param {String} seq - Contains the sequence provide by the user
	 * @param {Array} peakDataList - Contains the peak list provided by the user
	 * @param {Float} charge - Contains the chrage of the ion
	 * @param {String} pref_suffInd - Indicator to indiace prefix or suffix
	 */
	
	emass(seq,peakDataList,charge,pref_suffInd)
	{
		let AcidArray = seq ;
		let AcidArrayLen= AcidArray.length;
		let totDistributionList = [] ;
		
		for(let i = 0; i < AcidArrayLen ; i++)
		{
			let aminoAcidDist = getAminoAcidDistribution(AcidArray[i]);
			totDistributionList = this.getMassAndIntensity(totDistributionList,aminoAcidDist) ;
		}
		
		if(pref_suffInd == this.SUFFIX)
		{
			let waterDist = getAminoAcidDistribution("H2O");
			totDistributionList = this.getMassAndIntensity(totDistributionList,waterDist);
		}

		totDistributionList = this.getMZwithHighInte(totDistributionList,charge,peakDataList);
		totDistributionList = this.getNormalizedIntensity(totDistributionList,peakDataList);
		totDistributionList.sort(function(x, y){
		   return d3.ascending(x.mz, y.mz);
		})
		return totDistributionList ;
	}
	/**
	 * Logic to calculate distribution 
	 * @param {Array} totDistributionList - Array with current total distribution
	 * @param {Array} aminoAcidDist - Array with existing calculated distribution of amino acid 
	 */
	getMassAndIntensity(totDistributionList,aminoAcidDist)
	{
		let maxintensity = 0 ;
		if(totDistributionList.length == 0)
		{
			return aminoAcidDist ;
		}
		else
		{
			let len = totDistributionList.length + aminoAcidDist.length - 1;
	  		let completeDistributionList = new Array(len).fill(0) ;
			for(let i=0;i<totDistributionList.length;i++)
			{
				for(let j=0;j<aminoAcidDist.length;j++)
				{
					let intensity = 0;
					let index = i+j ;
					let mass = totDistributionList[i].mass+aminoAcidDist[j].mass ;
					intensity = totDistributionList[i].intensity * aminoAcidDist[j].intensity ;
					if(completeDistributionList[index] != 0) intensity = intensity + completeDistributionList[index].intensity ;
					completeDistributionList[index] = {mass:mass,intensity:intensity};
					if(intensity > maxintensity) maxintensity = intensity ;
			
				}
			}
			let completeDistributionList_temp = [];
			for(let i = 0; i<completeDistributionList.length; i++)
			{
				let tempintensityVal = (completeDistributionList[i].intensity/maxintensity)*100;
				//Taking float values till 6 decimals
				if( tempintensityVal > this.minintensity)
				{
					completeDistributionList[i].intensity = Math.round(tempintensityVal * 1000000) / 1000000;//tempintensityVal//Math.round(tempintensityVal * 1000000) / 1000000; //parseFloat(((tempDistributionList[i].intensity/maxintensity)*100).toFixed(6)) ;
					completeDistributionList_temp.push(completeDistributionList[i]);
				}
			}
			completeDistributionList = completeDistributionList_temp ;
			return completeDistributionList ;
		}
	}
	/**
	 * Code to remove the calculate MZ(mass/charge) value and remove low intensities
	 * @param {Array} totDistributionList - Total distribution calculated
	 * @param {Integer} charge - charge for the mass list by user
	 * @param {Array} peakDataList - peaklist entered by the user
	 */
	getMZwithHighInte(totDistributionList,charge,peakDataList)
	{
		let len = totDistributionList.length;
		let overaAllMaxIntensity = 0 ;
		let onePerctInte = 0;
		let peakListLen = peakDataList.length;
		for(let k=0;k<peakListLen ; k++)
		{
			if(peakDataList[k].intensity > overaAllMaxIntensity) overaAllMaxIntensity = peakDataList[k].intensity ;
		}
		onePerctInte = overaAllMaxIntensity/100 ;
		let newtotDistributionList = [];
		for(let i=0;i<len;i++)
		{
			let intensity = totDistributionList[i].intensity ;
			intensity = overaAllMaxIntensity * intensity/100 ;
			if(intensity > onePerctInte)
			{
				let mz = totDistributionList[i].mass/charge + this.protonMass;
				let intensity = totDistributionList[i].intensity ;
				let tempdistObj = {mz:mz,intensity:intensity};

				//if (mz >= 702.2700481840848 && mz <= 703.4556072446619) {
					//console.log("intensity ", intensity);
				//}


				newtotDistributionList.push(tempdistObj);
			}
		}
		return newtotDistributionList ;
	}
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
}