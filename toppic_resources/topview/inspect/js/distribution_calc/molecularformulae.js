/**
 * Code to calculate distribution based on molecular formula
 */
class MolecularFormulae{
	constructor(){
		this.averagine = "C4.9384H7.7583N1.3577O1.4773S0.0417" ;
		this.avrgMass = 111.0543055992;
		this.atomList = ["C","H","N","O","S"]
		this.carbonCnt = 4.9384;
		this.hydrogenCnt = 7.7583;
		this.nitrogenCnt = 1.3577;
		this.oxygenCnt = 1.4773
		this.sulphurCnt = 0.0417;
		this.minintensity = 0.0001;
		this.protonMass = 1.00727647  ;
		this.eachAtomCount;
		this.toleraceMassDiff = 0.01;
		this.mininte = 1;
		this.intensityTolerance = 1 ;
	}
	/**
	 * Function to calculate the emass distrubution fo the given sequence
	 * @param {Double} mass - Contains the mass at that position of the sequence
	 * @param {Array} peakDataList - Contains the peak list provided by the user
	 * @param {Float} charge - Contains the chrage of the ion
	 */
	emass(mass,charge,peakDataList){
		//Give the count of each atom in molecule
		this.eachAtomCount = mass/this.avrgMass ;
		let calculatedMass;
		let atomCountList ;
		[atomCountList,calculatedMass] = this.getMolecularFormula() ;
		let massError = mass - calculatedMass ;
		let len = atomCountList.length;
		let totDistributionList = [] ;
		let numOfAtoms = this.atomList.length ;
		for(let i=0;i<len;i++)
		{
			for(let j=0;j<atomCountList[i].count;j++)
			{
				let atomdist = getIsotopicMassRef(atomCountList[i].atom);
				totDistributionList = this.getMassAndIntensity(totDistributionList,atomdist);
			}
		}
		for(let k=0;k<totDistributionList.length;k++)
		{
			for(let i = 0; i < numOfAtoms ; i++)
			{
				let IsotopicMass = getIsotopicMassOfAtom(atomCountList[i].atom);
				totDistributionList[k].mass = totDistributionList[k].mass + (IsotopicMass[0].mass * atomCountList[i].count) ;
			}
		}
		console.log("totDistributionList 1 :", totDistributionList);
		totDistributionList = this.getMZwithHighInte(totDistributionList,charge,massError,peakDataList);
		totDistributionList = this.getNormalizedIntensity(totDistributionList,peakDataList);
		return totDistributionList ;
	}
	/**
	 * Logic to calculate distribution 
	 * @param {Array} totDistributionList - Array with current total distribution
	 * @param {Array} aminoAcidDist - Array with existing calculated distribution of amino acid 
	 */
	getMassAndIntensity(totDistributionList,atomdist){
		let maxintensity = 0 ;
		if(totDistributionList.length == 0)
		{
			return atomdist ;
		}
		else
		{
			let len = totDistributionList.length + atomdist.length - 1;
	  		let completeDistributionList = new Array(len).fill(0) ;
	  		for(let i=0;i<totDistributionList.length;i++)
			{
				for(let j=0;j<atomdist.length;j++)
				{
					let intensity = 0;
					let index = i+j ;
					let mass = totDistributionList[i].mass+atomdist[j].mass ;
					intensity = totDistributionList[i].intensity * atomdist[j].intensity ;
					if(completeDistributionList[index] != 0) intensity = intensity + completeDistributionList[index].intensity ;
					completeDistributionList[index] = {mass:mass,intensity:intensity};
					if(intensity > maxintensity) maxintensity = intensity ;
				}
			}
	  		let completeDistributionList_temp = [];
			for(let i = 0; i<completeDistributionList.length; i++)
			{
				let tempintensityVal = (completeDistributionList[i].intensity/maxintensity)*100;
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
	 * Code to get the molecular formular based on the mass.
	 * Number of atoms of each atom in the averagine * number of
	 * atoms of each atom calculated from mass gives the total number of atoms.
	 */
	getMolecularFormula(){
		let atomListwithCnt = new Array(5) ;
		let len = this.atomList.length;
		let mass = 0;
		for(let i=0;i<len;i++)
		{
			let count ;
			if(this.atomList[i] == "C"){
				count = parseInt(this.eachAtomCount * this.carbonCnt) ;
				mass = mass + count*getIsotopicMassOfAtom(this.atomList[i])[0].mass;
			}
			else if(this.atomList[i] == "H"){
				count = parseInt(this.eachAtomCount * this.hydrogenCnt) ;
				mass = mass + count*getIsotopicMassOfAtom(this.atomList[i])[0].mass;
			}
			else if(this.atomList[i] == "N"){
				count = parseInt(this.eachAtomCount * this.nitrogenCnt) ;
				mass = mass + count*getIsotopicMassOfAtom(this.atomList[i])[0].mass;
			}
			else if(this.atomList[i] == "O"){
				count = parseInt(this.eachAtomCount * this.oxygenCnt) ;
				mass = mass + count*getIsotopicMassOfAtom(this.atomList[i])[0].mass;
			}
			else if(this.atomList[i] == "S"){
				count = parseInt(this.eachAtomCount * this.sulphurCnt) ;
				mass = mass + count*getIsotopicMassOfAtom(this.atomList[i])[0].mass;
			}
			atomListwithCnt[i] = {atom:this.atomList[i],count:count};
		}
		return [atomListwithCnt,mass] ;
	}
	/**
	 * Code to calculate MZ(mass/charge) and remove low intensity values
	 * @param {Array} totDistributionList - Array with total distribution
	 * @param {Float} charge - Float from mass list 
	 * @param {Float} massError - Calculated massError needed to be added 
	 * as we taken int values of number of atoms eleminating mass from the float values. 
	 * @param {Array} peakDataList - Array of peak list from user
	 */
	getMZwithHighInte(totDistributionList,charge,massError,peakDataList)
	{
		let len = totDistributionList.length;
		let newtotDistributionList = [];
		let overaAllMaxIntensity = 0 ;
		let onePerctInte = 0;
		let peakListLen = peakDataList.length;
		for(let k=0;k<peakListLen ; k++)
		{
			if(peakDataList[k].intensity > overaAllMaxIntensity) overaAllMaxIntensity = peakDataList[k].intensity ;
		}
		onePerctInte = overaAllMaxIntensity/100 ;
		for(let i=0;i<len;i++)
		{
			let intensity = totDistributionList[i].intensity ;
			intensity = overaAllMaxIntensity * intensity/100 ;
			if(intensity > onePerctInte)
			{
				let mz = (totDistributionList[i].mass+massError)/charge + this.protonMass;
				let intensity = totDistributionList[i].intensity ;
				//Converting mass variable to mz(mass/charge)
				let tempdistObj = {mz:mz,intensity:intensity};
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
}