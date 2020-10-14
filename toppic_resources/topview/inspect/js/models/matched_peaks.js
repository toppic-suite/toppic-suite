/**
 * Class to calculate Matched Peaks and distribution of the Amoino Acids
 */
class MatchedPeaks {
	PREFIX = "PREFIX";
	SUFFIX = "SUFFIX";
	matchedUnMatchedList = [];
	matchedList = [];
	water = "H2O";
	
	CONST_A = "a";
	CONST_B = "b";
	CONST_C = "c";
	CONST_Y = "Y";
	colors = ["red","orange","blue","grey"];
	massShiftList = [] ;
	constructor(){
		this.matchedList = [];
		this.matchedUnMatchedList = [];
		this.million = 1000000;
	}
	/**
	 * Function to form object with necessary attributes at each amino acid position
	 * @param {Object} monoMassList_temp - Single element from mono Mass List
	 * @param {String} id - Unique Id for each peak
	 * @param {String} ion - Contains matched type(b/y) appended with charge
	 * @param {Integer} position - Position of the amino acid
	 * @param {Float} massDiff - Mass difference at particular position
	 * @param {Float} mass - Mass at the position of the acid 
	 * @param {Char} matchedInd - Indicator indicating matched/not matched. If matched containd "Y"
	 */
	matchedPeakAttributes(monoMassList_temp,id,ion,position,massDiff,mass,matchedInd)
	{
		let peak = monoMassList_temp;
		let thMass = "";
		let PPMerror = "";
		// Check if Matched Indicator is Yes "Y" else keep empty spaces in the new attributes
		if(matchedInd == this.CONST_Y)
		{
			thMass = Math.round(mass * 10000)/10000 ;
			massDiff = Math.round(massDiff * 10000)/10000;
			PPMerror = massDiff/thMass * this.million ;
			PPMerror = Math.round(PPMerror * 10000)/10000;
		}
		peak.peakId = id ;
		peak.ion = ion;
		peak.position = position;
		peak.massError = massDiff;
		peak.thMass = thMass ;
		peak.PPMerror = PPMerror;
		peak.matchedInd = matchedInd;
		// console.log("peak : ", peak);
		return peak;
	}
	/**
	 * Returns Complete list with matched Indication
	 * @param {Array} prefixOrSuffixMassList - List with either prefix or suffix mass list
	 * @param {Array} monoMassList - Contains Mono Mass list from html page
	 * @param {String} sequence - Protein sequence
	 * @param {Float} massErrorthVal - Threshhold value entered by user to determine the mass error and matched peaks
	 * @param {Float} ppmErrorthVal - Threshhold value entered by user to determine the mass error and matched peaks in ppm units
	 * @param {String} ionType - Contains either y/b
	 * @param {String} preOrSufInd - Contains information of either PREFIX/SUFFIX
	 */
	getMatchedPeakList(prefixOrSuffixMassList,monoMassList,sequence,massErrorthVal,ppmErrorthVal,ionType,preOrSufInd)
	{
		let matchedList = [];
		let monoMassList_temp = monoMassList.slice();
		let len = monoMassList_temp.length;
		let preOrSufListln = prefixOrSuffixMassList.length;
		let seqln = sequence.length;
		for(let i=0; i < len; i++)
		{
			let peakId = i+1;
			let peak = {};
			for(let j = 0; j<preOrSufListln; j++)
			{
				let massDiff = monoMassList_temp[i].mass - prefixOrSuffixMassList[j];
				if(massErrorthVal != null || isNaN(massErrorthVal))
				{
					if(Math.abs(massDiff) <= massErrorthVal)
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						let ion = ionType + j;
						let position = j;
						if(preOrSufInd == "suffix")
						{
							position = seqln - position;
						}
						let mass = prefixOrSuffixMassList[j];
						let matchedInd = "Y";
						peak = this.matchedPeakAttributes(monoMassList_temp[i],peakId,ion,position,
																massDiff,mass,matchedInd);
						matchedList.push(peak);		
					}
				}
				else{
					let prefDiff = Math.round(massDiff * 10000)/10000;
					let prefMassRounded = Math.round(prefixOrSuffixMassList[j] * 10000)/10000;
					let prePPMerror = prefDiff/prefMassRounded * this.million ;
					if(Math.abs(prePPMerror) <= ppmErrorthVal)
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						let ion = ionType + j;
						let position = j;
						if(preOrSufInd == "suffix")
						{
							position = seqln - position;
						}
						let mass = prefixOrSuffixMassList[j];
						let matchedInd = "Y";
						peak = this.matchedPeakAttributes(monoMassList_temp[i],peakId,ion,position,
																massDiff,mass,matchedInd);
						matchedList.push(peak);
					}
				}
			}
		}
		return matchedList;
	}
	/**
	 * Set matched indicator to complete data and Seperate matched and unmatched list of data with matchedInd attribute
	 * @param {Array} monoMassList - Contains mono Mass list data
	 * @param {Array} matchedList - Contains all the data of both matched and unmatched data
	 */
	getMatchedAndUnMatchedList(monoMassList,matchedList)
	{
		let matchedAndUnmatchedList = matchedList.map(x => ({...x}));
		//matchedAndUnmatchedList.concat(matchedList);
		let len = matchedList.length;
		// let MatchedPeaksObj = new MatchedPeaks();
		let self = this;
		monoMassList.forEach(function(eachElement,i){
			let matched = false;
			let peakId = i+1;
			for(let j=0; j < len; j++)
			{
				if(peakId == parseInt(matchedList[j].peakId))
				{
					matched = true;
					break;
				}
			}
			if(!matched)
			{
				let matchedInd = "N";
				let peak = self.matchedPeakAttributes(monoMassList[i],peakId,"","","","",matchedInd);
			  matchedAndUnmatchedList.push(peak);
			}
		})
		// completeCalData.matchedandunmatcheddata = matchedAndUnmatchedList;
		return matchedAndUnmatchedList ;
	}

	/**
	 * Get the envelopes/distribution for the sequence
	 * @param {Array} peakDataList - Contains peak list information
	 * @param {String} sequence - Contains protein sequence
	 * @param {Array} matchedUnMatchedList - Contains all matched and unmatched list 
	 */
	getDistribution(peakDataList,sequence,matchedUnMatchedList, completeMassShiftList)
	{
		let len = matchedUnMatchedList.length;
		let totalDistribution = [] ;
		let calEmassAndDistObj = new CalculateEmassAndDistribution();
		let molecularFormObj = new MolecularFormulae();
		let seqln = sequence.length;
		
		//matchedUnmatchedlist sort by intensity
		matchedUnMatchedList.sort(function(x,y){
			return d3.descending(x.intensity, y.intensity);
		})
		for(let i = 0; i < len; i++)
		{
			let distributionList = {};

			distributionList.mono_mass = matchedUnMatchedList[i].mass;
			distributionList.charge = matchedUnMatchedList[i].charge;
			distributionList.env_peaks = molecularFormObj.emass(distributionList.mono_mass,distributionList.charge,peakDataList);
			totalDistribution.push(distributionList);
      
			/*if(matchedUnMatchedList[i].matchedInd == "Y")
			{
				if(this.CONST_A == matchedUnMatchedList[i].ion[0] || this.CONST_B == matchedUnMatchedList[i].ion[0] 
																		|| this.CONST_C == matchedUnMatchedList[i].ion[0])
				{
					let matchedPos = matchedUnMatchedList[i].position ;
					let seq = sequence.slice(0,matchedPos) ;
					/*compare with completeMassShiftList to see if this seq includes mass shift
					* compare matchedPos with positons in completeMassShiftList.
					* if matchedPos is bigger, there is a mass shift inside seq. 
					* Send the position to emass so that the mass shift is reflected in the toDistribution list after the acid*/
					/*let massShiftList = [];
					for (let i = 0; i < completeMassShiftList.length; i++){
						if (matchedPos >= completeMassShiftList[i].position){
							massShiftList.push(completeMassShiftList[i]);
						}
					}
					distributionList.mono_mass = matchedUnMatchedList[i].mass;
					distributionList.charge = matchedUnMatchedList[i].charge;
					distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.PREFIX, massShiftList);
					totalDistribution.push(distributionList);

				}
				else
				{
					let matchedPos = matchedUnMatchedList[i].position ;
					let seq = sequence.slice(matchedPos,seqln) ;

					let massShiftList = [];
					for (let i = 0; i < completeMassShiftList.length; i++){
						if (matchedPos <= completeMassShiftList[i].position){
							//as the seq is a slice of original sequence, mass shift position should be adjusted
							let massData = {};
							massData["position"] = completeMassShiftList[i].position - matchedPos;
							massData["mass"] = completeMassShiftList[i].mass;
							massShiftList.push(massData);
						}
					}
					distributionList.mono_mass = matchedUnMatchedList[i].mass;
					distributionList.charge = matchedUnMatchedList[i].charge;
					distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.SUFFIX, massShiftList);
					totalDistribution.push(distributionList);
				}
			}
			else
			{
				distributionList.mono_mass = matchedUnMatchedList[i].mass;
				distributionList.charge = matchedUnMatchedList[i].charge;
				distributionList.env_peaks = molecularFormObj.emass(distributionList.mono_mass,distributionList.charge,peakDataList);
				totalDistribution.push(distributionList);
			}*/
		}
		if(totalDistribution.length !== 0)
		{	 
			totalDistribution.sort(function(x,y){
				return d3.ascending(x.env_peaks[0].mz, y.env_peaks[0].mz);
			})
		}
		// let envlength = totalDistribution.length;
		// let colorListsize = this.colors.length;
		
		// while(envlength--){
		// 	let index = envlength%colorListsize ;
		// 	totalDistribution[envlength].color = this.colors[index] ;
			
		// }

		return totalDistribution ;
	}
	/**
	 * Get All the matched positions and mass list seperated with a matched indicator
	 * @param {Array} prefixOrSuffixMassList - List of prefix and suffix masses 
	 * @param {Array} monoMassList - Contains Mono Mass list
	 * @param {Float} massErrorthVal - Mass Error Threshold value entered by user
	 * @param {Float} ppmErrorthVal - PPM Error Threshold value entered by user
	 * @param {String} prefixInd - contains if the complete data is of prefix/suffix mass list
	 */
	getMatchedAndUnmatchedPrefixAndSuffixMassList(prefixOrSuffixMassList, monoMassList,
																massErrorthVal,ppmErrorthVal,prefixInd)
	{
		// console.log("monoMassList : ", monoMassList);
		let MatchedAndUnMatchedList = [];
		let monoMassList_temp = monoMassList.slice();
		let len = monoMassList_temp.length;
		let preOrSufListln = prefixOrSuffixMassList.length;
		let MatchedAndUnMatchedListObj = {};
		// console.log("monoMassList_temp:", monoMassList_temp);
		for(let j = 0; j<preOrSufListln; j++)
		{
			let position = j + 1;
			let charge = null;
			if(prefixInd !== "prefix")
			{
				position = preOrSufListln-position+1;
			}
			let mass = prefixOrSuffixMassList[j];
			let matchedInd = "N";
			for(let i=0; i < len; i++)
			{
				let massDiff = monoMassList_temp[i].mass - prefixOrSuffixMassList[j] ;
				let prefDiff = Math.round(massDiff * 10000)/10000;
				let prefMassRounded = Math.round(prefixOrSuffixMassList[j] * 10000)/10000;
				let prePPMerror = prefDiff/prefMassRounded * this.million ;
				if((massErrorthVal != null || isNaN(massErrorthVal) )|| 
								(ppmErrorthVal != null || isNaN(ppmErrorthVal) ))
				{
					if((Math.abs(massDiff) <= massErrorthVal) || (Math.abs(prePPMerror) <= ppmErrorthVal))
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						position = j + 1;
						if(prefixInd !== "prefix")
						{
							position = preOrSufListln- position +1;
						}
						mass = prefixOrSuffixMassList[j];
						charge = monoMassList_temp[i].charge;
						matchedInd = "Y";
						break;
					}
				}
			}
			MatchedAndUnMatchedListObj = {position:position,mass:mass,matchedInd:matchedInd,charge:charge};
			MatchedAndUnMatchedList.push(MatchedAndUnMatchedListObj);
		}
		return MatchedAndUnMatchedList;
	}
  }