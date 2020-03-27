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

	matchedPeakAttributes(monoMassList_temp,id,ion,position,massDiff,mass,matchedInd)
	{
		let peak = monoMassList_temp;
		let thMass = "";
		let PPMerror = "";
		// Check if Matched Indicator is Yes "Y" eles keep empty spaces in the new attributes
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
		return peak;
	}

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
				let massDiff = monoMassList_temp[i].mass - prefixOrSuffixMassList[j].mass ;
				if(massErrorthVal != null || isNaN(massErrorthVal))
				{
					if(Math.abs(massDiff) <= massErrorthVal)
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						let ion = ionType + prefixOrSuffixMassList[j].position ;
						let position = prefixOrSuffixMassList[j].position ;
						if(preOrSufInd == "suffix")
						{
							position = seqln - position;
						}
						let mass = prefixOrSuffixMassList[j].mass;
						let matchedInd = "Y";
						peak = this.matchedPeakAttributes(monoMassList_temp[i],peakId,ion,position,
																massDiff,mass,matchedInd);
						matchedList.push(peak);		
					}
				}
				else{
					let prefDiff = Math.round(massDiff * 10000)/10000;
					let prefMassRounded = Math.round(prefixOrSuffixMassList[j].mass * 10000)/10000;
					let prePPMerror = prefDiff/prefMassRounded * this.million ;
					if(Math.abs(prePPMerror) <= ppmErrorthVal)
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						let ion = ionType + prefixOrSuffixMassList[j].position ;
						let position = prefixOrSuffixMassList[j].position ;
						if(preOrSufInd == "suffix")
						{
							position = seqln - position;
						}
						let mass = prefixOrSuffixMassList[j].mass;
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
	getMatchedAndUnMatchedList(monoMassList,matchedList)
	{
		let matchedAndUnmatchedList = matchedList.map(x => ({...x}));
		//matchedAndUnmatchedList.concat(matchedList);
		let len = matchedList.length;
		let MatchedPeaksObj = new MatchedPeaks();
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
				let peak = MatchedPeaksObj.matchedPeakAttributes(monoMassList[i],peakId,"","","","",matchedInd);
			  matchedAndUnmatchedList.push(peak);
			}
		})
		return matchedAndUnmatchedList ;
	}

	getDistribution(peakDataList,sequence,matchedUnMatchedList, completeMassShiftList)
	{
		let len = matchedUnMatchedList.length;
		let totalDistribution = [] ;
		let calEmassAndDistObj = new CalculateEmassAndDistribution();
		let molecularFormObj = new MolecularFormulae();
		let seqln = sequence.length;
		for(let i = 0; i < len; i++)
		{
			let distributionList = {};
      
			if(matchedUnMatchedList[i].matchedInd == "Y")
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
					let massShiftList = [];
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
			}
		}
		if(totalDistribution.length != 0)
		{
			totalDistribution.sort(function(x,y){
				return d3.ascending(x.env_peaks[0].mz, y.env_peaks[0].mz);
			})
		}
		
		let envlength = totalDistribution.length;
		let colorListsize = this.colors.length;
		while(envlength--){
			let index = envlength%colorListsize ;
			totalDistribution[envlength].color = this.colors[index] ;
		}

		let totalDistributionCopy = [];

		for(let i = 0; i < len; i++)
		{
			let distributionList = {};

			distributionList.mono_mass = matchedUnMatchedList[i].mass;
			distributionList.charge = matchedUnMatchedList[i].charge;
			distributionList.env_peaks = molecularFormObj.emass(distributionList.mono_mass,distributionList.charge,peakDataList);
			totalDistributionCopy.push(distributionList);
		}
		if(totalDistributionCopy.length != 0)
		{
			totalDistributionCopy.sort(function(x,y){
				return d3.ascending(x.env_peaks[0].mz, y.env_peaks[0].mz);
			})
		}
		
		envlength = totalDistributionCopy.length;
		colorListsize = this.colors.length;
		while(envlength--){
			let index = envlength%colorListsize ;
			totalDistributionCopy[envlength].color = this.colors[index] ;
		}
		this.compareTwoComputation(totalDistribution, totalDistributionCopy)
		return totalDistribution ;
	}
	compareTwoComputation(amino, molecular){
		let orderMatched = 0;
		let matchFound = false;
		let monomassMismatch = 0;
		let missingPeak = 0;
		let mzCount = {"tole_0.01":0, "tole_0.1":0, "tole_1":0, "tole_10":0, "tole_bigger_than_10":0};
		//console.log("amino envelopes: ", amino.length, ", molecular envelopes: ", molecular.length);
		for(let i = 0; i < molecular.length; i++){
			for (let j = 0; j < amino.length; j++){
				if (molecular[i].mono_mass == amino[j].mono_mass && molecular[i].charge == amino[j].charge){
					//compare env peaks
					if (molecular[i].env_peaks.length == amino[j].env_peaks.length){
						orderMatched++;
						for (let e = 0; e < molecular[i].env_peaks.length; e ++){
							let mzDiff = Math.abs(molecular[i].env_peaks[e].mz - amino[j].env_peaks[e].mz);
							let inteDiff = Math.abs(molecular[i].env_peaks[e].intensity - amino[j].env_peaks[e].intensity);
							
							if (mzDiff <= 0.01){
								mzCount["tole_0.01"]++;
							}
							else if (mzDiff <= 0.1){
								mzCount["tole_0.1"]++;
							}
							else if(mzDiff <= 1){
								mzCount["tole_1"]++;
							}
							else if (mzDiff <= 10){
								mzCount["tole_10"]++;
							}
							else{
								console.log("big difference in mz: ", molecular[i].env_peaks[e].mz, amino[j].env_peaks[e].mz)
								mzCount["tole_bigger_than_10"]++;
							}	
						}
						matchFound = true;
						break;
					}
					else{
						//console.log("env peaks length do not match")
						//console.log("amino env peak: ", molecular[i].env_peaks)
						//console.log("molecular env peak: ", amino[j].env_peaks)
						missingPeak++;
					}
				}
			}
			if (matchFound == false){
				//same mono mass and charge does not exist in second result
				console.log("envelope not found")
				console.log("molecular env peak", molecular[i])
				console.log("corresponding amino peak: ", amino[j]);
				
				monomassMismatch++;
			}
		}
		//result print
		console.log("peaks order match in ", orderMatched, " of ", molecular.length);
		console.log(monomassMismatch, " envelopes not matching")
		console.log(missingPeak, " envelopes have diff num of peaks")
		console.log("matchced mz by tolerannce level", mzCount);
	}
	getMatchedAndUnmatchedPrefixAndSuffixMassList(prefixOrSuffixMassList, monoMassList,
																massErrorthVal,ppmErrorthVal,prefixInd)
	{
		//console.log("monoMassList : ", monoMassList);
		let MatchedAndUnMatchedList = [];
		let monoMassList_temp = monoMassList.slice();
		let len = monoMassList_temp.length;
		let preOrSufListln = prefixOrSuffixMassList.length;
		let MatchedAndUnMatchedListObj = {};

		for(let j = 0; j<preOrSufListln; j++)
		{
			let position = prefixOrSuffixMassList[j].position;
			if(prefixInd != "prefix")
			{
				position = preOrSufListln-position+1;
			}
			let mass = prefixOrSuffixMassList[j].mass;
			let matchedInd = "N";
			let charge = 1;// Setting Default Value
			for(let i=0; i < len; i++)
			{
				let massDiff = monoMassList_temp[i].mass - prefixOrSuffixMassList[j].mass ;
				let prefDiff = Math.round(massDiff * 10000)/10000;
				let prefMassRounded = Math.round(prefixOrSuffixMassList[j].mass * 10000)/10000;
				let prePPMerror = prefDiff/prefMassRounded * this.million ;
				if((massErrorthVal != null || isNaN(massErrorthVal) )|| 
								(ppmErrorthVal != null || isNaN(ppmErrorthVal) ))
				{
					if((Math.abs(massDiff) <= massErrorthVal) || (Math.abs(prePPMerror) <= ppmErrorthVal))
					{
						//increment the position as the seq starts from 0 but the data in table starts from 1
						position = prefixOrSuffixMassList[j].position ;
						if(prefixInd != "prefix")
						{
							position = preOrSufListln-prefixOrSuffixMassList[j].position +1;
						}
						mass = prefixOrSuffixMassList[j].mass;
						matchedInd = "Y";
						charge = monoMassList_temp[i].charge;
						break;
					}
				}
			}
			MatchedAndUnMatchedListObj = {position:position,mass:mass,charge:charge,matchedInd:matchedInd};
			MatchedAndUnMatchedList.push(MatchedAndUnMatchedListObj);
		}
		return MatchedAndUnMatchedList;
	}
  }
