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

	getDistribution(peakDataList,sequence,matchedUnMatchedList)
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
					distributionList.mono_mass = matchedUnMatchedList[i].mass;
					distributionList.charge = matchedUnMatchedList[i].charge;
					distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.PREFIX);
					totalDistribution.push(distributionList);
				}
				else
				{
					let matchedPos = matchedUnMatchedList[i].position ;
					let seq = sequence.slice(matchedPos,seqln) ;
					distributionList.mono_mass = matchedUnMatchedList[i].mass;
					distributionList.charge = matchedUnMatchedList[i].charge;
					distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.SUFFIX);
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
		return totalDistribution ;
	}
	getMatchedAndUnmatchedPrefixAndSuffixMassList(prefixOrSuffixMassList, monoMassList,
																massErrorthVal,ppmErrorthVal,prefixInd)
	{
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
						break;
					}
				}
			}
			MatchedAndUnMatchedListObj = {position:position,mass:mass,matchedInd:matchedInd};
			MatchedAndUnMatchedList.push(MatchedAndUnMatchedListObj);
		}
		return MatchedAndUnMatchedList;
	}
  }