/*	Spectrum start point */
addSpectrum = function(id,peakList,envelopeList,monoMZ){
	let ratio = 1.000684;
	let specParameters = new SpectrumParameters();
	peakList.sort(function(x,y){
		return d3.ascending(x.mz, y.mz);
	})
	let listSize = peakList.length;
	let minMzData = peakList[0].mz;
	let maxMzData = peakList[listSize-1].mz;
	
	peakList.sort(function(x,y){
		return d3.ascending(x.intensity, y.intensity);
	})
	let maxIntensity = peakList[listSize-1].intensity;
	let minIntensity = peakList[0].intensity;
	let currminMz = minMzData ;
	let currmaxMz = maxMzData ;
	let currentMaxIntensity = maxIntensity ;
	if(monoMZ != null)
	{
		//Initializing with 1% of total intensity. If there exists no peaks in the 
		//current range the intensity can't be 0 or nothing which will produce an 
		//error in arithmatic calculations

		currentMaxIntensity = 0.01 * currentMaxIntensity ;
		monoMZ = monoMZ * ratio;
		currminMz = monoMZ - specParameters.onClickMassAdjacentRange ;
		currmaxMz = monoMZ + specParameters.onClickMassAdjacentRange ;
		for(let i = 0; i<listSize ; i++)
		{
			if(peakList[i].mz > currminMz && peakList[i].mz < currmaxMz)
			{
				if(peakList[i].intensity > currentMaxIntensity)
				{
					currentMaxIntensity = peakList[i].intensity;
				}
			}
		}
	}
	specParameters.initScale(currminMz,currmaxMz,maxIntensity,minIntensity,minMzData,maxMzData,currentMaxIntensity);
	let peakData = {};
	// Code must be included to get top 200 intensities at any given time
	peakList.sort(function(x,y){
		return d3.descending(x.intensity, y.intensity);
	})
	peakData.peak_list = peakList ;
	
	peakData.envelope_list = envelopeList ;
	id = "#"+id;
	let spectrumgraph = new SpectrumGraph(id,specParameters,peakData);
	return spectrumgraph;
}

