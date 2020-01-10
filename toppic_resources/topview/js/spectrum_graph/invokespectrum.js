/*	Spectrum start point */
addSpectrum = function(id,peakList,envelopeList,monoMZ){
  let specParameters = compSpectrumParameters(peakList, envelopeList, monoMZ);
	let peakData = {};
	peakData.peak_list = peakList ;
	peakData.envelope_list = sortEnvelopes(envelopeList) ;
	id = "#"+id;
	// console.log("peakData : ", peakData);
	let spectrumgraph = new SpectrumGraph(id,specParameters,peakData);
	return spectrumgraph;
}

compSpectrumParameters = function (peakList, envelopeList, monoMZ) {
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

  let maxEnvelope = -1;
  let minEnvelope = 1000000000;

	envelopeList.forEach(function(element){
		element.env_peaks.forEach(function(e){
			if (e.intensity > maxEnvelope){
				maxEnvelope = e.intensity;
			}
			else if (e.intensity < minEnvelope){
				minEnvelope = e.intensity;
			}
		})
	})

	let maxIntensity = peakList[listSize-1].intensity;
	let minIntensity = peakList[0].intensity;
	let currminMz = minMzData ;
	let currmaxMz = maxMzData ;
	//let currentMaxIntensity = maxIntensity ;
	let currentMaxIntensity;

	if (maxIntensity > maxEnvelope) {
		currentMaxIntensity = maxIntensity;
	}
	else {
		currentMaxIntensity = maxEnvelope;
	}

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

	//for specParameters, going to pass whichever value between peak max and envelope max that is bigger.

	if (maxIntensity < maxEnvelope){
		maxIntensity = maxEnvelope;
	}
	if (minIntensity > minEnvelope){
		minIntensity = minEnvelope;
	}
	specParameters.initScale(currminMz,currmaxMz,maxIntensity,minIntensity,minMzData,maxMzData,currentMaxIntensity);

	// Code must be included to get top 200 intensities at any given time
	peakList.sort(function(x,y){
		return d3.descending(x.intensity, y.intensity);
	})


  return specParameters;
}
/**
 * Sorting envelopes based on intensity to show top 200 envelops with high intensitites
 */
sortEnvelopes = function(envelopeList)
{
	envelopeList.forEach(function(envelope){
		// ...env_peak converts arguments into list as sort function is not allowed on arguments
		envelope.env_peaks.forEach(function(...env_peak){
			env_peak.sort(function(x,y){
				return d3.descending(x.intensity, y.intensity);
			})
		})
	})
	envelopeList.sort(function(x,y){
	 	return d3.descending(x.env_peaks[0].intensity, y.env_peaks[0].intensity);
	})
	return envelopeList ;
}
