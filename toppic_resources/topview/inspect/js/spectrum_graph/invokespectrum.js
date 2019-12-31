/*	Spectrum start point */
addSpectrum = function(peakList,envelopeList){
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
	specParameters.initScale(currminMz,currmaxMz,maxIntensity,minIntensity,minMzData,maxMzData,currentMaxIntensity);
	let peakData = {};
	peakData.peak_list = peakList ;
	peakData.envelope_list = envelopeList ;
	let spectrumgraph = new SpectrumGraph("#spectrum",specParameters,peakData);
	 return spectrumgraph;
}

