/**
 * Starting point of drawing spectrum graph
 * @param {String} id - gets the svg id of the graph 
 * @param {list} peakList - contains the list of data with mz and intensity used to draw lines on the graph 
 * @param {list} envelopeList - contains the list of data with actual mass, mz and intensity used to draw circles on graph
 * @param {float} monoMZ - Value to which tha graph as to point on click of mz value used to zoom the grpah to that location
 * @param {list} ionData - contains the list of data with mass and acid name
 * @param {object} graphFeatures - contains all the features needed for drawing the graphs
 */
addSpectrum = function(id,peakList,envelopeList,monoMZ, ionData, graphFeatures){
	let specParameters = compSpectrumParameters(peakList, envelopeList, monoMZ);
	specParameters.graphFeatures = graphFeatures;
	let peakData = {};
	peakData.peak_list = peakList ;
	if(envelopeList != null)
	{
		peakData.envelope_list = sortEnvelopes(envelopeList) ;
	}
	id = "#"+id;
	if(ionData == null)
	{
		specParameters.graphFeatures.showIons = false ;
	}
	else{
		ionData = groupBy(ionData,'ion');
	}
	specParameters.padding = graphFeatures.padding;
	specParameters.specWidth = graphFeatures.specWidth;
	specParameters.specHeight = graphFeatures.specHeight;
	specParameters.svgHeight = graphFeatures.svgHeight;
	specParameters.padding.head = graphFeatures.padding.head;
	let spectrumGraph = new SpectrumGraph(id,specParameters,peakData,ionData);
	return spectrumGraph;
}
 //  Function to set spectrum perameters based on the data
 function compSpectrumParameters(peakList, envelopeList, monoMZ){
	let ratio = 1; 
	let specParameters = new SpectrumParameters();
	// Sort by mz
	peakList.sort(function(x,y){
		return x.mz - y.mz;
	});
	let listSize = peakList.length;
	let minMzData = peakList[0].mz;
	let maxMzData = peakList[listSize-1].mz;
	// Sort by intensity
	peakList.sort(function(x,y){
		return x.intensity - y.intensity;
	});

	let maxEnvelope = -1;
	let minEnvelope = 1000000000; // A random onstant high value
	if(envelopeList != null)
	{
		ratio = 1.000684;
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
	}
	let maxIntensity = peakList[listSize-1].intensity;
	let minIntensity = peakList[0].intensity;
	let currminMz = minMzData ;
	let currmaxMz = maxMzData ;
	let currentMaxIntensity;
	if (maxIntensity > maxEnvelope) {
		currentMaxIntensity = maxIntensity;
	}
	else {
		currentMaxIntensity = maxEnvelope;
	}

	if(monoMZ != null)
	{
		//  Initializing with 1% of total intensity. If there exists no peaks in the 
		//  current range the intensity can't be 0 or nothing which will produce an 
		//  error in arithmatic calculations
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

	//  for specParameters, going to pass whichever value between peak max and envelope max that is bigger
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
	});
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

/**
 * 
 * @param {list} listData - contains list of data
 * @param {string} keyValue - contains keyword based on which new group is created
 * Function returns a map list with key and value 
 */
function groupBy(listData,keyValue){
	const map = new Map();
    listData.forEach((element)=>{
        const key = element[keyValue];
        const collection = map.get(key);
        if(!collection) map.set(key,[element]);
        else collection.push(element);
    });
    return map;
}
