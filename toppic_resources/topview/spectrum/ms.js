class Spectrum{
	constructor(folder,scanId,spectrumId,centerValue){
		this.folder = folder ;
		this.scanId = scanId ;
		this.spectrumId = spectrumId;
		this.centerValue = centerValue ;
	}
	getSpectrum(){
		console.log("this.scanId : "+ this.scanId);
		let scanId = this.scanId;
		let folder = this.folder;
		let spectrumId = this.spectrumId;
		$("."+this.spectrumId).show();
		var file_name = "../../topfd/"+folder+"/spectrum"+scanId+".js";
		let temp_script= document.createElement('script');
		temp_script.src = file_name;
		document.head.appendChild(temp_script);
		temp_script.onload = function()
		{
			let ms_data;
			if(folder == "ms1_json") ms_data = ms1_data ;
			else ms_data = ms2_data;
			let peakList = getPeakData(ms_data);
        	let envelopList = getEnvelopeData(ms_data);
			addSpectrum(spectrumId,peakList,envelopList,this.centerValue);
		}

	}
}

function getPeakData(json_data){
    let peakList = [];
	let i = json_data.peaks.length ;
	/*	Pass the data into peakmass and peakintensity ----------------------------*/
	while(i--)
	{
	    mz = parseFloat(json_data.peaks[i].mz);
	    intensity = parseFloat(json_data.peaks[i].intensity);
	    peak = {mz:mz, intensity:intensity}
	    peakList[i] = peak
	}
	peakList.sort(function(x,y){
		return d3.ascending(x.mz, y.mz);
	})
	return peakList;
}
function getEnvelopeData(json_data){
    circleColor_list = ["red","orange","blue","green"];
	let envelopList = [];
	json_data.envelopes.sort(function(x,y){
		return d3.ascending(x.env_peaks[0].mz, y.env_peaks[0].mz);
	})
	
	let i = json_data.envelopes.length ;
	let colorListsize = circleColor_list.length;
	
	while(i--)
	{
		let env_peaks = [];
		let mono_mass = parseFloat(json_data.envelopes[i].mono_mass);
		let charge = parseFloat(json_data.envelopes[i].charge);
		let color = circleColor_list[i%colorListsize];
		j = json_data.envelopes[i].env_peaks.length ;
		while(j--){
			let mz = parseFloat(json_data.envelopes[i].env_peaks[j].mz);
			let intensity = parseFloat(json_data.envelopes[i].env_peaks[j].intensity);
			let env_peak = {mz:mz,intensity:intensity}
			env_peaks[j] = env_peak ;
		}
		let envelope = {mono_mass:mono_mass,charge:charge,env_peaks:env_peaks,color:color};
		envelopList[i] = envelope;
    }
    return envelopList;
}