function getSpectrum(folderName)
{
	let scanId_1 = $("#scanid1").val();
	let scanId_2 = $("#scanid2").val();
	if(scanId_1 != null && scanId_2 == null)
	{
		scanId_2 = parseInt(scanId_1) + 1;
	}else if(scanId_1 == null && scanId_2 != null){
		scanId_1 = scanId_2 - 1; 
	}
	console.log("scanId 1: ",scanId_1);
	console.log("scanId 2: ",scanId_2);
	
	$("#ms1_spectrum").show();
	var file_name1 = "../data/"+folderName+"/ms1_json/spectrum"+scanId_1+".js";
	var body= document.getElementsByTagName('body')[0];
 	var script1= document.createElement('script');
 	script1.src = file_name1 ;
    body.append(script1);
    script1.onload = function(){
		console.log("spectrum_data : ", ms1_data);
        let peakList = getPeakData(ms1_data);
        let envelopList = getEnvelopeData(ms1_data);
		addSpectrum("ms1_spectrum",peakList,envelopList,null);
		$("#ms2_spectrum").show();
		var file_name2 = "../data/"+folderName+"/ms2_json/spectrum"+scanId_2+".js";
		console.log("file_name2 : ", file_name2);
		var body= document.getElementsByTagName('body')[0];
		var script2= document.createElement('script');
		script2.src = file_name2 ;
		body.append(script2);
		script2.onload = function(){
			console.log("spectrum_data : ", ms2_data);
			let peakList = getPeakData(ms2_data);
			let envelopList = getEnvelopeData(ms2_data);
			addSpectrum("ms2_spectrum",peakList,envelopList,null);
		}
	}
	
    
 	/* document.getElementById("spectrum_script").src = file_name; */
 	
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