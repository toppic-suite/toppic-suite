/*	On-click Actions to show spectrum graph	*/
function buttons(){
	let peakdata = new PeakData();
	let peakList = peakdata.getPeakData(ms2_data);
	let envelopList = peakdata.getEnvelopeData(ms2_data);
	let ms1peakList = peakdata.getPeakData(ms1_data);
	let ms1envelopList = peakdata.getEnvelopeData(ms1_data);
	/*	Invocation on click of Mono M/z at the table	*/
	$( ".peakRows" ).click(function() {
		document.getElementById("ms1svg").style.display = "block";
		document.getElementById("show_spectrum").innerHTML = "Hide Spectrum";
		// document.getElementById("download_spectrum").style.display = "block";
		document.getElementById("spectrum_help").style.display = "block";
		/*	get Mono M/z value till 3 decimal values	*/
		this.peakValue = parseFloat(this.innerHTML).toFixed(3) ;
		var spectrumGraph = addSpectrum("ms1svg",peakList,envelopList,this.peakValue);
		//var spectrumGraph = addSpectrum(this.peakValue);
	});
	/*	Invocation on click of show spectrum button		*/
	$("#show_spectrum" ).click(function() {
		 if($.trim($(this).text()) === 'Show Spectrum')
		 {
			 document.getElementById("ms1svg").style.display = "block";
			//  document.getElementById("download_spectrum").style.display = "block";
			 document.getElementById("spectrum_help").style.display = "block";
			 document.getElementById("a_show_spectrum").href = "#"; 
			 $(this).text('Hide Spectrum');
		 }
		 else
		 {
			 $(this).text('Show Spectrum'); 
			 document.getElementById("ms1svg").style.display = "none";
			//  document.getElementById("a_show_spectrum").href = "#";
			 document.getElementById("download_spectrum").style.display = "none";
			document.getElementById("spectrum_help").style.display = "none";
		 }
		 this.peakValue = null ;
		 var spectrumGraph = addSpectrum("ms1svg",peakList,envelopList,null);
		 //var spectrumGraph = addSpectrum(this.peakValue);
	});
	$("#download_spectrum").click(function(){
		let svg_element = d3.selectAll("#ms1svg").node();
			svg2svg(svg_element);
	});
	
	$( "#precursormz" ).click(function(){
		let value = $(this).html();
		var spectrumGraph = addSpectrum("popupspectrum",ms1peakList,ms1envelopList,value);
	});
}