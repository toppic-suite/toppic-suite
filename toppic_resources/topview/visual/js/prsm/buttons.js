/*	On-click Actions to show spectrum graph	*/
function buttons(){
	/*	Invocation on click of Mono M/z at the table	*/
	$( ".peakRows" ).click(function() {
		document.getElementById("ms1svg").style.display = "block";
		document.getElementById("spectrum_help").style.display = "block";
		/*	get Mono M/z value till 3 decimal values	*/
		let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
    ms2_graph.redraw(peak_value);
	});
	/*	Invocation on click of show spectrum button		*/
  /*
	$("#show_ms1_spectrum" ).click(function() {
		$("#show_ms2_spectrum" ).text("Show Ms2");

		let peak_data = new PeakData();
    let ms1_peak_list = peak_data.getPeakData(ms1_data);
		let ms1_envelope_list = peak_data.getEnvelopeData(ms1_data);
		ms1_graph = addSpectrum("ms1svg", ms1_peak_list,
		  								ms1_envelope_list,null);
		  
		 if($.trim($(this).text()) === 'Show Ms1')
		 {
			 document.getElementById("ms1svg").style.display = "block";
			 document.getElementById("spectrum_help").style.display = "block";
			 //document.getElementById("a_show_spectrum").href = "#"; 
			 $(this).text('Hide Ms1');
		 }
		 else
		 {
			 $(this).text('Show Ms1'); 
			 document.getElementById("ms1svg").style.display = "none";
			 document.getElementById("spectrum_help").style.display = "none";
			 //document.getElementById("download_spectrum").style.display = "none";
			 //document.getElementById("spectrum_help").style.display = "none";
		 }
	});
  */
	$("#show_ms2_spectrum" ).click(function() {
		//$("#show_ms1_spectrum" ).text("Show Ms1");

    /*
		let peak_data = new PeakData();
    let ms1_peak_list = peak_data.getPeakData(ms2_data);
		let ms1_envelope_list = peak_data.getEnvelopeData(ms2_data);
		ms1_graph = addSpectrum("ms1svg", ms1_peak_list,
		  								ms1_envelope_list,null);
                      */
		 
		 if($.trim($(this).text()) === 'Show Spectrum')
		 {
			 document.getElementById("ms1svg").style.display = "block";
			 document.getElementById("spectrum_help").style.display = "block";
			 //document.getElementById("a_show_spectrum").href = "#"; 
			 $(this).text('Hide Spectrum');
		 }
		 else
		 {
			 $(this).text('Show Spectrum'); 
			 document.getElementById("ms1svg").style.display = "none";
			 document.getElementById("spectrum_help").style.display = "none";
			 //document.getElementById("download_spectrum").style.display = "none";
			 //document.getElementById("spectrum_help").style.display = "none";
		 }
	});
	

	$("#precursormz").click(function(){
		let prec_mz = $("#precursormz").html();
    ms1_graph.redraw(prec_mz);
    $("#spectrumpop").draggable({
      appendTo: "body"
    });
	});

	$("#spectrum_help").click(function(){
    $("#spectrumHelp").draggable({
      appendTo: "body"
    });
	});

}
