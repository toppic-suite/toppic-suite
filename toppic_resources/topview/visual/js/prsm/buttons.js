/*	On-click Actions to show spectrum graph	*/
function buttons(){
	/*	Invocation on click of Mono M/z at the table	*/
	$("#show_ms2_spectrum" ).click(function() {
		 
		 if($.trim($(this).text()) === 'Show Spectrum')
		 {
			 showSpectrun(this);
		 }
		 else
		 {
			hideSpectrum();
		 }
	});
	
	// let specParameters = new SpectrumParameters();
	// let spectrumDownload = new SpectrumDownload();
	let x,y;
	d3.select("#download_popupms2_png").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		//function in prsmtohtml
		popupnamewindow("png", "popup_ms2_spectrum",x,y)
	})
	d3.select("#download_popupms2_svg").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("svg", "popup_ms2_spectrum",x,y)
	})
	d3.select("#download_popup_png").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("png", "popupspectrum",x,y)
	})
	d3.select("#download_popup_svg").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("svg", "popupspectrum",x,y)
	})

	$("#precursormz").click(function(){
		let prec_mz = $("#precursormz").html();
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
function showSpectrun(){
	$("#show_ms2_spectrum").text('Hide Spectrum');
	//document.getElementById("ms2svg").style.display = "block";
	document.getElementById("spectrum_help").style.display = "block";
	document.getElementById("graph_download").style.display = "block";
	document.getElementById("ms2svg_div").style.display = "block";
	document.getElementById("monoMassGraph").style.display = "block";
}
function hideSpectrum(){
	$("#show_ms2_spectrum").text('Show Spectrum'); 
	//document.getElementById("ms2svg").style.display = "none";
	document.getElementById("spectrum_help").style.display = "none";
	document.getElementById("graph_download").style.display = "none";
	document.getElementById("ms2svg_div").style.display = "none";
	document.getElementById("monoMassGraph").style.display = "none";
}
