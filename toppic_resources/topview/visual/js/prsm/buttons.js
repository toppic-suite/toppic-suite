/**
 * On-click Actions to show spectrum graph and invoke downloads and other button actions
 */
function buttons(){
	//	Invocation on click of Mono M/z at the table
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
	
	// Coordinates at which a pop window to be launched to give name for the image to be downloaded
	// All the actions for different image buttons on different spectrums
	let x,y;
	d3.select("#download_popupms2_png").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		// function inside prsmtohtml.js to pop up a window
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
	// On click of help button, allow it to be dragged
	$("#spectrum_help").click(function(){
		$("#spectrumHelp").draggable({
		appendTo: "body"
		});
	});
}
/**
 * Show spectrum graphs on click of Show Spectrum
 */
function showSpectrun(){
	$("#show_ms2_spectrum").text('Hide Spectrum');
	document.getElementById("spectrum_help").style.display = "block";
	document.getElementById("graph_download").style.display = "block";
	document.getElementById("ms2svg_div").style.display = "block";
	document.getElementById("monoMassGraph").style.display = "block";
}
/**
 * Hide spectrum graphs on click of hide spectrum
 */
function hideSpectrum(){
	$("#show_ms2_spectrum").text('Show Spectrum'); 
	document.getElementById("spectrum_help").style.display = "none";
	document.getElementById("graph_download").style.display = "none";
	document.getElementById("ms2svg_div").style.display = "none";
	document.getElementById("monoMassGraph").style.display = "none";
}
