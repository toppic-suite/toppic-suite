/*	On-click Actions to show spectrum graph	*/
function buttons(){
	/*	Invocation on click of Mono M/z at the table	*/
	$( ".peakRows" ).click(function() {
		let parent_id  = $(this).parent().parent().prop('id');
		let CurrentScanVal = document.getElementById(parent_id).firstChild.innerHTML;
		/*	get Mono M/z value till 3 decimal values	*/
		let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
		console.log("ms2_ScansWithData in buttons: ", ms2_ScansWithData);
		console.log("CurrentScanVal : ", CurrentScanVal);
		let [currentData,specId] = getCurrentData(ms2_ScansWithData,CurrentScanVal);
		console.log("currentData : ", currentData);
		generateCorrespondingGraph(currentData,"ms2svg",peak_value,specId);
		activateCurrentnavbar("ms2_graph_nav",CurrentScanVal)
		showSpectrun();
	});
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
	
	let specParameters = new SpectrumParameters();
	let spectrumDownload = new SpectrumDownload();
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
	$("#graph_download").click(function(){
		//let ms2svg = d3.select("#ms2svg");
		$("#ms2_graph_popup_nav li").remove();
		let scanIdList = prsm_data.prsm.ms.ms_header.scans.split(" ");
		let MultiScanObj = new MultiScan();
		MultiScanObj.createMs2NavEements(scanIdList,"ms2_graph_popup_nav");
		console.log("ms2_ScansWithData : ", ms2_ScansWithData);
		let [current_data,specId] = getCurrentData(ms2_ScansWithData,scanIdList[0]);
		//ms2_graph.redraw(null,"popup_ms2_spectrum");
		generateCorrespondingGraph(current_data,"popup_ms2_spectrum",null,specId);
		$("#ms2spectrumpop").draggable({
			appendTo: "body"
		});
		//scanbuttons();
	})
	
}
function showSpectrun(){
	$("#show_ms2_spectrum").text('Hide Spectrum');
	document.getElementById("ms2svg").style.display = "block";
	document.getElementById("spectrum_help").style.display = "block";
	document.getElementById("graph_download").style.display = "block";
	document.getElementById("ms2svg_div").style.display = "block";
}
function hideSpectrum(){
	$("#show_ms2_spectrum").text('Show Spectrum'); 
	document.getElementById("ms2svg").style.display = "none";
	document.getElementById("spectrum_help").style.display = "none";
	document.getElementById("graph_download").style.display = "none";
	document.getElementById("ms2svg_div").style.display = "none";
}
