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
		let currentData = getCurrentData(ms2_ScansWithData,CurrentScanVal);
		console.log("currentData : ", currentData);
		generateCorrespondingGraph(currentData,"ms2svg",peak_value);
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
	d3.select("#graph_download_png").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY + 80;
		//function in prsmtohtml
		popupnamewindow("png", "ms2svg",x,y)
	})
	d3.select("#graph_download_svg").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY + 40;
		popupnamewindow("svg", "ms2svg",x,y)
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
function scanbuttons(){
	//ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms2_ScansWithData.length;
		let currentData = getCurrentData(ms2_ScansWithData,value);
		generateCorrespondingGraph(currentData,"ms2svg",null);
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	$(".ms1_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms1_ScansWithData.length;
		let currentData = getCurrentData(ms1_ScansWithData,value);
		generateCorrespondingGraph(currentData,"popupspectrum",null);
		$("#ms1_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
}