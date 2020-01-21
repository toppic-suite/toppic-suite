/*	On-click Actions to show spectrum graph	*/
function buttons(){
	/*	Invocation on click of Mono M/z at the table	*/
	console.log("dummmmmmm");
	$( ".peakRows" ).click(function() {
		document.getElementById("ms2svg").style.display = "block";
		document.getElementById("spectrum_help").style.display = "block";
		/*	get Mono M/z value till 3 decimal values	*/
		let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
		console.log("this : ", this);
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
			 console.log("In Show Spectrum");
			 document.getElementById("ms2svg").style.display = "block";
			 document.getElementById("spectrum_help").style.display = "block";
			 document.getElementById("graph_download").style.display = "block";
			 document.getElementById("ms2svg_div").style.display = "block";
			 //document.getElementById("a_show_spectrum").href = "#"; 
			 let current_data = ms2_ScansWithData[0].value;
			 generateCorrespondingGraph(current_data,"ms2svg",null);
			 $(this).text('Hide Spectrum');
		 }
		 else
		 {
			 $(this).text('Show Spectrum'); 
			 document.getElementById("ms2svg").style.display = "none";
			 document.getElementById("spectrum_help").style.display = "none";
			 document.getElementById("graph_download").style.display = "none";
			 document.getElementById("ms2svg_div").style.display = "none";
			 //document.getElementById("download_spectrum").style.display = "none";
			 //document.getElementById("spectrum_help").style.display = "none";
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
		popupnamewindow("png", "ms2svg",x,y)
	})
	d3.select("#download_popup_png").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("png", "popupspectrum",x,y)
	})
	d3.select("#download_popup_svg").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("png", "popupspectrum",x,y)
	})

	$("#precursormz").click(function(){
		let prec_mz = $("#precursormz").html();
		let current_data = ms1_ScansWithData[0].value;
		generateCorrespondingGraph(current_data,"popupspectrum",prec_mz) 
		$("#spectrumpop").draggable({
		appendTo: "body"
    });
	});

	$("#spectrum_help").click(function(){
    $("#spectrumHelp").draggable({
      appendTo: "body"
    });
	});
	//ms2_scanIds is the Id of the nav tabs for multiple navs
	$(".ms2_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms2_ScansWithData.length;
		for(let i=0;i<len;i++)
		{
			if(ms2_ScansWithData[i].key == parseInt(value))
			{
				console.log("Equal");
				let current_data = ms2_ScansWithData[i].value;
				generateCorrespondingGraph(current_data,"ms2svg");
				break;
			}
		}
		$("#ms2_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	$(".ms1_scanIds").click(function(){
		let value = this.getAttribute('value')
		let len = ms2_ScansWithData.length;
		for(let i=0;i<len;i++)
		{
			if(ms2_ScansWithData[i].key == parseInt(value))
			{
				console.log("Equal");
				let current_data = ms2_ScansWithData[i].value;
				generateCorrespondingGraph(current_data,"popupspectrum");
				break;
			}
		}
		$("#ms1_graph_nav .active").removeClass("active");
   		$(this).addClass("active");
	})
	
}
