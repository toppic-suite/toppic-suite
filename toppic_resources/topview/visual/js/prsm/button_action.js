/**
 * Button Actions 
 */
function addButtonActions() {
	//	MS1 graph popup window 
	$("#precursor_mz").click(function(){
		$("#ms1_graph_popup_window").draggable({
			appendTo: "body"
		});
	});

	d3.select("#download_ms1_png_btn").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupNameWindow("png", "ms1_svg",x,y)
	})
	d3.select("#download_ms1_svg_btn").on("click",function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupNameWindow("svg", "ms1_svg",x,y)
	})

	// Save PrSM popup window
  d3.select('#save_prsm_btn').on("click",function(){
    //	Id of the pop_up svg
    let svgId = "prsm_popup_svg" ;
    let svg = d3.select("body").select("#"+svgId);
    svg.selectAll("#"+svgId + "_g").remove();
    let para = new parameters();
    let prsm = prsm_data.prsm ;
    [para,svgId] = buildSvg(para,prsm,svgId);
    //	Get the amount of skipped acid and write the amount 
    // 	of skipped acid at the start and end of the sequence 
    skippedAcidNotification(para,prsm,svgId) ;
    if(para.show_num)
    {
      // Get the numerical count at the start enad end of 
      // each row of sequence
      getNumValues(para,prsm,svgId);
    }
    //	Determine the start and end position of the sequence
    drawAnnoOfStartEndPosition(para,prsm,svgId);
    //	Draw Annotations
    annotations(para,prsm,svgId);
    //	Get the position of the fixed ptms and color them to red
    addColorToFixedPtms(prsm,svgId);
    //	Color the background of occurence of mass shift 
    massShiftBackgroundColor(para,prsm,svgId);
    //	set the dimensions of popup svg to default values
    document.getElementById("row-size").value = para.row_length ;
    document.getElementById("letter-width").value = para.letter_width ;
    document.getElementById("row-height").value = para.row_height ;
    document.getElementById("block-width").value = para.gap_width ;
    document.getElementById("num-width").value = para.numerical_width ;
    document.getElementsByName("show-num")[0].checked = para.show_num ;
    document.getElementsByName("show-skipped-lines")[0].checked = para.show_skipped_lines ;

    // Allows to drag the pop up windows
		$("#save_prsm_popup_window").draggable({
			appendTo: "body"
		});
  });

	d3.select('#prsm_graph_redraw_btn').on("click", function(){
    let svgId = "prsm_popup_svg" ;
    let svg = d3.select("body").select("#"+svgId);
    svg.selectAll("#"+svgId + "_g").remove();
		//	Get dimension parameters of svg
		var para = new parameters() ;
		para.row_length = parseInt(document.getElementById("row-size").value);
		para.letter_width = parseInt(document.getElementById("letter-width").value) ;
		para.row_height = parseInt(document.getElementById("row-height").value) ;
		para.gap_width = parseInt(document.getElementById("block-width").value) ;
		para.numerical_width = parseInt(document.getElementById("num-width").value) ;

		//	Check whether show numbers is checked
		if(document.getElementsByName("show-num")[0].checked)
		{
			para.show_num = true ;
		}
		else
		{
			//	Reduce the left and right margins when numbers are 
			// 	no needed at start and end of each row in svg
			para.show_num = false ;
			para.left_margin = 20;
			para.right_margin = 20;
		}
		//	Check to show skipped lines 
		if(document.getElementsByName("show-skipped-lines")[0].checked)
		{
			para.show_skipped_lines = true ;
		}
		else
		{
			para.show_skipped_lines = false ;
		}
		//	Redraw the svg with new dimension parameters
		prsm = prsm_data.prsm ;
		[para,svgId] = buildSvg(para,prsm,svgId);
		//	Get the amount of skipped acid and write the amount 
		//	of skipped acid at the start and end of the sequence 
		skippedAcidNotification(para,prsm,svgId) ;
		if(para.show_num)
		{
			// Get the numerical count at the start and end of 
			// each row of sequence
			getNumValues(para,prsm,svgId);
		}
		//	Determine the start and end position of the sequence
		drawAnnoOfStartEndPosition(para,prsm,svgId) ;
		//	Draw Annotations
		annotations(para,prsm,svgId);
		//	Get the position of the fixed ptms and color them to red
		addColorToFixedPtms(prsm,svgId);
		//	Color the background of occurence of mass shift
		massShiftBackgroundColor(para,prsm,svgId);
	});	

	$("#prsm_popup_help_btn").click(function(){
		$("#prsm_help_popup_window").draggable({
			appendTo: "body"
		});
	});

  //	Download the svg as ".svg" image
  d3.select('#prsm_popup_svg_btn').on("click",function(){
    let svgId = "prsm_popup_svg" ;
    let x = d3.event.pageX;
    let y = d3.event.pageY;
    popupNameWindow("svg",svgId,x,y);
  });
  //	Download svg as PNG Image
  d3.select('#prsm_popup_png_btn').on("click", function(){
    let svgId = "prsm_popup_svg" ;
    let x = d3.event.pageX;
    let y = d3.event.pageY;
    popupNameWindow("png",svgId,x,y);
  });

	// Show MS2 graph button 
	$("#ms2_graph_show_btn" ).click(function() {
    if($.trim($(this).text()) === 'Show Spectrum') {
      showMs2Graph(this);
    }
    else {
      hideMs2Graph();
    }
	});

	// MS2 graph help button 
	$("#ms2_graph_help_btn").click(function(){
		$("#ms2_graph_help_popup_window").draggable({
		appendTo: "body"
		});
	});

  // MS2 graph save button 
  $("#ms2_graph_save_btn").click(function(){
    // Get the Scan number of the Current spectrum graph showing on the screen
    let ms2Id = $('.ms2_graph_list.active')[0].id;
    let ms2Split = ms2Id.split("_");
    let ms2Index = parseInt(ms2Split[ms2Split.length-1]);
    console.log(ms2Index);
    let svgId = "popup_ms2_svg";
    let peaks = ms2SpecList[ms2Index].peaks;
    let envelopes = ms2SpecList[ms2Index].envelopes;
    let ions = {};
    let spGraph = new SpectrumGraph(svgId,peaks,envelopes,ions);
    spGraph.redraw();
    // This allows pop up window to be moved
    $("#ms2_graph_popup_window").draggable({
      appendTo: "body"
    });
  })

  // By changing the graph features, on click of redraw invoke this and generate the pop spectrun repectively to be downloaded
  $("#ms2_popup_redraw_btn").click(function(){
    let ms2Id = $('.ms2_graph_list.active')[0].id;
    let ms2Split = ms2Id.split("_");
    let ms2Index = parseInt(ms2Split[ms2Split.length-1]);
    let svgId = "popup_ms2_svg";
    let peaks = ms2SpecList[ms2Index].peaks;
    let envelopes = ms2SpecList[ms2Index].envelopes;
    let ions = [];
    let spGraph = new SpectrumGraph(svgId,peaks,envelopes,ions);
    let para = spGraph.para; 
    para.showEnvelopes = document.getElementsByName("show_envelops")[0].checked ;
    para.showIons = document.getElementsByName("show_ions")[0].checked ;
    spGraph.redraw();
  })
    
	
	// Coordinates at which a pop window to be launched to give name for the image to be downloaded
	// All the actions for different image buttons on different spectrums
	d3.select("#download_ms2_graph_png_btn").on("click",function(){
		let x = d3.event.pageX;
		let y = d3.event.pageY;
		// function inside prsmtohtml.js to pop up a window
		popupNameWindow("png", "popup_ms2_svg",x,y)
	})
	d3.select("#download_ms2_graph_svg_btn").on("click",function(){
		let x = d3.event.pageX;
		let y = d3.event.pageY;
		popupNameWindow("svg", "popup_ms2_svg",x,y)
	})

  // On click of mono mass mz, zoom all the graph to the corresponding point
  $(".row_mono_mz").click(function() {
    let parentId  = $(this).parent().parent().prop('id');
    let scanNum = document.getElementById(parentId).firstChild.innerHTML;
    /*	get Mono M/z value till 3 decimal values	*/
    let monoMz = parseFloat(this.innerHTML).toFixed(3) ;
    let index = 0;
    for (let i = 0; i < ms2SpecList.length; i++) {
      if (scanNum == ms2SpecList[i].scan) {
        index = i; 
        //console.log("match index", i);
        break;
      }
    }
    showMs2Graph();
    activateCurrentNavbar("ms2_graph_nav", index);
    let spGraph = ms2GraphList[index]; 
    // set monoMz to do
    spGraph.para.updateMzRange(monoMz);
    spGraph.redraw();
  });
}
/**
 * Show spectrum graphs on click of Show Spectrum
 */
function showMs2Graph(){
	$("#ms2_graph_show_btn").text('Hide Spectrum');
	document.getElementById("ms2_graph_help_btn").style.display = "block";
	document.getElementById("ms2_graph_save_btn").style.display = "block";
	document.getElementById("ms2_svg_div").style.display = "block";
	//document.getElementById("monoMassGraph").style.display = "block";
}
/**
 * Hide spectrum graphs on click of hide spectrum
 */
function hideMs2Graph(){
	$("#ms2_graph_show_btn").text('Show Spectrum'); 
	document.getElementById("ms2_graph_help_btn").style.display = "none";
	document.getElementById("ms2_graph_save_btn").style.display = "none";
	document.getElementById("ms2_svg_div").style.display = "none";
	//document.getElementById("monoMassGraph").style.display = "none";
}

/**
 * Function to produce a pop up window to provide name and set name to 
 * the image while downloading the image of Graph SVG and Sequence SVG
 * @param {String} type - Provides if the image is downloaded as svg or png
 * @param {String} id - Provides the id of the svg to be downloaded
 * @param {Float} x - Provides coordinate on where to show a tooltip block to enter name of the image to be downloaded
 * @param {Float} y - Provides coordinate on where to show a tooltip block to enter name of the image to be downloaded
 */
function popupNameWindow(type,id,x,y){
	d3.selectAll("#tooltip_imagename").remove() ;
	var div = d3.select("body").append("div")
	.attr("class", "tooltip")
	.attr("id","tooltip_imagename")
	.style("opacity", 1);

	// Provides a tooltip to enter a name for the image to be downloaded
	div.transition()
	.duration(200)
	.style("opacity", .9);
	div.html( 
			'<input type="text" placeholder="Image Name" id="imagename" />'+
			'<button id="saveimage" style = "none" type="button">save</button>'
			)
	.style("left", (x - 30) + "px")  // x Coordinate of the position of the tooltip           
	.style("top", (y - 60) + "px")	// y Coordinate of the position of the tooltip 
	.attr("box-sizing","border")
	.attr("display","inline-block")
	.attr("min-width","1.5em")
	.attr("padding","2px")
	.attr("margin-left","0px")
	.attr("text-align","center")
	.attr("text-decoration","none")
	.attr("border","1px solid #111111")
	.attr("background-color","white");
	
	// On click action to save the image on click of download button
	$("#saveimage").click(function(){
		let imagename = $("#imagename").val();
		if( imagename == null || imagename == "")
		{
			imagename = "spectrum";
		}
		// Check if the image needs to be downloaded as svg
		if(type == "svg"){
			d3.selectAll("#tooltip_imagename").remove() ;
      console.log(id);
			let svgContainer = d3.select("#"+id);
			let svgElement = svgContainer.node();
			svg2svg(svgElement,imagename);
		}
		// Check if the image needs to be downloaded as png
		if(type == "png"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let l_svgContainer = d3.select("#"+id);
			let svgString = getSVGString(l_svgContainer.node());
			let specParams =  new SpectrumParameters();
			let width = specParams.svgWidth;
			let height = specParams.svgHeight ;
			svgString2Image( svgString, 2*width, 2*height, 'png', save ); 
			function save( dataBlob, filesize ){
				saveAs( dataBlob, imagename ); 
			}
		}
	})
}

/**
 * Function highlights the active spectrum on clcik of the scan number
 * @param {String} id - Contains id of the SVG tag  
 * @param {int} currentValue - Contains scan number
 * To modify
 */
function activateCurrentNavbar(id, value){
  let childs = $("#"+id).children();
  let len = childs.length;
  /*
  for(let i=0;i<len;i++) {
    let grandchildValue = $(childs[i]).children().attr('value');
    if(grandchildValue == currentValue) {
      $("#"+id+" .active").removeClass("active");
      $(childs[i]).children().addClass("active");
      break;
    }
  }
  */
}
