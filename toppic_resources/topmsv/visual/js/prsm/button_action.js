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

  // On click of mono mass mz, zoom all the graph to the corresponding point
  $(".row_mono_mz").click(function() {
    let parentId  = $(this).parent().parent().prop('id');
    let scanNum = document.getElementById(parentId).firstChild.innerHTML;
    /*	get Mono M/z value till 3 decimal values	*/
    let monoMz = parseFloat(this.innerHTML).toFixed(3) ;
    for (let i = 0; i < ms2SpecList.length; i++) {
      let listId = "ms2_svg_div_graphlist_" + i;
      let graphId = "ms2_svg_div_graph_" + i;
      let monolistId = "ms2_svg_div_monographlist_" + i;
      let monoGraphId = "ms2_svg_div_mono_graph_" + i;
      let listElement = document.getElementById(listId);
      let graphElement = document.getElementById(graphId);
      let monoListElement = document.getElementById(monolistId);
      let monoGraphElement = document.getElementById(monoGraphId);
      if (scanNum == ms2SpecList[i].scan) {
        listElement.classList.add("active");
        graphElement.style.display="";
        monoListElement.classList.remove("active");
        monoGraphElement.style.display="none";
        let spGraph = ms2GraphList[i]; 
        // set monoMz to do
        spGraph.para.updateMzRange(monoMz);
        spGraph.redraw();
      }
      else {
        listElement.classList.remove("active");
        graphElement.style.display="none";
        monoListElement.classList.remove("active");
        monoGraphElement.style.display="none";
      }
    }
    showMs2Graph();
  });

  // On click of break points
  $(".break_point").click(function() {
    let pos = $(this).attr("ion_pos"); 
    showIonPeaks(pos); 
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



