/*	Create title and url for "all proteins","protein" and "proteoform" of the prsm */
function BuildUrl(folderpath)
{
	document.title = "Protein-Spectrum-Match for Spectrum #"+prsm_data.prsm.ms.ms_header.ids;
	let l_allproteins_url = "proteins.html?folder="+folderpath;
	document.getElementById("allprotein_url").href = l_allproteins_url;
	document.getElementById("allprotein_url_end").href = l_allproteins_url;
	let l_protein_URL = prsm_data.prsm.annotated_protein.sequence_name + " " + prsm_data.prsm.annotated_protein.sequence_description;
	document.getElementById("protien_url").innerHTML = l_protein_URL;
	document.getElementById("protien_url").href = "protein.html?folder="+folderpath+"&protein_Id="+prsm_data.prsm.annotated_protein.sequence_id;
	document.getElementById("protien_url_end").innerHTML = l_protein_URL;
	document.getElementById("protien_url_end").href = "protein.html?folder="+folderpath+"&protein_Id="+prsm_data.prsm.annotated_protein.sequence_id;
	let l_protroform_URL = "Proteoform #"+prsm_data.prsm.annotated_protein.proteoform_id ;
	document.getElementById("proteoform_url").innerHTML = l_protroform_URL;
	document.getElementById("proteoform_url").href = "proteoform.html?folder="+folderpath+"&proteoform_Id="+prsm_data.prsm.annotated_protein.proteoform_id;
	document.getElementById("proteoform_url_end").innerHTML = l_protroform_URL;
	document.getElementById("proteoform_url_end").href ="proteoform.html?folder="+folderpath+"&proteoform_Id="+prsm_data.prsm.annotated_protein.proteoform_id;
	document.getElementById("Protein-Spectrum-Match-Id-SpecId").innerHTML ="Protein-Spectrum-Match #"+prsm_data.prsm.prsm_id+" for Spectrum #" + prsm_data.prsm.ms.ms_header.ids;
}
/*	Get the data of the prsm from data variable prsm_data and convert to Html */
function loadDatafromJson2Html(){
	document.getElementById("File_name").innerHTML = prsm_data.prsm.ms.ms_header.spectrum_file_name;
	document.getElementById("PrSM_ID").innerHTML = prsm_data.prsm.prsm_id;
	document.getElementById("Scan").innerHTML = prsm_data.prsm.ms.ms_header.scans;
	document.getElementById("Precursor_charge").innerHTML = prsm_data.prsm.ms.ms_header.precursor_charge;
	document.getElementById("precursormz").innerHTML = prsm_data.prsm.ms.ms_header.precursor_mz ;
	document.getElementById("Precursor_mass").innerHTML = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
	document.getElementById("Proteoform_mass").innerHTML = prsm_data.prsm.annotated_protein.proteoform_mass;
	document.getElementById("matched_peaks").innerHTML = prsm_data.prsm.matched_peak_number;
	document.getElementById("matched_fragment_ions").innerHTML = prsm_data.prsm.matched_fragment_number;
	document.getElementById("unexpected_modifications").innerHTML = prsm_data.prsm.annotated_protein.unexpected_shift_number;
	document.getElementById("E_value").innerHTML = prsm_data.prsm.e_value;
	document.getElementById("P_value").innerHTML = prsm_data.prsm.p_value;
	document.getElementById("Q_value").innerHTML = prsm_data.prsm.fdr;
}
/*	Convert peak data into table format with link to m/z value */
function createTableElements(){
	var table = document.getElementById('spectrum');
	var tbdy = document.createElement('tbody');
	var ionArray = []; //contains the ion types used in the data


	let l_scans = prsm_data.prsm.ms.ms_header.scans.split(" ") ;
	let l_specIds = prsm_data.prsm.ms.ms_header.ids.split(" ") ;
	let l_matched_peak_count = 0;
	prsm_data.prsm.ms.peaks.peak.forEach(function(peak,i){
		/*	Check if peak contain matched_ions_num attribute	*/

		//for each peak, get the ion type and store it in ionArray to determine which ion type to be checked in Inspect
		if (parseInt(peak.matched_ions_num)>0){
			let ion = peak.matched_ions.matched_ion.ion_type;
			if (ionArray.length < 1){
				ionArray.push(ion);
			}			
			else if (ionArray.indexOf(ion) < 0) {
				ionArray.push(ion);
			}
		}
		
		if(peak.hasOwnProperty('matched_ions_num') && parseInt(peak.matched_ions_num)>1)
		{
			peak.matched_ions.matched_ion.forEach(function(matched_ion,i){
				peak.matched_ions.matched_ion = matched_ion ;
				loop_matched_ions(peak,i) ;
			})
		}
		else
		{
			loop_matched_ions(peak,i) ;
		}
	})
	
	//after looping through the prsm files, store the ion type data to local storage
	window.localStorage.setItem('ionType', ionArray);
	
	function loop_matched_ions(peak,i){
		/*	Create row for each peak value object in the table	*/
		var tr = document.createElement('tr');
		id = peak.spec_id+"peak"+peak.peak_id;
		let l_scan;
		if((parseInt(peak.peak_id) + 1)%2 == 0)
		{
			l_class = "unmatched_peak even";
		}
		else
		{
			l_class = "unmatched_peak odd";
		}
		if(peak.hasOwnProperty('matched_ions_num'))
		{
			id = id + peak.matched_ions.matched_ion.ion_type;
			if((parseInt(peak.peak_id) + 1)%2 == 0)
			{
				l_class = "matched_peak even";
			}
			else
			{
				l_class = "matched_peak odd";
			}
			l_matched_peak_count++;
			/*	create a name for each row */
			tr.setAttribute("name",peak.matched_ions.matched_ion.ion_position);
		}
		/*	Set "id","class name" and "role" for each row	*/
		tr.setAttribute("id", id);
		tr.setAttribute("class",l_class);
		tr.setAttribute("role","row");
		for(let i = 0;i<11;i++){
			var td = document.createElement('td');
			td.setAttribute("align","center");
			if(i == 0)
			{
				if(peak.spec_id == l_specIds[0]) l_scan = l_scans[0];
				else l_scan = l_scans[1];
				td.setAttribute("class","row_scanIds");
				td.innerHTML = l_scan ;
			}
			if(i == 1)
			{
				td.innerHTML = parseInt(peak.peak_id) + 1 ;
				td.setAttribute("class","row_peakNum");
			}
			if(i == 2)
			{
				td.innerHTML = peak.monoisotopic_mass;
				td.setAttribute("class","row_monoMass");
			}
			if(i == 3)
			{
				/*	provide link to click on m/z value to view spectrum */
				let a = document.createElement('a');
				a.href="#!"
				a.className = "peakRows"
				a.innerHTML = peak.monoisotopic_mz;
				td.appendChild(a);
			}
			if(i == 4)
			{
				td.innerHTML = peak.intensity;
				td.setAttribute("class","row_intensity");
			}
			if(i == 5)
			{
				td.innerHTML = peak.charge;
				td.setAttribute("class","row_charge");
			}
			if(peak.hasOwnProperty('matched_ions_num'))
			{
				if(i == 6)
				{
					td.innerHTML = peak.matched_ions.matched_ion.theoretical_mass;
				}
				if(i == 7)
				{
					td.innerHTML = peak.matched_ions.matched_ion.ion_type+peak.matched_ions.matched_ion.ion_display_position;
				}
				if(i == 8)
				{
					td.innerHTML = peak.matched_ions.matched_ion.ion_position;
				}
				if(i == 9)
				{
					td.innerHTML = peak.matched_ions.matched_ion.mass_error;
				}
				if(i == 10)
				{
					td.innerHTML = peak.matched_ions.matched_ion.ppm;
				}
			}
			tr.appendChild(td);
		}
		tbdy.appendChild(tr);
	}
	let l_All_Peaks = prsm_data.prsm.ms.peaks.peak.length;
	let l_not_matched_peak_count = l_All_Peaks - l_matched_peak_count; 
	document.getElementById("all_peak_count").innerHTML = "All peaks (" + l_All_Peaks + ")" ;
	document.getElementById("matched_peak_count").innerHTML = "Matched peaks (" + l_matched_peak_count + ")" ;
	document.getElementById("not_matched_peak_count").innerHTML = "Not Matched peaks (" + l_not_matched_peak_count + ")" ;
	
	table.appendChild(tbdy);

}
/*	Get occurence of "Variable" and "Fixed", convert the data to HTML */
function occurence_ptm(prsm)
{
	let variable_ptm = "";
	let fixed_ptm = "" ;
	if(prsm.annotated_protein.annotation.hasOwnProperty('ptm'))
	{
		let ptm = prsm.annotated_protein.annotation.ptm ;
		if(Array.isArray(ptm))
		{
			ptm.forEach(function(ptm, index){
				
				if(ptm.ptm_type == "Variable")
				{
					variable_ptm_raw = getVariablePtm(ptm).replace("[", " [");
					variable_ptm = variable_ptm + variable_ptm_raw + "] " ;
				}
				if(ptm.ptm_type == "Fixed")
				{
					fixed_ptm_raw = getFixedPtm(ptm).replace("[", " [");
					fixed_ptm = fixed_ptm + fixed_ptm_raw + "] " ;
				}
			})
		}
		else
		{
			if(ptm.ptm_type == "Variable")
			{
				variable_ptm_raw = getVariablePtm(ptm).replace("[", " [");
				variable_ptm = variable_ptm + variable_ptm_raw  + "]" ;
			}
			if(ptm.ptm_type == "Fixed")
			{
				fixed_ptm_raw = getFixedPtm(ptm).replace("[", " [");
				fixed_ptm = fixed_ptm + fixed_ptm_raw + "]" ;
			}
		}
	}
	/*	Convert ptms to Html	*/
	if(fixed_ptm != "")
	{
		let div = document.getElementById("ptm_abbreviation") ;
		let text1 = document.createElement("text");
		let text2 = document.createElement("text");
		text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"+"Fixed PTMs: ";
		text2.innerHTML = fixed_ptm ;
		text2.target = "_blank";
		text2.style = "color:red" ;
		div.appendChild(text1);
		div.appendChild(text2);
	}
	if(variable_ptm != "")
	{
		let div = document.getElementById("ptm_abbreviation") ;
		let text1 = document.createElement("text");
		let text2 = document.createElement("text");
		text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"+"Variable PTMs: ";
		text2.innerHTML = variable_ptm ;
		text2.target = "_blank";
		text2.style = "color:red" ;
		div.appendChild(text1);
		div.appendChild(text2);
	}
}
/* Get all the unknown ptms information */
function getUnknownPtms(prsm)
{
	let data = " [";
	if(prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
	{
		let mass_shift = prsm.annotated_protein.annotation.mass_shift ;
		let UnknownExist =  false;
			if(Array.isArray(mass_shift))
			{
				let len = mass_shift.length;
				mass_shift.forEach(function(each_mass_shift, i){
					if(each_mass_shift.shift_type == "unexpected")
					{
						UnknownExist = true;
						if(i != len-1)
						{
							data = data + each_mass_shift.anno + "," ;
						}
						else{
							data = data + each_mass_shift.anno + "]" ;
						}
					}
				})
			}
			else if(mass_shift.shift_type == "unexpected")
			{
				UnknownExist = true;
				data = data + mass_shift.anno + "]" ;
			}
			if(UnknownExist)
			{
				let val = "Unknown" + data;
				document.getElementById("ptm_unexpectedmodification").style.display = "block";

				let div = document.getElementById("ptm_unexpectedmodification") ;
				let text1 = document.createElement("text");
				let text2 = document.createElement("text");
				text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"+"Unexpected modifications: ";
				text2.innerHTML = val ;
				text2.target = "_blank";
				text2.style = "color:red" ;
				div.appendChild(text1);
				div.appendChild(text2);
			}
	}
	
}
/*	Get all the variable ptms	*/
function getVariablePtm(ptm)
{
	let variable_ptm = "[" ;
	let abbrevation = ptm.ptm.abbreviation ;
	if(Array.isArray(ptm.occurence))
	{
		ptm.occurence.forEach(function(occurence,i){
      let left = parseInt(occurence.left_pos) + 1;
			variable_ptm = variable_ptm + left + "-" + occurence.right_pos ;
			if(ptm.occurence.length-1 > i )
			{
				variable_ptm = variable_ptm + ";" ;
			}
		})
	}
	else
	{
    let left = parseInt(ptm.occurence.left_pos) + 1;
		variable_ptm = variable_ptm + left + "-" + ptm.occurence.right_pos ;
	}
	variable_ptm = ptm.ptm.abbreviation + variable_ptm ;
	return variable_ptm ;
}
/*	Get all the "Fixed" ptms*/
function getFixedPtm(ptm)
{
	let fixed_ptm = "[" ;
	let abbrevation = ptm.ptm.abbreviation ;
	if(Array.isArray(ptm.occurence))
	{
		ptm.occurence.forEach(function(occurence,i){
			//fixed_ptm = fixed_ptm + occurence.left_pos + "-" + occurence.right_pos ;
      fixed_ptm = fixed_ptm + occurence.right_pos;
			if(ptm.occurence.length-1 > i )
			{
				fixed_ptm = fixed_ptm + ";" ;
			}
		})
	}
	else
	{
		//fixed_ptm = fixed_ptm + ptm.occurence.left_pos + "-" + ptm.occurence.right_pos ;
    	fixed_ptm = fixed_ptm + occurence.right_pos;
	}
	fixed_ptm = ptm.ptm.abbreviation + fixed_ptm ;
	return fixed_ptm ;
}
/*	Create buttons to save the svg as png/svg and to redraw the svg with given dimensions*/
function buttonsAndAlerts(para,prsm,id)
{
	let x,y;
	/*	Id of the pop_up svg	*/
	id = "l_popup_svg" ;
	/*	On click action to get pop_up window	*/
	d3.select('#saveImage').on("click",function(){
		d3.selectAll(".l_popup_svg_g").remove();
		para = new parameters();
		let prsm = prsm_data.prsm ;
		[para,id] = buildSvg(para,prsm,id);
		/*	Get the amount of skipped acid and write the amount 
		 * 	of skipped acid at the start and end of the sequence 
		 */
		skippedAcidNotification(para,prsm,id) ;
		if(para.show_num)
		{
			/*Get the numerical count at the start enad end of 
			 * each row of sequence */
			getNumValues(para,prsm,id);
		}
		/*	Determine the start and end position of the sequence */
		drawAnnoOfStartEndPosition(para,prsm,id);
		/*	Draw Annotations */
		annotations(para,prsm,id);
		/*	Get the position of the fixed ptms and color them to red */
		addColorToFixedPtms(para,prsm,id);
		/*	Color the background of occurence of mass shift */
		massShiftBackgroundColor(para,prsm,id);
		/*	set the dimensions of popup svg to default values*/
		document.getElementById("row-size").value = para.row_length ;
		document.getElementById("letter-width").value = para.letter_width ;
		document.getElementById("row-height").value = para.row_height ;
		document.getElementById("block-width").value = para.gap_width ;
		document.getElementById("num-width").value = para.numerical_width ;
		document.getElementsByName("show-num")[0].checked = para.show_num ;
		document.getElementsByName("show-skipped-lines")[0].checked = para.show_skipped_lines ;

    $("#myModal").draggable({
		appendTo: "body"
		});
	});

	/*	Download the svg as ".svg" image	*/
	d3.select('#download_SVG').on("click",function(){
			x = d3.event.pageX;
			y = d3.event.pageY;
			popupnamewindow("svg",id,x,y);
		});
	/*	Download svg as PNG Image	*/
	d3.select('#download_PNG').on("click", function(){
		x = d3.event.pageX;
		y = d3.event.pageY;
		popupnamewindow("png",id,x,y);
	})	;

  d3.select('#image_help').on("click",function(){
    $("#helpModal").draggable({
      appendTo: "#myModal"
    });
  });

	/*	On Click action to resize the svg with user dimensions	*/
	d3.select('#resize').on("click", function(){
		d3.selectAll("."+id+"_g").remove();
		/*	Get dimension parameters of svg	*/
		var para = new parameters() ;
		para.row_length = parseInt(document.getElementById("row-size").value);
		para.letter_width = parseInt(document.getElementById("letter-width").value) ;
		para.row_height = parseInt(document.getElementById("row-height").value) ;
		para.gap_width = parseInt(document.getElementById("block-width").value) ;
		para.numerical_width = parseInt(document.getElementById("num-width").value) ;

		/*	Check whether show numbers is checked	*/
		if(document.getElementsByName("show-num")[0].checked)
		{
			para.show_num = true ;
		}
		else
		{
			/*	Reduce the left and right margins when numbers are 
			 * 	no needed at start and end of each row in svg		*/
			para.show_num = false ;
			//para.numerical_width = 0 ;
			para.left_margin = 20;
			para.right_margin = 20;
		}
		/*	Check to show skipped lines */
		if(document.getElementsByName("show-skipped-lines")[0].checked)
		{
			para.show_skipped_lines = true ;
			//para.numerical_width = parseInt(document.getElementById("num-width").value);
		}
		else
		{
			para.show_skipped_lines = false ;
		}
		/*	Redraw the svg with new dimension parameters*/
		prsm = prsm_data.prsm ;
		[para,id] = buildSvg(para,prsm,id);
		/*	Get the amount of skipped acid and write the amount 
		 * 	of skipped acid at the start and end of the sequence 
		 */
		skippedAcidNotification(para,prsm,id) ;
		if(para.show_num)
		{
			/*Get the numerical count at the start and end of 
			 * each row of sequence */
			getNumValues(para,prsm,id);
		}
		/*	Determine the start and end position of the sequence */
		drawAnnoOfStartEndPosition(para,prsm,id) ;
		/*	Draw Annotations */
		annotations(para,prsm,id);
		/*	Get the position of the fixed ptms and color them to red */
		addColorToFixedPtms(para,prsm,id);
		/*	Color the background of occurence of mass shift */
		massShiftBackgroundColor(para,prsm,id);
	});	
}
/* Function to produce a pop up window to provide name and set name to 
*  the image while downloading the image of Graph SVG and Sequence SVG */
function popupnamewindow(type,id,x,y){
	d3.selectAll("#tooltip_imagename").remove() ;
	var div = d3.select("body").append("div")
	.attr("class", "tooltip")
	.attr("id","tooltip_imagename")
	.style("opacity", 1);

	div.transition()
	.duration(200)
	.style("opacity", .9);
	div.html(
			'<input type="text" placeholder="Image Name" id="imagename" />'+
			'<button id="saveimage" style = "none" type="button">save</button>'
			)
	.style("left", (x - 30) + "px")             
	.style("top", (y - 60) + "px")
	// .style("transform","translateX(-35%)!important")
	.attr("box-sizing","border")
	.attr("display","inline-block")
	.attr("min-width","1.5em")
	.attr("padding","2px")
	.attr("margin-left","0px")
	.attr("text-align","center")
	.attr("text-decoration","none")
	.attr("border","1px solid #111111")
	.attr("background-color","white");

	$("#saveimage").click(function(){
		let imagename = $("#imagename").val();
		if( imagename == null || imagename == "")
		{
			imagename = "spectrum";
		}
		if(type == "svg"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let svg_element = d3.selectAll("#"+id).node();
			svg2svg(svg_element,imagename);
		}
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
