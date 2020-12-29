/**
 * Create title and navigation urls to "all proteins","protein" and "proteoform" of the prsm 
 * @param {String} folderpath - provides folder path to the data and helps in building urls
 */
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
/**
 * Get the data of the prsm from global data variable prsm_data.
 * Build the data into html to show the information about the prsm 
 */
function loadDatafromJson2Html(){
	document.getElementById("File_name").innerHTML = prsm_data.prsm.ms.ms_header.spectrum_file_name;
	document.getElementById("PrSM_ID").innerHTML = prsm_data.prsm.prsm_id;
	document.getElementById("Scan").innerHTML = prsm_data.prsm.ms.ms_header.scans;
	document.getElementById("Precursor_charge").innerHTML = prsm_data.prsm.ms.ms_header.precursor_charge;
	document.getElementById("precursor_mz").innerHTML = prsm_data.prsm.ms.ms_header.precursor_mz ;
	document.getElementById("Precursor_mass").innerHTML = prsm_data.prsm.ms.ms_header.precursor_mono_mass;
	document.getElementById("Proteoform_mass").innerHTML = prsm_data.prsm.annotated_protein.proteoform_mass;
	document.getElementById("matched_peaks").innerHTML = prsm_data.prsm.matched_peak_number;
	document.getElementById("matched_fragment_ions").innerHTML = prsm_data.prsm.matched_fragment_number;
	document.getElementById("unexpected_modifications").innerHTML = prsm_data.prsm.annotated_protein.unexpected_shift_number;
	document.getElementById("E_value").innerHTML = prsm_data.prsm.e_value;
	document.getElementById("Q_value").innerHTML = prsm_data.prsm.fdr;
}
/**
 * Build the monomass table with all the data from the peak variable of the prsm_data
 * Provide a unique class name to the m/z values to provide on click action to xoom the
 * spectrum graph to that position
 */
function createTableElements(){
	var table = document.getElementById('spectrum');
	var tbdy = document.createElement('tbody');
	var ionArray = []; //contains the ion types used in the data
	let l_scans = prsm_data.prsm.ms.ms_header.scans.split(" ") ;
	let l_specIds = prsm_data.prsm.ms.ms_header.ids.split(" ") ;
	let l_matched_peak_count = 0;
	let duplicatePeaks = [];//Id of peaks that were matched by more than 1 ion
	let ion;

	prsm_data.prsm.ms.peaks.peak.forEach(function(peak,i){
		// Check if peak contain matched_ions_num attribute	
		// for each peak, get the ion type and store it in ionArray 
    // to determine which ion type to be checked in Inspect
		if (parseInt(peak.matched_ions_num)>0){
			if (Array.isArray(peak.matched_ions.matched_ion)){
				peak.matched_ions.matched_ion.forEach((ion) => {
					ion = ion.ion_type;
					if (ionArray.length < 1){
						ionArray.push(ion);
					}			
					else if (ionArray.indexOf(ion) < 0) {
						ionArray.push(ion);
					}
				})
			}else{
				ion = peak.matched_ions.matched_ion.ion_type;
				if (ionArray.length < 1){
					ionArray.push(ion);
				}			
				else if (ionArray.indexOf(ion) < 0) {
					ionArray.push(ion);
				}
			}
		}
		
		if(peak.hasOwnProperty('matched_ions_num') && parseInt(peak.matched_ions_num)>1)
		{
			let origMatchedIons = [];
			//because matched ion array gets overwritten in the code below
			//which leads to missing ion information in mass graph
			//so store the original value in this array and re-overwrite the ion information
			//after writing to the table is finished
			
			for (let i = 0; i < peak.matched_ions.matched_ion.length; i++){
				origMatchedIons.push(peak.matched_ions.matched_ion[i]);
			}

			peak.matched_ions.matched_ion.forEach((matched_ion,i) => {
				peak.matched_ions.matched_ion = matched_ion ;
				loop_matched_ions(peak,i) ;
			})
			peak.matched_ions.matched_ion = origMatchedIons;
		}
		else
		{
			loop_matched_ions(peak,i) ;
		}
	})

	//after looping through the prsm files, store the ion type data to local storage
	window.localStorage.setItem('ionType', ionArray);

	/**
	 * Inner function to create a rows and columns for monomass table
	 * @param {object} peak - contains information of each peak 
	 * @param {int} i - index of the peak
	 */
	function loop_matched_ions(peak,i){
		/*	Create row for each peak value object in the table	*/
		var tr = document.createElement('tr');
		id = peak.spec_id+"peak"+peak.peak_id;
		let l_scan;
		if((parseInt(peak.peak_id) + 1)%2 == 0)
		{
			// class name helps to get unmatched peaks when clicking unmatched peaks
			l_class = "unmatched_peak even"; 
		}
		else
		{
			// class name helps to get unmatched peaks when clicking unmatched peaks
			l_class = "unmatched_peak odd"; 
		}
		if(peak.hasOwnProperty('matched_ions_num'))
		{
			id = id + peak.matched_ions.matched_ion.ion_type;
			if((parseInt(peak.peak_id) + 1)%2 == 0)
			{
				// class name helps to get matched peaks when clicking matched peaks
				l_class = "matched_peak even";
			}
			else
			{
				// class name helps to get matched peaks when clicking matched peaks
				l_class = "matched_peak odd";
			}
			l_matched_peak_count++;
			//	create a name for each row
			tr.setAttribute("name",peak.matched_ions.matched_ion.ion_position);

			if (parseInt(peak.matched_ions_num) > 1){
				let peakId = peak.peak_id;
				if (duplicatePeaks.indexOf(peakId) < 0){
					duplicatePeaks.push(peakId);
				}
			}
		}
		//	Set "id","class name" and "role" for each row
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
				//	provide link to click on m/z value to view spectrum 
				let a = document.createElement('a');
				a.href="#!"
				a.className = "row_mono_mz"
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
					//if c-term ion, pos = pos + 1
					let ion = peak.matched_ions.matched_ion.ion_type;
					if (ion.indexOf("X") >= 0 || ion.indexOf("Y") >= 0 || ion.indexOf("Z") >= 0)
					{
						td.innerHTML = parseInt(peak.matched_ions.matched_ion.ion_position) + 1;
					}
					else
					{
						td.innerHTML = peak.matched_ions.matched_ion.ion_position;
					}
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
	let l_duc_peak_count = duplicatePeaks.length;
	l_matched_peak_count = l_matched_peak_count - l_duc_peak_count;
	let l_not_matched_peak_count = l_All_Peaks - l_matched_peak_count; 
	
	document.getElementById("all_peak_count").innerHTML = "All peaks (" + l_All_Peaks + ")" ;
	document.getElementById("matched_peak_count").innerHTML = "Matched peaks (" + l_matched_peak_count + ")" ;
	document.getElementById("not_matched_peak_count").innerHTML = "Not Matched peaks (" + l_not_matched_peak_count + ")" ;
	
	table.appendChild(tbdy);
}
/**
 * Get occurence of "Variable" and "Fixed", convert the data to HTML
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function occurence_ptm(prsm)
{
	let variable_ptm = "";
	let fixed_ptm = "" ;
	// Check if annotation has ptm attribute inside it
	if(prsm.annotated_protein.annotation.hasOwnProperty('ptm'))
	{
		let ptm = prsm.annotated_protein.annotation.ptm ;
		// Check if ptm attribute is an array
		if(Array.isArray(ptm))
		{
			ptm.forEach(function(ptm, index){
				// Check if there exist variable or fixed ptms
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
			// Check if there exist variable or fixed ptms
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
	// Add the information of fixed ptms to html at id - ptm_abbreviation
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
	// Add the information of varibale ptms to html at id - ptm_abbreviation
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
/**
 * Get the information about all the unknown ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function getUnknownPtms(prsm)
{
	let data = " [";
	// Check if annotation attribute has mass_shift attribute inside it
	if(prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
	{
		let mass_shift = prsm.annotated_protein.annotation.mass_shift ;
		let UnknownExist =  false;
			if(Array.isArray(mass_shift))
			{
				let len = mass_shift.length;
				mass_shift.forEach(function(each_mass_shift, i){
					// Check for unexpected mass shifts
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
			// If unexpected modifications exist add them to html at id - ptm_unexpectedmodifications
			if(UnknownExist)
			{
				let val = "Unknown" + data;
				document.getElementById("ptm_unexpectedmodification").style.display = "block";
				let div = document.getElementById("ptm_unexpectedmodification") ;
				let text1 = document.createElement("text");
				let text2 = document.createElement("text");
				text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"+"Unexpected modifications: "; // Adding space using &nbsp;
				text2.innerHTML = val ;
				text2.target = "_blank";
				text2.style = "color:red" ;
				div.appendChild(text1);
				div.appendChild(text2);
			}
	}
	
}
/**
 * Get all the variable ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
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
/**
 * Get all the "Fixed" ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function getFixedPtm(ptm)
{
	let fixed_ptm = "[" ;
	let abbrevation = ptm.ptm.abbreviation ;
	if(Array.isArray(ptm.occurence))
	{
		ptm.occurence.forEach(function(occurence,i){
      	fixed_ptm = fixed_ptm + occurence.right_pos;
			if(ptm.occurence.length-1 > i )
			{
				fixed_ptm = fixed_ptm + ";" ;
			}
		})
	}
	else
	{
    	fixed_ptm = fixed_ptm + ptm.occurence.right_pos;
	}
	fixed_ptm = ptm.ptm.abbreviation + fixed_ptm ;
	return fixed_ptm ;
}


