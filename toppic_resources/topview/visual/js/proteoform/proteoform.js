/**
 * Create naviagtion urls to go back to protein page and all protein page
 * @param {String} folderpath - Path to data folder 
 */
function proteoformUrl(folderpath)
{
	// prsm_data is a global variable containing complete data of the proteoform from the data file
	l_proteoform_Url = prsm_data.compatible_proteoform.sequence_name+" " +prsm_data.compatible_proteoform.sequence_description;
	// Set title for the page
	document.title = "Proteoform #"+prsm_data.compatible_proteoform.proteoform_id+" from " + l_proteoform_Url ;	
	let l_allproteins_url = "proteins.html?folder="+folderpath;
	// Set the naviagtion urls to all proteins page and protein page
	document.getElementById("allprotein_url_start").href = l_allproteins_url;
	document.getElementById("allprotein_url_end").href = l_allproteins_url;
	document.getElementById("protein_url_start").innerHTML = l_proteoform_Url;
	document.getElementById("protein_url_start").href = "protein.html?folder="+folderpath+"&protein_Id="+prsm_data.compatible_proteoform.sequence_id;
	document.getElementById("protein_url_end").innerHTML = l_proteoform_Url;
	document.getElementById("protein_url_end").href = "protein.html?folder="+folderpath+"&protein_Id="+prsm_data.compatible_proteoform.sequence_id;
	document.getElementById("proteoform_header").innerHTML ="Proteoform #"+prsm_data.compatible_proteoform.proteoform_id;
	// Get the count of number of proteoform 
	if(Array.isArray(prsm_data.compatible_proteoform.prsm))
	{
		document.getElementById("prsm_count").innerHTML = prsm_data.compatible_proteoform.prsm.length+" PrSMs are identified for the proteoform";
	}
	else
	{
		document.getElementById("prsm_count").innerHTML = "1 PrSM is identified for the proteoform";
	}
}
/**
 * Create a table with prsm data and URL to navigate to appropriate prsm
 * @param {String} folderpath - Provides path to the data folder
 */
function createTableData(folderpath){
	let table = document.getElementById('proteoform_data');
	let tbdy = document.createElement('tbody');
	let count = 0;
	let sequence_name = prsm_data.compatible_proteoform.sequence_name ;
	// Iterate through the number of prsms
	if(Array.isArray(prsm_data.compatible_proteoform.prsm))
	{
		prsm_data.compatible_proteoform.prsm.forEach(function(prsm,i){
			let All_Peak_count = prsm.ms.peaks.peak.length;
			let tr = document.createElement('tr');
			for(let i=0;i<7;i++){
				var td = document.createElement('td');
				td.setAttribute("align","center");
				if(i === 0){
					td.innerHTML = prsm.ms.ms_header.scans ;
				}
				if(i === 1){
					td.innerHTML = sequence_name ;
				}
				if(i === 2){
					td.innerHTML = prsm.e_value ;
				}
				if(i === 3){
					td.innerHTML = All_Peak_count ;
				}
				if(i === 4){
					td.innerHTML = prsm.matched_peak_number;
				}
				if(i === 5){
					td.innerHTML = prsm.matched_fragment_number;
				}
				if(i === 6){
					// Create URL to navigate
					let a = document.createElement('a')
					l_href = "prsm.html?folder="+folderpath+"&prsm_id="+prsm.prsm_id;
					l_link = "link" + (count +1) ;
					a.setAttribute("href", l_href);
					a.setAttribute("id", l_link);
					a.innerHTML = "See PrSM&gt;&gt;"; // Adding >> marks
					td.appendChild(a);
					count ++ ;
				}
				tr.appendChild(td);
			}
			tbdy.appendChild(tr);
		})
	}
	else
	{
		prsm = prsm_data.compatible_proteoform.prsm ;

		let All_Peak_count = prsm.ms.peaks.peak.length;
		let tr = document.createElement('tr');
		for(let i=0;i<7;i++){
			var td = document.createElement('td');
			td.setAttribute("align","center");
			if(i === 0){
				td.innerHTML = prsm.ms.ms_header.scans ;
			}
			if(i === 1){
				td.innerHTML = sequence_name ;
			}
			if(i === 2){
				td.innerHTML = prsm.e_value ;
			}
			if(i === 3){
				td.innerHTML = All_Peak_count ;
			}
			if(i === 4){
				td.innerHTML = prsm.matched_peak_number;
			}
			if(i === 5){
				td.innerHTML = prsm.matched_fragment_number;
			}
			if(i === 6){
				// Create link to navigate
				let a = document.createElement('a')
				l_href = "prsm.html?prsm_id="+prsm.prsm_id;
				l_link = "link" + (count +1) ;
				a.setAttribute("href", l_href);
				a.setAttribute("id", l_link);
				a.innerHTML = "See PrSM&gt;&gt;";
				td.appendChild(a);
				count ++ ;
			}
			tr.appendChild(td);
		}
		tbdy.appendChild(tr);
	}
	table.appendChild(tbdy);
}