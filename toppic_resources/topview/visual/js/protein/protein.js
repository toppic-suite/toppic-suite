/*	Build Title,Description,URL and get best prsm for each proteoform */
function protein(folderpath){
	document.title = "Proteoforms for protein " + prsm_data.protein.sequence_name + " "+ prsm_data.protein.sequence_description;
	document.getElementById('sequence_description').innerHTML = prsm_data.protein.compatible_proteoform_number+" proteoforms for protein "
																+ prsm_data.protein.sequence_name + " "+ prsm_data.protein.sequence_description; 
	if(Array.isArray(prsm_data.protein.compatible_proteoform))
	{
		prsm_data.protein.compatible_proteoform.forEach(function(compatible_proteoform,index){
			proteoformToHtml(compatible_proteoform,index,folderpath);
		})
	}
	else
	{
		proteoformToHtml(prsm_data.protein.compatible_proteoform,0,folderpath);
	}
}
/*	Get best prsm and build HTML tags */
function proteoformToHtml(compatible_proteoform,index,folderpath)
{
	var div_container = document.getElementsByClassName("proteoformcontainer")[0];
	var div = document.createElement('div');
	let id = "p"+ compatible_proteoform.proteoform_id;
	div.setAttribute("id",id);
	var h2 = document.createElement('h3');
	let p = document.createElement("p");
	let e_value;
	let BestPrSM ;
	if(compatible_proteoform.prsm.length > 0)
	{
		h2.innerHTML = "Proteoform #" + compatible_proteoform.proteoform_id + " Feature intensity: "
						+ compatible_proteoform.prsm[0].ms.ms_header.feature_inte;
		let precursor_mass;
		let prsm_id ;
		[e_value,precursor_mass,prsm_id] = getBestPrsm(compatible_proteoform.prsm);
		p = Build_BestPrSM(e_value,precursor_mass,prsm_id,compatible_proteoform.proteoform_id, compatible_proteoform.prsm.length,folderpath);
		for(let i = 0; i< compatible_proteoform.prsm.length ; i++)
		{
			if(prsm_id == compatible_proteoform.prsm[i].prsm_id)
			{
				BestPrSM = compatible_proteoform.prsm[i] ;
				break;
			}
		}
	}
	else
	{
		h2.innerHTML = "Proteoform #" + compatible_proteoform.proteoform_id + " Feature intensity: "
										+ compatible_proteoform.prsm.ms.ms_header.feature_inte;
		p.setAttribute("style","font-size:16px;");
		let text1 = document.createElement("text");
		text1.innerHTML = "There is only ";
		p.appendChild(text1);
		let a_prsm = document.createElement("a");
		a_prsm.href = "prsm.html?folder="+folderpath+"&prsm_id="+compatible_proteoform.prsm.prsm_id;
		a_prsm.innerHTML = " 1 PrSM ";
		p.appendChild(a_prsm);
		let text2 = document.createElement("text");
		text2.innerHTML = "with an E-value "+compatible_proteoform.prsm.e_value +" and a precursor mass "+
							compatible_proteoform.prsm.ms.ms_header.precursor_mono_mass +".";
		p.appendChild(text2);
		e_value = compatible_proteoform.prsm.e_value ;
		BestPrSM = compatible_proteoform.prsm ;
	}
	div_container.appendChild(div);
	div_container.appendChild(h2);
	div_container.appendChild(p);
	
	let Svg_id = "svgContainer" + index;
	let svgContainer = document.createElement("div");
	svgContainer.setAttribute("id",Svg_id);
	svgContainer.setAttribute("class","svgContainer");
	div_container.appendChild(svgContainer);
	id = "l_svg"+ index;
	l_svgContainer = d3.select("#"+Svg_id).append("svg")
									.attr("id",id);
	let para  = new parameters();
	[para,id] = buildSvg(para,BestPrSM,id);
	/*	Get the amount of skipped acid and write the amount 
	 * 	of skipped acid at the start and end of the sequence 
	 */
	skippedAcidNotification(para,BestPrSM,id) ;
	if(para.show_num)
	{
		/*Get the numerical count at the start enad end of 
		 * each row of sequence */
		getNumValues(para,BestPrSM,id);
	}
	/*	Determine the start and end position of the sequence */
	drawAnnoOfStartEndPosition(para,BestPrSM,id) ;
	/*	Get the position of the fixed ptms and color them to red */
	addColorToFixedPtms(para,BestPrSM,id);
	/*	Color the background of occurence of mass shift */
	massShiftBackgroundColor(para,BestPrSM,id);
}
/*	Get "precursor mass","prsm Id" and least "e value" for each proteoform  */
function getBestPrsm(prsm)
{
	let e_value = " " ;
	let precursor_mass = " " ;
	let prsm_id = "";
	let temp = parseFloat(prsm[0].e_value);
	e_value = prsm[0].e_value;
	precursor_mass = prsm[0].ms.ms_header.precursor_mono_mass;
	prsm_id = prsm[0].prsm_id;
	for(let i = 1 ; i < (prsm.length) ; i++)
	{
		if(temp >= parseFloat(prsm[i].e_value))
		{
			temp = parseFloat(prsm[i].e_value)
			e_value = prsm[i].e_value;
			precursor_mass = prsm[i].ms.ms_header.precursor_mono_mass;
			prsm_id = prsm[i].prsm_id;
		}
	}
	return [e_value,precursor_mass,prsm_id];
}

/*	Create HTML URL link to navigate to best prsm and to navigae to proteoform page */
function Build_BestPrSM(e_value,precursor_mass,prsm_id,proteoform_id, PrSM_Count,folderpath )
{
	let p = document.createElement("p");
	p.setAttribute("style","font-size:16px;");
	let text1 = document.createElement("text");
	text1.innerHTML = "The ";
	p.appendChild(text1);
	let a_prsm = document.createElement("a");
	a_prsm.href = "prsm.html?folder="+folderpath+"&prsm_id="+prsm_id;
	a_prsm.innerHTML = " best PrSM ";
	p.appendChild(a_prsm);
	let text2 = document.createElement("text");
	text2.innerHTML = "has an E-value "+e_value+" and a precursor mass "+precursor_mass+". There are ";
	p.appendChild(text2);
	let a_proteoform = document.createElement("a");
	a_proteoform.href = "proteoform.html?folder="+folderpath+"&proteoform_id=" + proteoform_id;
	a_proteoform.innerHTML = PrSM_Count + " PrSMs ";
	p.appendChild(a_proteoform);
	let text3 = document.createElement("text");
	text3.innerHTML = "in total.";
	p.appendChild(text3);
	return p ;
}

