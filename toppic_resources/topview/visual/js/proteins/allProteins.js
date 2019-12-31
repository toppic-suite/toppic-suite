/*	Get the data of all the proteins and build html */
function allProteins(folderName)
{
	var l_proteins = prsm_data ;
	var count = 1 ;
	
	/*get the div container*/
	let div = document.getElementsByClassName("container")[0];
	let h2 = document.createElement('h2');
	
	if(Array.isArray(l_proteins.protein_list.proteins.protein)){
		count = l_proteins.protein_list.proteins.protein.length ;
		document.title = count + " proteins are identified" ;
		h2.innerHTML = count + " proteins are identified." ;
	}
	else
	{
		document.title = count + " protein is identified" ;
		h2.innerHTML = count + " protein is identified." ;
	}
	
	let br = document.createElement('br');
	/*create header with protein count*/
	div.appendChild(h2);
	div.appendChild(br);
	/*get the best prsm for each protein and form unique links for all the proteins*/
	if(Array.isArray(l_proteins.protein_list.proteins.protein))
	{
			l_proteins.protein_list.proteins.protein.forEach(function(protein,index){
				let div_temp = proteinToHtml(protein,folderName);
				let p = getBestPrsm(protein,folderName);
				let br1 = document.createElement('br');
				div_temp.appendChild(p);
				div_temp.appendChild(br1);
				div.appendChild(div_temp);
			})
	}
	else
	{
			let protein = l_proteins.protein_list.proteins.protein ;
			let div_temp = proteinToHtml(protein,folderName);
			let p = getBestPrsm(protein,folderName);
			let br1 = document.createElement('br');
			div_temp.appendChild(p);
			div_temp.appendChild(br1);
			div.appendChild(div_temp);
	}
	
}
/*convert the json protein data into HTML and create links for each protein to navigate*/
function proteinToHtml(protein,folderName)
{
	let div  = document.createElement('div');
	let id = "p"+ protein.sequence_id ;
	div.setAttribute("id",id);
	let p = document.createElement('p');
	p.setAttribute("style", "font-size:16px;");
	let a  = document.createElement('a');
	a.href = "protein"+".html" + "?folder="+folderName+"&protein="+ protein.sequence_id;
	a.innerHTML = protein.sequence_name +" "+protein.sequence_description ;
	p.appendChild(a);
	div.appendChild(p);
	return div ;
}
/*	Get the beat PrSM based on the PrSM e value and create link------ 
 * 	to navigate for the best PrSM ------------------------------------*/
function getBestPrsm(protein,folderName)
{
	let best_e_value = 100;
	let prsm_id = "" ;
	let proteoform_count = protein.compatible_proteoform.length;
	if(proteoform_count > 1)
	{
		protein.compatible_proteoform.forEach(function(proteoform,index){
			[best_e_value,prsm_id] = proteoformMultirow(proteoform,best_e_value,prsm_id);
		})
	}
	else
	{
		[best_e_value,prsm_id] = proteoformMultirow(protein.compatible_proteoform,best_e_value,prsm_id);
	}
	let p = document.createElement('p');
	p.setAttribute("style", "font-size:16px;");
	let text1 = document.createElement("text");
	text1.innerHTML = "The ";
	p.appendChild(text1) ;
	let a  = document.createElement('a');
	a.href = "prsm.html?folder="+folderName+"&"+"prsm_id="+prsm_id ;
	a.innerHTML = "best PrSM ";
	p.appendChild(a);
	let text2 = document.createElement("text");
	let val = "has an E-value "+best_e_value+". There";
	if(proteoform_count > 1)
	{
		val = val + " are "+proteoform_count + " proteoforms."
	}
	else
	{
		val = val + " is 1 proteoform."
	}
	text2.innerHTML = val ;
	p.appendChild(text2);
	return p;
}
/*Get the best prsm e value and prsm Id by looping through the proteoform array*/
function proteoformMultirow(proteoform, best_e_value, prsm_id)
{
	let l_best_e_value = best_e_value ;
	let l_prsm_id = prsm_id ;
	if(proteoform.prsm.length > 1)
	{
		proteoform.prsm.forEach(function(prsm,index){
			if(parseFloat(l_best_e_value) < parseFloat(prsm.e_value))
			{
				l_best_e_value = l_best_e_value ;
				l_prsm_id = l_prsm_id ;
			}
			else
			{
				l_best_e_value = prsm.e_value ;
				l_prsm_id = prsm.prsm_id ;
			}
		})
	}
	else
	{
		if(parseFloat(l_best_e_value) < parseFloat(proteoform.prsm.e_value))
		{
			l_best_e_value = l_best_e_value ;
			l_prsm_id = l_prsm_id ;
		}
		else
		{
			l_best_e_value = proteoform.prsm.e_value ;
			l_prsm_id = proteoform.prsm.prsm_id ;
		}
	}
	
	return [l_best_e_value,l_prsm_id] ;
}

