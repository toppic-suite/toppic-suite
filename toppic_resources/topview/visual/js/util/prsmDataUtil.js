/*	Get the cleavage positions from the prsm data	*/
function json2CleavagePositions(prsm)
{
	let matched_ion = [] ;
	let position ;
	let exist_n_ion;
	let exist_c_ion;
	prsm.annotated_protein.annotation.cleavage.forEach(function(cleavage,i){
		position = cleavage.position ;
		exist_n_ion = cleavage.exist_n_ion ;
		exist_c_ion = cleavage.exist_c_ion ;
		if(cleavage.matched_peaks != null)
		{
			if(cleavage.matched_peaks.matched_peak.length > 1)
			{
				cleavage.matched_peaks.matched_peak.forEach(function(matched_peak,i){
					inner_matched_peak(matched_peak) ;
				})
			}
			else
			{
				inner_matched_peak(cleavage.matched_peaks.matched_peak);
			}
		}
	});
	/*	Local function to get cleavage "position","ion type","ion display position","peak charge","ion position" */
	function inner_matched_peak(matched_peak)
	{
		let matched_ion_temp_array = {} ;
		
		matched_ion_temp_array.position = position ;
		matched_ion_temp_array.exist_n_ion = exist_n_ion ;
		matched_ion_temp_array.exist_c_ion = exist_c_ion ;
		// Ion Type => "Y"/"B"
		matched_ion_temp_array.ion_type = matched_peak.ion_type ;
		// Ion Display position
		matched_ion_temp_array.ion_display_position = matched_peak.ion_display_position ;
		// Ion Charge
		matched_ion_temp_array.peak_charge = matched_peak.peak_charge ;
		// ion_position
		matched_ion_temp_array.ion_position = matched_peak.ion_position ;
		matched_ion.push(matched_ion_temp_array) ;
	}
	return matched_ion ;
}
/*	Get occurence of fixed ptm positions*/
function json2FixedPtmOccurence(prsm){
	let occurence_list = [] ;
	if(prsm.annotated_protein.annotation.hasOwnProperty("ptm") )
	{
		if(Array.isArray(prsm.annotated_protein.annotation.ptm))
		{
			prsm.annotated_protein.annotation.ptm.forEach(function(ptm,index){
				if(ptm.ptm_type == "Fixed")
				{
					if(ptm.hasOwnProperty("occurence"))
					{
						if(Array.isArray(ptm.occurence))
						{
							ptm.occurence.forEach(function(occurence,i){
								occurence_list.push(occurence.left_pos);
							});
						}
						else
						{
							occurence_list.push(ptm.occurence.left_pos);
						}
					}
				}
			})
		}
		else
		{
			if(prsm.annotated_protein.annotation.ptm.hasOwnProperty("occurence"))
			{
				if(prsm.annotated_protein.annotation.ptm.ptm_type == "Fixed")
				{
					if(Array.isArray(prsm.annotated_protein.annotation.ptm.occurence))
					{
						prsm.annotated_protein.annotation.ptm.occurence.forEach(function(occurence,i){
							occurence_list.push(occurence.left_pos);
						});
					}
					else
					{
						occurence_list.push(prsm.annotated_protein.annotation.ptm.occurence.left_pos);
					}
				}
			}
		}
	}

	return occurence_list ;
}
/*	Get left and right positions of background color and mass shift value */
function json2BackgroundColorArray(prsm)
{
	let backgroundColorAndMassShift = [];
	if(prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
	{
		if(Array.isArray(prsm.annotated_protein.annotation.mass_shift)){
			prsm.annotated_protein.annotation.mass_shift.forEach(function(mass_shift,i){
				
				if(mass_shift.right_position != "0")
				{
					backgroundColorAndMassShift.push(mass_shift) ;
				}
			})
		}
		else
		{
			let mass_shift = prsm.annotated_protein.annotation.mass_shift ;
			
			if(mass_shift.right_position != "0")
			{
				backgroundColorAndMassShift.push(mass_shift) ;
			}
			
		}
	}
	if(prsm.annotated_protein.annotation.hasOwnProperty('ptm'))
	{
		let otherPtmList = json2OtherPtmOccurences(prsm);
		backgroundColorAndMassShift = backgroundColorAndMassShift.concat(otherPtmList);
	}
	return backgroundColorAndMassShift ;
}
/*	Get position and other ptm lists other than FIxed Ptms */
function json2OtherPtmOccurences(prsm)
{
	let backgroundColorAndMassShift = [];
	if(Array.isArray(prsm.annotated_protein.annotation.ptm))
	{
		prsm.annotated_protein.annotation.ptm.forEach(function(ptm,index){
			if(ptm.ptm_type != "Fixed")
			{
				if(ptm.hasOwnProperty("occurence"))
				{
					if(Array.isArray(ptm.occurence))
					{
						ptm.occurence.forEach(function(occurence,i){
							let tempObj = {};
							tempObj.anno = ptm.ptm.abbreviation;
							tempObj.left_position = occurence.left_pos;
							tempObj.right_position = occurence.right_pos;
							backgroundColorAndMassShift.push(tempObj);
						});
					}
					else
					{
						let tempObj = {};
						tempObj.anno = ptm.ptm.abbreviation;
						tempObj.left_position = ptm.occurence.left_pos;
						tempObj.right_position = ptm.occurence.right_pos;
						backgroundColorAndMassShift.push(tempObj);
					}
				}
			}
			
		})
	}
	else
	{
		if(prsm.annotated_protein.annotation.ptm.hasOwnProperty("occurence"))
		{
			let ptm = prsm.annotated_protein.annotation.ptm;
			if(ptm.ptm_type != "Fixed")
			{
				if(Array.isArray(prsm.annotated_protein.annotation.ptm.occurence))
				{
					prsm.annotated_protein.annotation.ptm.occurence.forEach(function(occurence,i){
						let tempObj = {};
						tempObj.anno = ptm.ptm.abbreviation;
						tempObj.left_position = occurence.left_pos;
						tempObj.right_position = occurence.right_pos;
						backgroundColorAndMassShift.push(tempObj);
					});
				}
				else
				{
					let tempObj = {};
					tempObj.anno = ptm.ptm.abbreviation;
					tempObj.left_position = ptm.occurence.left_pos;
					tempObj.right_position = ptm.occurence.right_pos;
					backgroundColorAndMassShift.push(tempObj);
				}
			}
		}
	}
	return backgroundColorAndMassShift;
}
function json2ErrorDataList(prsm){
	let errorDataList = [];
	prsm.ms.peaks.peak.forEach((peak) => {
		if(peak.hasOwnProperty('matched_ions_num'))
		{
			errorDataList.push(peak.matched_ions.matched_ion);
		}
	})
	return errorDataList;
}