var peakCalData = function(prsm,massShift,shiftPosition){
	this.massList = [] ;
	this.monoMassList = [];
	this.thresholdVal = 1 ;
	prsm.ms.peaks.peak.forEach(function(peak,i){
		let tempPeak = {};
		tempPeak.peak_id = peak.peak_id ;
		tempPeak.monoisotopic_mass = peak.monoisotopic_mass;
		this.monoMassList.push(tempPeak);
	})
	let residueList_len = prsm.annotated_protein.annotation.residue.length;
	for(i = 0; i<residueList_len ; i++)
	{
		this.massList[i] = 0 ;
	}
	console.log("massList : ", massList);
	this.thMassList = getTheoreticalMass(prsm,massList,massShift,shiftPosition);
	this.prefixMassList = getPrefixMass(prsm,massList,massShift,shiftPosition);
	this.suffixMassList = getSuffixMass(prsm,massList,massShift,shiftPosition);
	this.completeCalData = CalmassError(this.thMassList,this.monoMassList,this.thresholdVal);
	
	console.log("thMassList : ",this.thMassList);
	console.log("prefixMassList : ", prefixMassList);
	console.log("suffixMassList : ",suffixMassList);
	console.log("completeCalData : ",this.completeCalData);
} 
getTheoreticalMass = function(prsm,massList,massShift,shiftPosition)
{
	let thMassList = [];
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	let zero = 0;
	residue_list = prsm.annotated_protein.annotation.residue ;
	
	prsm.annotated_protein.annotation.residue.forEach(function(residue,i){
		if(parseInt(residue.position) >= parseInt(first_residue_position) && parseInt(residue.position) <= parseInt(last_residue_position))
		{
			let temp_thMass = {};
			temp_thMass.position = parseInt(residue.position);
			let curPosition = residue.position ;
			temp_thMass.mass = getThMass(prsm,curPosition,massList,massShift,shiftPosition) ;
			thMassList.push(temp_thMass) ;
		}
	})
	return thMassList ;
}
getThMass = function(prsm,curPosition,massList,massShift,shiftPosition){
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	residue_list = prsm.annotated_protein.annotation.residue ;
	let annotateMassShift = json2BackgroundColorArray(prsm);
	let mass = 0;
	for(let i = first_residue_position ; i <= curPosition ; i++ )
	{
		let acidMass = parseFloat(getAcidMass(residue_list[i].acid).mono);
		mass = mass + acidMass + massList[i] ;
	}
	if(curPosition >= shiftPosition)
	{
		mass = mass + massShift ;
	}
	for(let j=0;j<annotateMassShift.length;j++)
	{
		if(annotateMassShift[j].left_position <= curPosition)
		{
			mass = mass + parseFloat(annotateMassShift[j].anno);
		}
	}
	return mass ;
}
getPrefixMass = function(prsm,massList,massShift,shiftPosition)
{
	let prefixMassList = [];
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	let residue_list = prsm.annotated_protein.annotation.residue ;
	prsm.annotated_protein.annotation.residue.forEach(function(residue,i){
		if(parseInt(residue.position) >= parseInt(first_residue_position) && parseInt(residue.position) <= parseInt(last_residue_position))
		{
			let temp_prefixMass = {};
			temp_prefixMass.position = parseInt(residue.position);
			let curPosition = residue.position ;
			temp_prefixMass.mass = getPrMass(prsm,curPosition,massList,massShift,shiftPosition) ;
			prefixMassList.push(temp_prefixMass) ;
		}
	})
	return prefixMassList ;
}
getPrMass = function(prsm,curPosition,massList,massShift,shiftPosition){
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	residue_list = prsm.annotated_protein.annotation.residue ;
	let annotateMassShift = json2BackgroundColorArray(prsm);
	let mass = 0;
	for(let i = first_residue_position ; i <= curPosition ; i++ )
	{
		let acidMass = parseFloat(getAcidMass(residue_list[i].acid).mono);
		mass = mass + acidMass + massList[i] 
	}
	if(curPosition >= shiftPosition)
	{
		mass = mass + massShift ;
	}
	for(let j=0;j<annotateMassShift.length;j++)
	{
		if(annotateMassShift[j].left_position <= curPosition)
		{
			mass = mass + parseFloat(annotateMassShift[j].anno);
		}
	}
	return mass ;
}
/*massShift - amount of mass is shifted, shiftPosition - start of the shift position*/ 
getSuffixMass = function(prsm,massList,massShift,shiftPosition)
{
	let suffixMassList = [] ;
	let residueList =  prsm.annotated_protein.annotation.residue ;
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	i = last_residue_position + 1 ;
	while(i--)
	{
		if(i >=first_residue_position && i <= last_residue_position)
		{
			let temp_SuffixMass = {}
			temp_SuffixMass.position = i;
			let curPosition = i
			temp_SuffixMass.mass = getSufMass(prsm,curPosition,massList,massShift,shiftPosition);
			suffixMassList.push(temp_SuffixMass);
		}
	}
	return suffixMassList ;
}
getSufMass = function(prsm,curPosition,massList,massShift,shiftPosition)
{
	let first_residue_position = parseInt(prsm.annotated_protein.annotation.first_residue_position) ;
	let last_residue_position = parseInt(prsm.annotated_protein.annotation.last_residue_position) ;
	let mass = 0 ;
	let annotateMassShift = json2BackgroundColorArray(prsm);
	let residue_list = prsm.annotated_protein.annotation.residue ;
	for(let i = last_residue_position; i >= curPosition ;i--)
	{
		let acidMass = parseFloat(getAcidMass(residue_list[i].acid).mono);
		mass = mass + acidMass + massList[i] 
	}
	if(shiftPosition >= curPosition)
	{
		mass = mass + massShift ;
	}
	for(let j=0;j<annotateMassShift.length;j++)
	{
		if((annotateMassShift[j].right_position - 1) >= curPosition)
		{
			mass = mass + parseFloat(annotateMassShift[j].anno);
		}
	}
	return mass ;
}

function CalmassError(thMassList,monoMassList,thresholdVal){
	let calculatedDataList = [] ;
	for(let i = 0; i < monoMassList.length; i++)
	{
		for(let j=0; j<thMassList.length ; j++)
		{
			if((monoMassList[i].monoisotopic_mass - thMassList[j].mass) < thresholdVal && (monoMassList[i].monoisotopic_mass - thMassList[j].mass) > 0-thresholdVal)
			{
				let calculatedData = {} ;
				/*add 1 as the peak id start from one where as othe positions start from 0*/
				calculatedData.peak_rid = parseInt(monoMassList[i].peak_id) + 1 ;
				calculatedData.massEror = monoMassList[i].monoisotopic_mass - thMassList[j].mass;
				calculatedData.thMass = thMassList[j].mass ;
				calculatedData.position = thMassList[j].position ;
				calculatedData.monoisotopic_mass = monoMassList[i].monoisotopic_mass ;
				let PPMerror = getPPMerror(calculatedData.massError,thMassList[j].mass)
				calculatedData.PPMerror = PPMerror ;
				calculatedDataList.push(calculatedData);
			}
		}
	}
	return calculatedDataList ;
}
function getPPMerror(error,thMass)
{
	let million = 1000000;
	let PPMerror = error/thMass * million ;
	return PPMerror ;
}
