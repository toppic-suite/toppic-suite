class CalculatePrefixAndSuffixMass{
	constructor(){
		this.fixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},
								{name:"Carboxymethyl",acid:"C",mass:58.005479}];
	}

	getIonTypeMass(ionType){
		let ionTypeMassList={
			"A":-27.9949,
			"A-H2O":-46.0149,
			"A-NH3":-45.02542,
			"B":0,
			"B-H2O":-18.02,
			"B-NH3":-17.03052,
			"C":17.0266,
			"C-H2O":-0.9934,
			"C-NH3":-0.00392,
			"X":43.99,
			"X-H2O":25.97,
			"X-NH3":26.95948,
			"Y":18.0106,
			"Y-H2O":0,
			"Y-NH3":0.98008,
			"Z":0.984,
			"Z-H2O":-17.026,
			"Z-NH3":-16.04652,
			"Z0":1.9919,
			"Z0-H2O":-16.018664,
			"Z0-NH3":-15.03862
		};
		return ionTypeMassList[ionType.toUpperCase()];
	}
	getSequence(prsm_data){
		let sequence = [];
		let firstposition = prsm_data.prsm.annotated_protein.annotation.first_residue_position;
		let lastposition = prsm_data.prsm.annotated_protein.annotation.last_residue_position;
		prsm_data.prsm.annotated_protein.annotation.residue.forEach(function(eachrow,i){
			if(parseInt(eachrow.position) >= parseInt(firstposition)&&
				parseInt(eachrow.position) <= parseInt(lastposition))
			{
				sequence = sequence+eachrow.acid;
			}
		})
	   return sequence;
	}
	getUnknownMassList()
	{
		let unknownMassShiftList = [];
		let l_prsm = prsm_data;
		if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
		{
			let mass_shift = l_prsm.prsm.annotated_protein.annotation.mass_shift ;
				if(Array.isArray(mass_shift))
				{
					let len = mass_shift.length;
					mass_shift.forEach(function(each_mass_shift, i){
						// Removing -1 as the sequece in inspect elements takes from 0
						let position = parseInt(each_mass_shift.left_position) ;
						let mass = parseFloat(each_mass_shift.anno);
						unknownMassShiftList.push({"position":position,"mass":mass})
					})
				}
				else if(mass_shift.shift_type == "unexpected")
				{
					// Removing -1 as the sequece in inspect elements takes from 0
					let position = parseInt(mass_shift.left_position);
					let mass = parseFloat(mass_shift.anno);
					unknownMassShiftList.push({"position":position,"mass":mass})
				}
		}
		return unknownMassShiftList;
	}
	getPrefixMassList(sequence,massShiftList,massShift_in){
		let seqln = sequence.length;
		let emptyMassList = [] ;
		massShiftList = massShiftList.slice();
		massShiftList = this.getFixedPTMMassList(massShiftList,prsm_data.prsm);
		if(seqln != 0)
		{
			let prefixMassList = new Array(seqln);
			let shiftListln = massShiftList.length;
			for(let i=0;i<seqln;i++)
			{
				// Get Calculated AminoAcidDistribution from aminoAcidDistributuion.js 
				//evaluate the return value for the case when the user enters wrong character by mistake
				let AcidMass;
				
				try{
					AcidMass = getAminoAcidDistribution(sequence[i])[0].mass;

					if(i == 0)
					{
						// Add 1 to i as the seq start from 0 but the peak id starts from 1
						if(shiftListln != 0)
						{
							AcidMass = this.addMassShift(i,massShiftList,AcidMass);
						}
						let tempObj = {acid:sequence[i], position:(i+1),mass:AcidMass} ;
						prefixMassList[i] = tempObj;
					}
					else
					{
						let mass = prefixMassList[i-1].mass + AcidMass;
						if(shiftListln != 0)
						{
							mass = this.addMassShift(i,massShiftList,mass);
						}
						// Add 1 to i as the seq start from 0 but the peak id starts from 1
						let tempObj = {acid:sequence[i],position:(i+1),mass:mass};
						prefixMassList[i] = tempObj;
					}
				}
				catch(error){
					//invalid character entered for protein sequence
					//window.alert("Error! Invalid amino acid in the sequence.")
					$('#errorAlertModal').modal('show');
					break;
				}
			}
			// Adding Mass Shift based on selection of Ion
			// let utilFunctionsObj = new utilFunctions();
			// let ionMassShiftVal = utilFunctionsObj.getNTerminusMassShiftVal();
			if(massShift_in != 0)
			{
				for(let j=0;j<seqln;j++)
				{
					prefixMassList[j].mass = prefixMassList[j].mass + massShift_in;
				}
			}
			return prefixMassList;
		}
		return emptyMassList;
	}
	getSuffixMassList(sequence,massShiftList,massShift_in){
		let seqln = sequence.length;
		let emptyMassList = [];
		// As the elements starts from 0
		if(seqln != 0)
		{
			let temp_seqln = seqln-1;
			massShiftList = this.getFixedPTMMassList(massShiftList,prsm_data.prsm);;
			let shiftListln = massShiftList.length;
			let suffixMassList = new Array(seqln);
			for(let i=temp_seqln ; i >= 0 ; i--)
			{
				// Get Calculated AminoAcidDistribution from aminoAcidDistributuion.js
				let AcidMass = getAminoAcidDistribution(sequence[i])[0].mass;
				if(i == temp_seqln)
				{
					// Add 1 to i as the seq start from 0 but the peak id starts from 1
					let position = temp_seqln - i;
					// Adding water mass for suffix, indirectly this will add water mass to all the masses
					let mass = AcidMass;//
					if(suffixMassList != 0)
					{
						mass = this.addMassShift(i,massShiftList,mass)
					}
					let tempObj = {acid:sequence[i],position:(i+1),mass:mass};
					//let tempObj = {position:(position+1),mass:mass} ;
					suffixMassList[position] = tempObj;
				}
				else
				{
					//Don't add water here
					let position = temp_seqln - i;
					let mass = suffixMassList[(position-1)].mass + AcidMass;
					if(suffixMassList != 0)
					{
						mass = this.addMassShift(i,massShiftList,mass)
					}
					// Add 1 to i as the seq start from 0 but the peak id starts from 1
					let tempObj = {acid:sequence[i],position:(i+1),mass:mass};
					//let tempObj = {position:(position+1),mass:mass};
					suffixMassList[position] = tempObj;
				}
			}
			if(massShift_in != 0)
			{
				for(let j=0;j<seqln;j++)
				{
					suffixMassList[j].mass = suffixMassList[j].mass + massShift_in  ;
				}
			}
			return suffixMassList;
		}
		return emptyMassList ;
	}
	addMassShift(position,massShiftList,mass){
		let len = massShiftList.length;
		for(let i=0;i<len ; i++)
		{
			if(position == massShiftList[i].position)
			{
				mass = mass + massShiftList[i].mass ;
				return mass ;
			}
		}
		return mass ;
	}
	getFixedPTMMassList(massShiftList,prsm){
		let occurence_list = [] ;
		if(prsm.annotated_protein.annotation.hasOwnProperty("ptm") )
		{
			if(Array.isArray(prsm.annotated_protein.annotation.ptm))
			{
				prsm.annotated_protein.annotation.ptm.forEach(function(ptm,index){
					if(ptm.ptm_type == "Fixed")
					{
						let mass = this.getMassofFixedPtm(ptm.ptm.abbreviation);

						if(ptm.hasOwnProperty("occurence"))
						{
							if(Array.isArray(ptm.occurence))
							{
								ptm.occurence.forEach(function(occurence,i){
									let tempObj = {"position":occurence.left_pos,"mass":mass}
									massShiftList.push(tempObj);
								});
							}
							else
							{
								let tempObj = {"position":occurence.left_pos,"mass":mass}
								massShiftList.push(tempObj);
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
						let mass = this.getMassofFixedPtm(prsm.annotated_protein.annotation.ptm.ptm.abbreviation);
						if(Array.isArray(prsm.annotated_protein.annotation.ptm.occurence))
						{
							prsm.annotated_protein.annotation.ptm.occurence.forEach(function(occurence,i){
								let tempObj = {"position":occurence.left_pos,"mass":mass}
								massShiftList.push(tempObj);
							});
						}
						else
						{
							let tempObj = {"position":prsm.annotated_protein.annotation.ptm.occurence.left_pos,"mass":mass}
							massShiftList.push(tempObj);
						}
					}
				}
			}
		}
		return massShiftList ;
	}
	getMassofFixedPtm(abbrevation)
	{
		console.log(abbrevation);
		let len = this.fixedPtmList.length;
		for(let i=0;i<len;i++)
		{
			if(this.fixedPtmList[i].name == abbrevation)
			{
				return this.fixedPtmList[i].mass;
			}
		}
	}
	// Function to add water to SuffixMass List
	addWaterMass(){
		mass = getAminoAcidDistribution(WATER)[0].mass ;
		return mass ;
	}

}

// calculatePrefixAndSuffixMass = function(){
// 	const WATER = "H2O";
// 	let protSequence = '';
	
// 	this.getPrefixMassList = function(sequence,massShiftList,massShift_in){
// 		let seqln = sequence.length;
// 		let emptyMassList = [] ;
// 		massShiftList = massShiftList.slice();
// 		if(seqln != 0)
// 		{
// 			let prefixMassList = new Array(seqln);
// 			let shiftListln = massShiftList.length;
// 			for(let i=0;i<seqln;i++)
// 			{
// 				// Get Calculated AminoAcidDistribution from aminoAcidDistributuion.js 
// 				//evaluate the return value for the case when the user enters wrong character by mistake
// 				let AcidMass;
				
// 				try{
// 					AcidMass = getAminoAcidDistribution(sequence[i])[0].mass;

// 					if(i == 0)
// 					{
// 						// Add 1 to i as the seq start from 0 but the peak id starts from 1
// 						if(shiftListln != 0)
// 						{
// 							AcidMass = this.addMassShift(i,massShiftList,AcidMass);
// 						}
// 						let tempObj = {position:(i+1),mass:AcidMass} ;
// 						prefixMassList[i] = tempObj;
// 					}
// 					else
// 					{
// 						let mass = prefixMassList[i-1].mass + AcidMass;
// 						if(shiftListln != 0)
// 						{
// 							mass = this.addMassShift(i,massShiftList,mass);
// 						}
// 						// Add 1 to i as the seq start from 0 but the peak id starts from 1
// 						let tempObj = {position:(i+1),mass:mass};
// 						prefixMassList[i] = tempObj;
// 					}
// 				}
// 				catch(error){
// 					//invalid character entered for protein sequence
// 					//window.alert("Error! Invalid amino acid in the sequence.")
// 					$('#errorAlertModal').modal('show');
// 					break;
// 				}
// 			}
// 			// Adding Mass Shift based on selection of Ion
// 			// let utilFunctionsObj = new utilFunctions();
// 			// let ionMassShiftVal = utilFunctionsObj.getNTerminusMassShiftVal();
// 			for(let j=0;j<seqln;j++)
// 			{
// 				prefixMassList[j].mass = prefixMassList[j].mass + massShift_in;
// 			}
// 			return prefixMassList;
// 		}
// 		return emptyMassList;
// 	}
// 	this.getSuffixMassList = function(sequence,massShiftList,massShift_in){
// 		let seqln = sequence.length;
// 		let emptyMassList = [];
// 		// As the elements starts from 0
// 		if(seqln != 0)
// 		{
// 			temp_seqln = seqln-1;
// 			let shiftListln = massShiftList.length;
// 			let suffixMassList = new Array(seqln);
// 			for(let i=temp_seqln ; i >= 0 ; i--)
// 			{
// 				// Get Calculated AminoAcidDistribution from aminoAcidDistributuion.js
// 				let AcidMass = getAminoAcidDistribution(sequence[i])[0].mass;
// 				if(i == temp_seqln)
// 				{
// 					// Add 1 to i as the seq start from 0 but the peak id starts from 1
// 					let position = temp_seqln - i;
// 					// Adding water mass for suffix, indirectly this will add water mass to all the masses
// 					let mass = AcidMass;//
// 					if(suffixMassList != 0)
// 					{
// 						mass = this.addMassShift(i,massShiftList,mass)
// 					}
// 					let tempObj = {position:(position+1),mass:mass} ;
// 					suffixMassList[position] = tempObj;
// 				}
// 				else
// 				{
// 					//Don't add water here
// 					let position = temp_seqln - i;
// 					let mass = suffixMassList[(position-1)].mass + AcidMass;
// 					if(suffixMassList != 0)
// 					{
// 						mass = this.addMassShift(i,massShiftList,mass)
// 					}
// 					// Add 1 to i as the seq start from 0 but the peak id starts from 1
// 					let tempObj = {position:(position+1),mass:mass};
// 					suffixMassList[position] = tempObj;
// 				}
// 			}
// 			// let utilFunctionsObj = new utilFunctions();
// 			// let ionMassShift = utilFunctionsObj.getCTerminusMassShiftVal();
// 			for(let j=0;j<seqln;j++)
// 			{
// 				suffixMassList[j].mass = suffixMassList[j].mass + massShift_in  ;
// 			}
// 			return suffixMassList;
// 		}
// 		return emptyMassList ;
// 	}
// 	this.getTotalSeqMass = function(seq,massShiftList){
// 		let mass = 0 ;
// 		let len = seq.length;
// 		for(let i=0;i<len;i++)
// 		{
// 			mass = mass + getAminoAcidDistribution(seq[i])[0].mass;
// 		}
// 		let shiftlen = massShiftList.length;
// 		for(let j=0;j<shiftlen;j++)
// 		{
// 			mass = mass + massShiftList[j].mass;
// 		}
// 		mass = mass + this.addWaterMass();

// 		return mass ;
// 	}
// 	// Function to add mass shift
// 	this.addMassShift = function(position,massShiftList,mass){
// 		let len = massShiftList.length;
// 		for(let i=0;i<len ; i++)
// 		{
// 			if(position == massShiftList[i].position)
// 			{
// 				mass = mass + massShiftList[i].mass ;
// 				return mass ;
// 			}
// 		}
// 		return mass ;
// 	}
// 	// Function to add water to SuffixMass List
// 	this.addWaterMass = function(){
// 		mass = getAminoAcidDistribution(WATER)[0].mass ;
// 		return mass ;
// 	}
// }
