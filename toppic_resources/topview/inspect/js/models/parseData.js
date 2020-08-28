/**
 * receives file name to parse and replace "" and [] to make it same as the result from original JSON.parse function
 * @param {String} dataName - contains parsePeakMass in a sting format as we stringify while storing in local storage variable
 */
function parsePeakMass(dataName){
	let data = window.localStorage.getItem(dataName);
		
	if (data != "" && data != null){
		data = data.replace("[", "")
		data = data.replace("]", "")
		data = data.replace(/,/g, "\n");
		data = data.replace(/"/g, "");
	}
	return data;
}
/**
 * receives file name to parse and replace "" and [] to make it same as the result from original JSON.parse function
 * @param {String} dataName - contains sequence in a sting format as we stringify while storing in local storage variable
 */
function parseSeq(dataName){
	let data = window.localStorage.getItem(dataName);
	
	if (data != "" && data != null){
		data = data.replace("[", "")
		data = data.replace("]", "")
		data = data.replace(/"/g, "");
	}
	return data;
}
/**
 * Function returns a json formatted after getting fixed ptm data from local storage
 * @param {String} dataName - contains Fixed PTMs in a sting format as we stringify while storing in local storage variable
 */
function parsePTM(dataName){
	//ptm is an object inside array, so to preserve its structure, using JSON parse
	//it is relatively small to other lists so the performance should not deteriorate much
	let data = JSON.parse(window.localStorage.getItem(dataName));
	return data;
}
/**
 * Function returns a json formatted data after getting unknow mass list data from local storage
 * @param {String} dataName - contains unkwon mass lists in a sting format as we stringify while storing in local storage variable
 */
function parseUnknowmassList(dataName){
	let data = JSON.parse(window.localStorage.getItem(dataName));
	return data;
}
/**
 * Function returns a json formatted data after getting precursor mass data from local storage
 * @param {String} dataName - contains precursor mass in a sting format as we stringify while storing in local storage variable
 */
function parsePrecursorMass(dataName){
	let data = parseFloat(JSON.parse(window.localStorage.getItem(dataName)));
	return data;
}

/**
 * @param {string} seq - an argument with mass shift changes embeded in [] square bracket.
 * @return {string} parsedseq - sequence after removing the mass
 * Shifts. 
 * @returns {Array} massShiftList - Array with {position,mass} position-position at which 
 * mass shift occured, mass- mass shift value.
 */
function parseSequenceMassShift(seq){
	let massShiftList = [] ;
	let parsedseq = "";
	let splitStr = seq.split(/\[(.*?)\]/);
	let splitArraylen = splitStr.length;
	let position = 0;
	
	for(let i = 0 ; i<splitArraylen;i++)
	{
		if(isNaN(splitStr[i]))
		{
			parsedseq = parsedseq + splitStr[i] ;
			position = position + splitStr[i].length ;
		}
		else
		{
			let mass = parseFloat(splitStr[i]);
			/**
			 * remove 1 as the data starts from 0 and length starts from 1
			 */
			let tempPosition = position - 1;
			//Initially set the bg_color to null
			let shiftobj = {position:tempPosition, mass:mass, bg_color:null};
			/**
			 * when the split occur at the end we get an extra "" in 
			 * the list. This is to check if the mass is numeric.
			 */
			if(!isNaN(mass))
			{
				massShiftList.push(shiftobj);
			}
		}
	}
	return [parsedseq,massShiftList] ;
}

// form residues from sequence
let formResidues = (sequence) => {
	let residues = [];
	for (let i = 0; i < sequence.length; i++) {
		let tempObj = {
			position: i.toString(),
			acid: sequence.charAt(i).toUpperCase()
		}
		residues.push(tempObj);
	}
	return residues;
}

// form fixed ptms
let formFixedPtms = (fixedMassShiftList, fixedPtmNameList, sequence) => {
	let result = [];
	result.push(fixedPtmNameList[0]);
	let tempArray = [];
	fixedMassShiftList.forEach((element) => {
		let tempObj = {
			pos: element.position.toString(),
			acid: sequence.charAt(element.position)
		};
		tempArray.push(tempObj);
	});
	result[0].posList = tempArray;
	return result;
}

let formMassShifts = (unknownMassShiftList) => {
	let result = [];
	unknownMassShiftList.forEach((element)=> {
		let tempObj = {
			anno: element.mass.toString(),
			leftPos: (element.position).toString(),
			rightPos: (element.position + 1).toString()
		}
		result.push(tempObj);
	})
	return result;
}

/**
 * @function getTotalSeqMass
 * @description Returns total mass of the sequence
 * @param {String} seq - Contains Protein sequence
 * @param {Array} massShiftList - Contains list of mass shifts
 */
let getTotalSeqMass = (seq,massShiftList) => {
	let mass = 0 ;
	let len = seq.length;
	for(let i=0;i<len;i++)
	{
		mass = mass + getAminoAcidDistribution(seq[i])[0].mass;
	}
	let shiftlen = massShiftList.length;
	for(let j=0;j<shiftlen;j++)
	{
		mass = mass + massShiftList[j].mass;
	}
	mass = mass + getAminoAcidDistribution("H2O")[0].mass;

	return mass ;
}

let getIons = (monoMassList) => {
	let ionData = [];
	monoMassList.forEach((element) => {
		let tempIonData = {"mz":element.mass,"intensity":element.intensity,"text": element.ion, "error": element.massError};
        ionData.push(tempIonData);
	});
	return ionData;
}

let ifNIon = (ionType) => {
	if(ionType === "a" || ionType === "b" || ionType === "c") return true ;
	else return false ;
}

let formBreakPoints = (matchedPeakList, sequenceLength) => {
	let breakPointsMap = new Map();
	let result = [];
	matchedPeakList.forEach((element) => {
		if (breakPointsMap.has(element.position)) {
			let tempObj = breakPointsMap.get(element.position);
			tempObj.existNIon = (ifNIon(element.ion[0]) || tempObj.existNIon);
			tempObj.existCIon = (!ifNIon(element.ion[0]) || tempObj.existCIon);
			tempObj.anno = [tempObj.anno, [element.ion.toUpperCase(), element.charge + "+"].join(" ")].join(" ");
			breakPointsMap.set(element.position, tempObj);
		} else {
			let tempObj = {};
			tempObj.position = element.position;
			tempObj.anno = [element.ion.toUpperCase(), element.charge + "+"].join(" ");
			tempObj.existNIon = ifNIon(element.ion[0]);
			tempObj.existCIon = !ifNIon(element.ion[0]);
			breakPointsMap.set(element.position, tempObj);
		}
	});
	breakPointsMap.forEach((val, key) => {
		result.push(val);
	});
	return result;
}