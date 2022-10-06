/**
 * returns ion terminal based on ion name
 * @param {String} ionName - ion name like X, A, B, Z..
*/
function getIonTerminal(ionName: string): string | null{
	if (ionName.indexOf("A") > -1 || ionName.indexOf("B") > -1 || ionName.indexOf("C") > -1) {
		return "N";
	}
	else if (ionName.indexOf("X") > -1 || ionName.indexOf("Y") > -1 || ionName.indexOf("Z") > -1) {
		return "C";
	}
	return null;
}

/**
 * receives file name to parse and replace "" and [] to make it same as the result from original JSON.parse function
 * @param {String} dataName - contains parsePeakMass in a sting format as we stringify while storing in local storage variable
 */
function parsePeakMass(dataName: string): string | null {
	let data: string | null = window.localStorage.getItem(dataName);
		
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
function parseSeq(dataName: string): string | null {
	let data: string | null = window.localStorage.getItem(dataName);
	
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
function parsePTM(dataName: string): MassShift[] {
	//ptm is an object inside array, so to preserve its structure, using JSON parse
	//it is relatively small to other lists so the performance should not deteriorate much
	let massShiftData: string | null = window.localStorage.getItem(dataName);
	let massShiftList: MassShift[] = [];
	if (massShiftData) {
		let massShifts: any = JSON.parse(massShiftData);
        massShifts.forEach((d: any) => {
        	let massShift = new MassShift(d.leftPos_, d.rightPos_, d.massShift_, d.type_, d.annotation_, d.ptmList_);
            massShiftList.push(massShift);
        });
	}
	return massShiftList;
}
/**
 * Function returns a json formatted data after getting unknow mass list data from local storage
 * @param {String} dataName - contains unkwon mass lists in a sting format as we stringify while storing in local storage variable
 */
/*function parseUnknowmassList(dataName: string): MassShift[] | null {
	let massShiftData: string | null = (window.localStorage.getItem(dataName));
	let massShiftList: MassShift[] = [];
	
	if (massShiftData) {
		massShiftList = JSON.parse(massShiftData);
	}
	return massShiftList;
}*/
/**
 * Function returns a json formatted data after getting precursor mass data from local storage
 * @param {String} dataName - contains precursor mass in a sting format as we stringify while storing in local storage variable
 */
function parsePrecursorMass(dataName: string): number | null {
	let data: string | null = window.localStorage.getItem(dataName);
	let precMass: number | null = null;
	if (data) {
		precMass = parseFloat(data);
	}
	return precMass;
}

/**
 * @param {string} seq - an argument with mass shift changes embeded in [] square bracket.
 * @return {string} parsedseq - sequence after removing the mass
 * Shifts. 
 * @returns {Array} massShiftList - Array with {position,mass} position-position at which 
 * mass shift occured, mass- mass shift value.
 */
function parseSequenceMassShift(seq: string): [string,MassShift[], MassShift[],MassShift[]] | null {
	let unknownMassShiftList: MassShift[] = [];
	let variablePtmsList: MassShift[] = [];
	let protVarPtmsList: MassShift[] = [];
	let parsedseq: string = "";
	let position: number = 0;
	let massShift: string = "";
	let isMassShift: boolean = false;
	let isBracketClosed: boolean = true; //to detect invalid brackets; 
	
	for (let i: number = 0; i < seq.length; i++){
		/**
		 * remove 1 as the data starts from 0 and length starts from 1
		*/
		let tempPosition: number = position - 1;

		if (seq[i] == "["){//mass shift
			//make sure it is annotating a residue
			if (i > 0 && seq[i-1].charCodeAt(0) < 65 || seq[i-1].charCodeAt(0) > 90) {
				alert("Invalid mass shift annotation!");
				return null;
			}
			if (isBracketClosed) {
				isMassShift = true;
				isBracketClosed = false;
			}
			else {
				alert("Bracket [] is not closed properly in the sequence.");
				return null;
			}
		}
		else if (seq[i] == "]"){//mass shift
			if (isBracketClosed) {
				alert("Bracket [] is not closed properly in the sequence.");
				return null;
			}
			isMassShift = false;
			//check mass shift value
			//massShift = massShift.slice(1);//"[" is included in mass shift string
			if (isNaN(parseFloat(massShift))){
				if (massShift == '') {
					alert("Invalid PTM annotation in the sequence. There is an empty bracket [].");
					return null;
				}
				let isVariablePtm = false;
				//check if this is variable ptm
				let listOfPtm: string[] = massShift.split(";");

				for (let k = 0; k < listOfPtm.length; k++) {
					for (let j: number = 0; j < commonPtmList.length; j++) {
						if (listOfPtm[k] == commonPtmList[j].abbr.toUpperCase()) {
							let varPtm: MassShift = new MassShift(tempPosition, tempPosition + 1, commonPtmList[j].mass, "Variable", commonPtmList[j].abbr, new Mod(seq[i], commonPtmList[j].mass, commonPtmList[j].name));
							variablePtmsList.push(varPtm);
							isVariablePtm = true;
						}
					}
				}
				if (!isVariablePtm) {
					alert(massShift + " is not a valid name for variable PTM.");
					return null;
				}
			}
			else{
				let mass: number = parseFloat(massShift);
				/**
				 * when the split occur at the end we get an extra "" in 
				 * the list. This is to check if the mass is numeric.
				 */
				if(!isNaN(mass))
				{
					let unknownMassShift: MassShift = new MassShift(tempPosition, tempPosition + 1, mass, "unexpected", mass.toString());
					unknownMassShiftList.push(unknownMassShift);
				}
			}
			massShift = "";
			isBracketClosed = true;
		}
		else{
			if (isMassShift){
				massShift = massShift + seq[i];
			}
			else{
				if (seq[i].charCodeAt(0) >= 65 && seq[i].charCodeAt(0) <= 90){//if it is  alphabet
					parsedseq = parsedseq + seq[i];
					position++;
				}
				else{
					alert("Invalid character found in the sequence.");
				}
			}
		}
	}
	if (!isBracketClosed) {
		alert("Bracket [] is not closed properly in the sequence.");
	}
	return [parsedseq,unknownMassShiftList, protVarPtmsList, variablePtmsList] ;
}
/*function parseSequenceMassShift(seq){
	let unknownMassShiftList = [];
	let variablePtmsList = [];
	let protVarPtmsList = [];
	let parsedseq = "";
	let splitStr = seq.split(/\[(.*?)\]/);
	let splitArraylen = splitStr.length;
	let position = 0;

	for(let i = 0 ; i<splitArraylen;i++)
	{
		if(isNaN(splitStr[i]))
		{
			isPTM = false;
			/*if (VAR_PTM_LIST.length > 0){
				//console.log("var_ptm_list", VAR_PTM_LIST)
				for (let j = 0; j < VAR_PTM_LIST.length; j++){
					//because mass information is needed, which is not written in sequence
					if (VAR_PTM_LIST[j].name.toUpperCase() == splitStr[i]){
						//let tempPosition = position -1;
						let shiftobj = {posList:VAR_PTM_LIST[j].posList, name: VAR_PTM_LIST[j].name, mono_mass:VAR_PTM_LIST[j].mono_mass, bg_color:null};
						isPTM = true;
						if (VAR_PTM_LIST[j].type == "Protein variable")
						{
							protVarPtmsList.push(shiftobj);
						}
						else
						{
							variablePtmsList.push(shiftobj);
						}
						break;		
					}
				}
			}*/
			/*if (!isPTM) {
				parsedseq = parsedseq + splitStr[i] ;
				position = position + splitStr[i].length ;
			}
		}
		else
		{
			let mass = parseFloat(splitStr[i]);
			/**
			 * remove 1 as the data starts from 0 and length starts from 1
			 */
			/*let tempPosition = position - 1;
			//Initially set the bg_color to null
			let shiftobj = {leftPos:tempPosition, rightPos:tempPosition + 1, anno:mass, label:mass, bg_color:null};
			/**
			 * when the split occur at the end we get an extra "" in 
			 * the list. This is to check if the mass is numeric.
			 */
			/*if(!isNaN(mass))
			{
				unknownMassShiftList.push(shiftobj);
			}
		}
	}
	return [parsedseq,unknownMassShiftList, protVarPtmsList, variablePtmsList] ;
}*/
/** 
  * Get fixed mass selected by user
  */
function parseCheckedFixedPtm(seq: string) {
    let fixedPtmList: MassShift[] = [];
    let checkedFixedPtm: Mod[] = getFixedPtmCheckList();
    checkedFixedPtm.forEach(ptm => {
      let fixedPtm: Mod = {} as Mod;
	  let isCustomPtm: boolean = true;
      for(let j: number =0; j<COMMON_FIXED_PTM_LIST.length;j++) {
        if(ptm.getShift() == COMMON_FIXED_PTM_LIST[j].getShift()) {
        let name: string = COMMON_FIXED_PTM_LIST[j].getName();
          fixedPtm = new Mod(ptm.getResidue(), ptm.getShift(), name);
          break;
        }
      }
	  if (isCustomPtm) {
		for (let j: number = 0; j < USER_FIXED_PTM_LIST.length; j++) {
			if (ptm.getShift() == USER_FIXED_PTM_LIST[j].getShift()) {
				let name: string = USER_FIXED_PTM_LIST[j].getName();
				fixedPtm = new Mod(ptm.getResidue(), ptm.getShift(), name);
				break;
			}
		}
	}
      for(let i: number = 0 ; i < seq.length; i++) {
        if(seq[i] === ptm.getResidue()) {
          let massShift: MassShift = new MassShift(i, i + 1, fixedPtm.getShift(), "Fixed", fixedPtm.getName(), fixedPtm);
          fixedPtmList.push(massShift);
        }
      }
	})
	return fixedPtmList;
}
/*
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
}*/
/*
// form fixed ptms
const formFixedPtms = (fixedMassShiftList: MassShift[], sequence: string) => {
	let result = [];
	result.push({});
	// result.push(fixedPtmNameList[0]);
	let tempArray = [];
	fixedMassShiftList.forEach((element) => {
		for (let i: number = 0; i < element.posList.length; i++) {
			let tempObj = {
				pos: element.posList[i].pos.toString(),
				acid: sequence.charAt(element.posList[i].pos)
			};
			tempArray.push(tempObj);
		}
	});
	result[0].posList = tempArray;
	return result;
}

// form variable ptms
let formVariablePtms = (variablePtmsList, sequence) => {
	let result = [];
	result.push({});
	// result.push(fixedPtmNameList[0]);
	let tempArray = [];
	variablePtmsList.forEach((element) => {
		for (let i = 0; i < element.posList.length; i++) {
			let tempObj = {
				pos: element.posList[i].pos.toString(),
				acid: sequence.charAt(element.posList[i].pos)
			};
			tempArray.push(tempObj);
		}
	});
	result[0].posList = tempArray;
	return result;
}
let formMassShifts = (unknownMassShiftList) => {
	let result = [];
	unknownMassShiftList.forEach((element)=> {
		let tempObj = {
			anno: element.anno.toString(),
      shift: element.anno.toSring(),
			leftPos: (element.leftPos).toString(),
			rightPos: (element.rightPos).toString()
		}
		result.push(tempObj);
	})
	return result;
}
*/
/**
 * @function getTotalSeqMass
 * @description Returns total mass of the sequence
 * @param {String} seq - Contains Protein sequence
 * @param {Array} massShiftList - Contains list of mass shifts
 */
const getTotalSeqMass = (seq: string,massShiftList: MassShift[]): number => {
	let mass: number = 0 ;
	let len: number = seq.length;
	for(let i: number =0;i<len;i++)
	{
		let aminoMass: AminoDist[] | null = getAminoAcidDistribution(seq[i]);
		if (aminoMass) {
			mass = mass + aminoMass[0].mass;
		}
		else{
			console.error("invalid amino acid distribution");
		}
	}
	let shiftlen: number = massShiftList.length;
	for(let j: number=0;j<shiftlen;j++)
	{
		mass = mass + massShiftList[j].getShift();
	}
	let waterMass: AminoDist[] | null = getAminoAcidDistribution("H2O");
	if (waterMass) {
		mass = mass + waterMass[0].mass;
	}
	else{
		console.error("invalid water mass distribution");
	}
	return mass;
}
const ifNIon = (ionType: string): boolean => {
    if (ionType === "A" || ionType === "B" || ionType === "C")
        return true;
    else
        return false;
};
const formBreakPoints = (matchedPeakList: MatchedUnMatchedPeak[]): BreakPoints[] => {
    let breakPointsMap = new Map();
    let result: BreakPoints[] = [];
    matchedPeakList.forEach((element: MatchedUnMatchedPeak) => {
        if (breakPointsMap.has(element.position)) {
            let tempObj: BreakPoints = breakPointsMap.get(element.position);
            tempObj.existNIon = (ifNIon(element.ion[0]) || tempObj.existNIon);
            tempObj.existCIon = (!ifNIon(element.ion[0]) || tempObj.existCIon);
            tempObj.anno = [tempObj.anno, [element.ion.toUpperCase(), element.charge + "+"].join(" ")].join(" ");
            breakPointsMap.set(element.position, tempObj);
        }
        else {
            let tempObj: BreakPoints = {} as BreakPoints;
            tempObj.position = element.position.toString();
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
};
let convertIonName = (ionName: string): string => {
    //if Z_DOT, H2O, NH3 are included in the ion name, convert them
    if (ionName.indexOf("Z_DOT") > -1) {
        ionName = ionName.replace("Z_DOT", "Z˙");
    }
    if (ionName.indexOf("H2O") > -1) {
        ionName = ionName.replace("H2O", "H₂O");
    }
    if (ionName.indexOf("NH3") > -1) {
        ionName = ionName.replace("NH3", "NH₃");
    }
    return ionName;
};
let getIonsSpectrumGraph = (matchedPeakList: MatchedUnMatchedPeak[], envelopeList: Envelope[]): MatchedIon[] => {
    let ionData: MatchedIon[] = [];
    //generate ion list
    for (let i: number = 0; i < matchedPeakList.length; i++) {
        for (let j: number = 0; j < envelopeList.length; j++) {
            let env: Envelope = envelopeList[j];
            let peak: MatchedUnMatchedPeak = matchedPeakList[i];
            if (env.getMonoMass() == peak.mass) {
                env.getPeaks().sort(function (x, y) {
                    return d3.descending(x.getIntensity(), y.getIntensity());
                });
                let ion: {"env": Envelope, "error": number, "intensity": number,
				"mz": number, "text": string} = {
                    "env": env,
                    "error": peak.massError,
                    "intensity": env.getPeaks()[0].getIntensity(),
                    "mz": env.getPeaks()[0].getMonoMz(),
                    "text": peak.ion.toUpperCase()
                };
                //if it is z_dot ion, text should be converted 
                ion.text = convertIonName(ion.text);
                ionData.push(ion);
            }
        }
    }
    return ionData;
};
//modified version of addOneIon function from loadSpectra.js from visual/prsm. 
const addOneIonInspect = (ionList: MatchedIon[], ion: MatchedUnMatchedPeak) => {
    //console.log("ion", ion)
    let idx: number = -1;
    let ionId: string = ion.ion;
    if (ionId.indexOf("Z_DOT") > -1) { //if ion is z_dot, ionList[i].text is already converted from Z_DOT to Z˙
        ionId = ionId.replace("Z_DOT", "Z˙");
    }
    for (let i: number = 0; i < ionList.length; i++) {
        if (ionId == ionList[i].text) {
            idx = i;
            break;
        }
    }
    if (idx == -1) {
        let tempIonData:MatchedIon = { "mz": ion.mass, "intensity": ion.intensity, "text": ion.ion, "pos": ion.ionPos, "error": ion.massError };
        //if it is z_dot ion, text should be converted 
        tempIonData.text = convertIonName(tempIonData.text);
        ionList.push(tempIonData);
    }
    else {
        if (ion.intensity > ionList[idx].intensity) {
            ionList[idx].intensity = ion.intensity;
        }
    }
};
const getIonsMassGraph = (matchedPeakList: MatchedUnMatchedPeak[]): MatchedIon[] => {
    let ionData: MatchedIon[] = [];
    /*monoMassList.forEach((element) => {
        let tempIonData = {"mz":element.mass,"intensity":element.intensity,"text": element.ion, "error": element.massError};
        ionData.push(tempIonData);
    });*/
    matchedPeakList.forEach((element: MatchedUnMatchedPeak) => {
        addOneIonInspect(ionData, element);
    });
    return ionData;
};
