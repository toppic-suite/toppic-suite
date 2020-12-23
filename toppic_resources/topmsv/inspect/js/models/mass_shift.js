class MassShifts {
	sequence;
	massShiftList;
	
	constructor(sequence="", massShiftList = []){
		this.sequence = sequence;
		this.massShiftList = massShiftList;
	}

	/**
	 * Generate mass shift list by given unexpected mass shift list and fixed PTM shift list
	 * @param {Array} unexpectedMassShiftList
	 * @param {Array} fixedPtmShiftList
	 * @param {Array} variablePtmShiftList
	 */
	generateMassShiftList(unexpectedMassShiftList, protVarPtmsList, variablePtmsList, fixedPtmShiftList) {
		this.massShiftList = [];
		//console.log(unexpectedMassShiftList, protVarPtmsList, variablePtmsList, fixedPtmShiftList)
		const ifMatch = (position, ptmShiftList) => {
			for (let i = 0; i < ptmShiftList.length; i++) {
				for (let j = 0; j < ptmShiftList[i].posList.length; j++){
					if(position === ptmShiftList[i].posList[j].leftPos) {
						return true;
					}	
				}
			}
			return false;
		}
		//reformat fixed and variable ptm shift list and add to massShiftList
		const addToMassShiftList = (ptmList) => {
			for (let i = 0; i < ptmList.length; i++) {
				let tempObj = {position: ptmList[i].posList, mass:ptmList[i].mono_mass, bg_color: ptmList[i].bg_color};
				this.massShiftList.push(tempObj);
			}
		}
		addToMassShiftList(fixedPtmShiftList);
		addToMassShiftList(protVarPtmsList);
		addToMassShiftList(variablePtmsList);

		for (let i = 0; i < unexpectedMassShiftList.length; i++) {
			// add unexpected mass shift if its position does not match any fixed PTM shift or variable PTM shift
			if(!ifMatch(unexpectedMassShiftList[i].position, fixedPtmShiftList) && 
			   !ifMatch(unexpectedMassShiftList[i].position, protVarPtmsList) && 
			   !ifMatch(unexpectedMassShiftList[i].position, variablePtmsList)) {
				let tempObj = {position: unexpectedMassShiftList[i].leftPos, mass:unexpectedMassShiftList[i].anno, bg_color: unexpectedMassShiftList[i].bg_color};
				this.massShiftList.push(tempObj);
			}
		}
	}

	/**
	 * Get shift mass list based on selected fixed ptm list
	 * @param {Array} selectedFixedPtmList - an array which contains fixed mass shift selected by users, EX. [{acid: X, mass: 12}]
	 * @return {Array} - Returns an Array of fixed mass shift with positions. EX. [{position:i, mass:x, bg_color:null}]
	 */
	getFixedMassShiftList(selectedFixedPtmList){
		let result = [];
		let fixedPtmAcid = null;
		let fixedPtmMass = null;
		for(let k=0; k<selectedFixedPtmList.length; k++)
		{
			fixedPtmAcid = selectedFixedPtmList[k].acid;
			fixedPtmMass = selectedFixedPtmList[k].mass;

			let tempObj = {mono_mass:fixedPtmMass, name:fixedPtmAcid, posList:[]}
			for(let i = 0 ; i<this.sequence.length; i++)
			{
				if(this.sequence[i] === fixedPtmAcid)
				{
					tempObj.posList.push({leftPos:i, rightPos: i + 1, acid:fixedPtmAcid});
				}
			}
			result.push(tempObj);
		}
		return result;
	}

	/**
	 * Remove acid mass shift from mass shift list
	 * @param {Char} removeAcid - mass to be removed of a specific Acid
	 */
	removeFixedMassList(removeAcid)
	{
		removeAcid = removeAcid.toUpperCase()
		for(let i=0;i<this.massShiftList.length;i++)
		{
			let position = this.massShiftList[i].position;
			if(this.sequence[position] === removeAcid)
			{
				this.massShiftList.splice(i,1);
			}
		}
		return this.massShiftList;
	}

	/**
	 * @param {integer} massShiftPosition - contains the position of the new mass
	 * shift entered in text box on click of any amino acid.
	 * @param {integer} massShiftValue - contains value of the mass entered.
	 * @param {string} bg_color - backgroud color
	 * @return {Array} with the new mass Shift value and position entered
	 * or changes the existing mass shift value.
	 */
	appendtoMassShiftList(massShiftPosition,massShiftValue,bg_color){
		if(massShiftValue === 0) {
			return this.massShiftList;
		}
		let matchFound = false;
		for(let i=0; i<this.massShiftList.length;i++)
		{
			if(massShiftPosition === this.massShiftList[i].position)
			{
				this.massShiftList[i].mass = massShiftValue;
				this.massShiftList[i].bg_color = bg_color;
				matchFound = true ;
			}
		}
		if(!matchFound)
		{
			let tempShiftObj = {position:massShiftPosition,mass:massShiftValue,bg_color:bg_color};
			this.massShiftList.push(tempShiftObj);
		}
		return this.massShiftList;
	}
	
	
	/**
	 * forms the seq with all the mass lists and selected fixed ptms
	 * @return {string} result - sequence with mass shifts embedded in []
	 */
	formSequence(){
		let result = this.sequence;
		let count = 0;
		if(!this.massShiftList) {
			return result;
		}
		// sort mass shift list by position, ascending
		this.massShiftList.sort(function(x,y){
            return x.position - y.position;
		})
		for(let i=0; i<this.massShiftList.length; i++)
		{
			if(this.massShiftList[i].mass !== 0){
				if(i > 0)
				{
					// this is the previous added mass
					let tempString = "["+this.massShiftList[i-1].mass+"]";
					count = count + tempString.length;
				}
				// add +1 as the position need to be added after the position of the acid.
				let tempPosition = this.massShiftList[i].position + 1 + count;
				result = result.slice(0, tempPosition) + "["+ this.massShiftList[i].mass + "]" + result.slice(tempPosition);
			}
		}
		return result;
	}
}