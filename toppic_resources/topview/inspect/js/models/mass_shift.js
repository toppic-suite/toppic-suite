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
	 */
	generateMassShiftList(unexpectedMassShiftList, fixedPtmShiftList) {
		function ifMatch(position, fixedPtmShiftList) {
			for (let i = 0; i < fixedPtmShiftList.length; i++) {
				if(position === fixedPtmShiftList[i].position) {
					return true;
				}
			}
			return false;
		}

		this.massShiftList = [...fixedPtmShiftList];
		for (let i = 0; i < unexpectedMassShiftList.length; i++) {
			// add unexpected mass shift if its position does not match any fixed PTM shift
			if(!ifMatch(unexpectedMassShiftList[i].position, fixedPtmShiftList)) {
				let tempObj = {position: unexpectedMassShiftList[i].position, mass:unexpectedMassShiftList[i].mass, bg_color: unexpectedMassShiftList[i].bg_color};
				this.massShiftList.push(tempObj);
			}
		}
		return this.massShiftList;
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
		// let selectedFixedPtmList = getFixedPtmCheckList();
		for(let k=0; k<selectedFixedPtmList.length; k++)
		{
			fixedPtmAcid = selectedFixedPtmList[k].acid;
			fixedPtmMass = selectedFixedPtmList[k].mass;
			for(let i = 0 ; i<this.sequence.length; i++)
			{
				if(this.sequence[i] === fixedPtmAcid)
				{
					let tempObj = {position:i, mass:fixedPtmMass, bg_color:null}
					result.push(tempObj);
				}
			}
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
			if(sequence[position] === removeAcid)
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
		// sort mass shift list by position, ascending
		this.massShiftList.sort(function(x,y){
            return x.position - y.position;
			// return d3.ascending(x.position, y.position);
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