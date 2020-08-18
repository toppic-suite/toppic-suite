class MassShift {
	constructor(sequence){
		this.sequence = sequence;
	}

	/**
	 * @param {string} seq - an argument with mass shift changes embeded in [] square bracket.
	 * @return {string} parsedseq - sequence after removing the mass
	 * Shifts. 
	 * @returns {Array} massShiftList - Array with {position,mass} position-position at which 
	 * mass shift occured, mass- mass shift value.
	 */
	parseSequenceMassShift(seq){
		let massShiftList = [] ;
		let parsedseq = "";
		let splitStr = seq.split(/\[(.*?)\]/);
		let splitArraylen = splitStr.length;
		let position = 0 ;
		
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
				let shiftobj = {mass:mass,position:tempPosition,bg_color:null};
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


	/**
	 * Get mass shift list by given unexpected mass shift list and fixed PTM shift list
	 * @param {Array} unexpectedMassShiftList
	 * @param {Array} fixedPtmShiftList
	 */
	getMassShiftList(unexpectedMassShiftList, fixedPtmShiftList) {
		function ifMatch(position, fixedPtmShiftList) {
			for (let i = 0; i < fixedPtmShiftList.length; i++) {
				if(position === fixedPtmShiftList[i].position) {
					return true;
				}
			}
			return false;
		}

		let result = [...fixedPtmShiftList];
		for (let i = 0; i < unexpectedMassShiftList.length; i++) {
			// add unexpected mass shift if its position does not match any fixed PTM shift
			if(!ifMatch(unexpectedMassShiftList[i].position, fixedPtmShiftList)) {
				let tempObj = {mass:unexpectedMassShiftList[i].mass, position: unexpectedMassShiftList[i].position, bg_color: unexpectedMassShiftList[i].bg_color};
				result.push(tempObj);
			}
		}
		return result;
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
	 * Remove acid mass shift from fixed mass list
	 * @param {string} sequence - sequence without mass shifts
	 * @param {Array} fixedMassShiftList - Fixed mass shift list
	 * @param {Char} removeAcid - mass to be removed of a specific Acid
	 */
	removeFixedMassList(fixedMassShiftList,removeAcid)
	{
		let result  = [];
		removeAcid = removeAcid.toUpperCase()
		for(let i=0;i<fixedMassShiftList.length;i++)
		{
			let pos = fixedMassShiftList[i].position;
			if(this.sequence[pos] !== removeAcid)
			{
				result.push(fixedMassShiftList[i]);
			}
		}
		return result ;
	}

	/**
	 * @param {integer} shiftPosition - contains the position of the new mass
	 * shift entered in text box on click of any amino acid.
	 * @param {integer} massShiftVal - contains value of the mass entered.
	 * @param {Array} massShiftList - contains list of all existing mass 
	 * shifts with positions.
	 * @return {Array} with the new mass Shift value and position entered
	 * or changes the existing mass shift value.
	 */
	appendtoMassShiftList(shiftPosition,massShiftVal,massShiftList,bg_color){
		let newMassShiftList = [];
		let len = massShiftList.length;
		let matchFound = false ;
		for(let i=0; i<len;i++)
		{
			if(shiftPosition == massShiftList[i].position)
			{
				massShiftList[i].mass = massShiftVal;
				massShiftList[i].bg_color = bg_color;
				matchFound = true ;
			}
		}
		if(!matchFound)
		{
			let tempShiftObj = {mass:massShiftVal,position:shiftPosition,bg_color:bg_color};
			massShiftList.push(tempShiftObj);
			matchFound = false;
		}
		let newlen = massShiftList.length;
		for(let j=0; j<newlen; j++)
		{
			if(massShiftList[j].mass != 0)
			{
				newMassShiftList.push(massShiftList[j]);
			}
		}
		return newMassShiftList ;
	}
	
	
	/**
	 * forms the seq with all the mass lists and selected fixed ptms
	 * @param {string} seq - sequence with only aminoacids and without mass lists embedded
	 * @param {Array} massShiftList - List with all the combined mass shifts
	 * @return {string} newSeq - sequence with mass shifts embedded in []
	 */
	formSequence(seq,massShiftList){
		let newSeq = seq;
		let len = massShiftList.length; 
		let seqLen = seq.length ;
		let count = 0;
		/**
		 * sorting the lists with position
		 */
		massShiftList.sort(function(x,y){
            return x.position - y.position;
			// return d3.ascending(x.position, y.position);
		})
		for(let i=0;i<len;i++)
		{
			let newSeqlen = newSeq.length ;
			/**
			 * Dont show when the mass is 0 in the string
			 */
			if(massShiftList[i].mass !== 0){
				if(i === 0)
				{
					/**
					 * Add +1 as we need to append the mass after the current position
					 */
					let tempPosition = massShiftList[i].position+1 ;
					newSeq = newSeq.slice(0, tempPosition) + "["+ massShiftList[i].mass + "]"+newSeq.slice(tempPosition, newSeqlen);
				}
				else
				{
					/**
					 * Form the mass between []
					 */
					let tempString = "["+massShiftList[i-1].mass+"]";
					count = count + tempString.length;
					/**
					 * add +1 as the position need to be added after 
					 * the position of the acid.
					 */
					let tempPosition = massShiftList[i].position + 1 + count ;
					newSeq = newSeq.slice(0, tempPosition) + "["+ massShiftList[i].mass + "]" + newSeq.slice(tempPosition, newSeqlen);
				}
			}
		}
		return newSeq ;
	}
}