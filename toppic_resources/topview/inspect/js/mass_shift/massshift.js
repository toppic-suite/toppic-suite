class MassShifts {
	constructor(){
		this.SEQUENCEDATA_ID = "#sequencedata";
	}
	/**
	 * Get the sequence entered from the HTML.
	 */
	getSequenceFromUI(){
		var seq = $(this.SEQUENCEDATA_ID).val().trim();
		seq = seq.toUpperCase();//set the sequence to be upper case automatically -- for user convenience
		
		let massShiftList = [] ;
		[seq,massShiftList]= this.getMassShiftList(seq);
		/**
		 * Remove spaces if exists between sequences
		 */
		seq = seq.replace(/ +/g, "");
		return [seq,massShiftList] ;
	}
	/**
	 * @param {string} seq - an argument with mass shift changes embeded
	 * in [] square bracket.
	 * @return {string,Array} parsedseq - sequence after removing the mass
	 * Shifts. Array with {position,mass} position-position at which 
	 * mass shift occured, mass- mass shift value.
	 */
	getMassShiftList(seq){
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
	 * @param {integer} shiftPosition - contains the position of the new mass
	 * shift entered in text box on click of any amino acid.
	 * @param {integer} massShiftVal - contains value of the mass entered.
	 * @param {Array} massShiftList - contains list of all existing mass 
	 * shifts with positions.
	 * @return{Array} with the new mass Shift value and position entered
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
	 * Get the fixed ptm mass list with position
	 * @param {string} seq - plain sequence entered by the user.
	 * @return {Array} - Returns an Array of fixed mass shift with positions.
	 */
	getFixedMassList(seq){
		let fixedShiftList = [];
		let fixedPtmAcid = null;
		let fixedPtmMass = null;
		let fixedMassShiftList = this.getFixedPtmChecklist();
		let fixedMassLen = fixedMassShiftList.length;
		for(let k=0; k<fixedMassLen; k++)
		{
			fixedPtmAcid = fixedMassShiftList[k].acid;
			fixedPtmMass = fixedMassShiftList[k].mass;
			if( fixedPtmAcid != null && fixedPtmMass != null)
			{
				let seqln = seq.length ;
				for(let i = 0 ; i<seqln;i++)
				{
					if(seq[i] == fixedPtmAcid)
					{
						let tempObj = {position:i,mass:fixedPtmMass,bg_color:null}
						fixedShiftList.push(tempObj);
					}
				}
			}
		}
		return fixedShiftList ;
	}
	/**
	 * Remove the removed fixed mass from fixed mass list
	 * @param {string} sequence - sequence without mass shifts
	 * @param {Array} fixedMassShiftList - Fixed mass shift list
	 * @param {Char} removeAcid - mass to be removed of a specific Acid
	 */
	removeFixedMassList(sequence,fixedMassShiftList,removeAcid)
	{
		let newList  = [];
		let len = fixedMassShiftList.length;
		removeAcid = removeAcid.toUpperCase()
		for(let i=0;i<len;i++)
		{
			let pos = fixedMassShiftList[i].position ;
			if(sequence[pos] != removeAcid)
			{
				newList.push(fixedMassShiftList[i]);
			}
		}
		return newList ;
	}
	/**
	 * This returns combined List of both Fixed and user entered mass shifts
	 * @param {Array} massShiftList - List of all the mass Shifts from sequence
	 * @param {Array} fixedMassShiftList - List of all selected fixed masses
	 * @returns {Array} combinedMassShiftList - List of combined lists by checking 
	 * over lapping posiitons
	 */
	getCombinedMassShiftList(massShiftList,fixedMassShiftList){
		let combinedMassShiftList = massShiftList;
		let fixedMasslen = fixedMassShiftList.length;
		let massShiftlen = massShiftList.length;
		
		for(let i=0; i<fixedMasslen ; i++)
		{
			let matched = false;
			for(let j=0;j<massShiftlen;j++)
			{
				/**
				 * Check if both has mass shift at common position and over ride the mass shift with fixed mass shift
				 */
				if(combinedMassShiftList[j].position == fixedMassShiftList[i].position)
				{
					combinedMassShiftList[j].mass = fixedMassShiftList[i].mass ;
					combinedMassShiftList[j].bg_color = fixedMassShiftList[i].bg_color ;
					matched = true;
					break;
				}
			}
			/**
			 * If no match found then copy to the new combined list
			 */
			if(!matched)
			{
				combinedMassShiftList.push(fixedMassShiftList[i]);
			}
		}
		
		return combinedMassShiftList ;
	}
	/**
	 * @return {Array} FixedPtmList - return all the selected fixed ptms with acid and mass
	 */
	getFixedPtmChecklist()
	{
		let FixedPtmList = [];
		let divs = $( ".fixedptms").get();
		$( ".fixedptms" ).each(function( index ) {
		  let acid = $( this ).find('#fixedptmacid').val().toUpperCase();
		  let mass = parseFloat($( this ).find('#fixedptmmass').val());
		  if(acid.length !=0  && mass.length != 0 && !isNaN(mass))
		  	{
				let tempfixedptm = {acid:acid,mass:mass}
				FixedPtmList.push(tempfixedptm);
			}
		});
		return FixedPtmList;
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
			return d3.ascending(x.position, y.position);
		})
		for(let i=0;i<len;i++)
		{
			let newSeqlen = newSeq.length ;
			/**
			 * Dont show when the mass is 0 in the string
			 */
			if(massShiftList[i].mass != 0){
				if(i == 0)
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
	/**
	 * write the sequence with embedded mass in [] to the screen(sequence box)
	 * @param {string} seqToUI - sequence with mass shifts embedded in []
	 */
	writeSeqToTextBox(seqToUI){
		$(this.SEQUENCEDATA_ID).val(seqToUI);
	}

}