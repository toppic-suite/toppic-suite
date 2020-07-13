/**
 * @function utilFunctions
 * @description Utility function retrieve data from UI
 */
utilFunctions = function(){
	const N_TERMINUS = "NTERMINUS";
	const C_TERMINUS = "CTERMINUS";
	
	/**
	 * @function getTerminus
	 * @description function checks for ion type and returns whether it is C terminus or N terminus
	 * @param {String} ionType - Contains information of the ion type from UI
	 */
	this.getTerminus = function(ionType){
		if(ionType == "a" || ionType == "b" || ionType == "c") return N_TERMINUS ;
		else return C_TERMINUS ;
	}
	/**
	 * @function getMassShift
	 * @description Function checks all the ions selected and adds together to get the mass shift needed to add to suffix mass list
	 * @param {Array} ions - Contains all the ion types selected
	 */
	this.getMassShift = function(ions){
		let MassShift = 0;
		ions.forEach(function(mass){
			MassShift = MassShift + parseFloat(mass);
		})
		return MassShift;
	}
	/**
	 * @function getNTerminusMassShiftVal 
	 * @description Checks for all the ions selected and gets all the consolidated mass of all the N terminus ions selected by user
	 */
	this.getNTerminusMassShiftVal = function(){
		var ions = [];
		$.each($("input[name='nterminus']:checked"), function(){
				console.log("this : ", $(this).attr( "id" ));
                ions.push($(this).val());
           });
		console.log("checkedValues : ", ions);
		let massShift = this.getMassShift(ions);
		return massShift;
	}
	/**
	 * @function getCTerminusMassShiftVal
	 * @description Checks for all the ions selected and gets all the consolidated mass of all the C terminus ions selected by user
	 */
	this.getCTerminusMassShiftVal = function(){
		let ions = [];
		$.each($("input[name='cterminus']:checked"), function(){
                ions.push($(this).val());
           });
		console.log("ionType : ", ions);
		let massShift = this.getMassShift(ions);
		return massShift;
	}
} 