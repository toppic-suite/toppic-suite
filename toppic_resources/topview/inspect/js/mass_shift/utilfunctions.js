utilFunctions = function(){
	const N_TERMINUS = "NTERMINUS";
	const C_TERMINUS = "CTERMINUS";
	
	this.getTerminus = function(ionType){
		if(ionType == "a" || ionType == "b" || ionType == "c") return N_TERMINUS ;
		else return C_TERMINUS ;
	}
	this.getMassShift = function(ions){
		let MassShift = 0;
		ions.forEach(function(mass){
			MassShift = MassShift + parseFloat(mass);
		})
		console.log("MassShift : ", MassShift);
		return MassShift;
	}
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