class iontabledata{
	
	/**
	 * @function getActualIdvalues
	 * @description Dictionary to get correspoding actual heading for 
	 * all the Id's of Ion Fragmentation from UI.
	 * @param {String} ionType - Id of the ion from UI
	 */
	static getActualIdvalues(ionType){
		let dict = [];
		dict["a"] = "a";
		dict['a1'] = "a-H<sub>2</sub>O";
		dict['a2'] = "a-NH<sub>3</sub>";

		dict["b"] = "b";
		dict['b1'] = "b-H<sub>2</sub>O";
		dict['b2'] = "b-NH<sub>3</sub>";

		dict["c"] = "c";
		dict['c1'] = "c-H<sub>2</sub>O";
		dict['c2'] = "c-NH<sub>3</sub>";

		dict["x"] = "x";
		dict['x1'] = "x-H<sub>2</sub>O";
		dict['x2'] = "x-NH<sub>3</sub>";

		dict["y"] = "y";
		dict['y1'] = "y-H<sub>2</sub>O";
		dict['y2'] = "y-NH<sub>3</sub>";

		dict["z"] = "z";
		dict['z1'] = "z-H<sub>2</sub>O";
		dict['z2'] = "z-NH<sub>3</sub>";

		dict["z_"] = "z&deg;";
		dict['z_1'] = "z&deg;-H<sub>2</sub>O";
		dict['z_2'] = "z&deg;-NH<sub>3</sub>";

		return dict[ionType];
	}
	/**
	 * Get all the N terminus Ions from UI
	 */
	getNterminusCheckedList(){
		var ions = [];
		$.each($("input[name='nterminus']:checked"), function(){
				let id = $(this).attr( "id" );
				let value = $(this).val();
				let ionType = iontabledata.getActualIdvalues(id);
				let temp_Obj = {ionType:ionType,mass:value}
                ions.push(temp_Obj);
		   });

		return ions;
	}
	/**
	 * Get all the checked C ternimus ions fro UI
	 */
	getCterminusCheckedList(){
		var ions = [];
		$.each($("input[name='cterminus']:checked"), function(){
				let id = $(this).attr( "id" );
				let value = $(this).val();
				let ionType = iontabledata.getActualIdvalues(id);
				let temp_Obj = {ionType:ionType,mass:value}
                ions.push(temp_Obj);
		   });

		return ions;
	}
}