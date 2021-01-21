/* THIS FILE DOES NOT SEEM TO BE USED
/**class MassShiftList {
    prsm_data;
    fixedPtmList;

    constructor(prsm_data, fixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464}, {name:"Carboxymethyl",acid:"C",mass:58.005479}]) {
        this.prsm_data = prsm_data;
        this.fixedPtmList = fixedPtmList;
    }

    getMassShiftList(){
        let unknownMassShiftList = this.getUnknownMassList();
        return this.getFixedPTMMassList(unknownMassShiftList);
    }

    getUnknownMassList()
	{
		let unknownMassShiftList = [];
		let l_prsm = this.prsm_data;
		if(l_prsm.prsm.annotated_protein.annotation.hasOwnProperty('mass_shift'))
		{
			let mass_shift = l_prsm.prsm.annotated_protein.annotation.mass_shift ;
				if(Array.isArray(mass_shift))
				{
					mass_shift.forEach(function(each_mass_shift){
						let position = parseInt(each_mass_shift.left_position) ;
						let mass = parseFloat(each_mass_shift.anno);
						unknownMassShiftList.push({"position":position,"mass":mass})
					})
				}
				else if(mass_shift.shift_type == "unexpected")
				{
					let position = parseInt(mass_shift.left_position);
					let mass = parseFloat(mass_shift.anno);
					unknownMassShiftList.push({"position":position,"mass":mass})
				}
		}
		return unknownMassShiftList;
    }
    
	getFixedPTMMassList(massShiftList){
		if(this.prsm_data.prsm.annotated_protein.annotation.hasOwnProperty("ptm") )
		{
			if(Array.isArray(this.prsm_data.prsm.annotated_protein.annotation.ptm))
			{
				this.prsm_data.prsm.annotated_protein.annotation.ptm.forEach(function(ptm){
					if(ptm.ptm_type == "Fixed")
					{
						let mass = this.getMassofFixedPtm(ptm.ptm.abbreviation);

						if(ptm.hasOwnProperty("occurence"))
						{
							if(Array.isArray(ptm.occurence))
							{
								ptm.occurence.forEach(function(occurence){
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
				if(this.prsm_data.prsm.annotated_protein.annotation.ptm.hasOwnProperty("occurence"))
				{
					if(this.prsm_data.prsm.annotated_protein.annotation.ptm.ptm_type == "Fixed")
					{
						let mass = this.getMassofFixedPtm(this.prsm_data.prsm.annotated_protein.annotation.ptm.ptm.abbreviation);
						if(Array.isArray(this.prsm_data.prsm.annotated_protein.annotation.ptm.occurence))
						{
							this.prsm_data.prsm.annotated_protein.annotation.ptm.occurence.forEach(function(occurence){
								let tempObj = {"position":occurence.left_pos,"mass":mass}
								massShiftList.push(tempObj);
							});
						}
						else
						{
							let tempObj = {"position":this.prsm_data.prsm.annotated_protein.annotation.ptm.occurence.left_pos,"mass":mass}
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
		for(let i=0; i<this.fixedPtmList.length; i++)
		{
			if(this.fixedPtmList[i].name === abbrevation)
			{
				return this.fixedPtmList[i].mass;
			}
		}
	}
}*/
