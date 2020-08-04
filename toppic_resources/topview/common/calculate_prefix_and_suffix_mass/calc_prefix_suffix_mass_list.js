class CalcPrefixSuffixMassList {
    sequence;
    massShiftList;

    constructor(sequence, massShiftList) {
        this.sequence = sequence;
        this.massShiftList = massShiftList;
    }

    getPrefixMassList(ionMassShift) {
        if (this.sequence) {
            let prefixMassList = new Array(this.sequence.length);
            for (let i = 0; i < this.sequence.length; i++) {
                let acidMass = getAminoAcidDistribution(this.sequence[i])[0].mass;

                if(!acidMass) {
                    return undefined;
                }

                if (i === 0) {
                    if (this.massShiftList) {
                        acidMass = this.addMassShift(i, acidMass);
                    }
                    let tempObj = {acid: this.sequence[i], position: (i+1), mass: acidMass};
                    prefixMassList[i] = tempObj;
                } else {
                    let mass = prefixMassList[i-1].mass + acidMass;
                    if (this.massShiftList) {
                        mass = this.addMassShift(i, mass);
                    }
                    let tempObj = {acid: this.sequence[i], position: (i+1), mass: mass};
                    prefixMassList[i] = tempObj;
                }
            }

            if (ionMassShift) {
                for (let i = 0; i < this.sequence.length; i++) {
                    prefixMassList[i].mass += ionMassShift;
                }
            }
            return prefixMassList;
        } else {
            return [];
        }
    }

    getSuffixMassList(ionMassShift) {
        if (this.sequence) {
            let suffixMassList = new Array(this.sequence.length);
            for (let i = this.sequence.length - 1; i >= 0; i--) {
                let acidMass = getAminoAcidDistribution(this.sequence[i])[0].mass;

                if(!acidMass) {
                    return undefined;
                }

                if (i === this.sequence.length - 1) {
                    acidMass = this.addMassShift(i, acidMass);
                    let tempObj = {acid: this.sequence[i], position: (i+1), mass: acidMass};
                    suffixMassList[0] = tempObj;
                } else {
                    let position = this.sequence.length - 1 - i;
                    let mass = suffixMassList[position - 1].mass + acidMass;
                    mass = this.addMassShift(i, mass);
                    let tempObj = {acid: this.sequence[i], position: (i+1), mass: mass};
                    suffixMassList[position] = tempObj;
                }
            }
            if (ionMassShift) {
                for (let i = 0; i < this.sequence.length; i++) {
                    suffixMassList[i].mass += ionMassShift;
                }
            }
            return suffixMassList;
        }else {
            return [];
        }  
    }

	/**
	 * Add mass shifts to massShift list
	 * @param {int} position - position of the mass list
	 * @param {Float} mass - Mass shift to be added to the list
	 */
	addMassShift(position, mass){
		for(let i=0; i< this.massShiftList.length; i++)
		{
			if(position === this.massShiftList[i].position)
			{
				mass = mass + this.massShiftList[i].mass ;
				return mass ;
			}
		}
		return mass ;
	}
}