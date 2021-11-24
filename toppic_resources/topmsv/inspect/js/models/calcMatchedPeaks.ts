/**
 * Class to calculate Matched Peaks and distribution of the Amoino Acids
 */
class CalcMatchedPeaks {
	PREFIX: string = "PREFIX";
	SUFFIX: string = "SUFFIX";
	matchedUnMatchedList: MatchedUnMatchedPeak[] = [];
	matchedList: MatchedUnMatchedPeak[] = [];
	water: string = "H2O";
	million: number = 1000000;

	CONST_A: string = "a";
	CONST_B: string = "b";
	CONST_C: string = "c";
	CONST_Y: string = "Y";
	colors: string[] = ["red","orange","blue","grey"];
	massShiftList: MassShift[] = [];

	/**
	 * Function to form object with necessary attributes at each amino acid position
	 * @param {Object} monoMassList_temp - Single element from mono Mass List
	 * @param {String} id - Unique Id for each peak
	 * @param {String} ion - Contains matched type(b/y) appended with charge
	 * @param {Integer} position - Position of the amino acid
	 * @param {Float} massDiff - Mass difference at particular position
	 * @param {Float} mass - Mass at the position of the acid 
	 * @param {Char} matchedInd - Indicator indicating matched/not matched. If matched containd "Y"
	 */
	matchedPeakAttributes(monoMassList_temp: Peak, id: string, matchedInd: string, ion: string, ionPos: number,
		position: number, massDiff: number, mass: number): MatchedUnMatchedPeak {
		let peak: Peak = monoMassList_temp;
		let matchedPeak = {} as MatchedUnMatchedPeak;
		let thMass: number | undefined;
		let PPMerror: number | undefined;
        let charge: number | undefined = peak.getCharge();
        let monoMass: number | undefined = peak.getMonoMass();
        if (ion && ionPos && position && massDiff && mass) {
            // Check if Matched Indicator is Yes "Y" else keep empty spaces in the new attributes
            if (matchedInd == this.CONST_Y) {
                thMass = Math.round(mass * 10000) / 10000;
                massDiff = Math.round(massDiff * 10000) / 10000;
                PPMerror = massDiff / thMass * this.million;
                PPMerror = Math.round(PPMerror * 10000) / 10000;
            }
            matchedPeak.ion = ion;
            matchedPeak.ionPos = ionPos.toString();
            matchedPeak.position = position;
            matchedPeak.massError = massDiff;
        }
        if (charge) {
            matchedPeak.charge = charge;
        }
        if (monoMass) {
            matchedPeak.mass = monoMass;
        }
        if (PPMerror || PPMerror == 0) {
            matchedPeak.PPMerror = PPMerror;
        }
        if (thMass) {
            matchedPeak.thMass = thMass;
        }
        matchedPeak.peakId = id;
        matchedPeak.matchedInd = matchedInd;
        matchedPeak.intensity = peak.getIntensity();
        // console.log("peak : ", peak);
        return matchedPeak;
	}
	/**
	 * Returns Complete list with matched Indication
	 * @param {Array} prefixOrSuffixMassList - List with either prefix or suffix mass list
	 * @param {Array} monoMassList - Contains Mono Mass list from html page
	 * @param {String} sequence - Protein sequence
	 * @param {Float} massErrorthVal - Threshhold value entered by user to determine the mass error and matched peaks
	 * @param {Float} ppmErrorthVal - Threshhold value entered by user to determine the mass error and matched peaks in ppm units
	 * @param {String} ionObj - Ion object containing name, terminal, mass information
	 */
	getMatchedPeakList(prefixOrSuffixMassList: TheoMass[], monoMassList: Peak[],sequence: string, massErrorthVal: number, ppmErrorthVal: number, ionObj: Ion)
	{
        let matchedList: MatchedUnMatchedPeak[] = [];
        let monoMassList_temp: Peak[] = monoMassList.slice();
        let len: number = monoMassList_temp.length;
        let preOrSufListln: number = prefixOrSuffixMassList.length;
        let seqln: number = sequence.length;
        //console.log("monoMassList_temp", monoMassList_temp);
        for (let i: number = 0; i < len; i++) {
            let peakId: number = i + 1;
            for (let j: number = 0; j < preOrSufListln; j++) {
                let peakMass: number | undefined = monoMassList_temp[i].getMonoMass();
                if (!peakMass) {
                    console.log(monoMassList_temp[i]);
                    console.error("ERROR: peak has an invalid mass");
                    return matchedList;
                }
                let massDiff: number = peakMass - prefixOrSuffixMassList[j].getMass();
                if (massErrorthVal != null || isNaN(massErrorthVal)) {
                    if (Math.abs(massDiff) <= massErrorthVal) {
                        //increment the position as the seq starts from 0 but the data in table starts from 1
                        let ion: string = ionObj.getName() + j.toString();
                        let position: number = j;
                        if (ionObj.getTerminal() == "C") {
                            position = seqln - position;
                        }
                        let mass: number = prefixOrSuffixMassList[j].getMass();
                        prefixOrSuffixMassList[j].setIon(ionObj);
                        let matchedInd: string = "Y";
                        let peak: MatchedUnMatchedPeak = this.matchedPeakAttributes(monoMassList_temp[i], peakId.toString(), matchedInd, ion, j, position, massDiff, mass);
                        //console.log(peak)		
                        matchedList.push(peak);
                    }
                }
                else {
                    let prefDiff: number = Math.round(massDiff * 10000) / 10000;
                    let prefMassRounded: number = Math.round(prefixOrSuffixMassList[j].getMass() * 10000) / 10000;
                    let prePPMerror: number = prefDiff / prefMassRounded * this.million;
                    if (Math.abs(prePPMerror) <= ppmErrorthVal) {
                        //increment the position as the seq starts from 0 but the data in table starts from 1
                        let ion: string = ionObj.getName() + j;
                        let position: number = j;
                        if (ionObj.getTerminal() == "C") {
                            position = seqln - position;
                        }
                        let mass: number = prefixOrSuffixMassList[j].getMass();
                        prefixOrSuffixMassList[j].setIon(ionObj);
                        let matchedInd: string = "Y";
                        let peak: MatchedUnMatchedPeak = this.matchedPeakAttributes(monoMassList_temp[i], peakId.toString(), matchedInd, ion, j, position, massDiff, mass);
                        matchedList.push(peak);
                    }
                }
            }
        }
        return matchedList;
	}
    /**
     * Set matched indicator to complete data and Seperate matched and unmatched list of data with matchedInd attribute
     * @param {Array} monoMassList - Contains mono Mass list data
     * @param {Array} matchedList - Contains all the data of both matched and unmatched data
     */
	 getMatchedAndUnMatchedList(monoMassList: Peak[], matchedList: MatchedUnMatchedPeak[]) {
        let matchedAndUnmatchedList: MatchedUnMatchedPeak[] = matchedList.map(x => (Object.assign({}, x)));
        //matchedAndUnmatchedList.concat(matchedList);
        let len: number = matchedList.length;
        // let MatchedPeaksObj = new MatchedPeaks();
        let self: any = this;
        monoMassList.forEach(function (eachElement, i) {
            let matched: boolean = false;
            let peakId: number = i + 1;
            for (let j: number = 0; j < len; j++) {
				let id: string | undefined = matchedList[j].peakId;
                if (id) {
					if (peakId == parseInt(id)) {
						matched = true;
						break;
					}	
				}
            }
            if (!matched) {
                let matchedInd: string = "N";
                let peak: MatchedUnMatchedPeak = self.matchedPeakAttributes(monoMassList[i], peakId.toString(), matchedInd);
                matchedAndUnmatchedList.push(peak);
            }
        });
        // completeCalData.matchedandunmatcheddata = matchedAndUnmatchedList;
        return matchedAndUnmatchedList;
    }
    /**
     * Get the envelopes/distribution for the sequence
     * @param {Array} peakDataList - Contains peak list information
     * @param {String} sequence - Contains protein sequence
     * @param {Array} matchedUnMatchedList - Contains all matched and unmatched list
     */
    getDistribution(peakDataList: Peak[], matchedUnMatchedList: MatchedUnMatchedPeak[]) {
        let len: number = matchedUnMatchedList.length;
        let totalEnvelopes: Envelope[] = [];
        let molecularFormObj: MolecularFormulae = new MolecularFormulae();
        //matchedUnmatchedlist sort by intensity
        matchedUnMatchedList.sort(function (x, y) {
            return d3.descending(x.intensity, y.intensity);
        });
        for (let i: number = 0; i < len; i++) {
            let envelope: Envelope = new Envelope(matchedUnMatchedList[i].mass, matchedUnMatchedList[i].charge, matchedUnMatchedList[i].intensity);
            let theoPeaks: Peak[] = molecularFormObj.emass(envelope.getMonoMass(), envelope.getCharge(), peakDataList);
            for (let j: number = 0; j < theoPeaks.length; j++) {
                let peak: Peak = new Peak(j.toString(), theoPeaks[j].getPos(), theoPeaks[j].getMonoMz(), theoPeaks[j].getIntensity());
                envelope.addPeaks(peak);
            }
            //check if the envelope for the current peak has already been generated
            //push to the totalDistribution only if this envelope is unique
            if (totalEnvelopes.length < 1) {
                totalEnvelopes.push(envelope);
            }
            else {
                if (envelope.getMonoMass() != totalEnvelopes[totalEnvelopes.length - 1].getMonoMass()) {
                    totalEnvelopes.push(envelope);
                }
            }
            /*if(matchedUnMatchedList[i].matchedInd == "Y")
            {
                if(this.CONST_A == matchedUnMatchedList[i].ion[0] || this.CONST_B == matchedUnMatchedList[i].ion[0]
                                                                        || this.CONST_C == matchedUnMatchedList[i].ion[0])
                {
                    let matchedPos = matchedUnMatchedList[i].position ;
                    let seq = sequence.slice(0,matchedPos) ;
                    /*compare with completeMassShiftList to see if this seq includes mass shift
                    * compare matchedPos with positons in completeMassShiftList.
                    * if matchedPos is bigger, there is a mass shift inside seq.
                    * Send the position to emass so that the mass shift is reflected in the toDistribution list after the acid*/
            /*let massShiftList = [];
            for (let i = 0; i < completeMassShiftList.length; i++){
                if (matchedPos >= completeMassShiftList[i].position){
                    massShiftList.push(completeMassShiftList[i]);
                }
            }
            distributionList.mono_mass = matchedUnMatchedList[i].mass;
            distributionList.charge = matchedUnMatchedList[i].charge;
            distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.PREFIX, massShiftList);
            totalDistribution.push(distributionList);

        }
        else
        {
            let matchedPos = matchedUnMatchedList[i].position ;
            let seq = sequence.slice(matchedPos,seqln) ;

            let massShiftList = [];
            for (let i = 0; i < completeMassShiftList.length; i++){
                if (matchedPos <= completeMassShiftList[i].position){
                    //as the seq is a slice of original sequence, mass shift position should be adjusted
                    let massData = {};
                    massData["position"] = completeMassShiftList[i].position - matchedPos;
                    massData["mass"] = completeMassShiftList[i].mass;
                    massShiftList.push(massData);
                }
            }
            distributionList.mono_mass = matchedUnMatchedList[i].mass;
            distributionList.charge = matchedUnMatchedList[i].charge;
            distributionList.env_peaks = calEmassAndDistObj.emass(seq,peakDataList,matchedUnMatchedList[i].charge,this.SUFFIX, massShiftList);
            totalDistribution.push(distributionList);
        }
    }
    else
    {
        distributionList.mono_mass = matchedUnMatchedList[i].mass;
        distributionList.charge = matchedUnMatchedList[i].charge;
        distributionList.env_peaks = molecularFormObj.emass(distributionList.mono_mass,distributionList.charge,peakDataList);
        totalDistribution.push(distributionList);
    }*/
        }
        if (totalEnvelopes.length !== 0) {
            totalEnvelopes.sort(function (x, y) {
                return d3.ascending(x.getPeaks()[0].getPos(), y.getPeaks()[0].getPos());
            });
        }
        // let envlength = totalDistribution.length;
        // let colorListsize = this.colors.length;
        // while(envlength--){
        // 	let index = envlength%colorListsize ;
        // 	totalDistribution[envlength].color = this.colors[index] ;
        // }
        return totalEnvelopes;
    }
    /**
     * Get All the matched positions and mass list seperated with a matched indicator
     * @param {Array} prefixOrSuffixMassList - List of prefix and suffix masses
     * @param {Array} monoMassList - Contains Mono Mass list
     * @param {Float} massErrorthVal - Mass Error Threshold value entered by user
     * @param {Float} ppmErrorthVal - PPM Error Threshold value entered by user
     * @param {String} prefixInd - contains if the complete data is of prefix/suffix mass list
     */
    getMatchedAndUnmatchedPrefixAndSuffixMassList(prefixOrSuffixMassList: TheoMass[], monoMassList: Peak[], massErrorthVal: number, ppmErrorthVal: number, prefixInd: string): MatchedUnMatchedPeakSimplified[] {
        // console.log("monoMassList : ", monoMassList);
        let MatchedAndUnMatchedList: MatchedUnMatchedPeakSimplified[] = [];
        let monoMassList_temp: Peak[] = monoMassList.slice();
        let len: number = monoMassList_temp.length;
        let preOrSufListln: number = prefixOrSuffixMassList.length;
        let MatchedAndUnMatchedListObj: MatchedUnMatchedPeakSimplified = {} as MatchedUnMatchedPeakSimplified;
        // console.log("monoMassList_temp:", monoMassList_temp);
        for (let j: number = 0; j < preOrSufListln; j++) {
            let position: number = j + 1;
            let charge: number = -1;
            if (prefixInd !== "prefix") {
                position = preOrSufListln - position + 1;
            }
            let mass: number = prefixOrSuffixMassList[j].getMass();
            let matchedInd: string = "N";
            for (let i: number = 0; i < len; i++) {
                let monoMass: number | undefined = monoMassList_temp[i].getMonoMass();
                if (!monoMass) {
                    console.error("ERROR: invalid mono mass of a peak");
                    return MatchedAndUnMatchedList;
                }
                let massDiff: number = monoMass - prefixOrSuffixMassList[j].getMass();
                let prefDiff: number = Math.round(massDiff * 10000) / 10000;
                let prefMassRounded: number = Math.round(prefixOrSuffixMassList[j].getMass() * 10000) / 10000;
                let prePPMerror: number = prefDiff / prefMassRounded * this.million;
                if ((massErrorthVal != null || isNaN(massErrorthVal)) ||
                    (ppmErrorthVal != null || isNaN(ppmErrorthVal))) {
                    if ((Math.abs(massDiff) <= massErrorthVal) || (Math.abs(prePPMerror) <= ppmErrorthVal)) {
                        //increment the position as the seq starts from 0 but the data in table starts from 1
                        position = j + 1;
                        if (prefixInd !== "prefix") {
                            position = preOrSufListln - position + 1;
                        }
                        mass = prefixOrSuffixMassList[j].getMass();
                        let tmp: number | undefined = monoMassList_temp[i].getCharge();
                        if (tmp) {
                            charge = tmp;
                        }
                        matchedInd = "Y";
                        break;
                    }
                }
            }
            MatchedAndUnMatchedListObj = { position: position, mass: mass, matchedInd: matchedInd, charge: charge };
            MatchedAndUnMatchedList.push(MatchedAndUnMatchedListObj);
        }
        return MatchedAndUnMatchedList;
    }
}
