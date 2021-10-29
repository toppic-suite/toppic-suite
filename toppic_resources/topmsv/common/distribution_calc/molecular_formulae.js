"use strict";
/**
 * Code to calculate distribution based on molecular formula
 */
class MolecularFormulae {
    constructor() {
        this.averagine = "C4.9384H7.7583N1.3577O1.4773S0.0417";
        this.avrgMass = 111.0543055992;
        this.atomList = ["C", "H", "N", "O", "S"];
        this.carbonCnt = 4.9384;
        this.hydrogenCnt = 7.7583;
        this.nitrogenCnt = 1.3577;
        this.oxygenCnt = 1.4773;
        this.sulphurCnt = 0.0417;
        this.minintensity = 0.0001;
        this.protonMass = 1.00727647;
        this.eachAtomCount = 0;
        this.toleraceMassDiff = 0.01;
        this.mininte = 1;
        this.intensityTolerance = 1;
    }
    /**
     * Function to calculate the emass distrubution fo the given sequence
     * @param {Double} mass - Contains the mass at that position of the sequence
     * @param {Array} peakDataList - Contains the peak list provided by the user
     * @param {Float} charge - Contains the chrage of the ion
     */
    emass(mass, charge, modifiedPeakList) {
        //Give the count of each atom in molecule
        this.eachAtomCount = mass / this.avrgMass;
        let calculatedMass;
        let atomCountList;
        [atomCountList, calculatedMass] = this.getMolecularFormula();
        let massError = mass - calculatedMass;
        let atomListlen = atomCountList.length;
        let totDistributionList = [];
        for (let i = 0; i < atomListlen; i++) {
            for (let j = 0; j < atomCountList[i].count; j++) {
                let atomdist = getIsotopicMassRef(atomCountList[i].atom);
                totDistributionList = this.getMassAndIntensity(totDistributionList, atomdist);
            }
            for (let k = 0; k < totDistributionList.length; k++) {
                let IsotopicMass = getIsotopicMassOfAtom(atomCountList[i].atom);
                totDistributionList[k].mass = totDistributionList[k].mass + (IsotopicMass[0].mass * atomCountList[i].count);
            }
        }
        totDistributionList = this.getMZwithHighInte(totDistributionList, charge, massError, modifiedPeakList);
        let finalDistList = this.getNormalizedIntensityAndAdjustedEnvelopes(totDistributionList, modifiedPeakList);
        return finalDistList;
    }
    /**
     * Logic to calculate distribution
     * @param {Array} totDistributionList - Array with current total distribution
     * @param {Array} aminoAcidDist - Array with existing calculated distribution of amino acid
     */
    getMassAndIntensity(totDistributionList, atomdist) {
        let maxintensity = 0;
        if (totDistributionList.length == 0) {
            return atomdist;
        }
        else {
            let len = totDistributionList.length + atomdist.length - 1;
            let completeDistributionList = new Array(len).fill({ mass: 0, intensity: 0 });
            for (let i = 0; i < totDistributionList.length; i++) {
                for (let j = 0; j < atomdist.length; j++) {
                    let intensity = 0;
                    let index = i + j;
                    let mass = totDistributionList[i].mass + atomdist[j].mass;
                    intensity = totDistributionList[i].intensity * atomdist[j].intensity;
                    if (completeDistributionList[index].mass != 0)
                        intensity = intensity + completeDistributionList[index].intensity;
                    completeDistributionList[index] = { mass: mass, intensity: intensity };
                    if (intensity > maxintensity)
                        maxintensity = intensity;
                }
            }
            let completeDistributionList_temp = [];
            for (let i = 0; i < completeDistributionList.length; i++) {
                let tempintensityVal = (completeDistributionList[i].intensity / maxintensity) * 100;
                if (tempintensityVal > this.minintensity) {
                    completeDistributionList[i].intensity = Math.round(tempintensityVal * 1000000) / 1000000; //tempintensityVal//Math.round(tempintensityVal * 1000000) / 1000000; //parseFloat(((tempDistributionList[i].intensity/maxintensity)*100).toFixed(6)) ;
                    completeDistributionList_temp.push(completeDistributionList[i]);
                }
            }
            completeDistributionList = completeDistributionList_temp;
            return completeDistributionList;
        }
    }
    /**
     * Code to get the molecular formular based on the mass.
     * Number of atoms of each atom in the averagine * number of
     * atoms of each atom calculated from mass gives the total number of atoms.
     */
    getMolecularFormula() {
        let atomListwithCnt = new Array(5);
        let len = this.atomList.length;
        let mass = 0;
        for (let i = 0; i < len; i++) {
            let count = 0;
            if (this.atomList[i] == "C") {
                count = Math.floor(this.eachAtomCount * this.carbonCnt);
                mass = mass + count * getIsotopicMassOfAtom(this.atomList[i])[0].mass;
            }
            else if (this.atomList[i] == "H") {
                count = Math.floor(this.eachAtomCount * this.hydrogenCnt);
                mass = mass + count * getIsotopicMassOfAtom(this.atomList[i])[0].mass;
            }
            else if (this.atomList[i] == "N") {
                count = Math.floor(this.eachAtomCount * this.nitrogenCnt);
                mass = mass + count * getIsotopicMassOfAtom(this.atomList[i])[0].mass;
            }
            else if (this.atomList[i] == "O") {
                count = Math.floor(this.eachAtomCount * this.oxygenCnt);
                mass = mass + count * getIsotopicMassOfAtom(this.atomList[i])[0].mass;
            }
            else if (this.atomList[i] == "S") {
                count = Math.floor(this.eachAtomCount * this.sulphurCnt);
                mass = mass + count * getIsotopicMassOfAtom(this.atomList[i])[0].mass;
            }
            atomListwithCnt[i] = { atom: this.atomList[i], count: count };
        }
        return [atomListwithCnt, mass];
    }
    /**
     * Code to calculate MZ(mass/charge) and remove low intensity values
     * @param {Array} totDistributionList - Array with total distribution
     * @param {Float} charge - Float from mass list
     * @param {Float} massError - Calculated massError needed to be added
     * as we taken int values of number of atoms eleminating mass from the float values.
     * @param {Array} peakDataList - Array of peak list from user
     */
    getMZwithHighInte(totDistributionList, charge, massError, modifiedPeakList) {
        let len = totDistributionList.length;
        let newtotDistributionList = [];
        let overaAllMaxIntensity = 0;
        let onePerctInte = 0;
        let peakListLen = modifiedPeakList.length;
        for (let k = 0; k < peakListLen; k++) {
            if (modifiedPeakList[k].getIntensity() > overaAllMaxIntensity)
                overaAllMaxIntensity = modifiedPeakList[k].getIntensity();
        }
        onePerctInte = overaAllMaxIntensity / 100;
        for (let i = 0; i < len; i++) {
            let intensity = totDistributionList[i].intensity;
            intensity = overaAllMaxIntensity * intensity / 100;
            if (intensity > onePerctInte) {
                let mz = (totDistributionList[i].mass + massError) / charge + this.protonMass;
                let intensity = totDistributionList[i].intensity;
                //Converting mass variable to mz(mass/charge)
                let tempdistObj = { mass: mz, intensity: intensity };
                newtotDistributionList.push(tempdistObj);
            }
        }
        return newtotDistributionList;
    }
    aminoDistToPeaks(totDistributionList) {
        let newDistList = [];
        totDistributionList.forEach((peak, i) => {
            newDistList.push(new Peak(i.toString(), peak.mass, peak.mass, peak.intensity));
        });
        return newDistList;
    }
    //modification of getNormalizedIntensity to add envelope y pos adjustment 
    getNormalizedIntensityAndAdjustedEnvelopes(totDistributionList, modifiedPeakList) {
        let len = totDistributionList.length;
        let peakListLen = modifiedPeakList.length;
        let count = 0;
        let maxinte = 0;
        let mininte = 100;
        let peakMaxinte = 0;
        let peakMininte = 10000000;
        let maxMz = 0;
        let minMz = 10000000;
        let matchedPeakList = [];
        for (let i = 0; i < len; i++) { //iterating through actual peaks in this envelope
            for (let j = 0; j < peakListLen; j++) { //iterating through theo peaks in the data 	
                let mzDifference = Math.abs(totDistributionList[i].mass - modifiedPeakList[j].getPos());
                if (mzDifference <= this.toleraceMassDiff) {
                    if (maxMz < totDistributionList[i].mass) {
                        maxMz = totDistributionList[i].mass;
                    }
                    if (minMz > totDistributionList[i].mass) {
                        minMz = totDistributionList[i].mass;
                    }
                    matchedPeakList.push([i, j]); //i is env index, j is peak index
                    count = count + 1;
                }
            }
        }
        maxMz = maxMz + this.toleraceMassDiff;
        minMz = minMz - this.toleraceMassDiff;
        /*previous function skews the result if there are > 1 peaks in the mz range and later one has high intensity
make sure the max min value changed when it is the peak that is closest to the given env
for now, will check if the peak has any envelopes within error tolerance
*/
        for (let j = 0; j < peakListLen; j++) {
            if (modifiedPeakList[j].getPos() >= minMz && modifiedPeakList[j].getPos() <= maxMz) {
                if (peakMaxinte < modifiedPeakList[j].getIntensity()) {
                    //for all env dots, check if this peak really belongs to this envelope
                    for (let env = 0; env < totDistributionList.length; env++) {
                        if (Math.abs(totDistributionList[env].mass - modifiedPeakList[j].getPos()) <= this.toleraceMassDiff) {
                            peakMaxinte = modifiedPeakList[j].getIntensity();
                        }
                    }
                }
                if (peakMininte > modifiedPeakList[j].getIntensity()) {
                    for (let env = 0; env < totDistributionList.length; env++) {
                        if (Math.abs(totDistributionList[env].mass - modifiedPeakList[j].getPos()) <= this.toleraceMassDiff) {
                            peakMininte = modifiedPeakList[j].getIntensity();
                        }
                    }
                }
            }
        }
        for (let i = 0; i < len; i++) {
            if (minMz <= totDistributionList[i].mass && totDistributionList[i].mass <= maxMz) {
                if (maxinte < totDistributionList[i].intensity) {
                    maxinte = totDistributionList[i].intensity;
                }
                if (mininte > totDistributionList[i].intensity) {
                    mininte = totDistributionList[i].intensity;
                }
            }
        }
        /*when calculating new y pos, when a later envelope has higher y pos, evaluate again after reducing peak inte
        based on that higher envelope, not in the previous left-to-right order.
        
        keep marking peak as shared peak as moving from left to right, and when an envelope meets a shared peak,
        check if the previous envelope with which it shares the peak had higher intensity at that peak.
        In order to compare, save the original intensity value in each envelope property for envelopes with shared
        peak. Then compare, and rewrite the previous y pos if needed. */
        if (count != 0) {
            let avg;
            let distributionAvgInte;
            avg = (peakMaxinte + peakMininte) / 2;
            distributionAvgInte = (maxinte + mininte) / 2;
            for (let i = 0; i < totDistributionList.length; i++) {
                if (avg == 0) {
                    avg = 1;
                }
                totDistributionList[i].intensity = (avg * totDistributionList[i].intensity) / distributionAvgInte;
                if (totDistributionList[i].intensity < 0 || distributionAvgInte <= 0) {
                    totDistributionList[i].intensity = 0;
                }
                ;
                for (let j = 0; j < matchedPeakList.length; j++) {
                    if (matchedPeakList[j][0] == i) {
                        let p = matchedPeakList[j][1];
                        //store original peak intensity value only when it is evaluted for the first time
                        //not going to run when called from peakData because it is already evaluated ther
                        let newInte = modifiedPeakList[p].getIntensity() - totDistributionList[i].intensity;
                        modifiedPeakList[p].setIntensity(newInte);
                        if (modifiedPeakList[p].getIntensity() < 0) {
                            modifiedPeakList[p].setIntensity(0);
                        }
                        //modifiedPeakList[p]["isPeakShared"] = true;
                    }
                }
            }
        }
        //convert aminoDist to peaks
        let newDistList = this.aminoDistToPeaks(totDistributionList);
        //remove envelopes without matching peaks
        this.removeEnvelopes(newDistList, matchedPeakList);
        return newDistList;
    }
    /**
     * Code to normalize the Intensity.
     * Take the average of intensity from the peaks entered by the user.
     * Take the average of the calculated distribution for each Array element in the Array.
     * Make both of them equal and calculating the rest of the
     * distribution intensity based on the avg value from the peak list.
     * @param {Array} totDistributionList - Total distribution calculated
     * @param {Array} peakDataList - Peak Data entered by the user
     *
     * below is original function by Sreekanth
     */
    /*getNormalizedIntensity(totDistributionList: AminoDist[], peakDataList: Peak[])
    {
        let len: number = totDistributionList.length;
        let peakListLen: number = peakDataList.length;
        let intensity: number = 0;
        let count: number = 0 ;
        let distributionInte: number = 0;
        let maxinte: number = 0;
        let mininte: number = 100;

        let peakMaxinte: number = 0;
        let peakMininte: number = 10000000;

        let maxMz: number = 0;
        let minMz: number = 10000000;
        for(let i=0;i<len;i++)
        {
            for(let j=0;j<peakListLen;j++)
            {
                if((totDistributionList[i].mass - peakDataList[j].mz) <= this.toleraceMassDiff
                    && (totDistributionList[i].mass - peakDataList[j].mz) >= 0-this.toleraceMassDiff )
                {
                    if(maxMz < totDistributionList[i].mz){
                        maxMz = totDistributionList[i].mz;
                    }
                    if(minMz > totDistributionList[i].mz){
                        minMz = totDistributionList[i].mz;
                    }
                    count = count + 1;
                }
            }
        }
        maxMz = maxMz + this.toleraceMassDiff;
        minMz = minMz - this.toleraceMassDiff;
        for(let i=0;i<len;i++)
        {
            if(minMz <= totDistributionList[i].mz &&  totDistributionList[i].mz <= maxMz)
            {
                if(maxinte < totDistributionList[i].intensity){
                    maxinte = totDistributionList[i].intensity;
                }
                if(mininte > totDistributionList[i].intensity){
                    mininte = totDistributionList[i].intensity;
                }
            }
        }
        for(let j=0;j<peakListLen;j++)
        {
            if(peakDataList[j].mz >= minMz && peakDataList[j].mz <= maxMz)
            {
                if(peakMaxinte < peakDataList[j].intensity){
                    peakMaxinte = peakDataList[j].intensity;
                }
                if(peakMininte > peakDataList[j].intensity){
                    peakMininte = peakDataList[j].intensity;
                }
            }
        }
        if(count !=0 )
        {
            let avg ;
            let distributionAvgInte;

            avg = (peakMaxinte + peakMininte)/2;
            distributionAvgInte = (maxinte+mininte)/2;
            for(let i=0;i<len;i++)
            {
                totDistributionList[i].intensity = (avg * totDistributionList[i].intensity)/distributionAvgInte ;
            }
        }
        return totDistributionList ;
    }*/
    removeEnvelopes(envList, matchedPeakList) {
        let threshold = 0.95;
        let totalInte = 0;
        let remainInte = 0;
        /*keep removing dots with no matching peaks while
        (remainig peaks intensity/all peaks intensity) >= threshold
        remove from the right, but stop when prev dot is not removed -- no remove in the middle only*/
        for (let i = 0; i < envList.length; i++) {
            totalInte = totalInte + envList[i].getIntensity();
        }
        remainInte = totalInte;
        let start = 0;
        let end = envList.length - 1;
        let removeEndEnv = true;
        let removeStartEnv = true;
        let endEnvInte = Number.MAX_VALUE;
        let startEnvInte = Number.MAX_VALUE;
        ;
        /*removing dots from both side: start evaluating from start and end index.
        while end >= start, and keepRemove is true (at least one of the dots could be removed)
        keep removing the envelopes from the both sides in turn, removing whichever that has lower inte.
        
        first, iterate matchedPeakList to see if end index is there. If exists, no more removal from end.
        If not exists, remove the env, end = end - 1. The matchedPeakList entry can also be removed.
        Same with start index.

        The while loop exists when no removal from both end or no more env left to evaluate.
        
        */
        while (end >= start && (removeEndEnv || removeStartEnv)) {
            //keep removing peaks at idx if it is not in the matchedpeaklist
            for (let p = matchedPeakList.length - 1; p >= 0; p--) {
                if (matchedPeakList[p][0] == end) { //no removal
                    removeEndEnv = false;
                }
                else if (matchedPeakList[p][0] == start) { //no removal
                    removeStartEnv = false;
                }
            }
            if (removeEndEnv) {
                remainInte = remainInte - envList[end].getIntensity();
                if (remainInte / totalInte >= threshold) {
                    endEnvInte = envList[end].getIntensity();
                }
                else {
                    removeEndEnv = false;
                }
            }
            if (removeStartEnv) {
                remainInte = remainInte - envList[start].getIntensity();
                if (remainInte / totalInte >= threshold) {
                    startEnvInte = envList[start].getIntensity();
                }
                else {
                    removeStartEnv = false;
                }
            }
            //decide which one has lower intensity --and to be removed
            if (endEnvInte < startEnvInte && removeEndEnv) {
                envList.splice(end, 1);
                end--;
            }
            else if (startEnvInte <= endEnvInte && removeStartEnv) {
                envList.splice(start, 1);
                end--;
            }
        }
        return envList;
    }
}
