// @ts-nocheck
"use strict";
/**
 * Get the cleavage positions from the prsm data
 * @param {object} prsm - json obeject with complete prsm data
 */
function json2BreakPoints(prsm, firstPos) {
    let breakPoints = [];
    let dataBps = prsm.annotated_protein.annotation.cleavage;
    for (let i = 0; i < dataBps.length; i++) {
        let dataBp = dataBps[i];
        if (dataBp.exist_n_ion == 0 && dataBp.exist_c_ion == 0) {
            continue;
        }
        let bp = {};
        bp.position = dataBp.position;
        bp.existNIon = (dataBp.exist_n_ion == 1);
        bp.existCIon = (dataBp.exist_c_ion == 1);
        bp.anno = "";
        bp.masses = [];
        if (dataBp.matched_peaks != null) {
            let dataMasses = [];
            if (dataBp.matched_peaks.matched_peak.length > 1) {
                dataMasses = dataBp.matched_peaks.matched_peak;
            }
            else {
                dataMasses.push(dataBp.matched_peaks.matched_peak);
            }
            for (let j = 0; j < dataMasses.length; j++) {
                let dataMass = dataMasses[j];
                let mass = {};
                // Ion type
                mass.ionType = dataMass.ion_type;
                // Ion Display position
                mass.ionDispPos = parseInt(dataMass.ion_display_position);
                // Ion Charge
                mass.charge = parseInt(dataMass.peak_charge);
                // ion_position
                // mass.ionPos = parseInt(dataMass.ion_position);
                bp.masses.push(mass);
                if (bp.anno != "") {
                    bp.anno = bp.anno + " ";
                }
                bp.anno = bp.anno + mass.ionType + mass.ionDispPos + " " + mass.charge + "+";
            }
        }
        breakPoints.push(bp);
    }
    return breakPoints;
}
function getAminoAcidSequence(formFirstPos, formLastPos, residues) {
    let sequence = "";
    for (let i = formFirstPos; i <= formLastPos; i++) {
        sequence = sequence + residues[i].acid;
    }
    return sequence;
}
function getJsonList(item) {
    let valueList = [];
    if (Array.isArray(item)) {
        valueList = item;
    }
    else {
        valueList.push(item);
    }
    return valueList;
}
/**
 * Get occurence of fixed ptm positions
 * @param {object} prsm - json obeject with complete prsm data
 */
function json2Ptms(prsm) {
    let fixedPtmList = [];
    let protVarPtmList = [];
    let varPtmList = [];
    if (!prsm.annotated_protein.annotation.hasOwnProperty("ptm")) {
        return [fixedPtmList, protVarPtmList, varPtmList];
    }
    let dataPtmList = getJsonList(prsm.annotated_protein.annotation.ptm);
    for (let i = 0; i < dataPtmList.length; i++) {
        let dataPtm = dataPtmList[i];
        if (dataPtm.ptm_type == "Fixed" || dataPtm.ptm_type == "Protein variable"
            || dataPtm.ptm_type == "Variable") {
            if (dataPtm.hasOwnProperty("occurence")) {
                let occList = getJsonList(dataPtm.occurence);
                //console.log(occList);
                for (let j = 0; j < occList.length; j++) {
                    let occurence = occList[j];
                    let ptm = new Mod(occurence.anno, parseFloat(dataPtm.ptm.mono_mass), dataPtm.ptm.abbreviation);
                    let massShift = new MassShift(parseInt(occurence.left_pos), parseInt(occurence.right_pos), ptm.getShift(), dataPtm.ptm_type, ptm.getName(), ptm);
                    if (dataPtm.ptm_type == "Fixed") {
                        fixedPtmList.push(massShift);
                    }
                    else if (dataPtm.ptm_type == "Protein variable") {
                        protVarPtmList.push(massShift);
                    }
                    else {
                        varPtmList.push(massShift);
                    }
                }
            }
        }
    }
    return [fixedPtmList, protVarPtmList, varPtmList];
}
/**
 * Get left and right positions of background color and mass shift value
 * @param {object} prsm - json obeject with complete prsm data
 */
function json2MassShifts(prsm) {
    let massShifts = [];
    if (prsm.annotated_protein.annotation.hasOwnProperty('mass_shift')) {
        let dataMassShifts = getJsonList(prsm.annotated_protein.annotation.mass_shift);
        for (let i = 0; i < dataMassShifts.length; i++) {
            let dataShift = dataMassShifts[i];
            if (dataShift.shift_type == "unexpected" && dataShift.right_position != "0") {
                if (isNaN(parseFloat(dataShift.anno))) {
                    //then it is annotated with ptm name
                    let massShift = new MassShift(parseInt(dataShift.left_position), parseInt(dataShift.right_position), parseFloat(dataShift.shift), dataShift.shift_type, dataShift.anno);
                    massShifts.push(massShift);
                }
                else {
                    let massShift = new MassShift(parseInt(dataShift.left_position), parseInt(dataShift.right_position), parseFloat(dataShift.shift), dataShift.shift_type, dataShift.anno);
                    massShifts.push(massShift);
                }
            }
            else if (dataShift.right_position == 0) {
                console.error("Mass shift right position is 0!", dataShift);
            }
        }
    }
    return massShifts;
}
