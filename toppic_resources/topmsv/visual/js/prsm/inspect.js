"use strict";
function onclickTopView(e, prsmObj) {
    let topviewBtn = e.currentTarget;
    if (!topviewBtn) {
        console.error("ERROR: invalid event target");
        return;
    }
    let specID = topviewBtn.getAttribute('specid');
    let ms2Spec = prsmObj.getMs2Spectra();
    if (!ms2Spec) {
        console.error("ERROR: ms2 spectrum is empty");
        return;
    }
    let currentSpec = ms2Spec.find((spectrum) => {
        if (spectrum.getSpectrumId() == specID) {
            return spectrum;
        }
    });
    if (!currentSpec) {
        console.error("ERROR: ms2 spectrum is empty");
        return;
    }
    let massAndIntensityList = [];
    let peakAndIntensityList = getDataFromPRSMtoSpectralView(currentSpec, specID);
    let ionList = [];
    currentSpec.getNTerminalIon().forEach(ion => {
        ionList.push(ion.getName());
    });
    currentSpec.getCTerminalIon().forEach(ion => {
        ionList.push(ion.getName());
    });
    //console.log(ionList);
    massAndIntensityList = getMassAndIntensityData(currentSpec);
    //console.log(prsmGraph.data.proteoform);
    let proteoform = prsmObj.getProteoform();
    let sequence = proteoform.getSeq();
    //remove skipped residue
    sequence = sequence.slice(proteoform.getFirstPos());
    sequence = sequence.slice(0, proteoform.getLastPos() + 1 - proteoform.getFirstPos());
    let fixedPtmList = proteoform.getFixedPtm();
    //prsmGraph.data.proteoform.compMassShiftList();//recalculate mass shifts
    let unknownMassShiftList = proteoform.getUnknownMassShift();
    let protVarPtmsList = proteoform.getProtVarPtm();
    let variablePtmsList = proteoform.getVarPtm();
    let precursorMass = FormatUtil.formatFloat(currentSpec.getPrecMass(), "precMass");
    //if multiple variable ptm is in the same range, convert them to unknown mass shift to prevent incorrect annotation in the inspect page
    let varPtmListTmp = protVarPtmsList.concat(variablePtmsList);
    let varPtmFiltered = [];
    let protVarPtmFiltered = [];
    if (varPtmListTmp.length > 1) {
        varPtmListTmp.sort((x, y) => {
            return x.getLeftPos() - y.getLeftPos();
        });
        let i = 0;
        while (i < varPtmListTmp.length) {
            let k = i + 1;
            if (k >= varPtmListTmp.length) { //store data at [i] and finish since there is no more ptm to compare with
                if (varPtmListTmp[i].getType() == ModType.Variable) {
                    varPtmFiltered.push(varPtmListTmp[i]);
                }
                else if (varPtmListTmp[i].getType() == ModType.ProteinVariable) {
                    protVarPtmFiltered.push(varPtmListTmp[i]);
                }
            }
            else { //search for variable ptm with same range
                let totalShift = varPtmListTmp[i].getShift();
                while (varPtmListTmp[i].getLeftPos() == varPtmListTmp[k].getLeftPos() && varPtmListTmp[i].getRightPos() == varPtmListTmp[k].getRightPos()) {
                    totalShift = totalShift + varPtmListTmp[k].getShift();
                    k++;
                    if (k >= varPtmListTmp.length) {
                        break;
                    }
                }
                if (k > i + 1) { //mulitple variable ptm in the same range
                    unknownMassShiftList.push(new MassShift(varPtmListTmp[i].getLeftPos(), varPtmListTmp[i].getRightPos(), totalShift, "Unknown", totalShift.toString()));
                }
                else {
                    if (varPtmListTmp[i].getType() == ModType.Variable) {
                        varPtmFiltered.push(varPtmListTmp[i]);
                    }
                    else if (varPtmListTmp[i].getType() == ModType.ProteinVariable) {
                        protVarPtmFiltered.push(varPtmListTmp[i]);
                    }
                }
            }
            i = k;
        }
    }
    else {
        varPtmFiltered = variablePtmsList;
        protVarPtmFiltered = protVarPtmsList;
    }
    // Stores all the data in the variables respectively
    window.localStorage.setItem('peakAndIntensityList', JSON.stringify(peakAndIntensityList));
    window.localStorage.setItem('massAndIntensityList', JSON.stringify(massAndIntensityList));
    window.localStorage.setItem('ionType', ionList.toString());
    window.localStorage.setItem('sequence', JSON.stringify(sequence));
    window.localStorage.setItem('fixedPtmList', JSON.stringify(fixedPtmList));
    window.localStorage.setItem('protVarPtmsList', JSON.stringify(protVarPtmFiltered));
    window.localStorage.setItem('variablePtmsList', JSON.stringify(varPtmFiltered));
    window.localStorage.setItem('unknownMassShiftList', JSON.stringify(unknownMassShiftList));
    window.localStorage.setItem('precursorMass', precursorMass);
    //if some residues are going to be cut off in inspect page, adjust mod pos;
    if (proteoform.getFirstPos() > 0) {
        let newUnknownMassShifts = [];
        let newProtVarPtms = [];
        let newVarPtms = [];
        unknownMassShiftList.forEach((ptm) => {
            let newL = ptm.getLeftPos() - proteoform.getFirstPos();
            let newR = ptm.getRightPos() - proteoform.getFirstPos();
            let newPtm = new MassShift(newL, newR, ptm.getShift(), "unknown", ptm.getAnnotation());
            newPtm.setPtmList(ptm.getPtmList());
            newUnknownMassShifts.push(newPtm);
        });
        protVarPtmFiltered.forEach((ptm) => {
            let newL = ptm.getLeftPos() - proteoform.getFirstPos();
            let newR = ptm.getRightPos() - proteoform.getFirstPos();
            let newPtm = new MassShift(newL, newR, ptm.getShift(), "Protein variable", ptm.getAnnotation());
            newPtm.setPtmList(ptm.getPtmList());
            newProtVarPtms.push(newPtm);
        });
        varPtmFiltered.forEach((ptm) => {
            let newL = ptm.getLeftPos() - proteoform.getFirstPos();
            let newR = ptm.getRightPos() - proteoform.getFirstPos();
            let newPtm = new MassShift(newL, newR, ptm.getShift(), "Variable", ptm.getAnnotation());
            newPtm.setPtmList(ptm.getPtmList());
            newVarPtms.push(newPtm);
        });
        window.localStorage.setItem('protVarPtmsList', JSON.stringify(newProtVarPtms));
        window.localStorage.setItem('variablePtmsList', JSON.stringify(newVarPtms));
        window.localStorage.setItem('unknownMassShiftList', JSON.stringify(newUnknownMassShifts));
    }
    window.open("../inspect/spectrum.html");
}
/**
 * Get the peaklist from respective spectrum.js to set the data for inspect page
 * @param {object} ms2_data - json object with complete data spectrum for corresponding scan Id
 */
function getDataFromPRSMtoSpectralView(ms2Spec, specID) {
    let peakAndIntensity = [];
    if (!specID) {
        console.error("ERROR: spectrum id is invalid");
        return peakAndIntensity;
    }
    ms2Spec.getPeaks().forEach(peak => {
        let tempObj = peak.getMonoMz().toString() + " " + peak.getIntensity().toString();
        peakAndIntensity.push(tempObj);
    });
    return peakAndIntensity;
}
/**
 * Get the masslist from respective prsm.js to set the data for inspect page
 * @param {Integer} specId - Contians spec Id to get the data of corrsponding mass list
 */
function getMassAndIntensityData(ms2Spec) {
    let massAndIntensityList = [];
    let decovPeaks = ms2Spec.getDeconvPeaks();
    if (decovPeaks) {
        decovPeaks.forEach(peak => {
            let monoMass = peak.getMonoMass();
            let charge = peak.getCharge();
            if (monoMass && charge) {
                let tempObj = monoMass.toString() + " " + peak.getIntensity().toString() + " " + charge.toString();
                massAndIntensityList.push(tempObj);
            }
            else {
                console.error("Error: invalid mono mass or charge found in a peak");
            }
        });
    }
    return massAndIntensityList;
}
/**
 * Create HTML dropdown buttons based on the scan list
 * @param {Array} scanIdList - Contains Scan id numbers
 * @param {Array} specIdList - Contains Spec Id numbers
 */
function setDropDownItemsForInspectButton(scanIdList, specIdList) {
    let dropdown_menu = $(".dropdownscanlist .dropdown-menu");
    let len = scanIdList.length;
    for (let i = 0; i < len; i++) {
        let value = scanIdList[i];
        let specId = specIdList[i];
        let id = "scan_" + value;
        let a = document.createElement("a");
        a.setAttribute("class", "dropdown-item");
        a.setAttribute("href", "#!");
        a.setAttribute("id", id);
        a.setAttribute("value", value);
        a.setAttribute("specid", specId);
        a.innerHTML = "Scan " + value;
        dropdown_menu.append(a);
    }
}
/**
 * Onclick function, invoked on click of the inspect scn button
 */
function onClickToInspect(prsmObj) {
    $(".dropdownscanlist .dropdown-item ").click(function (e) {
        onclickTopView(e, prsmObj);
    });
}
