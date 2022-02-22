"use strict";
/**
 * Create title and navigation urls to "all proteins","protein" and "proteoform" of the prsm
 * @param {String} folderpath - provides folder path to the data and helps in building urls
 */
function BuildUrl(folderpath, prsmObj) {
    let proteoformObj = prsmObj.getProteoform();
    let specIds = [];
    let ms2Spec = prsmObj.getMs2Spectra();
    if (ms2Spec) {
        ms2Spec.forEach(spec => {
            specIds.push(spec.getSpectrumId().toString());
        });
    }
    document.title = "Protein-Spectrum-Match for Spectrum #" + specIds.toString();
    let l_allproteins_url = "proteins.html?folder=" + folderpath;
    document.getElementById("allprotein_url").href = l_allproteins_url;
    document.getElementById("allprotein_url_end").href = l_allproteins_url;
    let l_protein_URL = proteoformObj.getProtName() + " " + proteoformObj.getProtDesc();
    let l_protroform_URL = "Proteoform #" + proteoformObj.getId();
    let proteinUrl = document.getElementById("protien_url");
    let proteinUrlEnd = document.getElementById("protien_url_end");
    let protUrl = document.getElementById("proteoform_url");
    let protUrlEnd = document.getElementById("proteoform_url_end");
    let titleElement = document.getElementById("Protein-Spectrum-Match-Id-SpecId");
    proteinUrl.innerHTML = l_protein_URL;
    proteinUrl.href = "protein.html?folder=" + folderpath + "&protein_Id=" + proteoformObj.getSeqId();
    proteinUrlEnd.innerHTML = l_protein_URL;
    proteinUrlEnd.href = "protein.html?folder=" + folderpath + "&protein_Id=" + proteoformObj.getSeqId();
    protUrl.innerHTML = l_protroform_URL;
    protUrl.href = "proteoform.html?folder=" + folderpath + "&proteoform_Id=" + proteoformObj.getId();
    protUrlEnd.innerHTML = l_protroform_URL;
    protUrlEnd.href = "proteoform.html?folder=" + folderpath + "&proteoform_Id=" + proteoformObj.getId();
    titleElement.innerHTML = "Protein-Spectrum-Match #" + prsmObj.getId() + " for Spectrum #" + specIds.toString();
}
/**
 * Get the data of the prsm from global data variable prsm_data.
 * Build the data into html to show the information about the prsm
 */
function loadDatafromJson2Html(prsmObj) {
    let fileName = document.getElementById("File_name");
    let prsmId = document.getElementById("PrSM_ID");
    let scan = document.getElementById("Scan");
    let precCharge = document.getElementById("Precursor_charge");
    let precMz = document.getElementById("Precursor_mz");
    let precMass = document.getElementById("Precursor_mass");
    let protMass = document.getElementById("Proteoform_mass");
    let matchedPeak = document.getElementById("matched_peaks");
    let matchedFrag = document.getElementById("matched_fragment_ions");
    let unexpected = document.getElementById("unexpected_modifications");
    let eVal = document.getElementById("E_value");
    let qVal = document.getElementById("Q_value");
    let proteoformObj = prsmObj.getProteoform();
    let ms2Spectrum = prsmObj.getMs2Spectra();
    if (!ms2Spectrum) {
        console.error("ERROR: invalid ms2 spectrum");
        return;
    }
    if (fileName && prsmId && scan && precCharge && precMz && precMass && protMass && matchedPeak &&
        matchedFrag && unexpected && eVal && qVal) {
        fileName.innerHTML = prsmObj.getfileName();
        prsmId.innerHTML = prsmObj.getId();
        if (ms2Spectrum.length > 1) {
            let scanText = "";
            ms2Spectrum.forEach((spectra) => {
                scanText = scanText + spectra.getScanNum() + " ";
            });
            scan.innerHTML = scanText;
        }
        else {
            scan.innerHTML = ms2Spectrum[0].getScanNum();
        }
        precCharge.innerHTML = ms2Spectrum[0].getPrecCharge().toString();
        precMz.innerHTML = FormatUtil.formatFloat(ms2Spectrum[0].getPrecMz(), "precMz");
        precMass.innerHTML = FormatUtil.formatFloat(ms2Spectrum[0].getPrecMass(), "precMass");
        protMass.innerHTML = FormatUtil.formatFloat(proteoformObj.getMass(), "protMass");
        matchedPeak.innerHTML = prsmObj.getMatchedPeakCount().toString();
        unexpected.innerHTML = prsmObj.getUnexpectedModCount().toString();
        eVal.innerHTML = prsmObj.getEValue().toString();
        qVal.innerHTML = prsmObj.getQValue().toString();
        let ionCnt = prsmObj.getFragIonCount();
        if (ionCnt) {
            matchedFrag.innerHTML = ionCnt.toString();
        }
    }
}
/**
 * Get occurence of "Variable" and "Fixed", convert the data to HTML
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function occurence_ptm(prsmObj) {
    let varPtm = "";
    let fixedPtm = "";
    let fixedPtmObj = [];
    let varPtmObj = [];
    let variablePtms = (prsmObj.getProteoform().getVarPtm()).concat(prsmObj.getProteoform().getProtVarPtm());
    let fixedPtms = prsmObj.getProteoform().getFixedPtm();
    fixedPtms.sort(function (x, y) {
        return x.getLeftPos() - y.getLeftPos();
    });
    fixedPtms.forEach((ptm) => {
        let notInObj = true;
        for (let i = 0; i < fixedPtmObj.length; i++) {
            if (fixedPtmObj[i].name == ptm.getAnnotation()) {
                fixedPtmObj[i].pos.push(ptm.getLeftPos() + 1);
                notInObj = false;
                break;
            }
        }
        if (notInObj) {
            fixedPtmObj.push({ "name": ptm.getAnnotation(), "pos": [ptm.getLeftPos() + 1] });
        }
    });
    fixedPtmObj.forEach((p) => {
        let posString = "";
        p.pos.forEach((eachPos) => {
            posString = posString + eachPos.toString() + ";";
        });
        posString = posString.slice(0, posString.length - 1); //remove ";" at the end
        fixedPtm = fixedPtm + p.name + " [" + posString + "]";
    });
    variablePtms.sort(function (x, y) {
        return x.getLeftPos() - y.getLeftPos();
    });
    variablePtms.forEach(ptm => {
        let notInObj = true;
        for (let i = 0; i < varPtmObj.length; i++) {
            if (varPtmObj[i].name == ptm.getAnnotation()) {
                if (ptm.getLeftPos() + 1 == ptm.getRightPos()) {
                    varPtmObj[i].pos.push((ptm.getLeftPos() + 1).toString());
                }
                else {
                    varPtmObj[i].pos.push((ptm.getLeftPos() + 1).toString() + "-" + (ptm.getRightPos()).toString());
                }
                notInObj = false;
                break;
            }
        }
        if (notInObj) {
            if (ptm.getLeftPos() + 1 == ptm.getRightPos()) {
                varPtmObj.push({ "name": ptm.getAnnotation(), "pos": [(ptm.getLeftPos() + 1).toString()] });
            }
            else {
                varPtmObj.push({ "name": ptm.getAnnotation(), "pos": [(ptm.getLeftPos() + 1).toString() + "-" + (ptm.getRightPos()).toString()] });
            }
        }
    });
    varPtmObj.forEach((p) => {
        let posString = "";
        p.pos.forEach((eachPos) => {
            posString = posString + eachPos + ";";
        });
        posString = posString.slice(0, posString.length - 1); //remove ";" at the end
        varPtm = varPtm + p.name + " [" + posString + "] ";
    });
    // Add the information of fixed ptms to html at id - ptm_abbreviation
    if (fixedPtm != "") {
        let div = document.getElementById("ptm_abbreviation");
        let text1 = document.createElement("text");
        let text2 = document.createElement("text");
        text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + "Fixed PTMs: ";
        text2.innerHTML = fixedPtm;
        //text2.target = "_blank";
        //@ts-ignore
        text2.style = "color:red";
        if (!div) {
            console.error("ERROR: invalid div id for ptm annotation");
            return;
        }
        div.appendChild(text1);
        div.appendChild(text2);
    }
    // Add the information of varibale ptms to html at id - ptm_abbreviation
    if (varPtm != "") {
        let div = document.getElementById("ptm_abbreviation");
        let text1 = document.createElement("text");
        let text2 = document.createElement("text");
        text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + "Variable PTMs: ";
        text2.innerHTML = varPtm;
        //text2.target = "_blank";
        //@ts-ignore
        text2.style = "color:red";
        if (!div) {
            console.error("ERROR: invalid div id for ptm annotation");
            return;
        }
        div.appendChild(text1);
        div.appendChild(text2);
    }
}
/**
 * Get the information about all the unknown ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function getUnknownPtms(prsmObj) {
    let unknownShift = "";
    let shiftCount = 0;
    prsmObj.getProteoform().getUnknownMassShift().forEach(shift => {
        if (shiftCount > 0) {
            //if it is not the first ptm and it is not the only ptm
            unknownShift = unknownShift + ", ";
        }
        //unknownShift = unknownShift + "[" + shift.getAnnotation().toString() + "]";
        unknownShift = unknownShift + FormatUtil.formatFloat(shift.getAnnotation(), "massShift");
        shiftCount++;
    });
    // If unexpected modifications exist add them to html at id - ptm_unexpectedmodifications
    if (shiftCount > 0) {
        let val = "Unknown" + " [" + unknownShift + "]";
        let unexpectedModDiv = document.getElementById("ptm_unexpectedmodification");
        if (!unexpectedModDiv) {
            console.error("ERROR: invalid div id for unexpected modificiation");
            return;
        }
        unexpectedModDiv.style.display = "block";
        let div = document.getElementById("ptm_unexpectedmodification");
        let text1 = document.createElement("text");
        let text2 = document.createElement("text");
        text1.innerHTML = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + "Unexpected modifications: "; // Adding space using &nbsp;
        text2.innerHTML = val;
        //text2.target = "_blank";
        //@ts-ignore
        text2.style = "color:red";
        if (!div) {
            console.error("ERROR: invalid div id for ptm annotation");
            return;
        }
        div.appendChild(text1);
        div.appendChild(text2);
    }
}
/**
 * Get all the variable ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function getVariablePtm(ptm) {
    let variable_ptm = "[";
    if (Array.isArray(ptm.occurence)) {
        ptm.occurence.forEach(function (occurence, i) {
            variable_ptm = variable_ptm + occurence.right_pos;
            if (ptm.occurence.length - 1 > i) {
                variable_ptm = variable_ptm + ";";
            }
        });
    }
    else {
        variable_ptm = variable_ptm + ptm.occurence.right_pos;
    }
    variable_ptm = ptm.ptm.abbreviation + variable_ptm;
    return variable_ptm;
}
/**
 * Get all the "Fixed" ptms
 * @param {object} prsm - prsm is the data attribute inside global prsm_data variable
 */
function getFixedPtm(ptm) {
    let fixed_ptm = "[";
    if (Array.isArray(ptm.occurence)) {
        ptm.occurence.forEach(function (occurence, i) {
            fixed_ptm = fixed_ptm + occurence.right_pos;
            if (ptm.occurence.length - 1 > i) {
                fixed_ptm = fixed_ptm + ";";
            }
        });
    }
    else {
        fixed_ptm = fixed_ptm + ptm.occurence.right_pos;
    }
    fixed_ptm = ptm.ptm.abbreviation + fixed_ptm;
    return fixed_ptm;
}
