// @ts-nocheck
//parse data from js file
"use strict";
function geneProteinObj(callback) {
    let prsms = [];
    let prsmCnts = [];
    if (parseInt(prsm_data.protein.compatible_proteoform_number) > 1) {
        prsm_data.protein.compatible_proteoform.forEach((compatible_proteoform, index) => {
            let [e_value, precursor_mass, prsm_id] = getBestPrsm(compatible_proteoform.prsm);
            if (Array.isArray(compatible_proteoform.prsm)) {
                for (let p of compatible_proteoform.prsm) {
                    if (prsm_id == p.prsm_id) {
                        let prsm = p;
                        let prot = prsm.annotated_protein;
                        let [fixedPtms, protVarPtms, variablePtms] = json2Ptms(prsm);
                        let massShifts = json2MassShifts(prsm);
                        let sequence = getAminoAcidSequence(0, prot.annotation.residue.length - 1, prot.annotation.residue);
                        let breakPoints = json2BreakPoints(prsm, parseInt(prot.annotation.first_residue_position));
                        let proteoformObj = new Proteoform(prot.proteoform_id, prot.sequence_name, prot.sequence_description, sequence, prot.sequence_id, parseInt(prot.annotation.first_residue_position), parseInt(prot.annotation.last_residue_position), parseInt(prot.proteoform_mass), massShifts, fixedPtms, protVarPtms, variablePtms);
                        let bestPrSM = new Prsm(prsm.prsm_id, proteoformObj, "", "", breakPoints, prsm.matched_peak_number, prsm.matched_fragment_number, prsm.e_value, prsm.fdr, prsm.ms.ms_header.feature_inte, prsm.ms.ms_header.precursor_mono_mass);
                        //also have to save total number of proteoforms
                        prsmCnts.push(compatible_proteoform.prsm.length);
                        prsms.push(bestPrSM);
                        break;
                    }
                }
            }
            else {
                let prsm = compatible_proteoform.prsm;
                let prot = prsm.annotated_protein;
                let [fixedPtms, protVarPtms, variablePtms] = json2Ptms(prsm);
                let massShifts = json2MassShifts(prsm);
                let sequence = getAminoAcidSequence(0, prot.annotation.residue.length - 1, prot.annotation.residue);
                let breakPoints = json2BreakPoints(prsm, parseInt(prot.annotation.first_residue_position));
                let proteoformObj = new Proteoform(prot.proteoform_id, prot.sequence_name, prot.sequence_description, sequence, prot.sequence_id, parseInt(prot.annotation.first_residue_position), parseInt(prot.annotation.last_residue_position), parseInt(prot.proteoform_mass), massShifts, fixedPtms, protVarPtms, variablePtms);
                let bestPrSM = new Prsm(prsm.prsm_id, proteoformObj, "", "", breakPoints, prsm.matched_peak_number, prsm.matched_fragment_number, prsm.e_value, prsm.fdr, prsm.ms.ms_header.feature_inte, prsm.ms.ms_header.precursor_mono_mass);
                //also have to save total number of proteoforms
                prsmCnts.push(1);
                prsms.push(bestPrSM);
            }
        });
    }
    else {
        let [e_value, precursor_mass, prsm_id] = getBestPrsm(prsm_data.protein.compatible_proteoform.prsm);
        if (Array.isArray(prsm_data.protein.compatible_proteoform.prsm)) {
            for (let p of prsm_data.protein.compatible_proteoform.prsm) {
                if (prsm_id == p.prsm_id) {
                    let prsm = p;
                    let prot = prsm.annotated_protein;
                    let [fixedPtms, protVarPtms, variablePtms] = json2Ptms(prsm);
                    let massShifts = json2MassShifts(prsm);
                    let sequence = getAminoAcidSequence(0, prot.annotation.residue.length - 1, prot.annotation.residue);
                    let breakPoints = json2BreakPoints(prsm, parseInt(prot.annotation.first_residue_position));
                    let proteoformObj = new Proteoform(prot.proteoform_id, prot.sequence_name, prot.sequence_description, sequence, prot.sequence_id, parseInt(prot.annotation.first_residue_position), parseInt(prot.annotation.last_residue_position), parseInt(prot.proteoform_mass), massShifts, fixedPtms, protVarPtms, variablePtms);
                    let bestPrSM = new Prsm(prsm.prsm_id, proteoformObj, "", "", breakPoints, prsm.matched_peak_number, prsm.matched_fragment_number, prsm.e_value, prsm.fdr, prsm.ms.ms_header.feature_inte, prsm.ms.ms_header.precursor_mono_mass);
                    //also have to save total number of proteoforms
                    prsmCnts.push(prsm_data.protein.compatible_proteoform.prsm.length);
                    prsms.push(bestPrSM);
                    break;
                }
            }
        }
        else {
            let prsm = prsm_data.protein.compatible_proteoform.prsm;
            let prot = prsm.annotated_protein;
            let [fixedPtms, protVarPtms, variablePtms] = json2Ptms(prsm);
            let massShifts = json2MassShifts(prsm);
            let sequence = getAminoAcidSequence(0, prot.annotation.residue.length - 1, prot.annotation.residue);
            let breakPoints = json2BreakPoints(prsm, parseInt(prot.annotation.first_residue_position));
            let proteoformObj = new Proteoform(prot.proteoform_id, prot.sequence_name, prot.sequence_description, sequence, prot.sequence_id, parseInt(prot.annotation.first_residue_position), parseInt(prot.annotation.last_residue_position), parseInt(prot.proteoform_mass), massShifts, fixedPtms, protVarPtms, variablePtms);
            let bestPrSM = new Prsm(prsm.prsm_id, proteoformObj, "", "", breakPoints, prsm.matched_peak_number, prsm.matched_fragment_number, prsm.e_value, prsm.fdr, prsm.ms.ms_header.feature_inte, prsm.ms.ms_header.precursor_mono_mass);
            //also have to save total number of proteoforms
            prsmCnts.push(prsm_data.protein.compatible_proteoform.prsm.length);
            prsms.push(bestPrSM);
        }
    }
    callback(prsms, prsmCnts);
}
/**
 * Get "precursor mass","prsm Id" and least "e value" for each proteoform
 * @param {object} prsm - COntains information of the prsms of a proteoform
 */
function getBestPrsm(prsm) {
    let e_value = " ";
    let precursor_mass = " ";
    let prsm_id = "";
    if (Array.isArray(prsm)) {
        let temp = parseFloat(prsm[0].e_value);
        e_value = prsm[0].e_value;
        precursor_mass = prsm[0].ms.ms_header.precursor_mono_mass;
        prsm_id = prsm[0].prsm_id;
        for (let i = 1; i < (prsm.length); i++) {
            if (temp >= parseFloat(prsm[i].e_value)) {
                temp = parseFloat(prsm[i].e_value);
                e_value = prsm[i].e_value;
                precursor_mass = prsm[i].ms.ms_header.precursor_mono_mass;
                prsm_id = prsm[i].prsm_id;
            }
        }
    }
    else {
        e_value = prsm.e_value;
        precursor_mass = prsm.ms.ms_header.precursor_mono_mass;
        prsm_id = prsm.prsm_id;
    }
    return [e_value, precursor_mass, prsm_id];
}
