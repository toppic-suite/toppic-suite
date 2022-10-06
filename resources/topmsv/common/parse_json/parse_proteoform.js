// @ts-nocheck
//parse data from js file
"use strict";
function geneProteoformObj(callback) {
    let prsms = [];
    prsm_data.compatible_proteoform.prsm.forEach((prsm, index) => {
        let proteoformObj = new Proteoform(prsm_data.compatible_proteoform.proteoform_id, prsm_data.compatible_proteoform.sequence_name, prsm_data.compatible_proteoform.sequence_description, "", prsm_data.compatible_proteoform.sequence_id, "", "", prsm.ms.ms_header.precursor_mono_mass, [], [], [], []);
        let matchedPairList = [];
        let peakList = [];
        let peakId = 0;
        for (let peak of prsm.ms.peaks.peak) {
            //if it has matched ion property
            if ("matched_ions" in peak) {
                let peakObj = new Peak(peakId.toString(), -1, -1, -1);
                if (parseInt(peak.matched_ions_num) > 1) {
                    for (let matchedIon in peak.matched_ions.matched_ion) {
                        let ionObj = new Ion("", matchedIon.ion_type, "", matchedIon.mass_shift);
                        matchedPairList.push(new MatchedPeakEnvelopePair(-1, peakObj, ionObj));
                    }
                }
                else {
                    let ionObj = new Ion(peak.matched_ions.matched_ion.ion_type + peak.matched_ions.matched_ion.ion_position, peak.matched_ions.matched_ion.ion_type, "", peak.matched_ions.matched_ion.mass_shift);
                    matchedPairList.push(new MatchedPeakEnvelopePair(-1, peakObj, ionObj));
                }
                peakId++;
            }
            peakList.push(peak);
        }
        let spectrum = new Spectrum("", prsm.ms.ms_header.scans, -1, peakList, [], [], [], [], -1);
        let prsmObj = new Prsm(prsm.prsm_id, proteoformObj, null, [spectrum], [], matchedPairList, "", prsm.e_value, undefined, undefined, undefined, parseInt(prsm.matched_fragment_number));
        //prsmObj.setMatchedPeakEnvelopePairs(matchedPairList);
        prsms.push(prsmObj);
    });
    callback(prsms);
}
