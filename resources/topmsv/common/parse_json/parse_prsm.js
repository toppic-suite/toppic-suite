// @ts-nocheck
//parse data from js file
"use strict";
class ParsePrsm {
    constructor(drawMs1Spec = false, ms1SpecPath = null, drawMs2Spec = false, ms2SpecPath = null) {
        this.drawMs1Spec_ = drawMs1Spec;
        this.drawMs2Spec_ = drawMs2Spec;
        this.ms1SpecPath_ = ms1SpecPath;
        this.ms2SpecPath_ = ms2SpecPath;
    }
    geneDataObj(callback) {
        let protObj = this.geneProteoform();
        this.genePrsm(protObj, function (prsmObj) {
            callback(prsmObj);
        });
    }
    genePrsm(proteoformObj, callback) {
        let prsm = prsm_data.prsm;
        let prot = prsm_data.prsm.annotated_protein;
        let breakPoints = json2BreakPoints(prsm, parseInt(prot.annotation.first_residue_position));
        let [nIons, cIon] = this.parseIons();
        /*this.geneMs1Spectrum(nIons, cIon, (ms1Spec) => {
            this.geneMs2Spectrum(nIons, cIon, (ms2Spec) => {
                let machedPeakEnvPair = [];
                ms2Spec.forEach((spectra) => {
                    machedPeakEnvPair = machedPeakEnvPair.concat(this.geneMatchedPeakEnvelopePairs(spectra.getSpectrumId(), spectra.getEnvs()));
                });
                callback(new Prsm(prsm.prsm_id, proteoformObj, ms1Spec, ms2Spec, breakPoints, machedPeakEnvPair, prsm.ms.ms_header.spectrum_file_name, prsm.e_value, prsm.fdr, undefined, undefined, prsm.matched_fragment_number));
            });
        });*/
        this.geneMs2Spectrum(nIons, cIon, (ms2Spec, targetMz, maxMz, minMz) => {
            let machedPeakEnvPair = [];
            ms2Spec.forEach((spectra) => {
                machedPeakEnvPair = machedPeakEnvPair.concat(this.geneMatchedPeakEnvelopePairs(spectra.getSpectrumId(), spectra.getEnvs()));
            });
            this.geneMs1Spectrum(nIons, cIon, targetMz, maxMz, minMz, (ms1Spec) => {
                callback(new Prsm(prsm.prsm_id, proteoformObj, ms1Spec, ms2Spec, breakPoints, machedPeakEnvPair, prsm.ms.ms_header.spectrum_file_name, prsm.e_value, prsm.fdr, undefined, undefined, prsm.matched_fragment_number));
            });
        });
    }
    geneProteoform() {
        let prot = prsm_data.prsm.annotated_protein;
        let prsm = prsm_data.prsm;
        let fixedPtms = [];
        let protVarPtms = [];
        let variablePtms = [];
        [fixedPtms, protVarPtms, variablePtms] = json2Ptms(prsm);
        let massShifts = json2MassShifts(prsm);
        let sequence = getAminoAcidSequence(0, prot.annotation.residue.length - 1, prot.annotation.residue);
        return new Proteoform(prot.proteoform_id, prot.sequence_name, prot.sequence_description, sequence, prot.sequence_id, parseInt(prot.annotation.first_residue_position), parseInt(prot.annotation.last_residue_position), parseFloat(prot.proteoform_mass), massShifts, fixedPtms, protVarPtms, variablePtms);
    }
    geneMs1Spectrum(nIons, cIons, targetMz, minMz, maxMz, callback) {
        if (!this.drawMs1Spec_) {
            callback(null);
        }
        else if (this.ms1SpecPath_ == null) {
            console.error("invalid path for ms1 spectrum file");
            callback(null);
        }
        else {
            let prsm = prsm_data.prsm;
            let ms1SpecId = prsm.ms.ms_header.ms1_ids.split(" ")[0];
            let filename = this.ms1SpecPath_ + "spectrum" + ms1SpecId + ".js";
            let script = document.createElement('script');
            script.src = filename;
            document.head.appendChild(script);
            script.onload = function () {
                let peaks = [];
                let envelopes = [];
                for (let i = 0; i < ms1_data.peaks.length; i++) {
                    let peakObj = new Peak(i.toString(), parseFloat(ms1_data.peaks[i].mz), parseFloat(ms1_data.peaks[i].mz), parseFloat(ms1_data.peaks[i].intensity));
                    peaks.push(peakObj);
                }
                for (let i = 0; i < ms1_data.envelopes.length; i++) {
                    let env = ms1_data.envelopes[i];
                    let envObj = new Envelope(env.mono_mass, env.charge);
                    for (let j = 0; j < env.env_peaks.length; j++) {
                        let peak = new Peak(j.toString(), parseFloat(env.env_peaks[j].mz), parseFloat(env.env_peaks[j].mz), parseFloat(env.env_peaks[j].intensity));
                        envObj.addPeaks(peak);
                    }
                    envelopes.push(envObj);
                }
                callback(new Spectrum(ms1_data.id.toString(), ms1_data.scan.toString(), 1, peaks, null, envelopes, nIons, cIons, parseFloat(prsm.ms.ms_header.precursor_mono_mass), parseFloat(prsm.ms.ms_header.precursor_charge), targetMz, minMz, maxMz));
            };
        }
    }
    geneMs2Spectrum(nIons, cIons, callback) {
        if (!this.drawMs2Spec_) {
            callback([]);
        }
        else if (this.ms2SpecPath_ == null) {
            console.error("invalid path for ms1 spectrum file");
            callback([]);
        }
        else {
            let prsm = prsm_data.prsm;
            let specIdList = prsm.ms.ms_header.ids.split(" ");
            let scanIdList = prsm.ms.ms_header.scans.split(" ");
            let fileList = [];
            let specList = [];
            let specObjList = [];
            if (typeof (ms2ScanList) != "undefined") {
                scanIdList.forEach((scan) => {
                    ms2ScanList.push(scan);
                });
            }
            for (let i = 0; i < specIdList.length; i++) {
                let ms2Filename = this.ms2SpecPath_ + "spectrum" + specIdList[i] + ".js";
                fileList.push(ms2Filename);
            }
            let len = fileList.length;
            let cnt = 0;
            for (let i = 0; i < len; i++) {
                let filename = fileList[i];
                let script = document.createElement('script');
                script.src = filename;
                document.head.appendChild(script);
                script.onload = function () {
                    ms2_data.id = specIdList[i];
                    specList.push(ms2_data);
                    cnt = cnt + 1;
                    if (cnt == len) {
                        specList.sort(function (x, y) {
                            return parseInt(x.scan) - parseInt(y.scan);
                        });
                        for (let j = 0; j < specList.length; j++) {
                            let peaks = [];
                            for (let k = 0; k < specList[j].peaks.length; k++) {
                                let peak = specList[j].peaks[k];
                                let peakObj = new Peak(k.toString(), parseFloat(peak.mz), parseFloat(peak.mz), parseFloat(peak.intensity));
                                peaks.push(peakObj);
                            }
                            let envelopes = [];
                            for (let m = 0; m < specList[j].envelopes.length; m++) {
                                let env = specList[j].envelopes[m];
                                let envObj = new Envelope(parseFloat(env.mono_mass), parseFloat(env.charge));
                                for (let n = 0; n < env.env_peaks.length; n++) {
                                    let peak = new Peak(n.toString(), parseFloat(env.env_peaks[n].mz), parseFloat(env.env_peaks[n].mz), parseFloat(env.env_peaks[n].intensity));
                                    envObj.addPeaks(peak);
                                }
                                envelopes.push(envObj);
                            }
                            let deconvPeaks = [];
                            prsm_data.prsm.ms.peaks.peak.forEach((peak) => {
                                if (peak.spec_id == specList[j].id) {
                                    deconvPeaks.push(new Peak(peak.peak_id, parseFloat(peak.monoisotopic_mass), parseFloat(peak.monoisotopic_mz), parseFloat(peak.intensity), parseFloat(peak.monoisotopic_mass), parseInt(peak.charge), peak.spec_id));
                                }
                            });
                            //check if nIons and cIons were identified from prsm file
                            //if not, add information here
                            if (nIons.length < 1) {
                                nIons.push(new Ion(specList[j].n_ion_type, specList[j].n_ion_type, "N", -1));
                            }
                            if (cIons.length < 1) {
                                cIons.push(new Ion(specList[j].c_ion_type, specList[j].c_ion_type, "C", -1));
                            }
                            let ms2Spectrum = new Spectrum(specList[j].id.toString(), specList[j].scan.toString(), 1, peaks, deconvPeaks, envelopes, nIons, cIons, parseFloat(prsm.ms.ms_header.precursor_mono_mass), parseFloat(prsm.ms.ms_header.precursor_charge), parseFloat(prsm.ms.ms_header.precursor_mz));
                            specObjList.push(ms2Spectrum);
                        }
                        let targetMz = parseFloat(specList[0].target_mz);
                        let minMz = parseFloat(specList[0].min_mz);
                        let maxMz = parseFloat(specList[0].max_mz);
                        callback(specObjList, targetMz, minMz, maxMz);
                    }
                };
            }
        }
    }
    parseIons() {
        let nIons = [];
        let cIons = [];
        prsm_data.prsm.ms.peaks.peak.forEach(function (peak) {
            // Check if peak contain matched_ions_num attribute	
            // for each peak, get the ion type and store it in ionArray 
            // to determine which ion type to be checked in Inspect
            if (parseInt(peak.matched_ions_num) > 0) {
                if (Array.isArray(peak.matched_ions.matched_ion)) {
                    peak.matched_ions.matched_ion.forEach((ion) => {
                        let ionTerm;
                        if (ion.ion_type == "X" || ion.ion_type == "Y" || ion.ion_type == "Z") {
                            ionTerm = "C";
                            cIons.push(new Ion(ion.ion_type + ion.ion_display_position, ion.ion_type, ionTerm, parseFloat(ion.match_shift), parseFloat(ion.mass_error), parseFloat(ion.ppm)));
                        }
                        else if (ion.ion_type == "A" || ion.ion_type == "B" || ion.ion_type == "C") {
                            ionTerm = "N";
                            nIons.push(new Ion(ion.ion_type + ion.ion_display_position, ion.ion_type, ionTerm, parseFloat(ion.match_shift), parseFloat(ion.mass_error), parseFloat(ion.ppm)));
                        }
                    });
                }
                else {
                    let ion = peak.matched_ions.matched_ion;
                    let ionTerm;
                    if (ion.ion_type == "X" || ion.ion_type == "Y" || ion.ion_type == "Z") {
                        ionTerm = "C";
                        cIons.push(new Ion(ion.ion_type + ion.ion_display_position, ion.ion_type, ionTerm, parseFloat(ion.match_shift), parseFloat(ion.mass_error), parseFloat(ion.ppm)));
                    }
                    else if (ion.ion_type == "A" || ion.ion_type == "B" || ion.ion_type == "C") {
                        ionTerm = "N";
                        nIons.push(new Ion(ion.ion_type + ion.ion_display_position, ion.ion_type, ionTerm, parseFloat(ion.match_shift), parseFloat(ion.mass_error), parseFloat(ion.ppm)));
                    }
                }
            }
        });
        return [nIons, cIons];
    }
    geneMatchedPeakEnvelopePairs(specId, envelopes) {
        let matchedPairList = [];
        let deconvPeaks = prsm_data.prsm.ms.peaks.peak;
        deconvPeaks.forEach(function (element) {
            if (element.hasOwnProperty('matched_ions_num') && element.spec_id == specId) {
                if (element.hasOwnProperty('matched_ions_num')) {
                    let ionText = "";
                    let peakId = parseInt(element.peak_id);
                    let envPeaks = envelopes[peakId].getPeaks();
                    envPeaks.sort(function (x, y) {
                        return d3.descending(x.getIntensity(), y.getIntensity());
                    });
                    if (parseInt(element.matched_ions_num) == 1) {
                        let matchedIon = element.matched_ions.matched_ion;
                        let ionType = matchedIon.ion_type;
                        if (ionType == "Z_DOT") {
                            ionType = "Z\u02D9";
                        }
                        ionText = ionType + matchedIon.ion_display_position;
                        //create matched pair
                        let matchedPeakObj = new Peak(element.peak_id, parseFloat(element.monoisotopic_mz), parseFloat(element.monoisotopic_mz), parseFloat(element.intensity), parseFloat(element.monoisotopic_mass), parseInt(element.charge), element.spec_id);
                        let ionTerm = "";
                        if (element.matched_ions.matched_ion.ion_type[0] == "X" || element.matched_ions.matched_ion.ion_type[0] == "Y" ||
                            element.matched_ions.matched_ion.ion_type[0] == "Z") {
                            ionTerm = "C";
                        }
                        else if (element.matched_ions.matched_ion.ion_type[0] == "A" || element.matched_ions.matched_ion.ion_type[0] == "B" ||
                            element.matched_ions.matched_ion.ion_type[0] == "C") {
                            ionTerm = "N";
                        }
                        let matchedIonObj = new Ion(ionText, ionType, ionTerm, parseFloat(element.matched_ions.matched_ion.match_shift), parseFloat(element.matched_ions.matched_ion.mass_error), parseFloat(element.matched_ions.matched_ion.ppm));
                        let matchedPair = new MatchedPeakEnvelopePair(parseFloat(element.matched_ions.matched_ion.theoretical_mass), matchedPeakObj, matchedIonObj);
                        matchedPair.addEnvelope(envelopes[peakId]);
                        matchedPairList.push(matchedPair);
                    }
                    else {
                        for (let i = 0; i < parseInt(element.matched_ions_num); i++) {
                            let matchedIon = element.matched_ions.matched_ion[i];
                            let ionType = matchedIon.ion_type;
                            if (ionType == "Z_DOT") {
                                ionType = "Z\u02D9";
                            }
                            ionText = ionType + matchedIon.ion_display_position;
                            let matchedPeakObj = new Peak(element.peak_id, parseFloat(element.monoisotopic_mz), parseFloat(element.monoisotopic_mz), parseFloat(element.intensity), parseFloat(element.monoisotopic_mass), parseInt(element.charge), element.spec_id);
                            let ionTerm = "";
                            if (element.matched_ions.matched_ion[i].ion_type[0] == "X" || element.matched_ions.matched_ion[i].ion_type[0] == "Y" ||
                                element.matched_ions.matched_ion[i].ion_type[0] == "Z") {
                                ionTerm = "C";
                            }
                            else if (element.matched_ions.matched_ion[i].ion_type[0] == "A" || element.matched_ions.matched_ion[i].ion_type[0] == "B" ||
                                element.matched_ions.matched_ion[i].ion_type[0] == "C") {
                                ionTerm = "N";
                            }
                            let matchedIonObj = new Ion(ionText, ionType, ionTerm, parseFloat(element.matched_ions.matched_ion[i].match_shift), parseFloat(element.matched_ions.matched_ion[i].mass_error), parseFloat(element.matched_ions.matched_ion[i].ppm));
                            let matchedPair = new MatchedPeakEnvelopePair(parseFloat(element.matched_ions.matched_ion[i].theoretical_mass), matchedPeakObj, matchedIonObj);
                            matchedPair.addEnvelope(envelopes[peakId]);
                            matchedPairList.push(matchedPair);
                        }
                    }
                }
            }
        });
        return matchedPairList;
    }
}
