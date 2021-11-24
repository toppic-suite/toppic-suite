"use strict";
const domElements = {
    precursorMass: document.getElementById("precursormass"),
    customControlInput: document.getElementsByClassName("custom-control-input"),
    totalSeqMass: document.getElementById("totalseqmass_h6"),
    massVariation: document.getElementById("massvariation_h6"),
    dropDownMenuLink: document.getElementById("dropdownMenuLink"),
    fixedPtmList: document.getElementById("fixedptmslist"),
    matchedPeaks: document.getElementsByClassName('matched_peak'),
    unmatchedPeaks: document.getElementsByClassName('unmatched_peak'),
    ionTableContainer: document.getElementById("divselectediontablecontainer"),
    massDifference: document.getElementById("massvariation"),
    monoMassTableContainer: document.getElementById("divtableContainer"),
};
const jqueryElements = {
    peakData: $("#peakdata"),
    massData: $("#data"),
    sequenceData: $("#sequencedata"),
    errorDropdown: $('#error_dropdown'),
    errorValue: $("#errorval"),
    errorUnit: $("#errorunit"),
    submit: $("#submit"),
    hideTable: $('#hide_table'),
    monoMassTableContainer: $("#divtableContainer"),
    addFixedPtmRow: $('.addnewrow'),
    fixedPtms: $(".fixedptms"),
    ionTableContainer: $("#divselectediontablecontainer"),
    // matchedFragments : $(".matched_fragments"),
    totalMass: $("#totalmass"),
    monoMassTableBody: $("#tableContainer tbody"),
    monoMassTableColumns: $("#tableContainer thead tr th"),
    // peakRow : $(".peakRows"),
    allPeakCount: $("#all_peak_count"),
    matchedPeakCount: $("#matched_peak_count"),
    unmatchedPeakCount: $("#unmatched_peak_count"),
    peakCount: $("#peakCount"),
    precursorMassSubmit: $("#precursormass-submit"),
};
const COMMON_FIXED_PTM_LIST = [new Mod("C", 57.021464, "Carbamidomethylation"), new Mod("C", 58.005479, "Carboxymethyl")];
let USER_FIXED_PTM_LIST = []; //user-defined custom ptm list
let VAR_PTM_LIST = [];
//it is not const because the variable PTM in the original prsm 
//will be added to this list in index.js
/**
 * Class contains all the id's and class names used in HTML
 */
class Constants {
}
Constants.SEQSVGID = "seqsvg";
Constants.SVGDOWNLOADID = "save_prsm_btn";
Constants.SPECTRUMGRAPHID = "ms2_svg_div_graph_0";
Constants.MONOMASSGRAPHID = "ms2_svg_div_mono_graph_0";
Constants.PEAKCOUNTID = "peakCount";
Constants.MASSERROR = "masserror";
Constants.MASSSHIFT_CLASSID = "massshift_class";
Constants.TABLECONTAINERID = "tableContainer";
Constants.FRAGMENTIONTABLECONTAINER = "selectedIonTableContainer";
Constants.H_FRAGMENTEDTABLE = "h_fragmentedtable";
Constants.DIVTABLECONTAINER = "divtableContainer";
Constants.TABLEHEIGHT = "400px";
Constants.GRAPHDOWNLOAD = "ms2_graph_save_btn";
Constants.GRAPHTABDIV = "ms2_svg_div";
Constants.GRAPHTABNAV = "ms2_graph_nav";
