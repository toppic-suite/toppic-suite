const domElements = {
    precursorMass : document.getElementById("precursormass"),
    customControlInput : document.getElementsByClassName("custom-control-input"),
    totalSeqMass : document.getElementById("totalseqmass_h6"),
    massVariation : document.getElementById("massvariation_h6"),
    dropDownMenuLink : document.getElementById("dropdownMenuLink"),
    fixedPtmList : document.getElementById("fixedptmslist"),
    matchedPeaks : document.getElementsByClassName('matched_peak'),
    unmatchedPeaks : document.getElementsByClassName('unmatched_peak'),
    ionTableContainer : document.getElementById("divselectediontablecontainer"),
    massDifference : document.getElementById("massvariation"),
    monoMassTableContainer : document.getElementById("divtableContainer"),
};

const jqueryElements = {
    peakData : $("#peakdata"),
    massData : $("#data"),
    sequenceData : $("#sequencedata"),
    errorDropdown : $('#error_dropdown'),
    errorValue : $("#errorval"),
    errorUnit : $("#errorunit"),
    submit : $("#submit"),
    hideTable : $('#hide_table'),
    monoMassTableContainer : $("#divtableContainer"),
    addFixedPtmRow : $('.addnewrow'),
    fixedPtms : $(".fixedptms"),
    ionTableContainer : $("#divselectediontablecontainer"),
    // matchedFragments : $(".matched_fragments"),
    totalMass : $("#totalmass"),
    monoMassTableBody : $("#tableContainer tbody"),
    monoMassTableColumns : $("#tableContainer thead tr th"),
    // peakRow : $(".peakRows"),
    allPeakCount : $("#all_peak_count"),
    matchedPeakCount : $("#matched_peak_count"),
    unmatchedPeakCount : $("#unmatched_peak_count"),
    peakCount: $("#peakCount"),
    precursorMassSubmit : $("#precursormass-submit"),
}


const COMMON_FIXED_PTM_LIST = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];
// const commonFixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];

let VAR_PTM_LIST = [];
//it is not const because the variable PTM in the original prsm 
//will be added to this list in index.js

/**
 * Class contains all the id's and class names used in HTML
 */
class Constants{
    static SEQSVGID = "seqsvg";
    static SVGDOWNLOADID = "save_prsm_btn"
    static SPECTRUMGRAPHID = "spectrum";
    static MONOMASSGRAPHID = "monoMassGraph";
    static PEAKCOUNTID = "peakCount";
    static MASSERROR = "masserror";
    static MASSSHIFT_CLASSID = "massshift_class";
    static TABLECONTAINERID = "tableContainer";
    static FRAGMENTIONTABLECONTAINER = "selectedIonTableContainer";
    static H_FRAGMENTEDTABLE = "h_fragmentedtable";
    static DIVTABLECONTAINER = "divtableContainer";
    static TABLEHEIGHT = "400px";
    static GRAPHDOWNLOAD = "ms2_graph_save_btn";
    static GRAPHTABDIV = "ms2_svg_div";
    static GRAPHTABNAV = "ms2_graph_nav";
}