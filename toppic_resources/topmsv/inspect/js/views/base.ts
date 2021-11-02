const domElements = {
    precursorMass : <HTMLInputElement>document.getElementById("precursormass"),
    customControlInput : <HTMLCollectionOf<HTMLInputElement>>document.getElementsByClassName("custom-control-input"),
    totalSeqMass : <HTMLElement>document.getElementById("totalseqmass_h6"),
    massVariation : <HTMLElement>document.getElementById("massvariation_h6"),
    dropDownMenuLink : <HTMLElement>document.getElementById("dropdownMenuLink"),
    fixedPtmList : <HTMLElement>document.getElementById("fixedptmslist"),
    matchedPeaks : document.getElementsByClassName('matched_peak'),
    unmatchedPeaks : document.getElementsByClassName('unmatched_peak'),
    ionTableContainer : <HTMLElement>document.getElementById("divselectediontablecontainer"),
    massDifference : <HTMLElement>document.getElementById("massvariation"),
    monoMassTableContainer : <HTMLElement>document.getElementById("divtableContainer"),
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


const COMMON_FIXED_PTM_LIST: Mod[] = [new Mod("C", 57.021464, "Carbamidomethylation"), new Mod("C", 58.005479, "Carboxymethyl")]
let USER_FIXED_PTM_LIST: Mod[] = [];//user-defined custom ptm list

let VAR_PTM_LIST: Mod[] = [];
//it is not const because the variable PTM in the original prsm 
//will be added to this list in index.js

/**
 * Class contains all the id's and class names used in HTML
 */
class Constants{
    static SEQSVGID: string = "seqsvg";
    static SVGDOWNLOADID: string = "save_prsm_btn"
    static SPECTRUMGRAPHID: string = "ms2_svg_div_graph_0";
    static MONOMASSGRAPHID: string = "ms2_svg_div_mono_graph_0";
    static PEAKCOUNTID: string = "peakCount";
    static MASSERROR: string = "masserror";
    static MASSSHIFT_CLASSID: string = "massshift_class";
    static TABLECONTAINERID: string = "tableContainer";
    static FRAGMENTIONTABLECONTAINER: string = "selectedIonTableContainer";
    static H_FRAGMENTEDTABLE: string = "h_fragmentedtable";
    static DIVTABLECONTAINER: string = "divtableContainer";
    static TABLEHEIGHT: string = "400px";
    static GRAPHDOWNLOAD: string = "ms2_graph_save_btn";
    static GRAPHTABDIV: string = "ms2_svg_div";
    static GRAPHTABNAV: string = "ms2_graph_nav";
}