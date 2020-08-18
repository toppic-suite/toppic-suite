const domElements = {
    precursorMass = document.getElementById("precursormass"),
    customControlInput = document.getElementsByClassName("custom-control-input"),
    totalSeqMass = document.getElementById("totalseqmass_h6"),
    massVariation = document.getElementById("massvariation_h6")
};

const jqueryElements = {
    peakData = $("#peakdata"),
    massData = $("#mass"),
    sequenceData = $("#sequencedata"),
    errorDropdown = $('#error_dropdown'),
    errorValue = $("#errorval"),
    errorUnit = $("#errorunit"),
    submit = $("#submit"),
    hideTable = $('#hide_table'),
    divTableContainer = $("#divtableContainer")
}

const commonfixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];