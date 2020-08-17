const elements = {
    peakData = $("#peakdata"),
    massData = $("#mass"),
    sequenceData = $("#sequencedata"),
    precursorMass = document.getElementById("precursormass"),
    customControlInput = document.getElementsByClassName("custom-control-input")
};

const commonfixedPtmList = [{name:"Carbamidomethylation",acid:"C",mass:57.021464},{name:"Carboxymethyl",acid:"C",mass:58.005479}];