"use strict";
/**
 * Set Precursor mass on html
 */
function setPrecursorMass(precursorMass) {
    domElements.precursorMass.setAttribute("value", precursorMass.toString());
}
/**
 * get Precursor mass
 */
function getPrecursorMass() {
    let mass = domElements.precursorMass.getAttribute("value");
    if (mass) {
        return parseFloat(mass);
    }
    return null;
}
function setPrecursorMassEventHandler() {
    jqueryElements.precursorMassSubmit.click(function () {
        let precursorMass = domElements.precursorMass.getAttribute("value");
        let totalMass = jqueryElements.totalMass.val();
        if (precursorMass && typeof (totalMass) == "string") {
            setMassDifference(parseFloat(precursorMass), parseFloat(totalMass));
        }
        else {
            console.error("ERROR: precursor mass is null");
        }
    });
}
/**
 * Set Total mass on to the html
 * @param {*} totalMass
 */
function setTotalSeqMass(mass) {
    let totalMass = FormatUtil.formatFloat(mass, "protMass");
    jqueryElements.totalMass.html(totalMass);
    domElements.totalSeqMass.setAttribute("style", 'block');
}
/**
 * Set Mass difference on to the html
 * @param {Float} precursorMass - Contains Precursor mass
 * @param {Float} proteinMass - Contains calculated protein Mass
 */
function setMassDifference(precursorMass, proteinMass) {
    let diff = precursorMass - proteinMass;
    domElements.massDifference.innerHTML = FormatUtil.formatFloat(diff, "massDiff");
    domElements.massVariation.setAttribute("style", 'block');
}
