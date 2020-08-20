/**
 * Set Precursor mass on html 
 */
function setPrecursorMass(precursorMass)
{
    domElements.precursorMass.innerHTML = precursorMass;
}

/**
 * Set Total mass on to the html
 * @param {*} totalMass 
 */
function setTotalSeqMass(totalMass){
    totalMass = totalMass.toFixed(4);
    jqueryElements.totalMass.html(totalMass);
}

/**
 * Set Mass difference on to the html
 * @param {Float} precursorMass - Contains Precursor mass
 * @param {Float} proteinMass - Contains calculated protein Mass
 */
function setMassDifference(precursorMass, proteinMass){
    let diff = proteinMass - precursorMass ;
    domElements.massDifference.innerHTML = diff.toFixed(4);
}