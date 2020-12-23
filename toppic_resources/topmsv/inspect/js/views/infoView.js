/**
 * Set Precursor mass on html 
 */
function setPrecursorMass(precursorMass)
{
    domElements.precursorMass.value = precursorMass;
}
/**
 * get Precursor mass 
 */
function getPrecursorMass()
{
    return domElements.precursorMass.value;
}
function setPrecursorMassEventHandler()
{
    jqueryElements.precursorMassSubmit.click(function(){
        let precursorMass = domElements.precursorMass.value;
        let totalMass = jqueryElements.totalMass.html();
        setMassDifference(precursorMass,totalMass);
    })
}
/**
 * Set Total mass on to the html
 * @param {*} totalMass 
 */
function setTotalSeqMass(totalMass){
    totalMass = totalMass.toFixed(4);
    jqueryElements.totalMass.html(totalMass);
    domElements.totalSeqMass.style = 'block';
}

/**
 * Set Mass difference on to the html
 * @param {Float} precursorMass - Contains Precursor mass
 * @param {Float} proteinMass - Contains calculated protein Mass
 */
function setMassDifference(precursorMass, proteinMass){
    let diff = precursorMass - proteinMass;
    domElements.massDifference.innerHTML = diff.toFixed(4);
    domElements.massVariation.style = 'block';
}