/**
 * Set Precursor mass on html 
 */
function setPrecursorMass(precursorMass: number): void{
    domElements.precursorMass.setAttribute("value", precursorMass.toString());
}
/**
 * get Precursor mass 
 */
function getPrecursorMass(): number | null {
    let mass: string | null = domElements.precursorMass.getAttribute("value");
    if (mass) {
        return parseFloat(mass);
    }
    return null;
}
function setPrecursorMassEventHandler(): void {
    jqueryElements.precursorMassSubmit.click(function(){
        let precursorMass: string | null = domElements.precursorMass.getAttribute("value");
        let totalMass: string | number | string[] | undefined = jqueryElements.totalMass.val();
        if (precursorMass && typeof(totalMass) == "string") {
            setMassDifference(parseFloat(precursorMass), parseFloat(totalMass));
        }
        else {
            console.error("ERROR: precursor mass is null");
        }
    })
}
/**
 * Set Total mass on to the html
 * @param {*} totalMass 
 */
function setTotalSeqMass(mass: number): void{
    let totalMass: string = FormatUtil.formatFloat(mass, "protMass");
    jqueryElements.totalMass.html(totalMass);
    domElements.totalSeqMass.setAttribute("style", 'block');
}

/**
 * Set Mass difference on to the html
 * @param {Float} precursorMass - Contains Precursor mass
 * @param {Float} proteinMass - Contains calculated protein Mass
 */
function setMassDifference(precursorMass: number, proteinMass: number): void{
    let diff: number = precursorMass - proteinMass;
    domElements.massDifference.innerHTML = FormatUtil.formatFloat(diff, "massDiff");
    domElements.massVariation.setAttribute("style", 'block');
}