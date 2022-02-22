"use strict";
/**
 * Function to set default error value to HTML in Da units
 * @param {Float} massErrorthVal - Contains Mass Error Value in Da units
 */
function setMassErrorValue(massErrorthVal) {
    jqueryElements.errorValue.val(massErrorthVal);
    jqueryElements.errorUnit.html("Da&nbsp;&nbsp;");
}
/**
 * Function to Set Default error value to HTML in ppm units
 * @param {Float} ppmErrorthVal - Contains Mass Error in ppm units
 */
function setPPMErrorValue(ppmErrorthVal) {
    jqueryElements.errorValue.val(ppmErrorthVal);
    jqueryElements.errorUnit.html("ppm&nbsp;&nbsp;");
}
/**
 * Set Default Mass errors into html of both Da and ppm units
 * @param {Float} massErrorthVal - Contains Mass error in Da units
 * @param {Float} ppmErrorthVal - Contains ppm error in ppm units
 */
function writeMassErrorThreshholdValueToUI(massErrorthVal, ppmErrorthVal) {
    if (massErrorthVal)
        jqueryElements.errorValue.val(massErrorthVal);
    else
        jqueryElements.errorValue.val(ppmErrorthVal);
}
