"use strict";
/**
 * Set Sequence on to html
 */
function setDataToSequence(sequence, massShiftList, protVarPtmsList, variablePtmsList) {
    let modifiedSequence = formSequence(sequence, massShiftList, protVarPtmsList, variablePtmsList);
    if (protVarPtmsList || variablePtmsList) {
        modifiedSequence = addVariablePtm(modifiedSequence, protVarPtmsList, variablePtmsList);
    }
    jqueryElements.sequenceData.val(modifiedSequence);
}
/**
 * Below function adds variable PTM annotation as mass shifts values
 */
/*function addToSequence(sequence: string, protVarPtms: MassShift, variablePtms: MassShift){
    let tempSeq: string;
    let isResidue: boolean = true;
    let residuePos: number = 0;
    for (let i = 0; i < sequence.length; i++){
            if (residuePos == protVarPtms.getLeftPos()){
                let decimal = (protVarPtms.mono_mass).indexOf(".");
                let mass = protVarPtms.mono_mass.slice(0, decimal + 5);
                let tempString = "["+mass+"]";
                let leftString = sequence.slice(0, i + 1);
                let rightString = sequence.slice(i + 1);
                tempSeq = leftString + tempString + rightString;
                return tempSeq;
            }
        
        for (let j = 0; j < variablePtms.posList.length; j++){
            if (residuePos == variablePtms.posList[j].leftPos){
                let decimal = (protVarPtms.mono_mass).indexOf(".");
                let mass = protVarPtms.mono_mass.slice(0, decimal + 5);
                let tempString = "["+mass+"]";
                let leftString = sequence.slice(0, i + 1);
                let rightString = sequence.slice(i + 1);
                tempSeq = leftString + tempString + rightString;
                return tempSeq;
            }
        }
        if (sequence[i] == "["){
            isResidue = false;
        }
        if (sequence[i] == "]"){
            isResidue = true;
            continue;
        }
        if (isResidue){
            residuePos++;
        }
    }
    return sequence;
}*/
/**
 * Below function adds variable PTM annotation as texts
 */
function addToSequence(sequence, variablePtms) {
    let tempSeq;
    let isResidue = true;
    let residuePos = 0;
    for (let i = 0; i < sequence.length; i++) {
        if (residuePos == variablePtms.getLeftPos()) {
            let tempString = "[" + variablePtms.getAnnotation() + "]";
            let leftString = sequence.slice(0, i + 1);
            let rightString = sequence.slice(i + 1);
            tempSeq = leftString + tempString + rightString;
            return tempSeq;
        }
        if (sequence[i] == "[") {
            isResidue = false;
        }
        if (sequence[i] == "]") {
            isResidue = true;
            continue;
        }
        if (isResidue) {
            residuePos++;
        }
    }
    return sequence;
}
/*function addToSequence(sequence, protVarPtms, variablePtms){
    let tempSeq;
    let isResidue = true;
    let residuePos = 0;

    for (let i = 0; i < sequence.length; i++){
        for (let j = 0; j < protVarPtms.posList.length; j++){
            if (residuePos == protVarPtms.posList[j].leftPos){
                let tempString = "["+protVarPtms.name+"]";
                let leftString = sequence.slice(0, i + 1);
                let rightString = sequence.slice(i + 1);
                tempSeq = leftString + tempString + rightString;
                return tempSeq;
            }
        }
        for (let j = 0; j < variablePtms.posList.length; j++){
            if (residuePos == variablePtms.posList[j].leftPos){
                let tempString = "["+variablePtms.name+"]";
                let leftString = sequence.slice(0, i + 1);
                let rightString = sequence.slice(i + 1);
                tempSeq = leftString + tempString + rightString;
                return tempSeq;
            }
        }
        if (sequence[i] == "["){
            isResidue = false;
        }
        if (sequence[i] == "]"){
            isResidue = true;
            continue;
        }
        if (isResidue){
            residuePos++;
        }
    }
    return sequence;
}*/
function addVariablePtm(sequence, protVarPtmsList, variablePtmsList) {
    let newSeq = sequence;
    for (let i = 0; i < protVarPtmsList.length; i++) {
        newSeq = addToSequence(newSeq, protVarPtmsList[i]);
    }
    for (let j = 0; j < variablePtmsList.length; j++) {
        newSeq = addToSequence(newSeq, variablePtmsList[j]);
    }
    return newSeq;
}
/**
 * Get the sequence entered from the HTML.
 */
function getSequenceFromUI() {
    let seqData = jqueryElements.sequenceData.val();
    let seq = '';
    if (typeof (seqData) == 'string') {
        seq = seqData.trim();
        seq = seq.toUpperCase(); //set the sequence to be upper case automatically -- for user convenience
        // Remove spaces if exists between sequences
        seq = seq.replace(/ +/g, "");
    }
    return seq;
    // let massShiftList = [];
    // [seq,massShiftList]= parseSequenceMassShift(seq);
    // completeCalData.sequence = seq;
    // return [seq,massShiftList] ;
}
/**
 * write the sequence with embedded mass in [] to the screen(sequence box)
 * @param {string} seqToUI - sequence with mass shifts embedded in []
 */
function writeSeqToTextBox(seqToUI) {
    jqueryElements.sequenceData.val(seqToUI);
}
/**
* forms the seq with all the mass lists and selected fixed ptms
* @return {string} result - sequence with mass shifts embedded in []
*/
function formSequence(sequence, massShiftList, protVarPtmsList, variablePtmsList) {
    let result = sequence;
    let count = 0;
    if (!massShiftList) {
        return result;
    }
    massShiftList.sort(function (x, y) {
        return x.getLeftPos() - y.getLeftPos();
    });
    for (let i = 0; i < massShiftList.length; i++) {
        if (massShiftList[i].getShift() !== 0) {
            if (i > 0) {
                // this is the previous added mass
                let tempString = "[" + FormatUtil.formatFloat(massShiftList[i - 1].getShift(), "massShift") + "]";
                count = count + tempString.length;
            }
            // add +1 as the position need to be added after the position of the acid.
            let tempPosition = massShiftList[i].getLeftPos() + 1 + count;
            result = result.slice(0, tempPosition) + "[" + FormatUtil.formatFloat(massShiftList[i].getShift(), "massShift") + "]" + result.slice(tempPosition);
        }
    }
    return result;
}
