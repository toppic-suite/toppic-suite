/**
 * Set Sequence on to html
 */
function setDataToSequence(sequence, massShiftList, protVarPtmsList, variablePtmsList){
    let massShiftObj = new MassShifts(sequence, massShiftList);
    let modifiedSequence = massShiftObj.formSequence();
    if(protVarPtmsList || variablePtmsList){
        modifiedSequence = addVariablePtm(modifiedSequence, protVarPtmsList, variablePtmsList);
    }
    jqueryElements.sequenceData.val(modifiedSequence);
}
/**
 * Below function adds variable PTM annotation as mass shifts values
 */
function addToSequence(sequence, protVarPtms, variablePtms){
    let tempSeq;
    let isResidue = true;
    let residuePos = 0;
    for (let i = 0; i < sequence.length; i++){
        for (let j = 0; j < protVarPtms.posList.length; j++){
            if (residuePos == protVarPtms.posList[j].leftPos){
                let decimal = (protVarPtms.mono_mass).indexOf(".");
                let mass = protVarPtms.mono_mass.slice(0, decimal + 5);
                let tempString = "["+mass+"]";
                let leftString = sequence.slice(0, i + 1);
                let rightString = sequence.slice(i + 1);
                tempSeq = leftString + tempString + rightString;
                return tempSeq;
            }
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
}
/**
 * Below function adds variable PTM annotation as texts
 */
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
function addVariablePtm(sequence, protVarPtmsList, variablePtmsList){
    let newSeq = sequence;
    for (let i = 0; i < protVarPtmsList.length; i++){
        newSeq = addToSequence(newSeq, protVarPtmsList[i]);
    }
    for (let i = 0; i < variablePtmsList.length; i++){
        newSeq = addToSequence(newSeq, variablePtmsList[i]);
    }
    return newSeq;
}

/**
 * Get the sequence entered from the HTML.
 */
function getSequenceFromUI(){
    var seq = jqueryElements.sequenceData.val().trim();
    seq = seq.toUpperCase();//set the sequence to be upper case automatically -- for user convenience
    // Remove spaces if exists between sequences
    seq = seq.replace(/ +/g, "");
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
function writeSeqToTextBox(seqToUI){
    jqueryElements.sequenceData.val(seqToUI);
}