/**
 * Set Sequence on to html
 */
function setDataToSequence(sequence, massShiftList, l_variablePtmList){
    let massShiftObj = new MassShifts(sequence, massShiftList);
    let modifiedSequence = massShiftObj.formSequence();
    if(l_variablePtmList){
        modifiedSequence = addVariablePtm(modifiedSequence, l_variablePtmList);
    }
    jqueryElements.sequenceData.val(modifiedSequence);
}
function addToSequence(sequence, l_variablePtm){
    let tempSeq;
    let isResidue = true;
    let residuePos = 0;

    for (let i = 0; i < sequence.length; i++){
        if (residuePos == l_variablePtm.position){
            let tempString = "["+l_variablePtm.name+"]";
            let leftString = sequence.slice(0, i + 1);
            let rightString = sequence.slice(i + 1);
            tempSeq = leftString + tempString + rightString;
            return tempSeq;
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
function addVariablePtm(sequence, l_variablePtmList){
    let newSeq = sequence;

    for (let i = 0; i < l_variablePtmList.length; i++){
        newSeq = addToSequence(newSeq, l_variablePtmList[i]);
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