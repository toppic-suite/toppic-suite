/**
 * Set Sequence on to html
 */
function setDataToSequence(sequence, massShiftList){
    let massShiftObj = new MassShifts(sequence, massShiftList);
    let modifiedSequence = massShiftObj.formSequence();
    jqueryElements.sequenceData.val(modifiedSequence);
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