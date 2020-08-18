/**
 * Set Sequence on to html
 */
function setDataToSequence(sequence, unknownMassShiftList){
    let massShiftsObj = new MassShifts();
    let modSequence = massShiftsObj.formSequence(sequence,unknownMassShiftList);
    jqueryElements.sequenceData.val(modSequence);
}
/**
 * Get the sequence entered from the HTML.
 */
function getSequenceFromUI(){
    var seq = jqueryElements.sequenceData.val().trim();
    seq = seq.toUpperCase();//set the sequence to be upper case automatically -- for user convenience
    
    let massShiftList = [] ;
    [seq,massShiftList]= getMassShiftList(seq);
    /**
     * Remove spaces if exists between sequences
     */
    seq = seq.replace(/ +/g, "");
    completeCalData.sequence = seq;
    return [seq,massShiftList] ;
}

/**
 * write the sequence with embedded mass in [] to the screen(sequence box)
 * @param {string} seqToUI - sequence with mass shifts embedded in []
 */
function writeSeqToTextBox(seqToUI){
    jqueryElements.sequenceData.val(seqToUI);
}