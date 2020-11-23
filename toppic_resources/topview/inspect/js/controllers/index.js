// Gets executed once HTML is loaded
$(document).ready(function(){
	onLoadOfHTML();
  	localStorage.clear();
  
	//ion type data needs to be preserved because it is recorded only once when the prsm.html opens. 
	//Without this line it is going to break when the user tries to open inspect window again after closing.
    //localStorage.setItem('ionType', ionType); 
});

/**
 * @function onLoadOfHTML
 * @description Gets invoked immediatley after loading html
 */
onLoadOfHTML = function()
{
    // Get the data from local storage 
	let peakAndIntensityList = parsePeakMass('peakAndIntensityList');
	// console.log(peakAndIntensityList);
	let massAndIntensityList = parsePeakMass('massAndIntensityList');
	let sequence = parseSeq('sequence');
	let l_fixedPtmList = parsePTM('fixedPtmList');
    let l_variablePtmList = parsePTM('variablePtmList');
    let unknownMassShiftList = parseUnknowmassList('unknownMassShiftList');
    let precursorMass = parsePrecursorMass("precursorMass");
        
    if(peakAndIntensityList !== null && massAndIntensityList !== null){
        setDataToPeakAndIntensity(peakAndIntensityList);
		setDataToMassAndIntensity(massAndIntensityList);
    }
    if(sequence) {
        setDataToSequence(sequence, unknownMassShiftList, l_variablePtmList);
    }
    if(l_fixedPtmList) {
        setFixedMasses(l_fixedPtmList);
    }
    if(precursorMass) {
		setPrecursorMass(precursorMass);
    }
	// if(peakAndIntensityList !== null && massAndIntensityList !== null
	// 	&& sequence !== null && precursorMass !== null)
	// {	
	// 	setDataToPeakAndIntensity(peakAndIntensityList);
	// 	setDataToMassAndIntensity(massAndIntensityList);
	// 	setDataToSequence(sequence, unknownMassShiftList);
    //     setFixedMasses(l_fixedPtmList);
	// 	setPrecursorMass(precursorMass);
    // }
    //set the checkbox based on the ion type used in the data, which is stored in local storage
    let ionType = getIonType();
    setIonCheckbox(ionType);

    let massErrorthVal = 0.1;
    let ppmErrorthVal = 15;
    /**
     * set All the common fixed PTM's to Fixed Ptm dropdown menu
     */
    setFixedPtmListToUI(COMMON_FIXED_PTM_LIST);
    /**
     * Set Error Threshhold value to default {massErrorthVal}
     */
    setMassErrorValue(massErrorthVal);
    /**
     * On Change Event handler. Changes the thresholds values from
     * from {massErrorthVal} to {ppmErrorthVal}
     */
    jqueryElements.errorDropdown.change(() => {
        let errorType = jqueryElements.errorDropdown.val();
        if(errorType === "masserror")
        {
            jqueryElements.errorValue.val(massErrorthVal);
            jqueryElements.errorUnit.html("Da&nbsp;&nbsp;");
        }
        else
        {
            jqueryElements.errorValue.val(ppmErrorthVal);
            jqueryElements.errorUnit.html("ppm&nbsp;&nbsp;");
        }
    });
    /**
     * On Click Event handler. Gets invoked on click of submit button
     * in HTML
     */
    jqueryElements.submit.click(function(){
        let errorVal ;
        let errorType = jqueryElements.errorDropdown.val();
        if(errorType === "masserror") {
            errorVal = parseFloat(jqueryElements.errorValue.val().trim());
            massErrorthVal = errorVal ;
        }
        else {
            errorVal = parseFloat(jqueryElements.errorValue.val().trim());
            ppmErrorthVal = errorVal ;
        }
        let executionObj = new SeqOfExecution();
        executionObj.sequenceOfExecution(errorType,errorVal,"", precursorMass);
        // domElements.totalSeqMass.style.display = "block";
        // domElements.massVariation.style.display = "block";
    })
    /**
     * On Click action to hide and show the table of calculate theoretical
     * masses with matched and unmatched masses
     */
    jqueryElements.hideTable.click(function(){
        let text_val = jqueryElements.hideTable.html().trim();
        if( text_val === "Hide Table"){
            jqueryElements.hideTable.html("Show Table");
            jqueryElements.monoMassTableContainer.hide();
        }
        if(text_val === "Show Table")
        {
            jqueryElements.hideTable.html("Hide Table");
            jqueryElements.monoMassTableContainer.show();
        }
    })
}



/**
 * @function showAllPeaks
 * @description Function to display all peaks of data in table. This handles on click action
 * from html of show all peaks button.
 */
showAllPeaks = function() {
	var elems = domElements.matchedPeaks;
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = '';
	}
	elems = domElements.unmatchedPeaks;
	for(var i = 0; elems.length > i; i++) {
	    elems[i].style.display = '';
	}
	$('div.dataTables_scrollBody').height(400);
}
/**
 * @function showMatchedPeaks
 * @description Function to display only matched peaks in table. This handles on click action 
 * from html of show matched peaks button.
 */
showMatchedPeaks = function()
{
	var elems = domElements.matchedPeaks;
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "";
	}
	elems = domElements.unmatchedPeaks;
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "none";
	}
	$('div.dataTables_scrollBody').height(400);
}
/**
 * @function showNonMatchedPeaks
 * @description Function to display only un matched peaks in table. This handles on click action
 * from html of show un matched peaks button.
 */
showNonMatchedPeaks = function() 
{
	var elems = domElements.matchedPeaks;
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "none";
	}
	elems = domElements.unmatchedPeaks;
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "";
	}
	$('div.dataTables_scrollBody').height(400);
}