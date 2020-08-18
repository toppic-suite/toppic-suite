/**
 * @function onLoadOfHTML
 * @description Gets invoked immediatley after loading html
 * @param {Float} precursorMass - Contains precursorMass value
 */
onLoadOfHTML = function(precursorMass)
{
    let massErrorthVal = 0.1 ;
    let ppmErrorthVal = 15 ;
    /**
     * set All the common fixed PTM's to Fixed Ptm dropdown menu
     */
    let commonFixedPtmsObj = new commonFixedPtms();
    commonFixedPtmsObj.setFixedPtmListToUI();
    /**
     * Set Error Threshhold value to default {massErrorthVal}
     */
    let UIHelperObj = new UIHelper();
    UIHelperObj.setMassErrorValue(massErrorthVal);
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
        executionObj.sequenceOfExecution(errorType,errorVal,"");
        domElements.totalSeqMass.style.display = "block";
        domElements.massVariation.style.display = "block";
    })
    /**
     * On Click action to hide and show the table of calculate theoretical
     * masses with matched and unmatched masses
     */
    jqueryElements.hideTable.click(function(){
        let text_val = jqueryElements.hideTable.html().trim();
        if( text_val === "Hide Table"){
            jqueryElements.hideTable.html("Show Table");
            jqueryElements.divTableContainer.hide();
        }
        if(text_val === "Show Table")
        {
            jqueryElements.hideTable.html("Hide Table");
            jqueryElements.divTableContainer.show();
        }
    })
}
/**
 * @function showAllPeaks
 * @description Function to display all peaks of data in table. This handles on click action
 * from html of show all peaks button.
 */
showAllPeaks = function()
{
	var elems = document.getElementsByClassName('matched_peak');
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = '';
	}
	elems = document.getElementsByClassName('unmatched_peak');
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
	var elems = document.getElementsByClassName("matched_peak");
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "";
	}
	elems = document.getElementsByClassName("unmatched_peak");
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
	var elems = document.getElementsByClassName("matched_peak");
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "none";
	}
	elems = document.getElementsByClassName("unmatched_peak");
	for(var i = 0; elems.length > i; i++) {
		elems[i].style.display = "";
	}
	$('div.dataTables_scrollBody').height(400);
}
/**
 * Function to show only matched on unmatched peaks on click of matched or unmatched peak buttons
 * @param {String} ids - Contains Ids of respective matched peaks or un matched peaks
 */
function showIonPeaks(ids) {
	console.log("ids : ", ids);
	  var elems = document.getElementsByClassName('matched_peak');
	  for(var i = 0; elems.length > i; i++) {
	    elems[i].style.display = 'none';
	  }
	  elems = document.getElementsByClassName('unmatched_peak');
	  for(var i = 0; elems.length > i; i++) {
	    elems[i].style.display = 'none';
	  }

	 elems = document.getElementsByName(ids);
	    for(var j = 0; elems.length > j; j++) {
	      elems[j].style.display  =  "";
	      elems[j].style.background  =  "#BEECFF";
	    }
}

/**
 * Function to zoom the graph to the mass point on click of matched mass
 */
function onClickofMatchedPeaks(){
    $(".matched_fragments").click(function(){
        let charge = $(this).attr("charge");
        let mass = $(this).html();
        let mz = mass/charge;
        console.log(mz);
        let graphFeatures = new GraphFeatures();
        ms2_graph.redraw(mz,graphFeatures);
    })
}