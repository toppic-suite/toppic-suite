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
    $('#error_dropdown').change(function(){
        let errorType = $("#error_dropdown").val();
        if(errorType == "masserror")
        {
            $("#errorval").val(massErrorthVal);
            $("#errorunit").html("Da&nbsp;&nbsp;");
        }
        else
        {
            $("#errorval").val(ppmErrorthVal);
            $("#errorunit").html("ppm&nbsp;&nbsp;");
        }
    });
    /**
     * On Click Event handler. Gets invoked on click of submit button
     * in HTML
     */
    $("#submit").click(function(){
        let errorVal ;
        let errorType = $("#error_dropdown").val();
        if(errorType == "masserror") {
            errorVal = parseFloat($("#errorval").val().trim());
            massErrorthVal = errorVal ;
        }
        else {
            errorVal = parseFloat($("#errorval").val().trim());
            ppmErrorthVal = errorVal ;
        }
        let executionObj = new SeqOfExecution();
        executionObj.sequenceOfExecution(errorType,errorVal,"");
        document.getElementById("totalseqmass_h6").style.display = "block";
        document.getElementById("massvariation_h6").style.display = "block";
    })
    /**
     * On Click action to hide and show the table of calculate theoretical
     * masses with matched and unmatched masses
     */
    $("#hide_table").click(function(){
        let text_val = $("#hide_table").html().trim();
        if( text_val == "Hide Table"){
            $("#hide_table").html("Show Table");
            $("#divtableContainer").hide();
        }
        if(text_val == "Show Table")
        {
            $("#hide_table").html("Hide Table");
            $("#divtableContainer").show();
        }
    })
}
/**
 * Function to display all peaks of data in table. This handles on click action
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
 * Function to display only matched peaks in table. This handles on click action 
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
 * Function to display only un matched peaks in table. This handles on click action
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