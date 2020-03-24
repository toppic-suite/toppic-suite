/**
 * Code to set width of the monomass tables and sorting the table based on the first column
 */ 
 $(document).ready(function() {
  $('#spectrum').dataTable( {
    "scrollY":        "400px",
    "scrollCollapse": true,
    "paging":         false,
    "order": [[ 1, "asc" ]],
    "bSortClasses": false
  } );
} ); 	

/**
 * Function to show only matched peaks on click of matched peaks
 */
function showMatchedPeaks() {
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
 * Function to show only unmatched peaks on click of unmatched peaks
 */
function showNotMatchedPeaks() {
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
 * Function to show all peaks on click of All peaks
 */
function showAllPeaks() {
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
 * Function to add all the peaks on load of the page
 */
function showIonPeaks(ids) {
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