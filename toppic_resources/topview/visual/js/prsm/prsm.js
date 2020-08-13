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
 * Function to show only matched peaks on click of matched peaks button
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
  //$('div.dataTables_scrollBody').height(400);
}

/**
 * Function to show only unmatched peaks on click of unmatched peaks button
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
  //$('div.dataTables_scrollBody').height(400);
}

/**
 * Function to show all peaks on click of All peaks button
 */
function showAllPeaks() {
  var elems = document.getElementsByClassName('matched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "";
  }
  elems = document.getElementsByClassName('unmatched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "";
  }
  //$('div.dataTables_scrollBody').height(400);
}
/**
 * This gets invoked on click of annotation in the SVG of sequence at matched positions
 * Function to show only ions matched at a particular position
 * @param {String} ids - contains name of the tag
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
