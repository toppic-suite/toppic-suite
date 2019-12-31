$(document).ready(function() {
  $('#spectrum').dataTable( {
    "scrollY":        "400px",
    "scrollCollapse": true,
    "paging":         false,
    "order": [[ 1, "asc" ]],
    "bSortClasses": false
  } );
} ); 	

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