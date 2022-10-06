"use strict";
/**
 * This function waits till all the HTML tags are loaded.
 * Loads realted protoform based on the proteoform id
 */
$(document).ready(function () {
    // Get the URL from browser
    let x = location.href;
    let l_split = x.split("?")[1];
    let pathAndVal = l_split.split("&");
    let folderpath = pathAndVal[0].split("=")[1];
    // get the proteoform Id number by splitting url with "?" and "="
    let proteoform_id = pathAndVal[1].split("=")[1];
    let head = document.getElementsByTagName('head')[0];
    let script = document.createElement('script');
    // get the proteoform Id number by splitting url with ? and = 
    let file_name = "../data/" + folderpath + "/proteoforms/proteoform" + proteoform_id + ".js";
    script.type = 'text/javascript';
    script.src = file_name;
    // include the proteoform data_file  
    head.appendChild(script);
    // Wait till the proteoform data file is loaded
    script.onload = function () {
        geneProteoformObj(function (prsm) {
            // Function call to generate url to navigate to protein page and all proteins page
            proteoformUrl(folderpath, prsm);
            // Function call that generates table of content with prsm and navigation link to prsm
            createTableData(folderpath, prsm);
        });
    };
});
