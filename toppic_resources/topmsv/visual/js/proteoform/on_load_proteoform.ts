"use strict";
/**
 * This function waits till all the HTML tags are loaded.
 * Loads realted protoform based on the proteoform id
 */
$(document).ready(function () {
    // Get the URL from browser
    let x: string = location.href;
    let l_split: string = x.split("?")[1];
    let pathAndVal: string[] = l_split.split("&");
    let folderpath: string = pathAndVal[0].split("=")[1];
    // get the proteoform Id number by splitting url with "?" and "="
    let proteoform_id: string = pathAndVal[1].split("=")[1];
    let head: HTMLHeadElement = document.getElementsByTagName('head')[0];
    let script: HTMLScriptElement = document.createElement('script');
    // get the proteoform Id number by splitting url with ? and = 
    let file_name: string = "../data/" + folderpath + "/proteoforms/proteoform" + proteoform_id + ".js";
    script.type = 'text/javascript';
    script.src = file_name;
    // include the proteoform data_file  
    head.appendChild(script);
    // Wait till the proteoform data file is loaded
    script.onload = function () {
        geneProteoformObj(function (prsm: Prsm[]) {
            // Function call to generate url to navigate to protein page and all proteins page
            proteoformUrl(folderpath, prsm);
            // Function call that generates table of content with prsm and navigation link to prsm
            createTableData(folderpath, prsm);
        });
    };
});
