"use strict";
/**
 * This function waits till all the HTML tags are loaded.
 * Loads the related proteins.js data file
 */
$(document).ready(function () {
    // Get the URL from browser
    let x: string = location.href;
    let l_split: string = x.split("?")[1];
    let pathAndVal: string[] = l_split.split("&");
    let folderpath: string = pathAndVal[0].split("=")[1];
    // get the proteoform Id number by splitting url with "?" and "=" 
    let protein_seq_num: string = pathAndVal[1].split("=")[1];
    let head: HTMLHeadElement = document.getElementsByTagName('head')[0];
    let script: HTMLScriptElement = document.createElement('script');
    let file_name: string = "../data/" + folderpath + "/proteins/protein" + protein_seq_num + ".js";
    // Set All proteins URL path
    let allproteinUrl: string = "proteins.html?folder=" + folderpath;
    $("#allproteinurl_1").attr("href", allproteinUrl);
    $("#allproteinurl_2").attr("href", allproteinUrl);
    script.type = 'text/javascript';
    script.src = file_name;
    // include the proteoform json data_file as script
    document.head.appendChild(script);
    /*script.onload = function () {
      // invoke protein function to form HTML elements from protein.js
      //@ts-ignore
      protein(folderpath, prsm_data);
    };*/
    script.onload = function () {
        geneProteinObj(function (prsmObj: Prsm[], prsmCnts: number[]) {
            protein(folderpath, prsmObj, prsmCnts);
        });
    };
});
