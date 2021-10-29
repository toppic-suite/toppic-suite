"use strict";
$(document).ready(function () {
    //go to the prsm folder and parse ms1 and ms2 spectrum numbers from each prsm
    /*$(".card-row").css('display', 'none');
    let x: string = location.href;
    let l_split: string[] = x.split(/[?#]+/);
    let folderName: string = l_split[1].split("=")[1];
    let head: HTMLHeadElement = document.getElementsByTagName('head')[0];
    let script: HTMLScriptElement = document.createElement('script');
    let file_name: string = folderName + "/prsms.js";
    script.type = 'text/javascript';
    script.src = file_name;
    head.appendChild(script);
    // Wait till the data is loaded from prsms.js and start executing the code
    script.onload = function () {
      // Function builds the complete HTML
      allSpectrum(folderName);
    };*/
    $(".card-row").css('display', 'none');
    $(".search-box").css('display', 'none');
    let x: string = location.href;
    // get the folder path by splitting url with "?","=","#"
    let l_split: string[] = x.split(/[?#]+/);
    if (l_split.length != 2) {
        $(".card-row").css('display', 'block');
        $(".card").hover(function () {
            $(this).css('cursor', 'pointer');
        });
        $(".card-body").click(function () {
            let folder: string | undefined = $(this).attr("id");
            let url: string = "ms.html?folder=../../" + folder + "/data_js";
            window.open(url, "_self");
        });
    }
    else {
        $(".search-box").css('display', 'inline-block');
        let folderName: string = l_split[1].split("=")[1];
        let finalPath: string = folderName;
        let head: HTMLHeadElement = document.getElementsByTagName('head')[0];
        let script: HTMLScriptElement = document.createElement('script');
        let file_name: string = finalPath + "/prsms.js";
        script.type = 'text/javascript';
        script.src = file_name;
        head.appendChild(script);
        // Wait till the data is loaded from proteins.js and start executing the code
        script.onload = function () {
          // Function builds the complete HTML
          allSpectrum(folderName);
        };
    }
});
