"use strict";
/**
 * draw navigation bar while avoiding cross origin issue
 */
const drawNav = function () {
    //depending on where it is being called... is it index html?
    let x = location.href;
    let n = x.lastIndexOf("/");
    let htmlName = x.substring(n + 1, x.length);
    //the path to the other pages are different for index.html only
    let relPath;
    if (htmlName == "index.html") {
        relPath = " ";
    }
    else {
        relPath = "../";
    }
    let navCode = '<nav class="navbar navbar-expand-lg fixed-top navbar-dark bg-dark" id="nav-div"> \
  <div id="for-flex" style="display:flex; width: 100%">\
  <div class=" navcontainer">\
  <button class="navbar-toggler" type="button" data-toggle="collapse" data-target=".navbar-collapse" aria-controls="navbar-collapse" aria-expanded="false" aria-label="Toggle navigation">\
  <span class="navbar-toggler-icon"></span></button>\
  <div class="collapse navbar-collapse">\
  <a class="navbar-brand logo" href="' + relPath + 'index.html"><span class="headBlock"><h3><strong id="toppic_icon">T</strong>opMSV</h3></span></a>\
  <ul class="navbar-nav mr-auto mt-2 mt-lg-0"><li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'visual/proteins.html">Protein Identifications</a></li>\
  <li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'visual/ms.html">Spectrum Identifications</a></li>\
  <li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'inspect/spectrum.html">Visual Inspection</a></li></ul>\
  </div></div></div></nav>';
    $(function () {
        $("#nav-bar").html(navCode);
    });
}();
