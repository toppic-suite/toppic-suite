"use strict";
/**
 * Function to Create Navigation buttons to navigate between two types of graph
 * createMs2NavElement from visual/prsm/load_spectra.js
 * @param {Array} divId - id of the div containing nav bar
 * @param {String} navId - id of the nav bar
 * These functions need to be rewritten as a library later
 */
function clearMs2NavElement(navId) {
    $("#" + navId).empty();
}
function addCheckboxTabInspect(navId) {
    let ul = document.getElementById(navId);
    if (!ul) {
        console.error("ERROR: invalid navId");
        return;
    }
    let li = document.createElement("li");
    //let li_id: string = "checkbox-tab";
    //li.setAttribute("id", li_id);
    li.setAttribute("class", "nav-item");
    let div = document.createElement("div");
    div.setAttribute("class", "nav-link");
    div.setAttribute("id", "checkbox-tab");
    div.style.display = "none";
    let checkbox = document.createElement("input");
    checkbox.setAttribute("type", "checkbox");
    checkbox.setAttribute("id", "checkbox-anno-line");
    checkbox.setAttribute("name", "checkbox-anno-line");
    checkbox.setAttribute("checked", "true");
    checkbox.setAttribute("value", "true");
    let label = document.createElement("label");
    label.setAttribute("for", "checkbox-anno-line");
    let text = document.createTextNode("Show annotation lines");
    label.appendChild(text);
    div.appendChild(checkbox);
    div.appendChild(label);
    li.appendChild(div);
    ul.appendChild(li);
}
function createMs2NavElementInspect(i, divId, navId, specScan) {
    let ul = document.getElementById(navId);
    let li = document.createElement("li");
    let li_id = divId + "_graphlist_" + i;
    li.setAttribute("id", li_id);
    if (i == 0) {
        li.setAttribute("class", "nav-item ms2_graph_list active");
    }
    else {
        li.setAttribute("class", "nav-item ms2_graph_list");
    }
    let a = document.createElement("a");
    a.setAttribute("class", "nav-link");
    a.setAttribute("href", "#!");
    a.innerHTML = "Scan " + specScan;
    li.appendChild(a);
    ul.appendChild(li);
    li = document.createElement("li");
    li_id = divId + "_monographlist_" + i;
    li.setAttribute("id", li_id);
    li.setAttribute("class", "nav-item ms2_graph_list");
    a = document.createElement("a");
    a.setAttribute("class", "nav-link");
    a.setAttribute("href", "#!");
    a.innerHTML = "Mass scan " + specScan;
    li.appendChild(a);
    ul.appendChild(li);
}
function addEventNavBar(monoMassGraphObj) {
    // add action for nav bar
    $(".ms2_graph_list").click((e) => {
        let Id = e.currentTarget.id;
        let ms2GraphList = document.getElementsByClassName("ms2_graph_list");
        for (let i = 0; i < ms2GraphList.length; i++) {
            if (ms2GraphList[i].id == Id) {
                ms2GraphList[i].classList.add("active");
                let svgId = ms2GraphList[i].id;
                let svgIdSplit = svgId.split("_");
                let type = svgIdSplit[3];
                let spectrumTab = document.getElementById(Constants.SPECTRUMGRAPHID);
                let monoMassTab = document.getElementById(Constants.MONOMASSGRAPHID);
                let checkboxTab = document.getElementById("checkbox-tab");
                if (type == "graphlist") {
                    spectrumTab.style.display = "";
                    monoMassTab.style.display = "none";
                    if (checkboxTab != null) {
                        checkboxTab.style.display = "none";
                    }
                }
                else {
                    spectrumTab.style.display = "none";
                    monoMassTab.style.display = "";
                    if (checkboxTab != null) {
                        checkboxTab.style.display = "";
                    }
                }
            }
            else {
                ms2GraphList[i].classList.remove("active");
            }
        }
    });
    //add an event listner for checkbox
    $("#checkbox-anno-line").on("change", function () {
        monoMassGraphObj.redraw();
    });
}
function switchTab(graphType) {
    let ms2GraphList = document.getElementsByClassName("ms2_graph_list");
    for (let i = 0; i < ms2GraphList.length; i++) {
        let svgId = ms2GraphList[i].id;
        let svgIdSplit = svgId.split("_");
        let type = svgIdSplit[3];
        if (type == graphType) {
            ms2GraphList[i].classList.add("active");
        }
        else {
            ms2GraphList[i].classList.remove("active");
        }
    }
    if (graphType == "monographlist") {
        $("#" + Constants.SPECTRUMGRAPHID).hide();
        $("#" + Constants.MONOMASSGRAPHID).show();
    }
    else {
        $("#" + Constants.SPECTRUMGRAPHID).show();
        $("#" + Constants.MONOMASSGRAPHID).hide();
    }
}
