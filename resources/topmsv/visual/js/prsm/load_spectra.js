"use strict";
/**
 * @function loadMsOne
 * @description - This function load an MS1 spectrum.
 */
function loadMsOne(ms1Spec, ms1SvgId) {
    if (!ms1Spec) {
        console.error("ERROR: no ms1 spectra data");
        return;
    }
    let ions = [];
    let spectrumDataPeaks = new SpectrumFunction();
    let spectrumDataEnvs = new SpectrumFunction();
    spectrumDataPeaks.assignLevelPeaks(ms1Spec.getPeaks());
    spectrumDataEnvs.assignLevelEnvs(ms1Spec.getEnvs());
    let spGraph = new SpectrumView(ms1SvgId, ms1Spec.getPeaks());
    spGraph.addRawSpectrumAnno(ms1Spec.getEnvs(), ions);
    let precMonoMz = ms1Spec.getPrecMz();
    spGraph.getPara().updateMzRange(precMonoMz);
    spGraph.getPara().setHighlight(ms1Spec);
    spGraph.redraw();
}
function loadMsTwo(prsmObj, ms2GraphList, divId, navId) {
    let ms2Spec = prsmObj.getMs2Spectra();
    if (!ms2Spec) {
        console.error("ERROR: no ms2 spectra data");
        return;
    }
    let graphList = [];
    let monoGraphList = [];
    let deconvPeaks = [];
    ms2Spec.forEach((spectrum, index) => {
        let deconvPeakTemp = spectrum.getDeconvPeaks();
        if (deconvPeakTemp) {
            deconvPeaks = deconvPeaks.concat(deconvPeakTemp); //prepare deconvPeaks before graph drawing
        }
    });
    ms2Spec.forEach((spectrum, index) => {
        createMs2NavElement(index, divId, navId, spectrum.getScanNum());
        let show = false;
        if (index == 0) {
            show = true;
        }
        let svgId = divId + "_graph_" + index;
        createSvg(show, divId, svgId, "ms2_svg_graph_class");
        let [ions, monoIons] = getIons(prsmObj.getMatchedPeakEnvelopePairs(), spectrum.getSpectrumId());
        let spGraph = new SpectrumView(svgId, spectrum.getPeaks());
        spGraph.addRawSpectrumAnno(spectrum.getEnvs(), ions);
        let spectrumDataPeaks = new SpectrumFunction();
        let spectrumDataEnvs = new SpectrumFunction();
        spectrumDataPeaks.assignLevelPeaks(spectrum.getPeaks());
        spectrumDataEnvs.assignLevelEnvs(spectrum.getEnvs());
        spGraph.redraw();
        ms2GraphList.push(spGraph);
        //mono mass svg
        let monoSvgId = divId + "_mono_graph_" + index;
        show = false;
        createSvg(show, divId, monoSvgId, "ms2_svg_graph_class");
        if (!deconvPeaks) {
            console.error("ERROR: no deconvoluted peaks in ms2 spectrum");
            return;
        }
        let monoMasses = getMonoMasses(deconvPeaks);
        let nIonType = spectrum.getNTerminalIon()[0].getName();
        let cIonType = spectrum.getCTerminalIon()[0].getName();
        let spectrumDataMonoPeaks = new SpectrumFunction();
        spectrumDataMonoPeaks.assignLevelPeaks(monoMasses);
        let monoSpGraph = new SpectrumView(monoSvgId, monoMasses, prsmObj.getProteoform().getSeq().length);
        monoSpGraph.addMonoMassSpectrumAnno(monoIons, prsmObj.getProteoform(), nIonType, cIonType);
        monoSpGraph.getPara().setMonoMassGraph(true);
        monoSpGraph.redraw();
        monoGraphList.push(monoSpGraph);
    });
    //add a tab for checkbox for mono graph
    addCheckboxTab(navId);
    //add an event listner for checkbox
    $("#checkbox-anno-line").on("change", function () {
        for (let i = 0; i < ms2GraphList.length; i++) {
            let monoListId = "ms2_svg_div_monographlist_" + i;
            let monoGraphId = "ms2_svg_div_mono_graph_" + i;
            let monoListElement = document.getElementById(monoListId);
            let monoGraphElement = document.getElementById(monoGraphId);
            if (monoListElement) {
                if (monoListElement.classList.contains("active")) {
                    monoGraphList[i].redraw();
                }
                ;
            }
        }
    });
    // add action for nav bar
    $(".ms2_graph_list").on("click", function () {
        let ms2Id = this.id;
        //console.log("ms2id", ms2Id);
        let ms2Split = ms2Id.split("_");
        let ms2Index = parseInt(ms2Split[ms2Split.length - 1]);
        let type = ms2Split[ms2Split.length - 2];
        for (let i = 0; i < ms2GraphList.length; i++) {
            let listId = "ms2_svg_div_graphlist_" + i;
            let monoListId = "ms2_svg_div_monographlist_" + i;
            let graphId = "ms2_svg_div_graph_" + i;
            let monoGraphId = "ms2_svg_div_mono_graph_" + i;
            //console.log(listId, graphId);
            let listElement = document.getElementById(listId);
            let monoListElement = document.getElementById(monoListId);
            let graphElement = document.getElementById(graphId);
            let monoGraphElement = document.getElementById(monoGraphId);
            let checkboxTab = document.getElementById("checkbox-tab");
            if (i == ms2Index) {
                if (type == "graphlist") {
                    if (listElement && monoListElement && graphElement && monoGraphElement) {
                        listElement.classList.add("active");
                        monoListElement.classList.remove("active");
                        graphElement.style.display = "";
                        monoGraphElement.style.display = "none";
                        if (checkboxTab != null) {
                            checkboxTab.style.display = "none";
                        }
                    }
                    else {
                        console.error("ERROR: graph ID is invalid");
                    }
                }
                else {
                    if (listElement && monoListElement && graphElement && monoGraphElement) {
                        listElement.classList.remove("active");
                        monoListElement.classList.add("active");
                        graphElement.style.display = "none";
                        monoGraphElement.style.display = "";
                        if (checkboxTab != null) {
                            checkboxTab.style.display = "";
                        }
                    }
                    else {
                        console.error("ERROR: graph ID is invalid");
                    }
                }
            }
            else {
                if (listElement && monoListElement && graphElement && monoGraphElement) {
                    listElement.classList.remove("active");
                    monoListElement.classList.remove("active");
                    graphElement.style.display = "none";
                    monoGraphElement.style.display = "none";
                    if (checkboxTab) {
                        checkboxTab.style.display = "none";
                    }
                }
                else {
                    console.error("ERROR: graph ID is invalid");
                }
            }
        }
    });
    //if the below lines are outside this scope, they execute before
    //graphList and monoGraphList are returned with valid values
    let saveSpectrumObj = new SaveSpectrum(ms2GraphList, monoGraphList);
    saveSpectrumObj.main();
}
function addCheckboxTab(navId) {
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
/**
 * Function to Create Navigation buttons to navigate between spectrums
 * @param {Array} scanidList - Contains scan Id List
 * @param {String} id - Contains Id of the avg on which spectrum to be drawn
 */
function createMs2NavElement(i, divId, navId, specScan) {
    let ul = document.getElementById(navId);
    let li = document.createElement("li");
    let li_id = divId + "_graphlist_" + i;
    if (!ul) {
        console.error("ERROR: invalid navId");
        return;
    }
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
    a.innerHTML = "Scan " + specScan + " masses";
    li.appendChild(a);
    ul.appendChild(li);
}
/**
 * This generates spectrum for each spec Id
 * @param {String} divId - Contains Id of the div tag under which the monomass graphs are drawn
 * @param {String} svgId - Contains id as "monoMassSvg_" to which scan Id is added
 * @param {String} className - Contains class name to which the corresponding svg graphs are drawn
 */
function createSvg(show, divId, svgId, className) {
    let div = document.getElementById(divId);
    let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("id", svgId);
    svg.setAttribute("class", className);
    svg.style.backgroundColor = "#F8F8F8";
    if (show) {
        svg.style.display = "";
    }
    else {
        svg.style.display = "none";
    }
    if (!div) {
        console.error("ERROR: invalid divId");
        return;
    }
    div.appendChild(svg);
}
function getMonoMasses(peaks) {
    let masses = [];
    for (let i = 0; i < peaks.length; i++) {
        let monoMass = peaks[i].getMonoMass();
        let charge = peaks[i].getCharge();
        if (!monoMass || !charge) {
            console.error("ERROR: invalid mono mass or charge in deconvPeaks");
            return masses;
        }
        let peakObj = new Peak(peaks[i].getId(), monoMass, peaks[i].getMonoMz(), peaks[i].getIntensity(), monoMass, charge);
        masses.push(peakObj);
    }
    return masses;
}
/**
 * @function getIonData
 * @description gets ion data list with mz, intensity and ion name
 * This function gets matched ion data
 * @param {object} prsm_data - contains complete data of prsm
 * @param {int} specId - contains information of the spec Id
 * @param {object} json_data - contains complete data of spectrum
 */
function getIons(matchedPeakEnvPairs, specId) {
    let ions = [];
    let monoIons = [];
    matchedPeakEnvPairs.forEach((pair) => {
        if (specId == pair.getPeak().getSpecId()) {
            let ionData;
            let monoIonData;
            let ionText = "";
            let massError;
            let envPeaks = [];
            let env = pair.getEnvelope();
            if (env) {
                envPeaks = env.getPeaks();
            }
            envPeaks.sort(function (x, y) {
                return y.getIntensity() - x.getIntensity();
            });
            let x = envPeaks[0].getMonoMz();
            let y = envPeaks[0].getIntensity();
            let matchedIon = pair.getIon();
            let ionType = matchedIon.getName();
            if (ionType == "Z_DOT") {
                ionType = "Z\u02D9";
            }
            let pos = matchedIon.getId().match(/\d+/);
            if (!pos) {
                console.error("ERROR: invalid matched ion position");
                return;
            }
            ionText = ionType + pos[0];
            if (!matchedIon.getMassError() && matchedIon.getMassError() != 0) { //mass error can be 0
                console.error("ERROR: mass error is undefined");
                return;
            }
            let pairEnv = pair.getEnvelope();
            if (!pairEnv) {
                console.error("Error: invalid envelope");
                return [[], []];
            }
            //@ts-ignore
            massError = matchedIon.getMassError(); //undefined already checked above
            ionData = { "mz": x, "intensity": y, "text": ionText, "error": massError, "env": pairEnv };
            ions.push(ionData);
            let monoX = pairEnv.getMonoMass();
            let monoY = 0;
            envPeaks.forEach(element => monoY += element.getIntensity());
            monoIonData = { "mz": monoX, "intensity": monoY, "text": ionText, "pos": pos[0], "error": massError };
            addOneIon(monoIons, monoIonData);
        }
    });
    return [ions, monoIons];
}
function addOneIon(ionList, ion) {
    let idx = -1;
    for (let i = 0; i < ionList.length; i++) {
        if (ion.text == ionList[i].text) {
            idx = i;
            break;
        }
    }
    if (idx == -1) {
        ionList.push(ion);
    }
    else {
        if (ion.intensity > ionList[idx].intensity) {
            ionList[idx].intensity = ion.intensity;
        }
    }
}
