"use strict";
/**
 * @function loadMsOne
 * @description - This function load an MS1 spectrum.
 */

function loadMsOne(ms1Spec: Spectrum | null, ms1SvgId: string): void{
    if(!ms1Spec) {
        console.error("ERROR: no ms1 spectra data");
        return;
    }
    let ions: MatchedIon[] = [];
    let spectrumDataPeaks: SpectrumFunction = new SpectrumFunction();
    let spectrumDataEnvs: SpectrumFunction = new SpectrumFunction();
    spectrumDataPeaks.assignLevelPeaks(ms1Spec.getPeaks());
    spectrumDataEnvs.assignLevelEnvs(ms1Spec.getEnvs());
    let spGraph: SpectrumView = new SpectrumView(ms1SvgId, ms1Spec.getPeaks());
    spGraph.addRawSpectrumAnno(ms1Spec.getEnvs(), ions);
    let precMonoMz: number = ms1Spec.getPrecMz();
    spGraph.getPara().updateMzRange(precMonoMz);
    spGraph.getPara().setHighlight(ms1Spec);
    spGraph.redraw();
}
function loadMsTwo(prsmObj: Prsm, ms2GraphList: SpectrumView[], divId: string, navId: string): void {
    let ms2Spec: Spectrum[] | null = prsmObj.getMs2Spectra();
    if(!ms2Spec) {
        console.error("ERROR: no ms2 spectra data");
        return;
    }
    let graphList: SpectrumView[] = [];
    let monoGraphList: SpectrumView[] = [];
    let deconvPeaks: Peak[] = [];
    ms2Spec.forEach((spectrum: Spectrum, index) => {
        let deconvPeakTemp: Peak[] | null = spectrum.getDeconvPeaks();
        if (deconvPeakTemp) {
            deconvPeaks = deconvPeaks.concat(deconvPeakTemp); //prepare deconvPeaks before graph drawing
        }
    });
    ms2Spec.forEach((spectrum, index) => {
        createMs2NavElement(index, divId, navId, spectrum.getScanNum());
        let show: boolean = false;
        if (index == 0) {
            show = true;
        }
        let svgId: string = divId + "_graph_" + index;
        createSvg(show, divId, svgId, "ms2_svg_graph_class");
        let [ions, monoIons] = getIons(prsmObj.getMatchedPeakEnvelopePairs(), spectrum.getSpectrumId());
        let spGraph: SpectrumView = new SpectrumView(svgId, spectrum.getPeaks());
        spGraph.addRawSpectrumAnno(spectrum.getEnvs(), ions);
        let spectrumDataPeaks: SpectrumFunction = new SpectrumFunction();
        let spectrumDataEnvs: SpectrumFunction = new SpectrumFunction();
        spectrumDataPeaks.assignLevelPeaks(spectrum.getPeaks());
        spectrumDataEnvs.assignLevelEnvs(spectrum.getEnvs());
        spGraph.redraw();
        ms2GraphList.push(spGraph);

        //mono mass svg
        let monoSvgId: string = divId + "_mono_graph_" + index;
        show = false;
        createSvg(show, divId, monoSvgId, "ms2_svg_graph_class");

        if (!deconvPeaks) {
            console.error("ERROR: no deconvoluted peaks in ms2 spectrum");
            return;
          }
    
        let monoMasses: Peak[] = getMonoMasses(deconvPeaks);

        let nIonType: string = spectrum.getNTerminalIon()[0].getName();
        let cIonType: string = spectrum.getCTerminalIon()[0].getName();

        let spectrumDataMonoPeaks: SpectrumFunction = new SpectrumFunction();
        spectrumDataMonoPeaks.assignLevelPeaks(monoMasses);
        let monoSpGraph: SpectrumView = new SpectrumView(monoSvgId, monoMasses, prsmObj.getProteoform().getSeq().length);
        monoSpGraph.addMonoMassSpectrumAnno(monoIons, prsmObj.getProteoform(), nIonType, cIonType);
        monoSpGraph.getPara().setMonoMassGraph(true);
        monoSpGraph.redraw();
        monoGraphList.push(monoSpGraph);
    })
    //add a tab for checkbox for mono graph
    addCheckboxTab(navId);

    //add an event listner for checkbox
    $("#checkbox-anno-line").on("change", function() {
        for (let i = 0; i < ms2GraphList.length; i++) {
            let monoListId: string = "ms2_svg_div_monographlist_" + i;
            let monoGraphId: string = "ms2_svg_div_mono_graph_" + i;
            let monoListElement: HTMLElement | null = document.getElementById(monoListId);
            let monoGraphElement: HTMLElement | null = document.getElementById(monoGraphId);
            if (monoListElement) {
                if (monoListElement.classList.contains("active")) {
                    monoGraphList[i].redraw();
                };
            }
        }
    })

    // add action for nav bar
    $(".ms2_graph_list").on("click", function() {
        let ms2Id: string = this.id;
        //console.log("ms2id", ms2Id);
        let ms2Split: string[] = ms2Id.split("_");
        let ms2Index: number = parseInt(ms2Split[ms2Split.length - 1]);
        let type: string = ms2Split[ms2Split.length - 2];
        for (let i = 0; i < ms2GraphList.length; i++) {
            let listId: string = "ms2_svg_div_graphlist_" + i;
            let monoListId: string = "ms2_svg_div_monographlist_" + i;
            let graphId: string = "ms2_svg_div_graph_" + i;
            let monoGraphId: string = "ms2_svg_div_mono_graph_" + i;
            //console.log(listId, graphId);
            let listElement: HTMLElement | null = document.getElementById(listId);
            let monoListElement: HTMLElement | null = document.getElementById(monoListId);
            let graphElement: HTMLElement | null = document.getElementById(graphId);
            let monoGraphElement: HTMLElement | null = document.getElementById(monoGraphId);
            let checkboxTab: HTMLElement | null = document.getElementById("checkbox-tab");
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
                    else{
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
                    else{
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
                else{
                    console.error("ERROR: graph ID is invalid");
                }
            }
        }
    });
    //if the below lines are outside this scope, they execute before
    //graphList and monoGraphList are returned with valid values
    let saveSpectrumObj: SaveSpectrum = new SaveSpectrum(ms2GraphList, monoGraphList);
    saveSpectrumObj.main();
}
function addCheckboxTab(navId: string): void {
    let ul: HTMLElement | null = document.getElementById(navId);

    if (!ul) {
        console.error("ERROR: invalid navId");
        return;
    }

    let li: HTMLLIElement = document.createElement("li");
    //let li_id: string = "checkbox-tab";
    //li.setAttribute("id", li_id);
    li.setAttribute("class", "nav-item");

    let div: HTMLDivElement = document.createElement("div");
    div.setAttribute("class", "nav-link");
    div.setAttribute("id", "checkbox-tab");
    div.style.display = "none";

    let checkbox: HTMLInputElement = document.createElement("input");
    checkbox.setAttribute("type", "checkbox");
    checkbox.setAttribute("id", "checkbox-anno-line");
    checkbox.setAttribute("name", "checkbox-anno-line");
    checkbox.setAttribute("checked", "true");
    checkbox.setAttribute("value", "true");

    let label: HTMLLabelElement = document.createElement("label");
    label.setAttribute("for", "checkbox-anno-line");

    let text: Text = document.createTextNode("Show annotation lines");

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
function createMs2NavElement(i: number, divId: string, navId: string, specScan: string): void {
    let ul: HTMLElement | null = document.getElementById(navId);
    let li: HTMLLIElement = document.createElement("li");
    let li_id: string = divId + "_graphlist_" + i;
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
    let a: HTMLAnchorElement = document.createElement("a");
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
function createSvg(show: boolean, divId: string, svgId: string, className: string): void {
    let div: HTMLElement | null = document.getElementById(divId);
    let svg: SVGSVGElement = document.createElementNS("http://www.w3.org/2000/svg", "svg");
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
function getMonoMasses(peaks: Peak[]): Peak[] {
    let masses: Peak[] = [];
    
    for (let i = 0; i < peaks.length; i++) {
        let monoMass: number | undefined = peaks[i].getMonoMass();
        let charge: number | undefined = peaks[i].getCharge();
        if (!monoMass || !charge) {
            console.error("ERROR: invalid mono mass or charge in deconvPeaks");
            return masses;
        }
        let peakObj: Peak = new Peak(peaks[i].getId(), monoMass, peaks[i].getMonoMz(), peaks[i].getIntensity(), monoMass, charge);
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
function getIons(matchedPeakEnvPairs: MatchedPeakEnvelopePair[], specId: string): [MatchedIon[], MatchedIon[]] {
    let ions: MatchedIon[] = [];
    let monoIons: MatchedIon[] = [];
    matchedPeakEnvPairs.forEach((pair) => {
        if (specId == pair.getPeak().getSpecId()) {
            let ionData: MatchedIon;
            let monoIonData: MatchedIon;
            let ionText: string = "";
            let massError: number;
            let envPeaks: Peak[] = [];
            let env: Envelope | null = pair.getEnvelope();
            if (env) {
                envPeaks = env.getPeaks();
            }
            envPeaks.sort(function (x, y) {
                return y.getIntensity() - x.getIntensity();
            });
            let x: number = envPeaks[0].getMonoMz();
            let y: number = envPeaks[0].getIntensity();
            let matchedIon: Ion = pair.getIon();
            let ionType: string = matchedIon.getName();
            if (ionType == "Z_DOT") {
                ionType = "Z\u02D9";
            }
            let pos: RegExpMatchArray | null = matchedIon.getId().match(/\d+/);
            if (!pos) {
                console.error("ERROR: invalid matched ion position");
                return;
            }
            ionText = ionType + pos[0];
            if (!matchedIon.getMassError() && matchedIon.getMassError() != 0) { //mass error can be 0
                console.error("ERROR: mass error is undefined");
                return;
            }
            let pairEnv: Envelope | null = pair.getEnvelope();
            if (!pairEnv) {
                console.error("Error: invalid envelope");
                return [[], []];
            }
            //@ts-ignore
            massError = matchedIon.getMassError(); //undefined already checked above
            ionData = { "mz": x, "intensity": y, "text": ionText, "error": massError, "env": pairEnv };
            ions.push(ionData);
            let monoX: number = pairEnv.getMonoMass();
            let monoY: number = 0;
            envPeaks.forEach(element => monoY += element.getIntensity());
            monoIonData = { "mz": monoX, "intensity": monoY, "text": ionText, "pos": pos[0], "error": massError };
            addOneIon(monoIons, monoIonData);
        }
    });
    return [ions, monoIons];
}
function addOneIon(ionList: MatchedIon[], ion: MatchedIon): void {
    let idx: number = -1;
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
