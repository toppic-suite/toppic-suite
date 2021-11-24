"use strict";
/**
 * Build Title,Description,URL and get best prsm for each proteoform
 * @param {String} folderpath - Provides path to the data folder
 */
function protein(folderpath, prsm_data, prsmCnts) {
    // Generate title for the webpage from the data
    document.title = `Proteoforms for protein ${prsm_data[0].getProteoform().getProtName()} \ 
	${prsm_data[0].getProteoform().getProtDesc()}`;
    // Generate description of the protein from the data
    let seqDesc = document.getElementById('sequence_description');
    if (seqDesc) {
        seqDesc.innerHTML = `${prsm_data.length} proteoforms for protein \ 
		${prsm_data[0].getProteoform().getProtName()} ${prsm_data[0].getProteoform().getProtDesc()}`;
    }
    prsm_data.forEach((prsm, index) => {
        proteoformToHtml(prsm, index, folderpath, prsmCnts[index]);
    });
}
/**
 * This function builds a header for a proteoform, gets information of best prsm.
 * Forms SVG visualization of the best prsm.
 * @param {object} compatible_proteoform - Contains information of a single proteoform
 * @param {Int} index - Index of proteoform
 * @param {String} folderpath - Provides path to data files
 */
function proteoformToHtml(prsm, index, folderpath, prsmCnt) {
    // Get the div element with class name proteoformcontainer
    // All the details of a each proteoform is under the proteoformcontainer div
    let div_container = document.getElementsByClassName("proteoformcontainer")[0];
    let div = document.createElement('div');
    let id = "p" + prsm.getProteoform().getId();
    div.setAttribute("id", id);
    let h2 = document.createElement('h3');
    let p = document.createElement("p");
    let precMass = prsm.getPrecMass();
    if (prsmCnt > 1 && precMass) {
        // Forms header for a proteoform
        p = Build_BestPrSM(prsm.getEValue(), precMass, prsm.getId(), prsm.getProteoform().getId(), prsmCnt, folderpath);
        h2.innerHTML = `Proteoform #${prsm.getProteoform().getId()} Feature intensity: ${prsm.getFeatureInte()}`;
    }
    else {
        // Forms header for a proteoform
        h2.innerHTML = `Proteoform #${prsm.getProteoform().getId()} Feature intensity: ${prsm.getFeatureInte()}`;
        p.setAttribute("style", "font-size:16px;");
        let text1 = document.createElement("text");
        text1.innerHTML = "There is only ";
        p.appendChild(text1);
        let a_prsm = document.createElement("a");
        a_prsm.href = "prsm.html?folder=" + folderpath + "&prsm_id=" + prsm.getId();
        a_prsm.innerHTML = " 1 PrSM ";
        p.appendChild(a_prsm);
        let text2 = document.createElement("text");
        text2.innerHTML = `with an E-value ${prsm.getEValue()} and a precursor mass ${prsm.getPrecMass()}.`;
        p.appendChild(text2);
    }
    div_container.appendChild(div);
    div_container.appendChild(h2);
    div_container.appendChild(p);
    let containerId = "svg_container" + index;
    let svgContainer = document.createElement("div");
    svgContainer.setAttribute("id", containerId);
    svgContainer.setAttribute("class", "svg_container");
    div_container.appendChild(svgContainer);
    let svgId = "prsm_svg" + index;
    d3.select("#" + containerId).append("svg")
        .attr("id", svgId);
    //console.log(svgId, BestPrSM);
    let graph = new PrsmView(svgId, prsm);
    graph.redraw();
}
/**
 * Create HTML URL link to navigate to best prsm and to navigae to proteoform page
 * @param {Float} e_value - Contains e value of the best prsm
 * @param {Float} precursor_mass - Contains precursor mass of the best prsm
 * @param {Int} prsm_id - Contains the best prsm id
 * @param {Int} proteoform_id - Contains the proteoform Id
 * @param {Int} PrSM_Count - Contians numbe rof prsms for a proteoform
 * @param {String} folderpath - Contains path to the data folder
 */
function Build_BestPrSM(e_value, precursor_mass, prsm_id, proteoform_id, PrSM_Count, folderpath) {
    let p = document.createElement("p");
    p.setAttribute("style", "font-size:16px;");
    let text1 = document.createElement("text");
    text1.innerHTML = "The ";
    p.appendChild(text1);
    let a_prsm = document.createElement("a");
    a_prsm.href = "prsm.html?folder=" + folderpath + "&prsm_id=" + prsm_id;
    a_prsm.innerHTML = " best PrSM ";
    p.appendChild(a_prsm);
    let text2 = document.createElement("text");
    text2.innerHTML = "has an E-value " + e_value + " and a precursor mass " + precursor_mass + ". There are ";
    p.appendChild(text2);
    let a_proteoform = document.createElement("a");
    a_proteoform.href = "proteoform.html?folder=" + folderpath + "&proteoform_id=" + proteoform_id;
    a_proteoform.innerHTML = PrSM_Count + " PrSMs ";
    p.appendChild(a_proteoform);
    let text3 = document.createElement("text");
    text3.innerHTML = "in total.";
    p.appendChild(text3);
    return p;
}
