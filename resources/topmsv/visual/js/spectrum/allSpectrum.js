"use strict";
/**
 * Starting point of building spectrums page.
 * Gets the data of all spectrum and shows on the html
 * @param {String} folderName - Provides the path to the data folder.
 * Provides path which helps to navigate to protein and prsm pages.
 */
/**
 * generate a search box for ms1 and ms2 spectra
 */
function searchBoxEventListner() {
    //find scan using ID in the page
    let scanNum = document.getElementById("ms2_scanid");
    if (scanNum) {
        let searchResult = document.getElementById("ms2_" + scanNum.value);
        if (searchResult) {
            let navbar = document.getElementById("nav-div");
            if (navbar) {
                let topOfElement = searchResult.offsetTop - navbar.clientHeight;
                window.scroll({ top: topOfElement });
            }
        }
    }
}
function addSearchBoxEvent() {
    let submitBtn = document.getElementById("ms2_submit");
    if (submitBtn) {
        submitBtn.addEventListener("click", searchBoxEventListner);
    }
}
/**
 * convert the json prsm data into HTML and create links for each spectra to navigate
 * @param {object} prsm - Contains data of a single prsm
 * @param {String} folderName - Provides path to build navigation links
 */
function prsmToHtml(prsm, folderName) {
    let protein = prsm.annotated_protein;
    let div = document.createElement('li');
    let id = "p" + prsm.prsm_id;
    div.setAttribute("id", id);
    let childDiv = document.createElement('div');
    childDiv.style.display = "flex";
    let ms1Spec = document.createElement('p');
    ms1Spec.innerHTML = "MS1 Scan #" + prsm.ms.ms_header.ms1_scans;
    ms1Spec.classList.add("scan-title");
    ms1Spec.id = "ms1_" + prsm.ms.ms_header.ms1_scans;
    let ms2Spec = document.createElement('p');
    ms2Spec.innerHTML = "MS2 Scan #" + prsm.ms.ms_header.scans;
    ms2Spec.classList.add("scan-title");
    ms2Spec.id = "ms2_" + prsm.ms.ms_header.scans;
    let protName = document.createElement('p');
    protName.innerHTML = protein.sequence_name;
    protName.classList.add("scan-title");
    childDiv.appendChild(ms2Spec);
    childDiv.appendChild(ms1Spec);
    childDiv.appendChild(protName);
    div.appendChild(childDiv);
    // let residues: {"position": string, "acid": string}[] = protein.annotation.residue;
    /*for (let i = protein.annotation.first_residue_position; i < protein.annotation.last_residue_position; i++) {
      sequence = sequence.concat(residues[i].acid);
    }*/
    let sequence = protein.annotation.annotated_seq;
    let p = document.createElement('p');
    let a = document.createElement('a');
    a.href = "prsm" + ".html" + "?folder=" + folderName + "&protein=" + prsm.prsm_id;
    //a.innerHTML = protein.sequence_name + " " + "first residue = " + protein.annotation.first_residue_position + " " + "last residue = " + protein.annotation.last_residue_position;
    a.innerHTML = sequence;
    p.appendChild(a);
    div.appendChild(p);
    return div;
}
function allSpectrum(folderName) {
    //@ts-ignore
    let l_prsms = prsm_data.prsms.prsm;
    let count = 1;
    // get the div container 
    let div = document.getElementsByClassName("container")[0];
    //switch between spectrum identification and protein identification
    let p = document.createElement('p'); //need to assign class 
    let x = location.href;
    let l_split = x.split(/[?#]+/)[0];
    let idx = l_split.lastIndexOf('\\');
    if (idx < 0) {
        idx = l_split.lastIndexOf('/');
    }
    let newAddress = l_split.slice(0, idx + 1) + "proteins.html?folder=" + folderName;
    p.innerHTML = '<a href=' + newAddress + '>Switch to Protein Identification</a>';
    //div.prepend(p);
    let titleDiv = document.getElementById('title-div');
    let h2 = document.getElementById('id-title');
    if (!h2) {
        console.error("ERROR: Identification title header doesn't exist");
        return;
    }
    if (!titleDiv) {
        console.error("ERROR: Title div doesn't exist");
        return;
    }
    // Check to see if protein variable inside l_proteins is an array.
    // Checks for multiple proteins
    if (Array.isArray(l_prsms)) {
        count = l_prsms.length;
        document.title = count + " spectra are identified";
        h2.innerHTML = count + " spectra are identified.";
    }
    else {
        document.title = count + " spectrum is identified";
        h2.innerHTML = count + " spectrum is identified.";
    }
    let br = document.createElement('br');
    // create header with protein count 
    div.appendChild(br);
    // Creating ordered list
    let ol = document.createElement('ol');
    // get the best prsm for each protein and form unique links for all the proteins
    // Check to see if protein variable inside l_proteins is an array.
    if (Array.isArray(l_prsms)) {
        l_prsms.forEach(function (prsm, index) {
            let div_temp = prsmToHtml(prsm, folderName);
            let br1 = document.createElement('br');
            div_temp.appendChild(br1);
            ol.appendChild(div_temp);
        });
    }
    else {
        let prsm = l_prsms;
        let div_temp = proteinToHtml(prsm, folderName);
        let br1 = document.createElement('br');
        div_temp.appendChild(p);
        div_temp.appendChild(br1);
        ol.appendChild(div_temp);
    }
    div.appendChild(ol);
    addSearchBoxEvent();
}
