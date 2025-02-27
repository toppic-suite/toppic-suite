"use strict";
/**
 * Create naviagtion urls to go back to protein page and all protein page
 * @param {String} folderpath - Path to data folder
 */
function proteoformUrl(folderpath, prsm) {
    // prsm_data is a global variable containing complete data of the proteoform from the data file
    let l_proteoform_Url = `${prsm[0].getProteoform().getProtName()} ${prsm[0].getProteoform().getProtDesc()}`;
    // Set title for the page
    document.title = `Proteoform #${prsm[0].getProteoform().getId()} from ${l_proteoform_Url}`;
    let l_allproteins_url = "proteins.html?folder=" + folderpath;
    // Set the naviagtion urls to all proteins page and protein page
    let allProteinURLStart = document.getElementById("allprotein_url_start");
    let allProteinURLEnd = document.getElementById("allprotein_url_end");
    let proteinURLStart = document.getElementById("protein_url_start");
    let proteinURLEnd = document.getElementById("protein_url_end");
    let protHeader = document.getElementById("proteoform_header");
    let prsmCnt = document.getElementById("prsm_count");
    if (allProteinURLStart) {
        allProteinURLStart.href = l_allproteins_url;
    }
    if (allProteinURLEnd) {
        allProteinURLEnd.href = l_allproteins_url;
    }
    if (proteinURLStart) {
        proteinURLStart.innerHTML = l_proteoform_Url;
        proteinURLStart.href = "protein.html?folder=" + folderpath + "&protein_Id=" + prsm[0].getProteoform().getSeqId();
    }
    if (proteinURLEnd) {
        proteinURLEnd.innerHTML = l_proteoform_Url;
        proteinURLEnd.href = "protein.html?folder=" + folderpath + "&protein_Id=" + prsm[0].getProteoform().getSeqId();
    }
    if (protHeader) {
        protHeader.innerHTML = "Proteoform #" + prsm[0].getProteoform().getId();
    }
    // Get the count of number of proteoform 
    if (prsm.length > 1) {
        if (prsmCnt) {
            prsmCnt.innerHTML = `${prsm.length} PrSMs are identified for the proteoform`;
        }
    }
    else {
        if (prsmCnt) {
            prsmCnt.innerHTML = "1 PrSM is identified for the proteoform";
        }
    }
}
/**
 * Create a table with prsm data and URL to navigate to appropriate prsm
 * @param {String} folderpath - Provides path to the data folder
 */
function createTableData(folderpath, prsm) {
    let table = document.getElementById('proteoform_data');
    if (!table) {
        console.error("Error: Can't find proteoform table element");
        return;
    }
    let tbdy = document.createElement('tbody');
    let count = 0;
    let sequence_name = prsm[0].getProteoform().getProtName();
    // Iterate through the number of prsms
    if (prsm.length > 1) {
        prsm.forEach(function (onePrsm, i) {
            let spectra = onePrsm.getMs2Spectra();
            if (!spectra) {
                console.error("error: total peaks information is not correctly added");
                return;
            }
            let All_Peak_count = spectra[0].getPeaks().length;
            let tr = document.createElement('tr');
            for (let i = 0; i < 7; i++) {
                let td = document.createElement('td');
                td.setAttribute("align", "center");
                if (i === 0) {
                    td.innerHTML = spectra[0].getScanNum();
                    td.setAttribute("width", "5%");
                }
                if (i === 1) {
                    td.innerHTML = sequence_name;
                    td.setAttribute("width", "25%");
                }
                if (i === 2) {
                    td.innerHTML = onePrsm.getEValue().toString();
                    td.setAttribute("width", "15%");
                }
                if (i === 3) {
                    td.innerHTML = All_Peak_count.toString();
                    td.setAttribute("width", "12%");
                }
                if (i === 4) {
                    td.innerHTML = onePrsm.getMatchedPeakCount().toString();
                    td.setAttribute("width", "16%");
                }
                if (i === 5) {
                    let ionCnt = onePrsm.getFragIonCount();
                    if (ionCnt) {
                        td.innerHTML = ionCnt.toString();
                        td.setAttribute("width", "22%");
                    }
                }
                if (i === 6) {
                    // Create URL to navigate
                    let a = document.createElement('a');
                    let l_href = "prsm.html?folder=" + folderpath + "&prsm_id=" + onePrsm.getId();
                    let l_link = "link" + (count + 1);
                    a.setAttribute("href", l_href);
                    a.setAttribute("id", l_link);
                    //a.innerHTML = "See PrSM&gt;&gt;"; // Adding >> marks
                    a.innerHTML = "Click"; // Adding >> marks
                    td.appendChild(a);
                    td.setAttribute("width", "5%");
                    count++;
                }
                tr.appendChild(td);
            }
            tbdy.appendChild(tr);
        });
    }
    else {
        let spectra = prsm[0].getMs2Spectra();
        if (!spectra) {
            console.error("error: total peaks information is not correctly added");
            return;
        }
        let All_Peak_count = spectra[0].getPeaks().length;
        let tr = document.createElement('tr');
        for (let i = 0; i < 7; i++) {
            let td = document.createElement('td');
            td.setAttribute("align", "center");
            if (i === 0) {
                td.innerHTML = spectra[0].getScanNum();
            }
            if (i === 1) {
                td.innerHTML = sequence_name;
            }
            if (i === 2) {
                td.innerHTML = prsm[0].getEValue().toString();
            }
            if (i === 3) {
                td.innerHTML = All_Peak_count.toString();
            }
            if (i === 4) {
                td.innerHTML = prsm[0].getMatchedPeakCount().toString();
            }
            if (i === 5) {
                let ionCnt = prsm[0].getFragIonCount();
                if (ionCnt) {
                    td.innerHTML = ionCnt.toString();
                }
            }
            if (i === 6) {
                // Create link to navigate
                let a = document.createElement('a');
                let l_href = "prsm.html?prsm_id=" + prsm[0].getId();
                let l_link = "link" + (count + 1);
                a.setAttribute("href", l_href);
                a.setAttribute("id", l_link);
                a.innerHTML = "See PrSM&gt;&gt;";
                td.appendChild(a);
                count++;
            }
            tr.appendChild(td);
        }
        tbdy.appendChild(tr);
    }
    table.appendChild(tbdy);
}
