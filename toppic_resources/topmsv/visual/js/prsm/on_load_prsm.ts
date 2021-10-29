"use strict";

let ms2ScanList: string[] = [];
let ms2GraphList: SpectrumView[] = [];
/**
 * This function waits till all the HTML tags are loaded.
 * Invokes functions to loads the complete prsm page and visualization
 */
$(document).ready(function () {
    let x = location.href;
    let l_split = x.split(/[?#]+/)[1];
    let path_and_value = l_split.split("&");
    let folder_path = path_and_value[0].split("=")[1];
    // get the prsm Id number by splitting url with "?" and "=" 
    let prsm_seq_num = path_and_value[1].split("=")[1];
    let prsm_data_script = document.createElement('script');
    // get the prsm Id number by splitting url with "?","=","#"
    let prsm_data_file_name = "../data/" + folder_path + "/prsms/prsm" + prsm_seq_num + ".js";
    prsm_data_script.type = 'text/javascript';
    prsm_data_script.src = prsm_data_file_name;
    // Add data file to the script tag in the html
    let head = document.getElementsByTagName('head')[0];
    head.appendChild(prsm_data_script);
    // Wait till the data is losded before calling any functions
    prsm_data_script.onload = function () {
        let parsePrsm: ParsePrsm = new ParsePrsm(true, "../../topfd/ms1_json/", true, "../../topfd/ms2_json/");

        parsePrsm.geneDataObj(function(prsmObj: Prsm): void {
            //Build Urls to naviga back to proteoform page, proteins page and all protein page
            BuildUrl(folder_path, prsmObj);
            // Get the information of the PRSM to the HTML
            loadDatafromJson2Html(prsmObj);
            // Get occurence ptms in prsmtohtml.js
            occurence_ptm(prsmObj);
            // Get Unknown Ptms to show in the html in prsmtohtml.js
            getUnknownPtms(prsmObj);
            // Calling function with actions on click of buttons
            addButtonActions();
            // Get all the scanIds
            let ms2Spectrum: Spectrum[] | null = prsmObj.getMs2Spectra();
            if (!ms2Spectrum) {
                console.error("ERROR: ms2 spectrum is empty");
                return;
            }

            //let scanIds: string[] = ms2Spectrum[0].getScanNum().split(" ");
            // Get all the SpecIds
            //let specIds: string[] = ms2Spectrum[0].getSpectrumId().split(" ");
            let scanIds: string[] = [];
            let specIds: string[] = [];
            ms2Spectrum.forEach((spectra) => {
                scanIds.push(spectra.getScanNum());
                specIds.push(spectra.getSpectrumId());
            });
            // Add Buttong with dropdowns with Scan numbers to navigae to inspect page
            setDropDownItemsForInspectButton(scanIds, specIds);
            // Add all the data and set local storage variables
            onClickToInspect(prsmObj);
            // Using spectrum graph library
            // Get Ms1 Id to draw MS1 Spectrum
            loadMsOne(prsmObj.getMs1Spectra(), "ms1_svg");
            // Get Ms2 ids to draw MS2 Spectrum
            loadMsTwo(prsmObj, ms2GraphList, "ms2_svg_div", "ms2_graph_nav");
            
            let prsmView = new PrsmView("prsm_svg", prsmObj);
            prsmView.redraw();
            // add prsm graph to popup
            let savePrsmObj: SavePrsm = new SavePrsm(prsmView);
            savePrsmObj.main();  
            // Create peaks data into table content
            let dataTable = new DataTable(prsmObj, true, ms2GraphList);
            dataTable.setSpecSvgId("ms2_svg_div_graph_");
            dataTable.setMonoMassSvgId("ms2_svg_div_mono_graph_");
            dataTable.drawTable();
        })        
    };
});
