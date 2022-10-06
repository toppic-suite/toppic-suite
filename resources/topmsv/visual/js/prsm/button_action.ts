"use strict";
/**
 * Button Actions
 */

function addButtonActions(): void {
    //	MS1 graph popup window 
    $("#Precursor_mz").click(function () {
        // @ts-ignore
        $("#ms1_graph_popup_window").draggable({
            appendTo: "body"
        });
    });
    d3.select("#download_ms1_png_btn").on("click", function () {
        let x: number = d3.event.pageX;
        let y: number = d3.event.pageY;
        popupNameWindow("png", "ms1_svg", x, y);
    });
    d3.select("#download_ms1_svg_btn").on("click", function () {
        let x: number = d3.event.pageX;
        let y: number = d3.event.pageY;
        popupNameWindow("svg", "ms1_svg", x, y);
    });
    // Show MS2 graph button 
    $("#ms2_graph_show_btn").click(function () {
        if ($.trim($(this).text()) === 'Show Spectrum') {
            showMs2Graph();
        }
        else {
            hideMs2Graph();
        }
    });
    // MS2 graph help button 
    $("#ms2_graph_help_btn").click(function () {
        // @ts-ignore
        $("#ms2_graph_help_popup_window").draggable({
            appendTo: "body"
        });
    });
    // On click of mono mass mz, zoom all the graph to the corresponding point
    /*$(".row_mono_mz").click(function () {
        console.log("clicked");
        let parentId: string = $(this).parent().parent().prop('id');
        let scanNumDiv = <HTMLElement>(<HTMLElement>document.getElementById(parentId)).firstChild;
        let scanNum: string = scanNumDiv.innerHTML;
        console.log(parentId, scanNum);
        if (scanNum === null) {
            console.error("ERROR: scan number is null");
        }

        /*	get Mono M/z value till 3 decimal values	*/
    /*let monoMz: number = parseFloat(parseFloat(this.innerHTML).toFixed(3));
    for (let i = 0; i < ms2ScanList.length; i++) {
        let listId: string = "ms2_svg_div_graphlist_" + i;
        let graphId: string = "ms2_svg_div_graph_" + i;
        let monolistId: string = "ms2_svg_div_monographlist_" + i;
        let monoGraphId: string = "ms2_svg_div_mono_graph_" + i;
        let listElement: HTMLElement | null = document.getElementById(listId);
        let graphElement: HTMLElement | null = document.getElementById(graphId);
        let monoListElement: HTMLElement | null = document.getElementById(monolistId);
        let monoGraphElement: HTMLElement | null = document.getElementById(monoGraphId);
        if (scanNum == ms2ScanList[i]) {
            if (!listElement) {
                console.log("ERROR: listElement is null");
                return;
            }
            if (!graphElement) {
                console.log("ERROR: graphElement is null");
                return;
            }
            if (!monoListElement) {
                console.log("ERROR: monoListElement is null");
                return;
            }
            if (!monoGraphElement) {
                console.log("ERROR: monoGraphElement is null");
                return;
            }
            listElement.classList.add("active");
            graphElement.style.display = "";
            monoListElement.classList.remove("active");
            monoGraphElement.style.display = "none";
            
            let spGraph: SpectrumView = ms2GraphList[i];
            // set monoMz to do
            spGraph.getPara().updateMzRange(monoMz);
            spGraph.redraw();
        }
        else {
            if (!listElement) {
                console.log("ERROR: listElement is null");
                return;
            }
            if (!graphElement) {
                console.log("ERROR: graphElement is null");
                return;
            }
            if (!monoListElement) {
                console.log("ERROR: monoListElement is null");
                return;
            }
            if (!monoGraphElement) {
                console.log("ERROR: monoGraphElement is null");
                return;
            }
            listElement.classList.remove("active");
            graphElement.style.display = "none";
            monoListElement.classList.remove("active");
            monoGraphElement.style.display = "none";
        }
    }
    showMs2Graph();
});*/
}
/**
 * Show spectrum graphs on click of Show Spectrum
 */
function showMs2Graph(): void {
    $("#ms2_graph_show_btn").text('Hide Spectrum');
    let helpBtn: HTMLElement | null = document.getElementById("ms2_graph_help_btn");
    let saveBtn: HTMLElement | null = document.getElementById("ms2_graph_save_btn");
    let svgDiv: HTMLElement | null = document.getElementById("ms2_svg_div");

    if (!helpBtn || !saveBtn || !svgDiv) {
        console.error("ERROR: invalid button ID or SVG div ID");
        return;
    }
    helpBtn.style.display = "block";
    saveBtn.style.display = "block";
    svgDiv.style.display = "block";
}
/**
 * Hide spectrum graphs on click of hide spectrum
 */
function hideMs2Graph(): void {
    $("#ms2_graph_show_btn").text('Show Spectrum');
    let helpBtn: HTMLElement | null = document.getElementById("ms2_graph_help_btn");
    let saveBtn: HTMLElement | null = document.getElementById("ms2_graph_save_btn");
    let svgDiv: HTMLElement | null = document.getElementById("ms2_svg_div");

    if (!helpBtn || !saveBtn || !svgDiv) {
        console.error("ERROR: invalid button ID or SVG div ID");
        return;
    }
    helpBtn.style.display = "none";
    saveBtn.style.display = "none";
    svgDiv.style.display = "none";
}
