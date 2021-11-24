/**
 * Function to Create Navigation buttons to navigate between two types of graph
 * createMs2NavElement from visual/prsm/load_spectra.js
 * @param {Array} divId - id of the div containing nav bar
 * @param {String} navId - id of the nav bar
 * These functions need to be rewritten as a library later
 */
 function clearMs2NavElement(navId: string){
    $("#"+navId).empty();
}
function addCheckboxTabInspect(navId: string): void {
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
function createMs2NavElementInspect(i: number, divId: string, navId: string, specScan: string): void {
    let ul: HTMLElement = <HTMLElement>document.getElementById(navId);
    let li: HTMLLIElement = document.createElement("li");
    let li_id: string = divId+"_graphlist_"+ i;
    li.setAttribute("id",li_id);
    if(i == 0) {
        li.setAttribute("class","nav-item ms2_graph_list active");
    }
    else {
        li.setAttribute("class","nav-item ms2_graph_list");
    }
    let a: HTMLAnchorElement = document.createElement("a");
    a.setAttribute("class","nav-link");
    a.setAttribute("href","#!");
    a.innerHTML = "Scan "+ specScan;
    li.appendChild(a);
    ul.appendChild(li);

    li = document.createElement("li");
    li_id = divId+"_monographlist_"+ i;
    li.setAttribute("id",li_id);
    li.setAttribute("class","nav-item ms2_graph_list");
    a = document.createElement("a");
    a.setAttribute("class","nav-link");
    a.setAttribute("href","#!");
    a.innerHTML = "Mass scan "+ specScan ;
    li.appendChild(a);
    ul.appendChild(li);
}
function addEventNavBar(monoMassGraphObj: SpectrumView): void{
    // add action for nav bar
    $(".ms2_graph_list").click((e: JQuery.ClickEvent) => {
        let Id: string = e.currentTarget.id;
        let ms2GraphList: HTMLCollectionOf<Element> = document.getElementsByClassName("ms2_graph_list");
        for (let i = 0; i < ms2GraphList.length; i++){
            if (ms2GraphList[i].id == Id){
                ms2GraphList[i].classList.add("active");
                let svgId: string  = ms2GraphList[i].id;
                let svgIdSplit: string[] = svgId.split("_");
                let type: string = svgIdSplit[3];
                let spectrumTab = <HTMLElement>document.getElementById(Constants.SPECTRUMGRAPHID);
                let monoMassTab = <HTMLElement>document.getElementById(Constants.MONOMASSGRAPHID);
                let checkboxTab: HTMLElement | null = document.getElementById("checkbox-tab");

                if (type == "graphlist"){
                    spectrumTab.style.display="";
                    monoMassTab.style.display="none";
                    if (checkboxTab != null) {
                        checkboxTab.style.display = "none";
                    }
                }else{
                    spectrumTab.style.display="none";
                    monoMassTab.style.display="";
                    if (checkboxTab != null) {
                        checkboxTab.style.display = "";
                    }
                }
            }else{
                ms2GraphList[i].classList.remove("active");
            }
        }
    })
    //add an event listner for checkbox
    $("#checkbox-anno-line").on("change", function() {
        monoMassGraphObj.redraw();
    })
}
function switchTab(graphType: string): void{
    let ms2GraphList: HTMLCollectionOf<Element> = document.getElementsByClassName("ms2_graph_list");
    
    for (let i = 0; i < ms2GraphList.length; i++){
        let svgId: string = ms2GraphList[i].id;
        let svgIdSplit: string[] = svgId.split("_");
        let type: string = svgIdSplit[3];
        if (type == graphType){
            ms2GraphList[i].classList.add("active");
        }
        else{
            ms2GraphList[i].classList.remove("active");
        }
    }
    if (graphType == "monographlist"){
        $("#"+Constants.SPECTRUMGRAPHID).hide();
        $("#"+Constants.MONOMASSGRAPHID).show();    
    }
    else{
        $("#"+Constants.SPECTRUMGRAPHID).show();
        $("#"+Constants.MONOMASSGRAPHID).hide();    
    }
}