/**
 * Function to Create Navigation buttons to navigate between two types of graph
 * createMs2NavElement from visual/prsm/load_spectra.js
 * @param {Array} divId - id of the div containing nav bar
 * @param {String} navId - id of the nav bar
 * These functions need to be rewritten as a library later
 */
 function clearMs2NavElement(navId){
    $("#"+navId).empty();
}
function createMs2NavElement(i, divId, navId, specScan){
    let ul = document.getElementById(navId);
    let li = document.createElement("li");
    let li_id = divId+"_graphlist_"+ i;
    li.setAttribute("id",li_id);
    if(i == 0) {
        li.setAttribute("class","nav-item ms2_graph_list active");
    }
    else {
        li.setAttribute("class","nav-item ms2_graph_list");
    }
    let a = document.createElement("a");
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
function addEventNavBar(){
    // add action for nav bar
    $(".ms2_graph_list").click(function(){
        let Id = this.id;
        let ms2GraphList = document.getElementsByClassName("ms2_graph_list");
        for (let i = 0; i < ms2GraphList.length; i++){
            if (ms2GraphList[i].id == Id){
                ms2GraphList[i].classList.add("active");
                let svgId  = ms2GraphList[i].id;
                let svgIdSplit = svgId.split("_");
                let type = svgIdSplit[3];
                if (type == "graphlist"){
                    document.getElementById("spectrum").style.display="";
                    document.getElementById("monoMassGraph").style.display="none";
                }else{
                    document.getElementById("spectrum").style.display="none";
                    document.getElementById("monoMassGraph").style.display="";
                }
            }else{
                ms2GraphList[i].classList.remove("active");
            }
        }
    })
}
/*function addEventNavBar(){
    // add action for nav bar
    $(".ms2_graph_list").click(function(){
        let Id = this.id;
        let ms2GraphList = document.getElementsByClassName("ms2_graph_list");
        for (let i = 0; i < ms2GraphList.length; i++){
            if (ms2GraphList[i].id == Id){
                ms2GraphList[i].classList.add("active");
                let svgId  = ms2GraphList[i].id;
                let svgIdSplit = svgId.split("_");
                let type = svgIdSplit[3];
                //svgId = svgId.slice(svgId.lastIndexOf("_") + 1);
                
                document.getElementById(svgId).style.display="";
            }else{
                ms2GraphList[i].classList.remove("active");
                let svgId  = ms2GraphList[i].id;
                
                svgId = svgId.slice(svgId.lastIndexOf("_") + 1);
                console.log(ms2GraphList[i])
                document.getElementById(svgId).style.display="none";
            }
        }
    })
}*/
