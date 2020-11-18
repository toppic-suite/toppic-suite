/**
 * Function to Create Navigation buttons to navigate between two types of graph
 * Modifed version of createMs2NavElement from visual/prsm/load_spectra.js
 * @param {Array} divId - id of the div containing nav bar
 * @param {String} navId - id of the nav bar
 */
function createMs2NavElement(divId, navId){
    let ul = document.getElementById(navId);
    let li = document.createElement("li");
    let li_id = divId+"_spectrum";
    li.setAttribute("id",li_id);
    li.setAttribute("class","nav-item ms2_graph_list active");
    
    let a = document.createElement("a");
    a.setAttribute("class","nav-link");
    a.setAttribute("href","#!");
    a.innerHTML = "Scan ";
    li.appendChild(a);
    ul.appendChild(li);
  
    li = document.createElement("li");
    li_id = divId+"_monoMassGraph";
    li.setAttribute("id",li_id);
    li.setAttribute("class","nav-item ms2_graph_list");
    a = document.createElement("a");
    a.setAttribute("class","nav-link");
    a.setAttribute("href","#!");
    a.innerHTML = "Mass scan ";
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
                svgId = svgId.slice(svgId.lastIndexOf("_") + 1);
                document.getElementById(svgId).style.display="";
            }else{
                ms2GraphList[i].classList.remove("active");
                let svgId  = ms2GraphList[i].id;
                svgId = svgId.slice(svgId.lastIndexOf("_") + 1);
                document.getElementById(svgId).style.display="none";
            }
        }
    })
}
