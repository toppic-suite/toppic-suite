/**
 * @function loadMsOne
 * @description - This function load an MS1 spectrum. 
 */
function loadMsOne(filename, ms1SvgId){
  let script= document.createElement('script');
  script.src = filename;
  document.head.appendChild(script);
  script.onload = function(){
    let peaks = ms1_data.peaks;
    let envelopes = ms1_data.envelopes;
    let ions = [];
    spGraph = new SpectrumGraph(ms1SvgId,peaks,envelopes,ions);
    let precMonoMz = prsm_data.prsm.ms.ms_header.precursor_mz;
    spGraph.para.updateMzRange(precMonoMz);
    spGraph.para.setHighlight(precMonoMz);
    spGraph.redraw();
    return spGraph;
  }
}

function loadMsTwo(specIdList, fileList, divId, navId){
  let len = fileList.length;
  let cnt = 0;
  let specList = [];
  let graphList = [];
  for (let i = 0; i < len; i++) {
    let filename = fileList[i];
    let script= document.createElement('script');
    script.src = filename;
    document.head.appendChild(script);
    script.onload = function () {
      specList.push(ms2_data);
      cnt = cnt+1;
      // As data loading is an asynchronous process, 
      // we need to wait till all the data is loaded to execute the below functions
      if (cnt == len) {
        document.getElementById("dataLoading").remove();
        specList.sort(function(x,y){
            return d3.ascending(x.scan, y.scan);
        })
        for (let j = 0; j < specList.length; j++) {
          //console.log(specList[i]);
          let svgId = divId + "_graph_" + j;
          createMs2NavElement(j, divId, navId, specList[j].scan);
          createSvg(j, divId, svgId, "ms2_svg_graph_class");
          let peaks = specList[j].peaks;
          let envelopes = specList[j].envelopes;
          let ions = [];
          let spGraph = new SpectrumGraph(svgId,peaks,envelopes,ions);
          spGraph.redraw();
          graphList.push(spGraph);
        }
        // add action for nav bar
        $(".ms2_graph_list").click(function(){
          let ms2Id = this.id;
          //console.log("ms2id", ms2Id);
          let ms2Split = ms2Id.split("_");
          let ms2Index = parseInt(ms2Split[ms2Split.length-1]);
          for (let i = 0; i < ms2GraphList.length; i++) {
            let listId = "ms2_svg_div_list_" + i;
            let graphId = "ms2_svg_div_graph_" + i;
            //console.log(listId, graphId);
            let listElement = document.getElementById(listId);
            let graphElement = document.getElementById(graphId);
            if (i== ms2Index) {
              listElement.classList.add("active");
              graphElement.style.display="inline";
            }
            else {
              listElement.classList.remove("active");
              graphElement.style.display="none";
            }
          }
        })
      }
    }
  }
  return [specList, graphList];
}

/**
 * Function to Create Navigation buttons to navigate between spectrums
 * @param {Array} scanidList - Contains scan Id List
 * @param {String} id - Contains Id of the avg on which spectrum to be drawn
 */
function createMs2NavElement(i, divId, navId, specScan){
  let _ul = document.getElementById(navId);
  let li = document.createElement("li");
  let li_id = divId+"_list_"+ i;
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
  _ul.appendChild(li);
}

/**
 * This generates spectrum for each spec Id
 * @param {String} divId - Contains Id of the div tag under which the monomass graphs are drawn
 * @param {String} svgId - Contains id as "monoMassSvg_" to which scan Id is added
 * @param {String} className - Contains class name to which the corresponding svg graphs are drawn 
 */
function createSvg(i, divId, svgId, className){
  let div = document.getElementById(divId); 
  let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
  svg.setAttribute("id",svgId);
  svg.setAttribute("class",className);
  svg.style.backgroundColor = "#F8F8F8"; 
  if(i != 0) {
    svg.style.display = "none"; 
  }
  else {
    svg.style.display = "inline"; 
  }
  div.appendChild(svg);
}

