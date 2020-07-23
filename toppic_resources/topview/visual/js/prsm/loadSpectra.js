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
    let ions = {};
    spGraph = new SpectrumGraph(ms1SvgId,peaks,envelopes,ions);
    let prec_mono_mz = prsm_data.prsm.ms.ms_header.precursor_mz;
    spGraph.para.updateMzRange(prec_mono_mz);
    spGraph.para.setHighlight(prec_mono_mz);
    spGraph.redraw();
    return spGraph;
  }
}

function loadMsTwo(specIdList, fileList, divId, navId){
  let len = fileList.length;
  let specList = [];
  let graphList = [];
  for (let i = 0; i < len; i++) {
    let filename = fileList[i];
    let script= document.createElement('script');
    script.src = filename;
    document.head.appendChild(script);
    script.onload = function () {
      specList.push(ms2_data);
      // As data loading is an asynchronous process, 
      // we need to wait till all the data is loaded to execute the below functions
      if (i == len - 1) {
        document.getElementById("dataLoading").remove();
        for (let i = 0; i < specList.length; i++) {
          console.log(specList[i]);
          createMs2NavElement(i, divId, navId, specIdList[i], specList[i].scan);
          let svgId = divId + "_graph_" + specIdList[i];
          createSvg(i, divId, svgId, "ms2_svg_graph_class");
          let peaks = specList[i].peaks;
          let envelopes = specList[i].envelopes;
          let ions = {};
          let spGraph = new SpectrumGraph(svgId,peaks,envelopes,ions);
          spGraph.redraw();
          graphList.push(spGraph);
        }
      }
      // Set on click actions once tabs to naviage between spectrums are created
      //graphOnClickActions();
    }
  }
  return graphList;
}

/**
 * Function to Create Navigation buttons to navigate between spectrums
 * @param {Array} scanidList - Contains scan Id List
 * @param {String} id - Contains Id of the avg on which spectrum to be drawn
 */
function createMs2NavElement(i, divId, navId, specId, specScan){
  let _ul = document.getElementById(navId);
  let li = document.createElement("li");
  li.setAttribute("class","nav-item");
  let li_id = divId+"_list_"+ specId;
  li.setAttribute("id",li_id);
  let a = document.createElement("a");
  a.setAttribute("class","nav-link ms2_scanIds");
  if(i == 0) {
    a.setAttribute("class","nav-link ms2_scanIds active");
  }
  a.setAttribute("href","#!");
  a.setAttribute("value", divId+"_graph_" + specId);
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
  if(i != 0)
  {
    svg.style.display = "none"; 
  }
  div.appendChild(svg);
}

