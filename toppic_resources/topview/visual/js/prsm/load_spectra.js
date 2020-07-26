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
          //console.log(specList[i]);
          let svgId = divId + "_graph_" + i;
          createMs2NavElement(i, divId, navId, svgId, specList[i].scan);
          createSvg(i, divId, svgId, "ms2_svg_graph_class");
          let peaks = specList[i].peaks;
          let envelopes = specList[i].envelopes;
          let ions = [];
          let spGraph = new SpectrumGraph(svgId,peaks,envelopes,ions);
          spGraph.redraw();
          graphList.push(spGraph);
        }
      }
      // Set on click actions once tabs to naviage between spectrums are created
      //addGraphOnClickActions();
    }
  }
  return [specList, graphList];
}

/**
 * Function to Create Navigation buttons to navigate between spectrums
 * @param {Array} scanidList - Contains scan Id List
 * @param {String} id - Contains Id of the avg on which spectrum to be drawn
 */
function createMs2NavElement(i, divId, navId, svgId, specScan){
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
  a.setAttribute("value", svgId);
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


/**
 * Generating Jquery Onclick actions
 */
/*
function addMsTwoGraphOnClickActions(){
  // On click of mono mass mz, zoom all the graph to the corresponding point
  $(".fragMonoMz").click(function() {
    let parent_id  = $(this).parent().parent().prop('id');
    let CurrentScanVal = document.getElementById(parent_id).firstChild.innerHTML;
    //	get Mono M/z value till 3 decimal values	
    let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
    let [currentData,specId] = getCurrentData(ms2_ScansWithData,CurrentScanVal);
    let id = "ms2svg_"+CurrentScanVal;
    showCorrespondingGraph(id,".ms2_svg_graph_class");
    generateCorrespondingGraph(currentData,id,peak_value,specId);
    id = "monoMassSvg_"+CurrentScanVal;
    showCorrespondingGraph(id,".monoMass_svg_graph_class");
    [currentData,specId] = getCurrentData(monoMassDataList,CurrentScanVal);
    let CurrentMonoMassVal = $("#"+parent_id+" .row_monoMass").html();
    generateCorrespondingGraph(currentData,id,parseFloat(CurrentMonoMassVal),specId);
    activateCurrentnavbar("ms2_graph_nav",CurrentScanVal);
    activateCurrentnavbar("monoMass_nav",CurrentScanVal);

    showSpectrun();
  });
  // ms2_scanIds is the Id of the nav tabs for multiple navs.
  // On Click shows corresponding graph by hiding others.
  $(".ms2_scanIds").click(function(){
    let value = this.getAttribute('value');
    let [currentData,specId] = getCurrentData(ms2_ScansWithData,value);
    id = "ms2svg_"+value;
    // Hide all the graphs except the one clicked
    showCorrespondingGraph(id,".ms2_svg_graph_class");
    $("#ms2_graph_nav .active").removeClass("active");
    $(this).addClass("active");
  })
    
  // Hide all othe graphs of monomass spectrum other than the one clicked
  $(".monoMass_scanIds").click(function(){
    let value = this.getAttribute('value');
    let [currentData,specId] = getCurrentData(monoMassDataList,value);
    id = "monoMassSvg_"+value;
    // Hide all the graphs except the one clicked
    showCorrespondingGraph(id,".monoMass_svg_graph_class");
    $("#monoMass_nav .active").removeClass("active");
    $(this).addClass("active");
  })
}
*/
