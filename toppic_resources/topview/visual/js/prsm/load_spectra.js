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
    spGraph = new SpectrumGraph(ms1SvgId,peaks);
    spGraph.addRawSpectrumAnno(envelopes, ions);
    let precMonoMz = prsm_data.prsm.ms.ms_header.precursor_mz;
    spGraph.para.updateMzRange(precMonoMz);
    spGraph.para.setHighlight(precMonoMz);
    spGraph.redraw();
    return spGraph;
  }
}

function loadMsTwo(specIdList, fileList, proteoform, divId, navId){
  let len = fileList.length;
  let cnt = 0;
  let specList = [];
  let graphList = [];
  let monoGraphList = [];
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
          createMs2NavElement(j, divId, navId, specList[j].scan);
          let show = false;
          if (j == 0) {
            show = true;
          }
          let svgId = divId + "_graph_" + j;
          createSvg(show, divId, svgId, "ms2_svg_graph_class");
          let specId = specList[j].id;
          let peaks = specList[j].peaks;
          let envelopes = specList[j].envelopes;
          let deconvPeaks = prsm_data.prsm.ms.peaks.peak;
          let [ions, monoIons] = getIons(specId, deconvPeaks, envelopes);
          specList[j].ions = ions;
          let spGraph = new SpectrumGraph(svgId,peaks);
          spGraph.addRawSpectrumAnno(envelopes,ions); 
          spGraph.redraw();
          graphList.push(spGraph);
          //mono mass svg
          let monoSvgId = divId + "_mono_graph_" +j;
          show = false;
          createSvg(show, divId, monoSvgId, "ms2_svg_graph_class");
          let monoMasses = getMonoMasses(deconvPeaks);
          //console.log(deconvPeaks, monoMasses);
          specList[j].monoMasses = monoMasses;
          specList[j].monoIons = monoIons;
          let nIonType = specList[j].n_ion_type;
          let cIonType = specList[j].c_ion_type;
          let monoSpGraph = new SpectrumGraph(monoSvgId,monoMasses); 
          monoSpGraph.addMonoMassSpectrumAnno(monoIons,proteoform, nIonType, cIonType);
          monoSpGraph.para.setMonoMassGraph(true);
          monoSpGraph.redraw();
          monoGraphList.push(monoSpGraph);
        }
        // add action for nav bar
        $(".ms2_graph_list").click(function(){
          let ms2Id = this.id;
          //console.log("ms2id", ms2Id);
          let ms2Split = ms2Id.split("_");
          let ms2Index = parseInt(ms2Split[ms2Split.length-1]);
          let type = ms2Split[ms2Split.length - 2];
          for (let i = 0; i < ms2GraphList.length; i++) {
            let listId = "ms2_svg_div_graphlist_" + i;
            let monoListId = "ms2_svg_div_monographlist_" + i;
            let graphId = "ms2_svg_div_graph_" + i;
            let monoGraphId = "ms2_svg_div_mono_graph_" + i;
            //console.log(listId, graphId);
            let listElement = document.getElementById(listId);
            let monoListElement = document.getElementById(monoListId);
            let graphElement = document.getElementById(graphId);
            let monoGraphElement = document.getElementById(monoGraphId);
            if (i== ms2Index) {
              if (type == "graphlist") {
                listElement.classList.add("active");
                monoListElement.classList.remove("active");
                graphElement.style.display="";
                monoGraphElement.style.display="none";
              }
              else {
                listElement.classList.remove("active");
                monoListElement.classList.add("active");
                graphElement.style.display="none";
                monoGraphElement.style.display="";
              }
            }
            else {
              listElement.classList.remove("active");
              monoListElement.classList.remove("active");
              graphElement.style.display="none";
              monoGraphElement.style.display="none";
            }
          }
        })
      }
    }
  }
  //if the below lines are outside this scope, they execute before
  //graphList and monoGraphList are returned with valid values
  let saveSpectrumObj = new SaveSpectrum(graphList, monoGraphList);
  saveSpectrumObj.main();

  return [specList, graphList, monoGraphList];
}

/**
 * Function to Create Navigation buttons to navigate between spectrums
 * @param {Array} scanidList - Contains scan Id List
 * @param {String} id - Contains Id of the avg on which spectrum to be drawn
 */
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
  a.innerHTML = "Scan "+ specScan + " masses";
  li.appendChild(a);
  ul.appendChild(li);
}

/**
 * This generates spectrum for each spec Id
 * @param {String} divId - Contains Id of the div tag under which the monomass graphs are drawn
 * @param {String} svgId - Contains id as "monoMassSvg_" to which scan Id is added
 * @param {String} className - Contains class name to which the corresponding svg graphs are drawn 
 */
function createSvg(show, divId, svgId, className){
  let div = document.getElementById(divId); 
  let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
  svg.setAttribute("id",svgId);
  svg.setAttribute("class",className);
  svg.style.backgroundColor = "#F8F8F8"; 
  if(show) {
    svg.style.display = ""; 
  }
  else {
    svg.style.display = "none"; 
  }
  div.appendChild(svg);
}

function getMonoMasses(deconvPeaks) {
  let masses = [];
  for (let i = 0; i < deconvPeaks.length; i++) {
    let mass = {};
    mass.mz = deconvPeaks[i].monoisotopic_mass;
    mass.intensity = deconvPeaks[i].intensity;
    masses.push(mass);
  }
  return masses;
}

/**
 * @function getIonData
 * @description gets ion data list with mz, intensity and ion name 
 * This function gets matched ion data
 * @param {object} prsm_data - contains complete data of prsm 
 * @param {int} specId - contains information of the spec Id
 * @param {object} json_data - contains complete data of spectrum
 */
function getIons(specId, deconvPeaks, envelopes){
  envelopes.sort(function(x,y) {
    return d3.ascending(x.id, y.id);
  })
  //console.log(specId, deconvPeaks);
  let ions = [];
  let monoIons = [];
  //console.log(deconvPeaks);
  deconvPeaks.forEach(function(element) {
    if(element.hasOwnProperty('matched_ions_num') && 
      element.spec_id == specId) {
      let ionText = "";
      let massError = 0;
      if (parseInt(element.matched_ions_num) == 1) {
        let matchedIon = element.matched_ions.matched_ion;
        let ionType = matchedIon.ion_type;
        if (ionType == "Z_DOT") {
          ionType = "Z\u02D9";
        }
        ionText = ionType + matchedIon.ion_display_position;
        massError = parseFloat(matchedIon.mass_error);
      }
      else {
        //console.log(element.matched_ions.length);
        for (let i = 0; i < element.matched_ions.length; i++) {
          let matchedIon = element.matched_ions[i];
          let ionType = matchedIon.ion_type;
          if (ionType == "Z_DOT") {
            ionType = "Z\u02D9";
          }
          let curIonText = ionType + matchedIon.ion_display_position;
          ionText = ionText + curIonText + " ";
          if (parseFloat(matchedIon.mass_error) > 0) {
            massError = parseFloat(matchedIon.mass_error);
          }
        }
      }
      let peakId = element.peak_id;
      //console.log(peakId, envelopes.length, specId);
      //console.log(envelopes[peakId]);
      let envPeaks = envelopes[peakId].env_peaks;
      envPeaks.sort(function(x,y){
        return d3.descending(x.intensity, y.intensity);
      });
      let x = parseFloat(envPeaks[0].mz);
      let y = parseFloat(envPeaks[0].intensity);
      let ionData = {"mz": x, "intensity": y, "text": ionText, "error": massError};
      //console.log(ionData);
      ionData.env = envelopes[peakId]; 
      ions.push(ionData);
      let monoX = parseFloat(envelopes[peakId].mono_mass);
      let monoY = 0; 
      envPeaks.forEach(element => monoY += element.intensity); 
      let monoIonData = {"mz": monoX, "intensity": monoY, "text": ionText, "error": massError};
      addOneIon(monoIons, monoIonData);
    }
  });
  return [ions, monoIons];
}

function addOneIon(ionList, ion) {
  let idx = -1;
  for (let i = 0; i < ionList.length; i++) {
    if (ion.text == ionList[i].text) {
      idx = i;
      break;
    }
  }
  if (idx == -1) {
    ionList.push(ion);
  }
  else {
    if (ion.intensity > ionList[idx].intensity) {
      ionList[idx].intensity = ion.intensity;
    }
  }
}
