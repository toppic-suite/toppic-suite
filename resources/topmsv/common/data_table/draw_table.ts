/**
 * Build the monomass table with all the data from the peak variable of the prsm_data
 * Provide a unique class name to the m/z values to provide on click action to xoom the
 * spectrum graph to that position
 */
 class DataTable {
    private prsmObj_: Prsm;
    private showScanNum_: boolean;
    private ms2GraphList_: SpectrumView[];
    private addOneToPeakNum_: boolean;
    private specSvgId_: string = "";
    private monoMassSvgId_: string = "";

    constructor(prsmObj: Prsm, showScanNum: boolean, ms2GraphList: SpectrumView[]) {
        this.prsmObj_ = prsmObj;
        this.showScanNum_ = showScanNum;
        this.ms2GraphList_ = ms2GraphList;
        if (this.showScanNum_) {
            this.addOneToPeakNum_ = true;
        }
        else {
            this.addOneToPeakNum_ = false;
        }
    }
    setSpecSvgId(svgId: string): void {
      this.specSvgId_ = svgId;
    }
    setMonoMassSvgId(svgId: string): void {
      this.monoMassSvgId_ = svgId;
    }
    addEventHandlers(): void {
        var _a: HTMLElement | null, _b: HTMLElement | null, _c: HTMLElement | null;
        (_a = document.getElementById("all_peak_count")) === null || _a === void 0 ? void 0 : _a.addEventListener("click", this.showAllPeaks);
        (_b = document.getElementById("matched_peak_count")) === null || _b === void 0 ? void 0 : _b.addEventListener("click", this.showMatchedPeaks);
        (_c = document.getElementById("not_matched_peak_count")) === null || _c === void 0 ? void 0 : _c.addEventListener("click", this.showNotMatchedPeaks);
        // On click of break points
        $(".break_point").click((e) => {
            let pos: string | undefined = $(e.target).attr("ion_pos");
            if (!pos) {
                console.error("ERROR: invalid ion position");
                return;
            }
            this.showIonPeaks(pos);
        });
        //mono mz click
        $(".row_mono_mz").click((e) => {
            /*	get Mono M/z value till 3 decimal values	*/
            let monoMz: number = parseFloat(parseFloat(e.currentTarget.innerHTML).toFixed(3));

            if (typeof ms2ScanList != "undefined") {
              let parentId: string = $(e.currentTarget).parent().parent().prop('id');//ion type

              let parentDiv = document.getElementById(parentId);
              if (!parentDiv) {
                  return;
              }
              let scanNumDiv: HTMLElement | null = <HTMLElement>parentDiv.firstChild;
              if (!scanNumDiv) {
                  return;
              }
              let scanNum: string = scanNumDiv.innerHTML;
              if (scanNum === null) {
                  console.error("ERROR: scan number is null");
              }
              for (let i = 0; i < ms2ScanList.length; i++) {
                  let listId: string = "ms2_svg_div_graphlist_" + i;
                  let graphId: string = this.specSvgId_ + i;
                  let monolistId: string = "ms2_svg_div_monographlist_" + i;
                  let monoGraphId: string = this.monoMassSvgId_ + i;
                  let listElement: HTMLElement | null = document.getElementById(listId);
                  let graphElement: HTMLElement | null = document.getElementById(graphId);
                  let monoListElement: HTMLElement | null = document.getElementById(monolistId);
                  let monoGraphElement: HTMLElement | null = document.getElementById(monoGraphId);
                  if (scanNum == ms2ScanList[i]) {
                    if (listElement) {
                      listElement.classList.add("active");
                    }
                    if (graphElement) {
                      graphElement.style.display = "";
                    }
                    if (monoListElement) {
                      monoListElement.classList.remove("active");
                    }
                    if (monoGraphElement) {
                      monoGraphElement.style.display = "none";
                    }
                    let spGraph = this.ms2GraphList_[i];
                    // set monoMz to do
                    spGraph.getPara().updateMzRange(monoMz);
                    spGraph.redraw();
                  }
                  else {
                    if (listElement) {
                      listElement.classList.remove("active");
                    }
                    if (graphElement) {
                      graphElement.style.display = "none";
                    }
                    if (monoListElement) {
                      monoListElement.classList.remove("active");
                    }
                    if (monoGraphElement) {
                      monoGraphElement.style.display = "none";
                    }
                  }
              }
              showMs2Graph();
            }
            else {
              let listId = "ms2_svg_div_graphlist_0";
              let graphId = this.specSvgId_
              let monolistId = "ms2_svg_div_monographlist_0";
              let monoGraphId = this.monoMassSvgId_;
              let listElement = document.getElementById(listId);
              let graphElement = document.getElementById(graphId);
              let monoListElement = document.getElementById(monolistId);
              let monoGraphElement = document.getElementById(monoGraphId);
              //console.log("graphId", graphId);
              //console.log("monolistId", monolistId);
              switchTab("graphlist");

             /* if (listElement) {
                listElement.classList.add("active");
              }
              if (graphElement) {
                graphElement.style.display = "";
              }
              if (monoListElement) {
                monoListElement.classList.remove("active");
              }
              if (monoGraphElement) {
                monoGraphElement.style.display = "none";
              }*/
              let spGraph = this.ms2GraphList_[0];
              // set monoMz to do
              spGraph.getPara().updateMzRange(monoMz);
              spGraph.redraw();
            }
        });
    }
    addTable() {
        if (this.showScanNum_) {
            //@ts-ignore
            $('#spectrum').DataTable({
                "scrollY": "400px",
                "scrollCollapse": true,
                "paging": false,
                "order": [[1, "asc"]],
                "bSortClasses": false,
                "columns": [
                    { "type": "num" },
                    { "type": "num" },
                    { "type": "num" },
                    null,
                    { "type": "num" },
                    { "type": "num" },
                    { "type": "num" },
                    null,
                    { "type": "num" },
                    { "type": "num" },
                    { "type": "num" }
                ]
            });
        }
        else {
            //@ts-ignore
            $('#spectrum').DataTable({
                "scrollY": "400px",
                "scrollCollapse": true,
                "paging": false,
                "destroy": true,
                "order": [[1, "asc"]],
                "bSortClasses": false,
                "columns": [
                    { "type": "num", "visible": false },
                    { "type": "num" },
                    { "type": "num" },
                    null,
                    { "type": "num" },
                    { "type": "num" },
                    { "type": "num" },
                    null,
                    { "type": "num" },
                    { "type": "num" },
                    { "type": "num" }
                ]
            });
        }
    }
    drawTable(): void {
      var _a: HTMLElement | null;
      let self: any = this;
      let ms2Spectrum: Spectrum[] | null = this.prsmObj_.getMs2Spectra();
      let seqLength = this.prsmObj_.getProteoform().getLastPos() - this.prsmObj_.getProteoform().getFirstPos() + 1; 
      if (!ms2Spectrum) {
        console.error("ERROR: invalid ms2 spectrum");
        return;
      }
      let table: HTMLElement = <HTMLElement>document.getElementById('spectrum');
      if (document.getElementById("peaksTableBody")) {
        (_a = document.getElementById("peaksTableBody")) === null || _a === void 0 ? void 0 : _a.remove();
      }
      let tbdy: HTMLTableSectionElement = document.createElement('tbody');
      tbdy.setAttribute("id", "peaksTableBody");
      let l_scans: string[] = [];
      let l_specIds: string[] = [];
      ms2Spectrum.forEach((spectra) => {
        l_scans.push(spectra.getScanNum());
        l_specIds.push(spectra.getSpectrumId());
      });
      let l_matched_peak_count: number = 0;
      let l_duc_peak_count: number = 0;
      let recordedPeaks: string[] = []; //Id of peaks that were matched by more than 1 ion
          
      ms2Spectrum.forEach((spectra, index) => {
        // let specId: string = spectra.getSpectrumId();
        let decovPeaks: Peak[] | null = spectra.getDeconvPeaks();
        let peakMz: number[] = [];
        if (!decovPeaks) {
          console.error("ERROR: no deconvoluted peaks in ms2 spectrum");
          return;
        }  
        this.prsmObj_.getMatchedPeakEnvelopePairs().forEach((matchedPeakPair) => {
          peakMz.push(matchedPeakPair.getPeak().getMonoMz());
          if (index < 1) { //avoid same matched peak added multiple times
              loop_matched_ions(matchedPeakPair.getPeak(), matchedPeakPair.getPeak().getSpecId(), true, matchedPeakPair, seqLength);
          }
        })
        //use monoMz to see if the peak has already been added
        decovPeaks.forEach((deconvPeak) => {
          if (peakMz.indexOf(deconvPeak.getMonoMz()) < 0) {
              loop_matched_ions(deconvPeak, deconvPeak.getSpecId(), false);
          }
      });
      })
      //after looping through the prsm files, store the ion type data to local storage
      //window.localStorage.setItem('ionType', ionArray.toString());
      /**
       * Inner function to create a rows and columns for monomass table
       * @param {object} peak - contains information of each peak
       * @param {int} i - index of the peak
       */
      function loop_matched_ions(peak: Peak, specId: string | undefined, matchedPeaks: boolean, 
        matchedPeakPair?: MatchedPeakEnvelopePair, seqLen?: number) {
        /*Create row for each peak value object in the table*/
        let tr: HTMLTableRowElement = document.createElement('tr');
        let id: string = "";
        if (specId) {
          id = specId + "peak" + peak.getId();
        }
        let l_scan: string;
        let l_class: string;
          
        if((parseInt(peak.getId()) + 1)%2 == 0){
          // class name helps to get unmatched peaks when clicking unmatched peaks
          l_class = "unmatched_peak even"; 
        }
        else{
          // class name helps to get unmatched peaks when clicking unmatched peaks
          l_class = "unmatched_peak odd"; 
        }
        if (matchedPeaks && matchedPeakPair) {
          id = id + matchedPeakPair.getIon().getName();
          if ((parseInt(peak.getId()) + 1) % 2 == 0) {
            // class name helps to get matched peaks when clicking matched peaks
            l_class = "matched_peak even";
          }
          else {
            // class name helps to get matched peaks when clicking matched peaks
            l_class = "matched_peak odd";
          }
          l_matched_peak_count++;
          //	create a name for each row
          let ion: string = matchedPeakPair.getIon().getName();
          let ionPos: number = parseInt((matchedPeakPair.getIon().getId()).slice(1));
          if (ion.indexOf("X") >= 0 || ion.indexOf("Y") >= 0 || ion.indexOf("Z") >= 0) {
              if (seqLen) {
                  ionPos = seqLen - (ionPos);
              }
          }
          tr.setAttribute("name", ionPos.toString());
          let peakId: string = peak.getId();
          if (recordedPeaks.indexOf(peakId) < 0) {
            recordedPeaks.push(peakId);
          }
          else{
            l_duc_peak_count++;
          }

        }
        //	Set "id","class name" and "role" for each row
        tr.setAttribute("id", id);
        tr.setAttribute("class", l_class);
        tr.setAttribute("role", "row");
        for (let i = 0; i < 11; i++) {
          var td: HTMLTableDataCellElement = document.createElement('td');
          td.setAttribute("align", "center");
          if (i == 0) {
            if (specId == l_specIds[0])
              l_scan = l_scans[0];
            else
              l_scan = l_scans[1];
              td.setAttribute("class", "row_scanIds");
              td.innerHTML = l_scan;
          }
          if (i == 1) {
            if (self.addOneToPeakNum_) {
                td.innerHTML = (parseInt(peak.getId()) + 1).toString();
            }
            else {
                td.innerHTML = peak.getId().toString();
            }
            td.setAttribute("class", "row_peakNum");
        }
          if (i == 2) {
            let monoMass: number | undefined = peak.getMonoMass();
            if (!monoMass) {
              console.error("ERROR: mono peak does not have mono mass");
              return;
            }
            td.innerHTML = FormatUtil.formatFloat(monoMass.toString(), "dataTable");
            td.setAttribute("class", "row_monoMass");
          }
          if (i == 3) {
            //	provide link to click on m/z value to view spectrum 
            let a: HTMLAnchorElement = document.createElement('a');
            a.href = "#!";
            a.className = "row_mono_mz";
            a.innerHTML = FormatUtil.formatFloat(peak.getMonoMz().toString(), "dataTable");
            td.appendChild(a);
          }
          if (i == 4) {
            //td.innerHTML = peak.getIntensity().toString();
            let intensity: string = (peak.getIntensity()).toExponential();
            td.innerHTML = Number.parseFloat(intensity).toPrecision(3);
            td.setAttribute("class", "row_intensity");
          }
          if (i == 5) {
            let charge: number | undefined = peak.getCharge();
            if (!charge) {
              console.error("ERROR: mono peak does not have charge");
              return;
            }
            td.innerHTML = charge.toString();
            td.setAttribute("class", "row_charge");
          }
          if (matchedPeaks && matchedPeakPair) {
            if (i == 6) {
              td.innerHTML = FormatUtil.formatFloat(matchedPeakPair.getTheoMass().toString(), "dataTable");
            }
            if (i == 7) {
              let ionPos: string = matchedPeakPair.getIon().getId();
              td.innerHTML = matchedPeakPair.getIon().getName() + ionPos.slice(1);
            }
            if (i == 8) {
            //if c-term ion, pos = pos + 1
              let ion: string = matchedPeakPair.getIon().getName();
              let ionPos: number = parseInt((matchedPeakPair.getIon().getId()).slice(1));
              if (ion.indexOf("X") >= 0 || ion.indexOf("Y") >= 0 || ion.indexOf("Z") >= 0) {
                if (seqLen){
                  td.innerHTML = (seqLen - ionPos).toString();
                }
                else{
                  console.error("ERROR: seqLen is not provided");
                }
              }
              else {
                td.innerHTML = ionPos.toString();
              }
            }
            if (i == 9) {
              let massError: number | undefined = matchedPeakPair.getIon().getMassError();
              if (massError == undefined) {
                console.error("ERROR: massError is not provided");
              }
              else {
                  td.innerHTML = FormatUtil.formatFloat(massError.toString(), "dataTable");
              }
            }
            if (i == 10) {
              let ppmError: number | undefined = matchedPeakPair.getIon().getPpmError();
              if (ppmError == undefined) {
                console.log(matchedPeakPair.getIon());
                console.error("ERROR: ppmError is not provided");
              }
              else {
                  td.innerHTML = FormatUtil.formatFloat(ppmError.toString(), "ppmError");
              }
            }
          }
            tr.appendChild(td);
          }
          tbdy.appendChild(tr);
      }
      //let l_All_Peaks: number = peakCnt;
      let l_All_Peaks: number = 0;
      ms2Spectrum.forEach((spectra) => {
          let decovPeaks: Peak[] | null = spectra.getDeconvPeaks();
          if (!decovPeaks) {
              console.error("ERROR: no deconvoluted peaks in ms2 spectrum");
              return;
          }
          l_All_Peaks = l_All_Peaks + decovPeaks.length;
      });
      l_matched_peak_count = l_matched_peak_count - l_duc_peak_count;
      let l_not_matched_peak_count: number = l_All_Peaks - l_matched_peak_count;
      let allPeakCnt: HTMLElement | null = document.getElementById("all_peak_count");
      let matchedPeakCnt: HTMLElement | null = document.getElementById("matched_peak_count");
      let notMatchedPeakCnt: HTMLElement | null = document.getElementById("not_matched_peak_count");
      if (allPeakCnt) {
          allPeakCnt.innerHTML = "All peaks (" + l_All_Peaks.toString() + ")";
      }
      if (matchedPeakCnt) {
          matchedPeakCnt.innerHTML = "Matched peaks (" + l_matched_peak_count.toString() + ")";
      }
      if (notMatchedPeakCnt) {
          notMatchedPeakCnt.innerHTML = "Not Matched peaks (" + l_not_matched_peak_count.toString() + ")";
      }
      if (!table) {
        console.error("ERROR: table element is not created correctly");
        return;
      }
      table.style.display = "";
      table.appendChild(tbdy);
      this.addTable();
      this.addEventHandlers();
      if (document.getElementById("peakCount")) {
          document.getElementById("peakCount")!.style.display = "inline-block";
      }    
    }
    showMatchedPeaks() {
        var elems: HTMLCollectionOf<HTMLElement> = <HTMLCollectionOf<HTMLElement>> document.getElementsByClassName("matched_peak");
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "";
        }
        elems = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName("unmatched_peak");
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "none";
        }
        //$('div.dataTables_scrollBody').height(400);
    }
    /**
     * Function to show only unmatched peaks on click of unmatched peaks button
     */
    showNotMatchedPeaks() {
        var elems: HTMLCollectionOf<HTMLElement> = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName("matched_peak");
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "none";
        }
        elems = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName("unmatched_peak");
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "";
        }
        //$('div.dataTables_scrollBody').height(400);
    }
    /**
     * Function to show all peaks on click of All peaks button
     */
    showAllPeaks() {
        var elems: HTMLCollectionOf<HTMLElement> = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName('matched_peak');
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "";
        }
        elems = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName('unmatched_peak');
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = "";
        }
        //$('div.dataTables_scrollBody').height(400);
    }
    /**
     * This gets invoked on click of annotation in the SVG of sequence at matched positions
     * Function to show only ions matched at a particular position
     * @param {String} ids - contains name of the tag
     */
    showIonPeaks(ids: string): void {
        var elems: HTMLCollectionOf<HTMLElement> = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName('matched_peak');
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = 'none';
        }
        elems = <HTMLCollectionOf<HTMLElement>>document.getElementsByClassName('unmatched_peak');
        for (var i: number = 0; elems.length > i; i++) {
            elems[i].style.display = 'none';
        }
        let elemsTemp: NodeListOf<HTMLElement> = document.getElementsByName(ids);
        for (var j: number = 0; elemsTemp.length > j; j++) {
            elemsTemp[j].style.display = "";
            //elems[j].style.background  =  "#BEECFF";
        }
    }
  }
