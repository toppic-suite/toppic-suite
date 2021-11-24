/**
 * Function to form table for all selected fragmented ions
 * @param {string} sequence - user entered sequence without mass shifts embedded
 * @param {Array} matchedUnMatchedPeaks - list of all the calculated masses
 */
function createTableForSelectedFragmentIons(sequence: string, matchedUnMatchedPeaks: MatchedUnMatchedObj[], spectrumGraph: SpectrumView) {
  /**
  * Remove if table already exist and rebuild the table
  */
  $("#selectedIonTableContainer").remove();
  $("#divselectediontablecontainer #selectedIonTableContainer_wrapper").remove();
  // jqueryElements.ionTableContainer.remove();
  let div: HTMLElement = <HTMLElement>document.getElementById("divselectediontablecontainer");
  let table: HTMLTableElement = document.createElement("table");
  table.setAttribute("id","selectedIonTableContainer");
  table.setAttribute("class","table table-striped display");

  let thead: HTMLTableSectionElement = document.createElement("thead");

  let tbody: HTMLTableSectionElement = document.createElement("tbody");
  tbody.setAttribute("id","selectedIonTableContainertbody");

  let tr: HTMLTableRowElement = document.createElement("tr");
  tr.setAttribute("role","row");
    
  let len: number = matchedUnMatchedPeaks.length;
  let seqln: number = sequence.length;
  /**
  * Create ID column Header
  */
  let th1: HTMLTableHeaderCellElement = document.createElement("th");
  let th0: HTMLTableHeaderCellElement = document.createElement("th");
  th0.setAttribute("class","th-sm");
  th0.innerHTML = "Id";
  /**
  * Create Acid Column Header
  */
  th1.setAttribute("class","th-sm");
  th1.innerHTML = "Amino acid";
  tr.appendChild(th0);
  tr.appendChild(th1);
  /**
  * Create other header of fragmented ion by taking from list
  */
  for(let i=0;i<len;i++)
  {
    let th: HTMLTableHeaderCellElement = document.createElement("th");
    let ionName: string = matchedUnMatchedPeaks[i].ionFragment;
    ionName = convertIonName(ionName); 
        
    th.setAttribute("class","th-sm");
    th.innerHTML = ionName;
    tr.appendChild(th);
  }
  /**
  * Create columns from the input list of matched and unmatched peaks
  */
  for(let j=0;j<seqln;j++)
  {
    let tr1: HTMLTableRowElement = document.createElement("tr");
    tr1.setAttribute("role","row");

    let td: HTMLTableDataCellElement = document.createElement("td");
    td.setAttribute("class","td_fragments");
    td.innerHTML = sequence[j];

    let td0: HTMLTableDataCellElement = document.createElement("td");
    td0.setAttribute("class","td_fragments");
    /**
    * position starts from 1
    */
    td0.innerHTML = (j+1).toString();
    tr1.appendChild(td0);
    tr1.appendChild(td);
        
    /**
    * Add mass data to respected columns 
    */
    for(let k=0;k<len;k++)
    {
      let td1: HTMLTableDataCellElement = document.createElement("td");
      let index: number = j;
      let isNull: boolean = false;//write Null if true

      if(matchedUnMatchedPeaks[k].ionFragment[0] === "X" || matchedUnMatchedPeaks[k].ionFragment[0] === "Y"
      || matchedUnMatchedPeaks[k].ionFragment[0] === "Z")
      {
        if (index === 0)//at the first row in the c-term column
        {
          isNull = true;
        }
        /**
        * position when suffix mass list is written to table
        */
        index = seqln-j-1;
      }
      else if(matchedUnMatchedPeaks[k].ionFragment[0] === "A" || matchedUnMatchedPeaks[k].ionFragment[0] === "B"
      || matchedUnMatchedPeaks[k].ionFragment[0] === "C")
      {
        if (index === seqln - 1)//at the last row in the n-term column
        {
          isNull = true;
        }
      }
      if(isNull)
      {
        td1.setAttribute("class","td_fragments");
        td1.setAttribute("charge", "");
        td1.innerHTML = "-";
      }
      else
      {
        td1.setAttribute("class","td_fragments");
        if(matchedUnMatchedPeaks[k].massList[index].matchedInd === "Y")
        {
          td1.setAttribute("class","td_fragments matched_fragments");
        }
        td1.setAttribute("charge",matchedUnMatchedPeaks[k].massList[index].charge.toString());
        td1.innerHTML = matchedUnMatchedPeaks[k].massList[index].mass.toFixed(4);
      }
      tr1.appendChild(td1);
    }
    tbody.appendChild(tr1);
  }
  thead.appendChild(tr);
  table.appendChild(thead);
  table.appendChild(tbody);
  div.appendChild(table);
  onClickofMatchedPeaks(spectrumGraph);
}

/**
 * Function to zoom the graph to the mass point on click of matched mass
 */
function onClickofMatchedPeaks(spectrumGraph: SpectrumView): void {
  $(".matched_fragments").click((cell: JQuery.ClickEvent) => {
    let charge: string | null = cell.currentTarget.getAttribute("charge");
    let mass: number = parseFloat(cell.currentTarget.innerHTML);

    if (!charge) {
      console.error("ERROR: charge of selected peak is null");
      return;
    }
    // console.log("mass:",mass);
    // console.log("mz:", mz);
    spectrumGraph.redraw(mass);

    //switch the graph tab to mass graph 
    switchTab("monographlist");
  })
}

/**
 * Get all the N terminus Ions from UI
 */
 function getNterminusCheckedList(): Ion[] {
  let ions: Ion[] = [];
  let cnt: number = 0;
  $.each($("input[name='nterminus']:checked"), function() {
    let id: string | undefined = $(this).attr("id");
    let value: string | number | string[] | undefined = $(this).val();
    if (typeof (id) != 'string' || typeof (value) != 'string') {
      console.error("ERROR: Nterminus ion value is not parsed correctly");
      return ions;
    }
    let ionType: string = getActualIdvalues(id);
    let ion = new Ion(cnt.toString(), ionType, "N", parseFloat(value));
    ions.push(ion);
    cnt++;
  });
  return ions;
}
/**
 * Get all the checked C ternimus ions from UI
 */
 function getCterminusCheckedList(): Ion[] {
  let ions: Ion[] = [];
  let cnt: number = 0;
  $.each($("input[name='cterminus']:checked"), function() {
    let id: string | undefined = $(this).attr("id");
    let value: string | number | string[] | undefined = $(this).val();
    if (typeof (id) != 'string' || typeof (value) != 'string') {
      console.error("ERROR: Cterminus ion value is not parsed correctly");
      return ions;
    }
    let ionType: string = getActualIdvalues(id);
    let ion = new Ion(cnt.toString(), ionType, "C", parseFloat(value));
    ions.push(ion);
    cnt++;
  });
  return ions;
}

/**
 * @function getActualIdvalues
 * @description Dictionary to get correspoding actual heading for 
 * all the Id's of Ion Fragmentation from UI.
 * @param {String} ionType - Id of the ion from UI
 */
function getActualIdvalues(ionType: string): string {
  let dict: {[key:string]:string} = {};
  dict["a"] = "A";
  dict['a1'] = "A-H2O";
  dict['a2'] = "A-NH3";

  dict["b"] = "B";
  dict['b1'] = "B-H2O";
  dict['b2'] = "B-NH3";

  dict["c"] = "C";
  dict['c1'] = "C-H2O";
  dict['c2'] = "C-NH3";

  dict["x"] = "X";
  dict['x1'] = "X-H2O";
  dict['x2'] = "X-NH3";

  dict["y"] = "Y";
  dict['y1'] = "Y-H2O";
  dict['y2'] = "Y-NH3";

  dict["z"] = "Z";
  dict['z1'] = "Z-H2O";
  dict['z2'] = "Z-NH3";

  dict["z_"] = "Z_DOT";
  dict['z_1'] = "Z_DOT-H2O";
  dict['z_2'] = "Z_DOT-NH3";
  return dict[ionType];
}