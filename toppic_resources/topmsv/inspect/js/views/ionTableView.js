/**
 * Function to form table for all selected fragmented ions
 * @param {string} sequence - user entered sequence without mass shifts embedded
 * @param {Array} matchedUnMatchedPeaks - list of all the calculated masses
 */
function createTableForSelectedFragmentIons(sequence,matchedUnMatchedPeaks,spectrumGraph){
    /**
     * Remove if table already exist and rebuild the table
     */
    $("#selectedIonTableContainer").remove();
    $("#divselectediontablecontainer #selectedIonTableContainer_wrapper").remove();
    // jqueryElements.ionTableContainer.remove();
    let div = document.getElementById("divselectediontablecontainer");
    let table = document.createElement("table");
    table.setAttribute("id","selectedIonTableContainer");
    table.setAttribute("class","table table-striped display");

    let thead = document.createElement("thead");

    let tbody = document.createElement("tbody");
    tbody.setAttribute("id","selectedIonTableContainertbody");

    let tr = document.createElement("tr");
    tr.setAttribute("role","row");
    
    let len = matchedUnMatchedPeaks.length;
    let seqln = sequence.length;
    /**
     * Create ID column Header
     */
    let th1 = document.createElement("th");
    let th0 = document.createElement("th");
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
        let th = document.createElement("th");
        let ionName = matchedUnMatchedPeaks[i].ionFragment;
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
        let tr1 = document.createElement("tr");
        tr1.setAttribute("role","row");

        let td = document.createElement("td");
        td.setAttribute("class","td_fragments");
        td.innerHTML = sequence[j];

        let td0 = document.createElement("td");
        td0.setAttribute("class","td_fragments");
        /**
         * position starts from 1
         */
        td0.innerHTML = j+1;
        tr1.appendChild(td0);
        tr1.appendChild(td);
        
        /**
         * Add mass data to respected columns 
         */
        for(let k=0;k<len;k++)
        {
            let td1 = document.createElement("td");
            let index = j;
            let isNull = false;//write Null if true

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
                td1.setAttribute("charge",null);
                td1.innerHTML = "Null";
            }
            else
            {
                td1.setAttribute("class","td_fragments");
                if(matchedUnMatchedPeaks[k].massList[index].matchedInd === "Y")
                {
                    td1.setAttribute("class","td_fragments matched_fragments");
                }
                td1.setAttribute("charge",matchedUnMatchedPeaks[k].massList[index].charge);
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
function onClickofMatchedPeaks(spectrumGraph){
    $(".matched_fragments").click(function() {
        let charge = $(this).attr("charge");
        let mass = parseFloat($(this).html());
        let mz = mass/charge;
        // console.log("mass:",mass);
        // console.log("mz:", mz);
        spectrumGraph.redraw(mass);
    })
}

/**
 * Get all the N terminus Ions from UI
 */
let getNterminusCheckedList = () => {
    let ions = [];
    $.each($("input[name='nterminus']:checked"), function(){
            let id = $(this).attr("id");
            let value = $(this).val();
            let ionType = getActualIdvalues(id);
            let temp_Obj = {ionType:ionType,mass:value}
            ions.push(temp_Obj);
        });
    return ions;
}
/**
 * Get all the checked C ternimus ions from UI
 */
let getCterminusCheckedList = () => {
    let ions = [];
    $.each($("input[name='cterminus']:checked"), function(){
            let id = $(this).attr("id");
            let value = $(this).val();
            let ionType = getActualIdvalues(id);
            let temp_Obj = {ionType:ionType,mass:value}
            ions.push(temp_Obj);
        });
    return ions;
}

/**
 * @function getActualIdvalues
 * @description Dictionary to get correspoding actual heading for 
 * all the Id's of Ion Fragmentation from UI.
 * @param {String} ionType - Id of the ion from UI
 */
let getActualIdvalues = (ionType) => {
    let dict = [];
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