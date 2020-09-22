/**
 * Function to form table for all selected fragmented ions
 * @param {string} sequence - user entered sequence without mass shifts embedded
 * @param {Array} matchedUnMatchedPeaks - list of all the calculated masses
 */
function createTableForSelectedFragmentIons(sequence,matchedUnMatchedPeaks,spectrumGraph){
    // console.log("matchedUnMatchedPeaks : ", matchedUnMatchedPeaks);
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
        th.setAttribute("class","th-sm");
        th.innerHTML = matchedUnMatchedPeaks[i].ionFragment;
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
        if (j === seqln - 1) {
            for (let k =0; k < len; k++) {
                let td1 = document.createElement("td");
                td1.setAttribute("class","td_fragments");
                td1.setAttribute("charge",null);
                td1.innerHTML = "Null";
                tr1.appendChild(td1);
            }
            tbody.appendChild(tr1);
            break;
        }
        /**
         * Add mass data to respected columns 
         */
        for(let k=0;k<len;k++)
        {
            let td1 = document.createElement("td");
            let index = j;
            if(matchedUnMatchedPeaks[k].ionFragment[0] === "x" || matchedUnMatchedPeaks[k].ionFragment[0] === "y"
                || matchedUnMatchedPeaks[k].ionFragment[0] === "z")
            {
                /**
                 * position when suffix mass list is written to table
                 */
                index = seqln-j-2;
            }
            td1.setAttribute("class","td_fragments");
            if(matchedUnMatchedPeaks[k].massList[index].matchedInd === "Y")
            {
                td1.setAttribute("class","td_fragments matched_fragments");
            }
            td1.setAttribute("charge",matchedUnMatchedPeaks[k].massList[index].charge);
            td1.innerHTML = matchedUnMatchedPeaks[k].massList[index].mass.toFixed(4);
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
    dict["a"] = "a";
    dict['a1'] = "a-H<sub>2</sub>O";
    dict['a2'] = "a-NH<sub>3</sub>";

    dict["b"] = "b";
    dict['b1'] = "b-H<sub>2</sub>O";
    dict['b2'] = "b-NH<sub>3</sub>";

    dict["c"] = "c";
    dict['c1'] = "c-H<sub>2</sub>O";
    dict['c2'] = "c-NH<sub>3</sub>";

    dict["x"] = "x";
    dict['x1'] = "x-H<sub>2</sub>O";
    dict['x2'] = "x-NH<sub>3</sub>";

    dict["y"] = "y";
    dict['y1'] = "y-H<sub>2</sub>O";
    dict['y2'] = "y-NH<sub>3</sub>";

    dict["z"] = "z";
    dict['z1'] = "z-H<sub>2</sub>O";
    dict['z2'] = "z-NH<sub>3</sub>";

    dict["z_"] = "z&deg;";
    dict['z_1'] = "z&deg;-H<sub>2</sub>O";
    dict['z_2'] = "z&deg;-NH<sub>3</sub>";

    return dict[ionType];
}