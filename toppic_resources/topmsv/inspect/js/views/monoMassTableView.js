/**
 * Function to create table
 */
function createMonoMassTable(){
    // Remove if table already exist and rebuild the table
    $("#tableContainer").remove();
    $("#divtableContainer #tableContainer_wrapper").remove();
    // Remove if tableContainer_wrapper already exist and rebuild the table, wrapper is an inbuild class of bootstrap table
    // jqueryElements.monoMassTableContainer.remove();
    let div = domElements.monoMassTableContainer;
    let table = document.createElement("table");
    table.setAttribute("id","tableContainer");
    table.setAttribute("class","table table-striped display");
    let thead = document.createElement("thead");
    let tbody = document.createElement("tbody");
    tbody.setAttribute("id","tableContainertbody");
    let tr = document.createElement("tr");
    tr.setAttribute("role","row");
    let colCount = 10;
    // Create table header
    for(let i = 0;i < colCount;i++)
    {
        let th = document.createElement("th");
        th.setAttribute("class","th-sm");
        if(i === 0) th.innerHTML = "Peak";
        if(i === 1) th.innerHTML = "Mono mass";
        if(i === 2) th.innerHTML = "Mono m/z";
        if(i === 3) th.innerHTML = "Intensity";
        if(i === 4) th.innerHTML = "Charge";
        if(i === 5) th.innerHTML = "Theoretical mass";
        if(i === 6) th.innerHTML = "Ion";
        if(i === 7) th.innerHTML = "Pos";
        if(i === 8) th.innerHTML = "Mass error";
        if(i === 9) th.innerHTML = "PPM error";
        tr.appendChild(th);
    }
    thead.appendChild(tr);
    table.appendChild(thead);
    table.appendChild(tbody);
    div.appendChild(table);

    addPrsmGraphClickHandler();
}
/** add event handler to the table
 * (filter rows when breakpoint in clicked in prsm graph)
*/
function addPrsmGraphClickHandler(){
    $(".break_point").click(function() {
        let pos = $(this).attr("ion_pos"); 
        showIonPeaks(pos); 
    });
}
/**
 * Add Data to the table created
 * @param {Array} matchedPeaks - Contains List of Matched and unmatched peaks
 */
function addMassDataToTable(matchedPeaks, graphObj)
{
    let dataContainer_tbody = $("#tableContainer tbody");
    const totalColCount = $("#tableContainer thead tr th").length;
    let len = matchedPeaks.length;
    for(let i=0; i<len; i++)
    {
        // let rowSize = $("#tableContainer tbody tr").length;
        let tr = document.createElement("tr");
        tr.setAttribute("id",i+"_row");
        tr.setAttribute("name",matchedPeaks[i].position);
        let id = 0;
        for(let j = 0; j < totalColCount; j++)
        {
            let td = document.createElement("td");
            if(j === 0)
            {
                id = matchedPeaks[i].peakId;
                td.innerHTML = matchedPeaks[i].peakId;
                td.style.fontWeight = "bold";
            }else if(j === 1) td.innerHTML = matchedPeaks[i].mass;
            else if(j === 4) td.innerHTML = matchedPeaks[i].charge;
            else if(j === 2) {
                let mz = matchedPeaks[i].mass / matchedPeaks[i].charge + 1.007276466879;
                let a = document.createElement('a');
                a.href="#!"
                a.className = "peakRows"
                a.innerHTML = mz.toFixed(4); 
                td.appendChild(a);
            }
            else if(j === 3) td.innerHTML = matchedPeaks[i].intensity ;
            else if(j === 5){
                td.className = "th_mass";
                td.innerHTML = matchedPeaks[i].thMass;
            } 
            else if(j === 6) {
                let ionName = matchedPeaks[i].ion;
                ionName = convertIonName(ionName);
                td.innerHTML = ionName;
            }
            else if(j === 7) {
                //if c-term ion, pos = pos + 1
                let ion = matchedPeaks[i].ion;
                if (ion.indexOf("X") >= 0 || ion.indexOf("Y") >= 0 || ion.indexOf("Z") >= 0)
                {
                    td.innerHTML = parseInt(matchedPeaks[i].position) + 1;
                }
                else
                {
                    td.innerHTML = matchedPeaks[i].position;
                }
            }
            else if(j === 8) td.innerHTML = matchedPeaks[i].massError;
            else if(j === 9) td.innerHTML = matchedPeaks[i].PPMerror;
            tr.appendChild(td);
        }
        tr.setAttribute("role","row");
        let classname = "" ;
        // Creating classes with matched_peak even and odd, this will help to show only matched peaks on click of matched peaks
        if(matchedPeaks[i].matchedInd === "Y")
        {
            if(id%2 === 0) classname = "matched_peak even";
            else classname = "matched_peak odd";
        }
        else
        {
            if(id%2 === 0) classname = "unmatched_peak even";
            else classname = "unmatched_peak odd";
        }
        tr.setAttribute("class",classname);
        dataContainer_tbody.append(tr);
    }

    $(".peakRows").click(function() {
        /*	get Mono M/z value till 3 decimal values	*/
        let peak_value = parseFloat(this.innerHTML).toFixed(3);
        let parent_id  = $(this).parent().parent().prop('id');
        let th_mass_val = $("#"+parent_id+" .th_mass").text();
        // console.log("th_mass_val : ",th_mass_val);
        graphObj.redraw(peak_value);
    });
}

/**
 * Function to display matched count and un-matched count
 */
function showPeakCounts(monoMassList, matchedIonList)
{
    let totalCount = monoMassList.length;
    let peakId = [];//array containing unique peak ids
    for (let i = 0; i < matchedIonList.length; i++){
        if (peakId.indexOf(matchedIonList[i].peakId) < 0){
            peakId.push(matchedIonList[i].peakId);
        }
    }
    let matchedCount = peakId.length;
    let unMatchedCount = totalCount - matchedCount;
    jqueryElements.allPeakCount.html("All Peaks ("+totalCount+")");
    jqueryElements.matchedPeakCount.html("Matched Peaks ("+ matchedCount +")");
    jqueryElements.unmatchedPeakCount.html("Non Matched Peaks ("+unMatchedCount +")");
}

/**
 * Function to show only matched on unmatched peaks on click of matched or unmatched peak buttons
 * @param {String} id - Contains Ids of respective matched peaks or un matched peaks
 */
function showIonPeaks(id) {
	//console.log("ids : ", ids);
    let elems = domElements.matchedPeaks;
    for(let i = 0; i < elems.length; i++) {
        elems[i].style.display = 'none';
    }
    elems = domElements.unmatchedPeaks;
    for(let i = 0; i < elems.length; i++) {
        elems[i].style.display = 'none';
    }

    elems = document.getElementsByName(id);
    for(let j = 0; j < elems.length; j++) {
        elems[j].style.display = "";
        /*elems[j].style.background = "#BEECFF";*/
    }
}