/**
 * Function to create table
 */
function createMonoMassTable(){
    // Remove if table already exist and rebuild the table
    // $("#tableContainer").remove();
    // Remove if tableContainer_wrapper already exist and rebuild the table, wrapper is an inbuild class of bootstrap table
    jqueryElements.monoMassTableContainer.remove();
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
        if(i === 0) th.innerHTML = "Id";
        if(i === 1) th.innerHTML = "Mono mass";
        if(i === 2) th.innerHTML = "Charge";
        if(i === 3) th.innerHTML = "Mono m/z";
        if(i === 4) th.innerHTML = "Intensity";
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
}

/**
 * Add Data to the table created
 * @param {Array} matchedPeaks - Contains List of Matched and unmatched peaks
 */
function addMassDataToTable(matchedPeaks)
{
    let dataContainer_tbody = jqueryElements.monoMassTableBody;
    const totalColCount = jqueryElements.monoMassTableColumns.length;
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
            else if(j === 2) td.innerHTML = matchedPeaks[i].charge;
            else if(j === 3) {
                let mz = matchedPeaks[i].mass / matchedPeaks[i].charge + 1.007276466879;
                let a = document.createElement('a');
                a.href="#!"
                a.className = "peakRows"
                a.innerHTML = mz.toFixed(4); 
                td.appendChild(a);
            }
            else if(j === 4) td.innerHTML = matchedPeaks[i].intensity ;
            else if(j === 5){
                td.className = "th_mass";
                td.innerHTML = matchedPeaks[i].thMass;
            } 
            else if(j === 6) td.innerHTML = matchedPeaks[i].ion;
            else if(j === 7) td.innerHTML = matchedPeaks[i].position;
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

    jqueryElements.peakRow.click(() => {
        /*	get Mono M/z value till 3 decimal values	*/
        let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
        let graphFeatures = new GraphFeatures();
        ms2_graph.redraw(peak_value,graphFeatures);
        // console.log("completeCalData : ", completeCalData);
        let parent_id  = $(this).parent().parent().prop('id');
        // console.log("parent_id : ",parent_id);
        let th_mass_val = $("#"+parent_id+" .th_mass").text();
        // console.log("th_mass_val : ",th_mass_val);
        let monoMassList = completeCalData.monomasslist;
        generateMonoMassGraph(monoMassList,th_mass_val);
    });
}

/**
 * Function to display matched count and un-matched count
 */
function showPeakCounts()
{
    let matched_elems = domElements.matchedPeaks;
    let unmatche_elems = domElements.unmatchedPeaks;
    let totalCount = matched_elems.length + unmatche_elems.length;
    let matchedCount = matched_elems.length;
    let unMatchedCount = unmatche_elems.length;
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
        elems[j].style.background = "#BEECFF";
    }
}