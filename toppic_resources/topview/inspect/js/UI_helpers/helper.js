/**
 * Class to help the get and put data into html file.
 * This also helps to create HTML tags like creating table
 */
class UIHelper{
    /**
     * Function to set default error value to HTML in Da units
     * @param {Float} massErrorthVal - Contains Mass Error Value in Da units
     */
    setMassErrorValue(massErrorthVal){
        $("#errorval").val(massErrorthVal);
        $("#errorunit").html("Da&nbsp;&nbsp;");
    }
    /**
     * Function to Set Default error value to HTML in ppm units
     * @param {Float} ppmErrorthVal - Contains Mass Error in ppm units
     */
    setPPMErrorValue(ppmErrorthVal){
        $("#errorval").val(ppmErrorthVal);
        $("#errorunit").html("ppm&nbsp;&nbsp;");
    }
    /**
     * Set Total mass on to the html
     * @param {*} totalMass 
     */
    setTotalSeqMass(totalMass){
        totalMass = totalMass.toFixed(4);
        $("#totalmass").html(totalMass);
    }
    /**
     * Set Mass difference on to the html
     * @param {Float} precursorMass - Contains Precursor mass
     * @param {Float} proteinMass - Contains calculated protein Mass
     */
    setMassDifference(precursorMass, proteinMass){
        let diff = proteinMass - precursorMass ;
        document.getElementById("massvariation").innerHTML = diff.toFixed(4);
        return (proteinMass - precursorMass);
    }
    /**
     * Set Default Mass errors into html of both Da and ppm units 
     * @param {Float} massErrorthVal - Contains Mass error in Da units
     * @param {Float} ppmErrorthVal - Contains ppm error in ppm units
     */
    writeMassErrorThreshholdValueToUI(massErrorthVal,ppmErrorthVal){
        if(massErrorthVal != "") $("#errorval").val(massErrorthVal);
        else $("#errorval").val(ppmErrorthVal);
    }
    /**
     * getlist of Mass,Intensity and charge from UI
     */
    getMassListFromUI()
    {
        let spectrumDataList = [];
        // Read data line by line from the mass and intensity box
        var lines = $('#data').val().split('\n');
        for(var i = 0; i < lines.length;i++){
            let massAndInte = lines[i].trim() ;
            if(massAndInte.length !=  0 )
            {
                let spectrumData = {};
                // Get Mass,intensity and charge either by space seperated or tab seperated
                let massInte = massAndInte.split(/[\s]+/);
                if(massInte[0] != undefined && massInte[1] != undefined 
                        && !isNaN(massInte[0]) && !isNaN(massInte[1]))
                {
                    spectrumData.mass = parseFloat(massInte[0]);
                    spectrumData.intensity = parseFloat(massInte[1]);
                    spectrumData.charge = parseFloat(massInte[2]);
                }
                if(!jQuery.isEmptyObject(spectrumData))
                {
                    if((spectrumData.mass !== undefined) &&
                        (spectrumData.intensity !== undefined))
                    {
                        spectrumDataList.push(spectrumData) ;
                    }
                }
            }
        }
        completeCalData.monomasslist = spectrumDataList;
        return spectrumDataList ;
    }
    /**
     * Function to get data of peaks and intensity from UI
     */
    getPeakListFromUI()
    {
        let spectrumDataList = [];
        // Read data line by line from peak and intensity box
        var lines = $('#peakdata').val().split('\n');
        for(var i = 0; i < lines.length;i++){
            let peakAndInte = lines[i].trim();
            if(peakAndInte.length !=  0 )
            {
                let spectrumData = {};
                let peakInte = peakAndInte.split(/[\s]+/);
                if(peakInte[0] != undefined && peakInte[1] != undefined 
                        && !isNaN(peakInte[0]) && !isNaN(peakInte[1]))
                {
                    spectrumData.mz = parseFloat(peakInte[0])///spectrumData.charge ;
                    spectrumData.intensity = parseFloat(peakInte[1]);
                }
                if(!jQuery.isEmptyObject(spectrumData))
                {
                    if((spectrumData.mz !== undefined) &&
                        (spectrumData.intensity !== undefined))
                    {
                        spectrumDataList.push(spectrumData) ;
                    }
                }
            }
        }	
        completeCalData.peakdatalist = spectrumDataList;
	    return spectrumDataList ;
    }
    /**
     * Function to create table
     */
    createMonoMassTable(){
        // Remove if table already exist and rebuild the table
        $("#tableContainer").remove();
        // Remove if tableContainer_wrapper already exist and rebuild the table, wrapper is an inbuild class of bootstrap table
        $("#divtableContainer #tableContainer_wrapper").remove();
        let div = document.getElementById("divtableContainer");
        let table = document.createElement("table");
        table.setAttribute("id","tableContainer");
        table.setAttribute("class","table table-striped display");
        let thead = document.createElement("thead");
        let tbody = document.createElement("tbody");
        tbody.setAttribute("id","tableContainertbody");
        let tr = document.createElement("tr");
        tr.setAttribute("role","row");
        let ColCount = 10;
        // Create table header
        for(let i = 0;i < ColCount;i++)
        {
            let th = document.createElement("th");
            th.setAttribute("class","th-sm");
            if(i == 0) th.innerHTML = "Id";
            if(i == 1) th.innerHTML = "Mono mass";
            if(i == 2) th.innerHTML = "Charge";
            if(i == 3) th.innerHTML = "Mono m/z";
            if(i == 4) th.innerHTML = "Intensity";
            if(i == 5) th.innerHTML = "Theoretical mass";
            if(i == 6) th.innerHTML = "Ion";
            if(i == 7) th.innerHTML = "Pos";
            if(i == 8) th.innerHTML = "Mass error";
            if(i == 9) th.innerHTML = "PPM error";
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
    addMassDataToTable(matchedPeaks)
    {
        let dataContainer_tbody = $("#tableContainer tbody");
        const totalColCount = $("#tableContainer thead tr th").length;
        let len = matchedPeaks.length;
        for(let i=0; i<len ; i++)
        {
            let rowSize = $("#tableContainer tbody tr").length;
            let tr = document.createElement("tr");
            tr.setAttribute("id",i+"_row");
            tr.setAttribute("name",matchedPeaks[i].position);
            let id = 0;
            for(let j = 0; j < totalColCount ; j++)
            {
                let td = document.createElement("td");
                if(j == 0)
                {
                    id = matchedPeaks[i].peakId;
                    td.innerHTML = matchedPeaks[i].peakId;
                    td.style.fontWeight = "bold";
                }else if(j == 1) td.innerHTML = matchedPeaks[i].mass ;
                else if(j == 2) td.innerHTML = matchedPeaks[i].charge;
                else if(j == 3) {
                  let mz = matchedPeaks[i].mass / matchedPeaks[i].charge + 1.007276466879;
                  let a = document.createElement('a');
                  a.href="#!"
                  a.className = "peakRows"
                  a.innerHTML = mz.toFixed(4); 
                  td.appendChild(a);
                }
                else if(j == 4) td.innerHTML = matchedPeaks[i].intensity ;
                else if(j == 5){
                    td.className = "th_mass";
                    td.innerHTML = matchedPeaks[i].thMass;
                } 
                else if(j == 6) td.innerHTML = matchedPeaks[i].ion;
                else if(j == 7) td.innerHTML = matchedPeaks[i].position;
                else if(j == 8) td.innerHTML = matchedPeaks[i].massError;
                else if(j == 9) td.innerHTML = matchedPeaks[i].PPMerror;
                tr.appendChild(td);
            }
            tr.setAttribute("role","row");
            let classname = "" ;
            // Creating classes with matched_peak even and odd, this will help to show only matched peaks on click of matched peaks
            if(matchedPeaks[i].matchedInd == "Y")
            {
                if(id%2 == 0) classname = "matched_peak even";
                else classname = "matched_peak odd";
            }
            else
            {
                if(id%2 == 0) classname = "unmatched_peak even";
                else classname = "unmatched_peak odd";
            }
            tr.setAttribute("class",classname);
            dataContainer_tbody.append(tr);
        }

      $(".peakRows").click(function() {
        /*	get Mono M/z value till 3 decimal values	*/
        let peak_value = parseFloat(this.innerHTML).toFixed(3) ;
        let graphFeatures = new GraphFeatures();
        ms2_graph.redraw(peak_value,graphFeatures);
        console.log("completeCalData : ", completeCalData);
        let parent_id  = $(this).parent().parent().prop('id');
        console.log("parent_id : ",parent_id);
        let th_mass_val = $("#"+parent_id+" .th_mass").text();
        console.log("th_mass_val : ",th_mass_val);
        let monoMassList = completeCalData.monomasslist;
        generateMonoMassGraph(monoMassList,th_mass_val);
      });
    }
    /**
     * Function to display matched count and un-matched count
     */
    showPeakCounts()
    {
        var matched_elems = document.getElementsByClassName("matched_peak");
        var unmatche_elems = document.getElementsByClassName("unmatched_peak");
        let totalCount = matched_elems.length + unmatche_elems.length;
        let matchedCount = matched_elems.length ;
        let unMatchedCount = unmatche_elems.length ;
        $("#all_peak_count").html(function(){return "All Peaks ("+totalCount+")";}) 
        $("#matched_peak_count").html(function(){return "Matched Peaks ("+ matchedCount +")";})
        $("#unmatched_peak_count").html(function(){return "Non Matched Peaks ("+unMatchedCount +")";}) 
    }
    /**
     * Function to show only matched or unmatched peaks
     */
    showIonPeaks(ids) 
    {
        var elems = document.getElementsByClassName('matched_peak');
        for(var i = 0; elems.length > i; i++) {
            elems[i].style.display = 'none';
        }
        elems = document.getElementsByClassName('unmatched_peak');
        for(var i = 0; elems.length > i; i++) {
            elems[i].style.display = 'none';
        }

        elems = document.getElementsByName(ids);
            for(var j = 0; elems.length > j; j++) {
            elems[j].style.display  =  "";
            elems[j].style.background  =  "#BEECFF";
            }
    }
}
