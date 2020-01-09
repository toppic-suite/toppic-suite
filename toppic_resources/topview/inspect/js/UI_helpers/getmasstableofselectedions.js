/**
 * Class contains {createTableForSelectedFragmentIons}
 */
class GetMassTableOfSelectedIons{

    /**
     * Function to form table for all selected fragmented ions
     * @param {string} sequence - user entered sequence without mass shifts embedded
     * @param {Array} matchedUnMatchedPeaks - list of all the calculated masses
     */
    createTableForSelectedFragmentIons(sequence,matchedUnMatchedPeaks){
        /**
         * Remove if table already exist and rebuild the table
         */
        $("#selectedIonTableContainer").remove();
        $("#divselectediontablecontainer #selectedIonTableContainer_wrapper").remove();
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
            /**
             * Add mass data to respected columns 
             */
            for(let k=0;k<len;k++)
            {
                let td1 = document.createElement("td");
                let index = j;
                if(matchedUnMatchedPeaks[k].ionFragment[0] == "x" || matchedUnMatchedPeaks[k].ionFragment[0] == "y"
                    || matchedUnMatchedPeaks[k].ionFragment[0] == "z")
                {
                    /**
                     * position when suffix mass list is writtent o table
                     */
                    index = seqln-j-1;
                }
                td1.setAttribute("class","td_fragments");
                if(matchedUnMatchedPeaks[k].massList[index].matchedInd == "Y")
                {
                    td1.setAttribute("id","matched_fragments");
                }
                td1.innerHTML = matchedUnMatchedPeaks[k].massList[index].mass.toFixed(4);
                tr1.appendChild(td1);
            }
            tbody.appendChild(tr1);
        }
        
        thead.appendChild(tr);
        table.appendChild(thead);
        table.appendChild(tbody);
        div.appendChild(table);
    }
}
