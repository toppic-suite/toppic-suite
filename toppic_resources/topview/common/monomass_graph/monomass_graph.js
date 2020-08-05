class MonomassGraph {
    // prsm_data;
    navID;
    divID;
    specIDList;
    scanIDList;
    peakList;

    constructor (prsm_data, navID = "monoMass_nav", divID = "monomass_div") {
        this.prsm_data = prsm_data;
        this.navID = navID;
        this.divID = divID;
        this.specIDList = prsm_data.prsm.ms.ms_header.ids.split(" ");
        this.scanIDList = prsm_data.prsm.ms.ms_header.scans.split(" ");
        this.peakList = prsm_data.prsm.ms.peaks.peak;
    }

    drawGraph() {
        let specIDList = this.getUniqueList(this.specIDList);
        let scanIDList = this.getUniqueList(this.scanIDList);

        this.createMonoMassNavElements(scanIDList, this.navID);

        let dataWithScanIdList = [];

        specIDList.forEach((element,i) => {
            let tempObj = {"specId":element, "scanId":scanIDList[i], value:[]};
            dataWithScanIdList.push(tempObj);
        });

        this.peakList.forEach(function(element){
            let ion = "";
            if(element.hasOwnProperty('matched_ions_num'))
            {   
                ion = element.matched_ions.matched_ion.ion_type + element.matched_ions.matched_ion.ion_display_position;
            }
            let tempObj = {"mz":parseFloat(element.monoisotopic_mass),"intensity":parseFloat(element.intensity),"ion":ion};
            for(let i=0; i < specIDList.length; i++)
            {
                if(dataWithScanIdList[i].specId === element.spec_id || dataWithScanIdList[i].specId === parseInt(element.spec_id, 10))
                {
                    dataWithScanIdList[i].value.push(tempObj);
                }
            }
        })

        this.createMultipleSvgs(this.divID, "monoMassSvg_","monoMass_svg_graph_class",dataWithScanIdList);

        // add onclick events to switch nav tabs
        let navID = this.navID;
        let self = this;
        $(".monoMass_scanIds").click(function(){
            let value = this.getAttribute('value');
            let id = "monoMassSvg_"+value;
            // Hide all the graphs except the one clicked
            self.showCorrespondingGraph(id,".monoMass_svg_graph_class");
            $("#"+navID+" .active").removeClass("active");
            $(this).addClass("active");
        })
    }

    getUniqueList(multiScanList){
        let uniqueList = [];
        let uniqueIdSet = new Set();
        multiScanList.forEach(function(element){
            uniqueIdSet.add(element);
        });
        //spread operator converts array to list
        uniqueList = [...uniqueIdSet];
        return uniqueList ;
    }

    createMonoMassNavElements(scanIDList,id){
        let _ul = document.getElementById(id);
        scanIDList.forEach(function(element,i){
            let li = document.createElement("li");
            li.setAttribute("class","nav-item");
            let li_id = id+"_"+element;
            li.setAttribute("id",li_id);
            let a = document.createElement("a");
            a.setAttribute("class","nav-link monoMass_scanIds");
            if(i === 0)
            {
                a.setAttribute("class","nav-link monoMass_scanIds active");
            }
            a.setAttribute("href","#!");
            a.setAttribute("value",element);
            a.innerHTML = "Scan "+ element;
            li.appendChild(a);
            _ul.appendChild(li);
         })
    }

    generateCorrespondingGraph(current_data,id,prec_mz){
        // Gets the svg id of the spectrum
        let startOfId = id.split("_")[0];
        
        if(startOfId == "monoMassSvg")
        {
            let graphFeatures = new GraphFeatures();
    
            let sequenceObj = new Sequence(this.prsm_data);
            let sequence = sequenceObj.getSequence();
    
            let ionMassShiftObj = new IonMassShift("B");
            let ionShift = ionMassShiftObj.getIonTypeMass();
    
            let massShiftListObj = new MassShiftList(this.prsm_data);
            let massShiftList = massShiftListObj.getMassShiftList();
    
            let cpsmObj = new CalcPrefixSuffixMassList(sequence, massShiftList);
            let prefixMassList = cpsmObj.getPrefixMassList(ionShift);
    
            ionMassShiftObj.ionType = "Y";
            ionShift = ionMassShiftObj.getIonTypeMass();
    
            let suffixMassList = cpsmObj.getSuffixMassList(ionShift);
    
            // Setting the graphFeatures object with all the features needed for the graph
            graphFeatures.showSequence = true;
            graphFeatures.addErrorPlot = true;
            graphFeatures.prefixSequenceData = prefixMassList;
            graphFeatures.suffixSequeceData = suffixMassList;
            graphFeatures.svgHeight = graphFeatures.svgHeight + graphFeatures.adjustableHeightVal + graphFeatures.heightForErrorPlot;
            graphFeatures.padding.head = graphFeatures.padding.head + graphFeatures.adjustableHeightVal;
            graphFeatures.padding.bottom = graphFeatures.padding.bottom + graphFeatures.heightForErrorPlot;
            graphFeatures.adjustableIonPosition = 10; // Random tested value for alignment
            // Gets the data list with mass error to plot in the monomass spectrum 
            let prsmDataUtilObj = new PrsmDataUtil(prsm_data);
            graphFeatures.errorListData = prsmDataUtilObj.json2ErrorDataList(); 
            // Gets the absolute max and minimum value for upper bound and lower bound of y axis to draw the error plot
            graphFeatures.errorThreshHoldVal = this.getAbsoluteMaxValfromList(graphFeatures.errorListData);
            // Invoking spectrum function to draw the spectrum
            let spectrumgraph = new addSpectrum(id, current_data, null, prec_mz, current_data,graphFeatures);
        }
    }

    createMultipleSvgs(divId,svgId,className,dataWithScanIdList){
        let div = document.getElementById(divId); 
        let self = this;
        dataWithScanIdList.forEach(function(element,i){
            let id = svgId+element.scanId;
            let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");;
            svg.setAttribute("id",id);
            svg.setAttribute("class",className);
            svg.style.backgroundColor = "#F8F8F8"; 
            if(i != 0)
            {
                svg.style.display = "none"; 
            }
            div.appendChild(svg);
            self.generateCorrespondingGraph(element.value, id, null);
        });
    }

    showCorrespondingGraph(id,className)
    {
        $(className).hide();
        document.getElementById(id).style = "block";
    }
    
    getAbsoluteMaxValfromList(errorDataList){
        let max = 0;
        errorDataList.forEach((element)=>{
            let val = Math.abs(element.mass_error);
            if(max < val) max = val; 
        })
        //Getting the round off fraction value
        max = max * 100;
        max = Math.ceil(max)/100;
        return max;
    }
}