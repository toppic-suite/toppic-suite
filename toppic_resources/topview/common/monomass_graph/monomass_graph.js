class MonomassGraph {
    prsm_data;
    nav_id;
    div_id;

    constructor (prsm_data, nav_id = "monoMass_nav", div_id = "monomass_div") {
        this.prsm_data = prsm_data;
        this.nav_id = nav_id;
        this.div_id = div_id;
    }

    drawGraph() {
        let MultiScanObj = new MultiScan();
        let ms2_ids = this.prsm_data.prsm.ms.ms_header.ids;
        let ms2_id_list = ms2_ids.split(" ");
        // Get scan Ids of MS2 Spectrum
        let scanIdList = this.prsm_data.prsm.ms.ms_header.scans.split(" ");
        scanIdList = MultiScanObj.getUniqueScanIdList(scanIdList);
        let ms2_uniqueList = MultiScanObj.getUniqueScanIdList(ms2_id_list);
        
        // Add tabs of scan Ids for Mono Mass Spectrum graph
        MultiScanObj.createMonoMassNavEements(scanIdList,this.nav_id);
        let dataWithScanIdList = [];
        ms2_uniqueList.forEach((element,i) => {
            let tempObj = {"specId":element, "scanId":scanIdList[i], value:[]};
            dataWithScanIdList.push(tempObj);
        });

        console.log("dataWithScanIdList:", dataWithScanIdList);

        let len = ms2_uniqueList.length;
        this.prsm_data.prsm.ms.peaks.peak.forEach(function(element,i){
            let ion = "";
            if(element.hasOwnProperty('matched_ions_num'))
            {   
                ion = element.matched_ions.matched_ion.ion_type + element.matched_ions.matched_ion.ion_display_position;
            }
            let tempObj = {"mz":parseFloat(element.monoisotopic_mass),"intensity":parseFloat(element.intensity),"ion":ion};
            for(let i=0;i<len;i++)
            {
                if(dataWithScanIdList[i].specId == element.spec_id)
                {
                    dataWithScanIdList[i].value.push(tempObj);
                }
            }
        })

        createMultipleSvgs(this.div_id, "monoMassSvg_","monoMass_svg_graph_class",dataWithScanIdList);

        // add onclick events to switch nav tabs
        let nav_id = this.nav_id;
        $(".monoMass_scanIds").click(function(){
            let value = this.getAttribute('value');
            let id = "monoMassSvg_"+value;
            // Hide all the graphs except the one clicked
            showCorrespondingGraph(id,".monoMass_svg_graph_class");
            $("#"+nav_id+" .active").removeClass("active");
            $(this).addClass("active");
        })
    }
}