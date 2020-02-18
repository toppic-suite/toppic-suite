function drawMonoMassSpectrum(specidList,scanIdList){
    let MultiScanObj = new MultiScan();
    MultiScanObj.createMonoMassNavEements(scanIdList,"monoMass_nav");
    dataWithScanIdList = [];
    specidList.forEach((element,i) => {
        let tempObj = {"specId":element,"scanId":scanIdList[i],value:[]};
        dataWithScanIdList.push(tempObj);
    });
    let len = specidList.length;
    prsm_data.prsm.ms.peaks.peak.forEach(function(element,i){
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
    document.getElementById("monoMassDataLoading").remove();
    createMultipleSvgs("monomass_div","monoMassSvg_","monoMass_svg_graph_class",dataWithScanIdList);
    monoMassDataList = [...dataWithScanIdList];
}