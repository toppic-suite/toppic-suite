function monomassGraph(prsm_data) {
    let MultiScanObj = new MultiScan();
    let ms2_ids = prsm_data.prsm.ms.ms_header.ids;
    let ms2_id_list = ms2_ids.split(" ");
    // Get scan Ids of MS2 Spectrum
    let scanIdList = prsm_data.prsm.ms.ms_header.scans.split(" ");
    scanIdList = MultiScanObj.getUniqueScanIdList(scanIdList);
    let ms2_uniqueList = MultiScanObj.getUniqueScanIdList(ms2_id_list);
    // Add tabs of scan Ids for Mono Mass Spectrum graph
    // This function also adds data into monoMassDataList 
    console.log("ms2_uniqueList", ms2_uniqueList);
    console.log("scanIdList", scanIdList);
    getMonoMassDataList(ms2_uniqueList,scanIdList);
}