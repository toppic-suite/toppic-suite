function generateCorrespondingGraph(peakDataList,distributionList,prec_mz){
    let graphParams = graphOptions();
    let ionData = null;
    ms2_graph = addSpectrum("spectrum",peakDataList,distributionList,prec_mz,ionData,graphParams);
}