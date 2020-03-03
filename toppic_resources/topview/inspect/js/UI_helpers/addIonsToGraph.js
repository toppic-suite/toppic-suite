function generateCorrespondingGraph(peakDataList,distributionList,prec_mz){
    let ionData = null;
    let graphFeatures = new GraphFeatures();
    ms2_graph = addSpectrum("spectrum",peakDataList,distributionList,prec_mz,ionData,graphFeatures);
}