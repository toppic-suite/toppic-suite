//
// Created by abbash on 8/30/22.
//

#include "write_feature.hpp"

namespace toppic {
namespace write_feature {

  void writeHeader(std::ofstream &of) {
    of.precision(16);
    of << "FeatureID" << "\t"
       << "MinScan" << "\t"
       << "MaxScan" << "\t"
       << "MinCharge" << "\t"
       << "MaxCharge" << "\t"
       << "MonoMass" << "\t"
       << "RepCharge" << "\t"
       << "RepMz" << "\t"
       << "Abundance" << "\t"
       << "MinElutionTime" << "\t"
       << "MaxElutionTime" << "\t"
       << "ApexElutionTime" << "\t"
       << "ElutionLength" << "\t"
       << "EnvCNNScore" << "\t"
       << "PercentMatchedPeaks" << "\t"
       << "IntensityCorrelation" << "\t"
       << "Top3Correlation" << "\t"
       << "EvenOddPeakRatios" << "\t"
       << "PercentConsecPeaks" << "\t"
       << "NumTheoPeaks" << "\t"
       << "MzErrorSum" << "\t"
       << "MzErrorSumBase" << "\t"
       << "Score" << "\t"
       << "Label"
       << std::endl;
  }

  void writeOneFeature(std::ofstream &of, Feature feature) {
    of << feature.getFeatureId() << "\t"
       << feature.getMinScan() << "\t"
       << feature.getMaxScan() << "\t"
       << feature.getMinCharge() << "\t"
       << feature.getMaxCharge() << "\t"
       << feature.getMonoMass() << "\t"
       << feature.getRepCharge() << "\t"
       << feature.getRepMz() << "\t"
       << feature.getAbundance() << "\t"
       << feature.getMinElutionTime() << "\t"
       << feature.getMaxElutionTime() << "\t"
       << feature.getApexElutionTime() << "\t"
       << feature.getElutionLength() << "\t"
      << feature.getEnvcnnScore() << "\t"
      << feature.getPercentMatchedPeaks() << "\t"
      << feature.getIntensityCorrelation() << "\t"
      << feature.getTop3Correlation() << "\t"
      << feature.getEvenOddPeakRatios() << "\t"
      << feature.getPercentConsecPeaks() << "\t"
      << feature.getNumTheoPeaks() << "\t"
      << feature.getMzErrorSum() << "\t"
      << feature.getScore() << "\t"
       << feature.getLabel()
       << std::endl;
  }

  void writeFeatures(const std::string &output_file_name, const std::vector<Feature> &features) {
    std::ofstream of(output_file_name);
    writeHeader(of);
    for (auto feature : features)
      writeOneFeature(of, feature);
    of.close();
  }
}
}