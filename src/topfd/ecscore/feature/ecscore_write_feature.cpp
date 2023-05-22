//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <fstream>

#include "topfd/ecscore/feature/ecscore_write_feature.hpp"

namespace toppic {

namespace ecscore_write_feature {

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
    << "MapMaxElutionTime" << "\t"
    << "ElutionLength" << "\t"
    << "EnvCNNScore" << "\t"
    << "PercentMatchedPeaks" << "\t"
    << "IntensityCorrelation" << "\t"
    << "Top3Correlation" << "\t"
    << "EvenOddPeakRatios" << "\t"
    << "PercentConsecPeaks" << "\t"
    << "NumTheoPeaks" << "\t"
    << "MzErrorSum" << "\t"
    << "Score" << "\t"
    << "Label"
    << std::endl;
}

void writeOneFeature(std::ofstream &of, FeaturePtr feature) {
  of << feature->getFeatureId() << "\t"
    << feature->getMinScan() << "\t"
    << feature->getMaxScan() << "\t"
    << feature->getMinCharge() << "\t"
    << feature->getMaxCharge() << "\t"
    << feature->getMonoMass() << "\t"
    << feature->getRepCharge() << "\t"
    << feature->getRepMz() << "\t"
    << feature->getAbundance() << "\t"
    << feature->getMinElutionTime() << "\t"
    << feature->getMaxElutionTime() << "\t"
    << feature->getApexElutionTime() << "\t"
    << feature->getElutionLength() << "\t"
    << feature->getMapMaxElutionTime() << "\t"
    << feature->getEnvcnnScore() << "\t"
    << feature->getPercentMatchedPeaks() << "\t"
    << feature->getIntensityCorrelation() << "\t"
    << feature->getTop3Correlation() << "\t"
    << feature->getEvenOddPeakRatio() << "\t"
    << feature->getPercentConsecPeaks() << "\t"
    << feature->getNumTheoPeaks() << "\t"
    << feature->getMzErrorSum() << "\t"
    << feature->getScore() << "\t"
    << feature->getLabel()
    << std::endl;
}

void writeFeatures(const std::string &output_file_name, 
                   FeaturePtrVec features) {
  std::ofstream of(output_file_name);
  writeHeader(of);
  for (auto feature: features) {
    if (feature->getElutionLength() > 0)
      writeOneFeature(of, feature);
  }
  of.close();
}

}
}
