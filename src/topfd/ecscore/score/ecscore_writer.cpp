//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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

#include "topfd/ecscore/score/ecscore_writer.hpp"

namespace toppic {

namespace ecscore_writer {

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
    << "MapMaxElutionTime" << "\t"
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

void writeOneScore(std::ofstream &of, ECScorePtr score) {
  of << score->getMinScan() << "\t"
    << score->getMaxScan() << "\t"
    << score->getMinCharge() << "\t"
    << score->getMaxCharge() << "\t"
    << score->getMonoMass() << "\t"
    << score->getRepCharge() << "\t"
    << score->getRepMz() << "\t"
    << score->getAbundance() << "\t"
    << score->getMinElutionTime() << "\t"
    << score->getMaxElutionTime() << "\t"
    << score->getApexElutionTime() << "\t"
    << score->getElutionLength() << "\t"
    << score->getMapMaxElutionTime() << "\t"
    << score->getEnvcnnScore() << "\t"
    << score->getPercentMatchedPeaks() << "\t"
    << score->getIntensityCorrelation() << "\t"
    << score->getTop3Correlation() << "\t"
    << score->getEvenOddPeakRatio() << "\t"
    << score->getPercentConsecPeaks() << "\t"
    << score->getNumTheoPeaks() << "\t"
    << score->getMzErrorSum() << "\t"
    << score->getScore() << "\t"
    << score->getLabel()
    << std::endl;
}

void writeScores(const std::string &output_file_name, 
                 ECScorePtrVec scores) {
  std::ofstream of(output_file_name);
  writeHeader(of);
  for (auto score: scores) {
    if (score->getElutionLength() > 0)
      writeOneScore(of, score);
  }
  of.close();
}

}
}
