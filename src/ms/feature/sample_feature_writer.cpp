//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <set>
#include <algorithm>
#include <iomanip>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/sample_feature_writer.hpp"

namespace toppic {

namespace sample_feature_writer {

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "Sample_ID" << "\t"
      << "Feature_ID" << "\t"
      << "Monoisotopic_mass" << "\t"
      << "Intensity" << "\t"
      << "Min_time" << "\t"
      << "Max_time" << "\t"
      << "Min_scan" << "\t"
      << "Max_scan" << "\t"
      << "Min_charge" << "\t"
      << "Max_charge" << "\t"
      << "Min_fraction_ID" << "\t"
      << "Max_fraction_ID" << "\t"
      << "Rep_fraction_ID" << "\t"
      << "Rep_Apex_time" << "\t"
      << "Rep_Apex_scan" << "\t"
      << "Rep_Apex_intensity" << "\t"
      << "Rep_charge" << "\t"
      << "Rep_average_mz" << "\t"
      << "Envelope_number" << "\t"
      << "Rep_EC_score" << "\t"
      << "Elution_length" 
      << std::endl;
}

void writeOneFeature(std::ofstream &of, SampleFeaturePtr feature) {
  of << feature->getSampleId() << "\t"
      << feature->getFeatId() << "\t"
      << std::fixed << std::setprecision(6) 
      << feature->getMonoMass() << "\t"
      << feature->getIntensity() << "\t"
      << std::setprecision(2) 
      << feature->getTimeBegin()/60 << "\t"
      << feature->getTimeEnd()/60 << "\t"
      << std::setprecision(0) 
      << feature->getMinScan() << "\t"
      << feature->getMaxScan() << "\t"
      << feature->getMinCharge() << "\t"
      << feature->getMaxCharge() << "\t"
      << feature->getMinFracId() << "\t"
      << feature->getMaxFracId() << "\t" 
      << feature->getRepFracId() << "\t"
      << std::setprecision(2) 
      << feature->getRepApexTime()/60 << "\t"
      << std::setprecision(0) 
      << feature->getRepApexScan() << "\t"
      << std::setprecision(6)
      << feature->getRepApexInte() << "\t"
      << std::setprecision(0) 
      << feature->getRepCharge() << "\t"
      << std::setprecision(6)
      << feature->getRepAvgMz() << "\t";
  of << feature->getEnvNum() << "\t"
      << feature->getRepEcScore() << "\t"
      << std::setprecision(4) 
      << feature->getElutionLen()/60
      << std::setprecision(6) 
      << std::endl;
}

void writeFeatures(const std::string &output_file_name,
                   const SampleFeaturePtrVec &features) {
  std::ofstream of(output_file_name);
  writeHeader(of);

  for (size_t i = 0; i < features.size(); i++) {
    SampleFeaturePtr feature = features[i];
    writeOneFeature(of, feature);
  }
  of.close();
}

}

} /* namespace toppic */
