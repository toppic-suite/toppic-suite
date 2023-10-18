//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/sample_feature_writer.hpp"

namespace toppic {

namespace sample_feature_writer {

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "Sample ID" << "\t"
      << "Feature ID" << "\t"
      << "Monoisotopic mass" << "\t"
      << "Intensity" << "\t"
      << "Min time (min)" << "\t"
      << "Max time (min)" << "\t"
      << "Min scan" << "\t"
      << "Max scan" << "\t"
      << "Min charge" << "\t"
      << "Max charge" << "\t"
      << "Apex time (min)" << "\t"
      << "Apex scan" << "\t"
      << "Apex intensity" << "\t"
      << "Rep charge" << "\t"
      << "Rep average mz" << "\t"
      << "Envelope number" << "\t"
      << "EC score" << "\t"
      << "Min fraction ID" << "\t"
      << "Max fraction ID" << "\t"
      << "Elution length (sec)" 
      << std::endl;
}

void writeOneFeature(std::ofstream &of, SampleFeaturePtr feature) {
  of << feature->getSampleId() << "\t"
      << feature->getId() << "\t"
      << feature->getMonoMass() << "\t"
      << feature->getIntensity() << "\t"
      << feature->getTimeBegin()/60 << "\t"
      << feature->getTimeEnd()/60 << "\t"
      << feature->getMinScan() << "\t"
      << feature->getMaxScan() << "\t"
      << feature->getMinCharge() << "\t"
      << feature->getMaxCharge() << "\t"
      << feature->getApexTime()/60 << "\t"
      << feature->getApexScan() << "\t"
      << feature->getApexInte() << "\t"
      << feature->getRepCharge() << "\t"
      << feature->getRepAvgMz() << "\t"
      << feature->getEnvNum() << "\t"
      << feature->getEcScore() << "\t"
      << feature->getMinFracId() << "\t"
      << feature->getMaxFracId() << "\t" 
      << feature->getElutionLen()
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
