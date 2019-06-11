//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include "feature/sample_feature_writer.hpp"

namespace toppic {

namespace sample_feature_writer {

void writeHeader(std::ofstream &of) {
  of.precision(16);
  of << "Sample_ID" << "\t"
      << "ID" << "\t"
      << "File_name" << "\t"
      << "Mass" << "\t"
      << "Intensity" << "\t"
      << "Time_begin" << "\t"
      << "Time_end" << "\t"
      << "Minimum_charge_state" << "\t"
      << "Maximum_charge_state" << "\t"
      << "Minimum_fraction_id" << "\t"
      << "Maximum_fraction_id"
      << std::endl;
}

void writeOneFeature(std::ofstream &of, SampleFeaturePtr feature) {
  of << feature->getSampleId() << "\t"
      << feature->getId() << "\t"
      << feature->getMonoMass() << "\t"
      << feature->getIntensity() << "\t"
      << feature->getTimeBegin() << "\t"
      << feature->getTimeEnd() << "\t"
      << feature->getMinCharge() << "\t"
      << feature->getMaxCharge() << "\t"
      << feature->getMinFracId() << "\t"
      << feature->getMaxFracId() 
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
